#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""Class to represent catalogue of SVs."""


from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pysam


class SVCatalogue:
    """Represents a catalogue of structural variants from a VCF, including
    metadata and normalized variant information.

    Because SV VCFs are not super standardized like SNP VCFs, we caution that
    the implementation of this class is ad hoc and if users plan to co-opt for
    their own reference catalogues, proceed with caution. PySam does not always
    play friendly with the VCFs, so some amount of wrangling will be required.

    Attributes:
      vcf_path: Path to the VCF file.
      source_name: Identifier for callset source.
      priority: Priotiry value for merging (lower is higher priority).
      read_type: Type of sequencing read (short or long).
      variants: List of normalized variant dictionaries.
      remove_contigs: Flag to remove contig names from the VCF.

    Examples::
        # Create an SVCatalogue instance and load variants from a VCF file
        >>> catalogue = SVCatalogue(
        ...     vcf_path="path/to/vcf_file.vcf",
        ...     source_name="example_source",
        ...     priority=priority,
        ...     read_type="short_read"
        ...     source_rules=source_rules,
        ... )
    """

    def __init__(
        self,
        vcf_path: str,
        source_name: str,
        priority: int,
        read_type: str,
        source_rules: Dict[Any, Any],
        remove_contigs: bool = True,
    ) -> None:
        """Initialize the SVCatalogue class."""
        self.vcf_path = vcf_path
        self.source_name = source_name
        self.priority = priority
        self.read_type = read_type

        self.variants: List[Dict] = []
        self.config = source_rules.get(self.source_name, {})

        self._load_variants()
        if remove_contigs:
            self.variants = self._remove_contigs_and_m()

    def _load_variants(self) -> None:
        """Load and normalize variants from the VCF."""
        path = Path(self.vcf_path)
        if not path.exists():
            print(f"[Warning] {self.vcf_path} not found; skipping.")
            return

        with pysam.VariantFile(str(path), "r") as vcf:
            for rec in vcf.fetch():
                parsed = self._parse_record(rec)
                if parsed is not None:
                    self.variants.append(parsed)

    def _parse_record(self, rec: pysam.VariantRecord) -> Optional[Dict]:
        """Parse one record using flags from config. Ensure that insertions
        parse sequence.
        """
        svtype = self._parse_sv_type(rec)
        end_val, size_val = self._compute_end_and_size(rec, svtype)

        new_var = {
            "chrom": rec.chrom,
            "pos": rec.pos,
            "end": end_val,
            "size": size_val,
            "svtype": svtype,
            "source": self.source_name,
            "priority": self.priority,
        }

        if svtype == "INS":
            insertion_seq = self._parse_insertion_seq(rec)
            if insertion_seq is None:
                return None
            new_var["seq"] = insertion_seq

        if self.read_type == "short_read" and not self._filter_short_read_support(rec):
            return None

        return new_var

    def _parse_sv_type(self, rec: pysam.VariantRecord) -> str:
        """Normalize SV type for the record."""
        if "SVTYPE" in rec.info:
            val = rec.info["SVTYPE"]
            if isinstance(val, tuple):
                val = val[0]
            return self._normalize_svtype(str(val))

        if rec.id:
            id_upper = rec.id.upper()
            for candidate in ("DEL->", "INS->", "DUP->", "INV->"):
                if candidate in id_upper:
                    return self._normalize_svtype(candidate.replace("->", ""))

        raise ValueError(
            f"[Error] Cannot parse SVTYPE for {rec.id} at {rec.chrom}:{rec.pos} "
            f"in source {self.source_name}"
        )

    def _get_evidence_flag(self, svtype: str) -> Tuple[bool, bool]:
        """Get trust_end and trust_svlen from config."""
        trust_end = self.config.get("trust_end", False)
        trust_svlen = self.config.get("trust_svlen", False)

        if "exceptions" in self.config and svtype in self.config["exceptions"]:
            override = self.config["exceptions"][svtype]
            trust_end = override.get("trust_end", trust_end)
            trust_svlen = override.get("trust_svlen", trust_svlen)
        return trust_end, trust_svlen

    def _compute_end_and_size(
        self, rec: pysam.VariantRecord, svtype: str
    ) -> Tuple[int, int]:
        """Compute the end coordinate and size of the variant."""
        pos = rec.pos
        trust_end, trust_svlen = self._get_evidence_flag(svtype)

        # Use the END position
        if trust_end:
            if result := self._compute_end_and_size_from_end(rec, pos):
                return result

        # Use SVLEN instead
        if trust_svlen:
            if result := self._compute_end_and_size_from_svlen(rec, pos):
                return result

        # Fallback calculation from pos and REF/ALT
        if svtype == "DEL":
            return self._compute_end_and_size_for_del(rec, pos)

        # Use the ALT allele length
        if result := self._compute_end_and_size_from_alt(rec, pos):
            return result

        raise ValueError(
            f"Cannot determine END for record {rec.id} at {rec.chrom}:{rec.pos} "
            f"with svtype={svtype} from source {self.source_name}."
        )

    def _parse_insertion_seq(self, rec: pysam.VariantRecord) -> Optional[str]:
        """Parse insertion sequence from the record.

        - "alt": means read from the first ALT if it's not symbolic <INS>
        - "info_inseq": means read rec.info["INSEQ"]
        - "none": means we do not store an insertion sequence
                  (or the source doesn't have INS)
        """
        ins_mode = self.config.get("ins_seq_mode", "none")

        if ins_mode == "none":
            return None

        elif ins_mode == "alt":
            if not rec.alts or len(rec.alts) == 0:
                return ""
            alt_allele = rec.alts[0]
            return None if alt_allele.upper().startswith("<INS") else alt_allele

        elif ins_mode == "info_inseq":
            if "INSEQ" not in rec.info:
                return None
            return rec.info["INSEQ"]

        raise ValueError(
            f"[Error] Unknown ins_seq_mode '{ins_mode}' for source='{self.source_name}'."
        )

    def _remove_contigs_and_m(self) -> List[Dict]:
        """Remove contig names from the record by filtering for variants with
        "_" in chrom.
        """
        if not self.variants:
            raise ValueError(
                "[Error] Cannot remove contigs from empty variant list. Run _load_variants first."
            )

        filtered_variants = [
            var
            for var in self.variants
            if "_" not in var["chrom"] and var["chrom"] != "chrM"
        ]
        if len(filtered_variants) == len(self.variants):
            print(f"[Warning] No contigs found in {self.source_name} variants. ")
            return self.variants
        print(
            f"[Info] Removed contigs from {self.source_name} variants. "
            f"Original count: {len(self.variants)}, "
            f"Filtered count: {len(filtered_variants)}."
        )
        return filtered_variants

    @staticmethod
    def _normalize_svtype(raw_svtype: str) -> str:
        """Normalize e.g. <DEL> => 'DEL', <INV> => 'INV'."""
        raw_svtype = raw_svtype.upper().replace("<", "").replace(">", "")
        if raw_svtype.startswith("DEL"):
            return "DEL"
        if raw_svtype.startswith("DUP"):
            return "DUP"
        if raw_svtype.startswith("INS"):
            return "INS"
        return "INV" if raw_svtype.startswith("INV") else raw_svtype

    @staticmethod
    def _filter_short_read_support(
        rec: pysam.VariantRecord, min_support: int = 10
    ) -> bool:
        """Filter short-read variants based on support."""
        if "SU" in rec.info:
            support = rec.info["SU"]
            if isinstance(support, tuple):
                support = support[0]
            return support >= min_support
        return False

    @staticmethod
    def _compute_end_and_size_from_end(
        rec: pysam.VariantRecord, pos: int
    ) -> Optional[Tuple[int, int]]:
        """Compute end and size using END info field."""
        if "END" in rec.info:
            endpos = rec.info["END"]
            return endpos, abs(endpos - pos)

        # try rec.stop
        return (rec.stop, abs(rec.stop - pos)) if rec.stop > pos else None

    @staticmethod
    def _compute_end_and_size_from_svlen(
        rec: pysam.VariantRecord, pos: int
    ) -> Optional[Tuple[int, int]]:
        """Compute end and size using SVLEN info field."""
        if "SVLEN" in rec.info:
            raw_svlen = rec.info["SVLEN"]
            if isinstance(raw_svlen, tuple):
                raw_svlen = raw_svlen[0]
            svlen = abs(int(raw_svlen))
            return (pos + svlen, svlen)
        return None

    @staticmethod
    def _compute_end_and_size_for_del(
        rec: pysam.VariantRecord, pos: int
    ) -> Tuple[int, int]:
        """Compute end and size for DEL SVTYPE."""
        ref_len = len(rec.ref) if rec.ref else 0
        end_pos = pos + ref_len
        return (end_pos, ref_len)

    @staticmethod
    def _compute_end_and_size_from_alt(
        rec: pysam.VariantRecord, pos: int
    ) -> Optional[Tuple[int, int]]:
        """Compute end and size from ALT allele length."""
        if rec.alts and len(rec.alts) == 1:
            alt_allele = rec.alts[0]
            alt_len = len(alt_allele)
            end_pos = pos + alt_len
            return (end_pos, alt_len)
        return None
