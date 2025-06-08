#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""Merge SV callsets to produce a genotype reference panel."""


from typing import Callable, Dict, List, Optional, Union

from intervaltree import IntervalTree  # type: ignore
import pysam  # type: ignore
from tqdm import tqdm  # type: ignore


class SVReferenceMerger:
    """Merges short-read (SR) and long-read (LR) SV callsets into a unified,
    non-redundant reference set using IntervalTrees for efficient lookups.

    This class implements a multi-step pipeline to combine structural variant
    callsets from different technologies (short-read and long-read sequencing).
    The core strategy prioritizes calls derived from long reads and then merges
    variants within each category (SR and LR) based on type-specific overlap and
    proximity criteria.

    The merging pipeline consists of the following key stages:

    1.  Input VCF files are separated into short-read and long-read groups
        (handled externally and passed into the constructor).
    2.  Short-read variants that significantly overlap with long-read variants
        are identified and removed to ensure that high-confidence long-read
        calls take precedence.
    3.  The remaining (non-removed) short-read variants are merged amongst
        themselves using IntervalTrees. Within clusters of overlapping variants
        of the same type, only the highest-priority call (based on source) is
        retained.
    4.  Long-read variants are similarly merged amongst themselves using
        IntervalTrees with their specific overlap criteria and priority rules.
    5.  The resulting set combines the merged long-read variants and the
        filtered, merged short-read variants into a single, comprehensive
        callset.

    Overlap criteria: Filtering short-read SVs that overlap long-read SVs:
     - DEL/DUP/INS (< 5 kbp): >= 10% Reciprocal Overlap (RO) with LR of same
       type.
     - DEL/DUP/INS (>= 5 kbp): >= 50% RO with LR of same type.
     - INV: >= 20% RO with LR INV.
     - INS vs INS special case: SR INS breakpoint within 100 bp of LR INS
     - INS vs DUP special Case: SR INS breakpoint within 100 bp of LR DUP start
       coordinate.

    Merging short-read vs short-read variants OR long-read vs long-read
    variants:
     - INV: >= 20% RO.
     - DEL/DUP/INS (< 5 kbp): Breakpoints within 200 bp AND >= 80% size
       similarity.
     - DEL/DUP/INS (>= 5 kbp): >= 50% RO.

    Attributes:
      short_read_svs: List of short-read variants.
      long_read_svs: List of long-read variants.
      output_vcf: Path to the final merged VCF file.

    Examples::
      # Instantiate the SVReferenceMerger class and run the merge pipeline
      >>> merger = SVReferenceMerger(
          short_read_svs=short_read_svs,
          long_read_svs=long_read_svs,
          ref_fasta=args.reference,
          )
      >>> merger.merge_callsets()
    """

    def __init__(
        self,
        short_read_svs: List[Dict],
        long_read_svs: List[Dict],
        ref_fasta: str,
    ) -> None:
        """Initialize the SVReferenceMerger class.

        Args:
          short_read_svs: List of structural variants derived from short-read
          sequencing.
          long_read_svs: List of structural variants derived from long-read
          sequencing.
          output_vcf: Path to which the final merged set of variants should be
          written.
        """
        self.merged_variants: List[Dict] = []

        self.short_read_svs = short_read_svs
        self.long_read_svs = long_read_svs
        self.ref = pysam.FastaFile(ref_fasta)
        self.longread_trees = self._build_interval_trees(self.long_read_svs)

    def merge_callsets(self) -> None:
        """Run the multi-step merge pipeline, producing a final non-redundant
        reference set of structural variants.
        """
        short_filtered = self._filter_short_read_svs(self.short_read_svs)
        print(
            f"[RefMerger] Removed {len(self.short_read_svs) - len(short_filtered)} "
            f"short-read calls that overlap with long-read calls"
        )

        short_merged = self._merge_svs(short_filtered, self._should_merge_svs)
        print(
            f"[RefMerger] Merged {len(short_filtered) - len(short_merged)} "
            "short-read calls."
        )

        long_merged = self._merge_svs(self.long_read_svs, self._should_merge_svs)
        print(
            f"[RefMerger] Merged {len(self.long_read_svs) - len(long_merged)} "
            "long-read calls."
        )

        self.merged_variants = long_merged + short_merged
        self.merged_variants.sort(key=lambda v: (v["chrom"], v["pos"], v["priority"]))
        print(f"[RefMerger] Final merged panel has {len(self.merged_variants)} calls.")

    def _build_interval_trees(self, variants: List[Dict]) -> Dict[str, IntervalTree]:
        """Build an IntervalTree for each chromosome for the provided variants.
        Each interval is stored as [pos, end).

        Args:
          variants: List of variant dictionaries.

        Returns:
          Dictionary mapping chromosome -> IntervalTree of variants.
        """
        trees: Dict[str, IntervalTree] = {}
        for var in variants:
            chrom = var["chrom"]
            if chrom not in trees:
                trees[chrom] = IntervalTree()
            start, end = var["pos"], var["end"]
            trees[chrom].addi(start, end, var)
        return trees

    def _filter_short_read_svs(self, short_read_svs: List[Dict]) -> List[Dict]:
        """Remove short-read calls that overlap same-type long reads.

        Removal criteria:
         - DEL/DUP/INS (< 5 kbp): >= 10% Reciprocal Overlap (RO) with LR of same
           type.
         - DEL/DUP/INS (>= 5 kbp): >= 50% RO with LR of same type.
         - INV: >= 20% RO with LR INV.
         - INS vs DUP Special Case: SR INS breakpoint within 100 bp of LR DUP
           start coordinate.

        Args:
          short_read_svs: List of short-read variant dictionaries.

        Returns:
          A list of short-read variants that did NOT overlap with long-read
          variants under the above criteria.
        """
        filtered_short: List[Dict] = []
        for variant in tqdm(short_read_svs, desc="Filtering short overlaps"):
            variant_sv_type = variant["svtype"]
            chrom = variant["chrom"]
            start = variant["pos"]
            end = variant["end"]

            remove_it = False
            if chrom in self.longread_trees:
                overlapping_lr_svs = self.longread_trees[chrom].overlap(start, end)
            else:
                overlapping_lr_svs = []

            for interval in overlapping_lr_svs:
                lr_var = interval.data

                if (
                    variant_sv_type == lr_var["svtype"]
                    and variant_sv_type
                    in (
                        "DEL",
                        "DUP",
                        "INS",
                    )
                    and self._check_del_dup_ins_overlap(variant, lr_var)
                ):
                    remove_it = True
                    break

                elif variant_sv_type == "INV":
                    if self._check_inv_overlap(variant, lr_var):
                        remove_it = True
                        break

                # Special case: short INS vs long INS => remove if pos < 100 bp
                # difference and end < 100 bp difference
                elif variant_sv_type == "INS" and lr_var["svtype"] == "INS":
                    if self._check_ins_vs_ins_overlap(variant, lr_var):
                        remove_it = True
                        break

                # Special case: short INS vs. long DUP => remove if pos <
                # 100 bp difference
                elif variant_sv_type == "INS" and lr_var["svtype"] == "DUP":
                    if self._check_ins_vs_dup_overlap(variant, lr_var):
                        remove_it = True
                        break

            if not remove_it:
                filtered_short.append(variant)

        return filtered_short

    def _merge_svs(
        self, variants: List[Dict], merge_func: Callable[[Dict, Dict], bool]
    ) -> List[Dict]:
        """Merge a list of variant dictionaries. Each variant is sorted by
        ascending priority. For each new variant:
         - Look up only those intervals that overlap its [pos, end) coordinates.
         - Check if the new variant can merge with an existing one.
         - If merge criteria are met, keep only the variant with higher priority
           (numerically lower 'priority' value).
         - Otherwise, add it to the merged set.

        Args:
          variants: The list of variant dictionaries to merge.
          merge_func: A function that decides whether two variants should merge.

        Returns:
            A list of merged variant dictionaries.
        """
        variants_sorted = sorted(variants, key=lambda v: v["priority"])
        merged: List[Dict] = []
        merged_trees: Dict[str, IntervalTree] = {}

        for var in tqdm(variants_sorted, desc="Merging (IntervalTree)"):
            found_merge = False
            chrom = var["chrom"]
            start, end = var["pos"], var["end"]

            if chrom not in merged_trees:
                merged_trees[chrom] = IntervalTree()

            overlapping_intervals = merged_trees[chrom].overlap(start, end)
            for interval in overlapping_intervals:
                existing_var = interval.data
                if merge_func(var, existing_var):
                    if var["priority"] < existing_var["priority"]:
                        merged_trees[chrom].remove(interval)
                        merged.remove(existing_var)
                    else:
                        found_merge = True
                    break

            if not found_merge:
                merged.append(var)
                merged_trees[chrom].addi(start, end, var)

        return merged

    @staticmethod
    def _should_merge_svs(
        var_a: Dict[str, Union[str, int]],
        var_b: Dict[str, Union[str, int]],
    ) -> bool:
        """Merge two SVs according to merging critera.

        Merging criteria:
        - INV: >= 20% reciprocal overlap (RO).
        - DEL/DUP/INS:
            1) If >= 50% reciprocal overlap => merge
            2) Otherwise, if breakpoints within 200 bp *and* >= 80% size ratio
            => merge
            3) If breakpoints within 3x the LARGER SV size AND >=80% size ratio
            => merge

        Args:
            var_a: First variant dictionary.
            var_b: Second variant dictionary.

        Returns:
            True if the two variants can be merged, False otherwise.
        """
        if not SVReferenceMerger._verify_sv_pair(var_a, var_b):
            return False

        pos_a, end_a, size_a = int(var_a["pos"]), int(var_a["end"]), int(var_a["size"])
        pos_b, end_b, size_b = int(var_b["pos"]), int(var_b["end"]), int(var_b["size"])
        sv_type = var_a["svtype"]

        if sv_type == "INV":
            return SVReferenceMerger._merge_inv_svs(pos_a, end_a, pos_b, end_b)
        elif sv_type in ("DEL", "DUP", "INS"):
            return SVReferenceMerger._merge_del_dup_ins_svs(
                pos_a, end_a, size_a, pos_b, end_b, size_b
            )
        return False

    def _get_ref_base(self, variant: Dict) -> Union[str, None]:
        """Get a valid reference base for the variant."""
        var_ref = variant.get("ref")
        if var_ref and var_ref != "N":
            ref_base = var_ref
        else:
            ref_base = self.ref.fetch(
                variant["chrom"], variant["pos"] - 1, variant["pos"]
            ).upper()

        if ref_base not in ("A", "C", "G", "T"):
            print(
                f"[Warning] Reference base for {variant} not recognized. Likely out of bounds."
            )
            return None

        return ref_base

    def _write_vcf_header(self) -> pysam.VariantHeader:
        """Write a minimal VCF header to the specified file."""
        contigs = {var["chrom"] for var in self.merged_variants}

        header = pysam.VariantHeader()
        header.add_meta("source", "RefMerger SV genotype reference panel.")
        header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">')
        header.add_line(
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of SV">'
        )
        header.add_line(
            '##INFO=<ID=SIZE,Number=1,Type=Integer,Description="Calculated size of the SV (end - pos)">'
        )
        header.add_line(
            '##INFO=<ID=SEQ,Number=1,Type=String,Description="Insertion sequence">'
        )
        header.add_line(
            '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Original data source">'
        )
        header.add_line(
            '##INFO=<ID=PRIORITY,Number=1,Type=Integer,Description="Priority value">'
        )
        for contig in contigs:
            header.add_line(f"##contig=<ID={contig}>")

        return header

    def write_merged(self, outfile: str) -> None:
        """Write final merged variants to an output VCF file."""
        if not self.merged_variants:
            print("[RefMerger] No merged variants to write. Run the pipeline first.")
            return

        header = self._write_vcf_header()

        with pysam.VariantFile(outfile, mode="w", header=header) as out_vcf:
            for var in self.merged_variants:
                try:
                    if record := self._create_vcf_record(var, out_vcf):
                        out_vcf.write(record)
                except Exception as e:
                    self._log_record_error(var, e)
                    continue

        print(f"[RefMerger] Wrote final merged calls to {outfile}")

    def _create_vcf_record(
        self, var: Dict[str, Union[str, int]], out_vcf: pysam.VariantFile
    ) -> Optional[pysam.VariantRecord]:
        """Create a VCF record from a variant dictionary."""
        svtype = var["svtype"]
        ref_base = self._get_ref_base(var)
        if not ref_base:
            return None

        rec = out_vcf.new_record()
        rec.chrom = var["chrom"]
        rec.pos = var["pos"]
        rec.stop = var["end"]
        size = var["end"] - var["pos"]  # type: ignore
        rec.id = f"{var['source']}_{svtype}_{var['pos']}_{size}"
        rec.ref = ref_base

        self._set_alt_and_info(
            rec=rec,
            var=var,
            svtype=svtype,  # type: ignore
            size=size,
        )
        return rec

    def _set_alt_and_info(
        self,
        rec: pysam.VariantRecord,
        var: Dict[str, Union[str, int]],
        svtype: str,
        size: int,
    ) -> None:
        """Set ALT and INFO fields based on SV type."""
        if svtype == "DEL":
            rec.alts = ("<DEL>",)
        elif svtype == "INS":
            rec.alts = ("<INS>",)
            insertion_seq = var.get("seq", "")
            rec.info["SEQ"] = insertion_seq
        elif svtype == "DUP":
            rec.alts = ("<DUP>",)
        elif svtype == "INV":
            rec.alts = ("<INV>",)
        else:
            print(
                f"[Warning] SV type {svtype} not recognized, writing as symbolic <{svtype}>"
            )
            rec.alts = (f"<{svtype}>",)

        rec.info["SVTYPE"] = svtype
        rec.info["SOURCE"] = var["source"]
        rec.info["PRIORITY"] = var["priority"]
        rec.info["SIZE"] = size

    def _log_record_error(
        self, var: Dict[str, Union[str, int]], error: Exception
    ) -> None:
        """Log an error when writing a record fails."""
        print(f"[Error] Failed to write record for variant={var}: {error}")

    def _check_del_dup_ins_overlap(self, variant: Dict, lr_var: Dict) -> bool:
        """Check overlap for DEL/DUP/INS SV types based on size."""
        variant_size = variant["size"]
        ro = self.reciprocal_overlap(
            variant["pos"], variant["end"], lr_var["pos"], lr_var["end"]
        )
        return (variant_size >= 5000 and ro >= 0.50) or (
            variant_size < 5000 and ro >= 0.10
        )

    def _check_inv_overlap(self, variant: Dict, lr_var: Dict) -> bool:
        """Check overlap for INV SV types."""
        ro_inv = self.reciprocal_overlap(
            variant["pos"], variant["end"], lr_var["pos"], lr_var["end"]
        )
        return ro_inv >= 0.20

    @staticmethod
    def _check_ins_vs_ins_overlap(variant: Dict, lr_var: Dict) -> bool:
        """Check overlap for INS vs INS special case."""
        return (
            abs(variant["pos"] - lr_var["pos"]) <= 100
            and abs(variant["end"] - lr_var["end"]) <= 100
        )

    @staticmethod
    def _check_ins_vs_dup_overlap(variant: Dict, lr_var: Dict) -> bool:
        """Check overlap for INS vs DUP special case."""
        return abs(variant["pos"] - lr_var["pos"]) <= 100

    @staticmethod
    def reciprocal_overlap(
        start_a: int,
        end_a: int,
        start_b: int,
        end_b: int,
    ) -> float:
        """Compute the reciprocal overlap for two intervals on the same
        chromosome.

        Returns:
          Float reciprocal overlap between the two intervals.
        """
        if end_a <= start_a or end_b <= start_b:
            return 0.0

        overlap = max(0, min(end_a, end_b) - max(start_a, start_b))
        return 0.0 if overlap == 0 else overlap / max(end_a - start_a, end_b - start_b)

    @staticmethod
    def _merge_inv_svs(pos_a: int, end_a: int, pos_b: int, end_b: int) -> bool:
        """Determine if two INV SVs should be merged based on reciprocal
        overlap.
        """
        ro = SVReferenceMerger.reciprocal_overlap(pos_a, end_a, pos_b, end_b)
        return ro >= 0.20

    @staticmethod
    def _merge_del_dup_ins_svs(
        pos_a: int,
        end_a: int,
        size_a: int,
        pos_b: int,
        end_b: int,
        size_b: int,
    ) -> bool:
        """Determine if two DEL/DUP/INS SVs should be merged based on
        criteria.
        """
        # Check for >= 50% reciprocal overlap
        ro = SVReferenceMerger.reciprocal_overlap(pos_a, end_a, pos_b, end_b)
        if ro >= 0.50:
            return True

        # Check for breakpoints within 200 bp and size ratio >= 80%
        start_dist = abs(pos_a - pos_b)
        end_dist = abs(end_a - end_b)
        size_ratio = SVReferenceMerger._calculate_size_ratio(size_a, size_b)
        if (start_dist <= 200) and (end_dist <= 200) and (size_ratio >= 0.80):
            return True

        # breakpoints <= 3× LARGER SV size AND size ratio >=80%
        larger_size = max(size_a, size_b)
        return larger_size > 0 and (
            (start_dist <= 3 * larger_size)
            and (end_dist <= 3 * larger_size)
            and (size_ratio >= 0.80)
        )

    @staticmethod
    def _calculate_size_ratio(size_a: int, size_b: int) -> float:
        """Calculate the size ratio (or size overlap) of two SVs"""
        larger_size = max(size_a, size_b)
        return min(size_a, size_b) / float(larger_size) if larger_size > 0 else 0.0

    @staticmethod
    def _verify_sv_pair(
        var_a: Dict[str, Union[str, int]],
        var_b: Dict[str, Union[str, int]],
    ) -> bool:
        """Ensure two SVs are on the same chromosome, of the same type, and
        possess required keys.
        """
        if var_a["chrom"] != var_b["chrom"]:
            return False

        if var_a["svtype"] != var_b["svtype"]:
            return False

        required_keys = ["pos", "end", "size"]
        return not any(key not in var_a or key not in var_b for key in required_keys)
