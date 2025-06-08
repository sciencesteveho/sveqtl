# test_sv_reference_merger.py
# -*- coding: utf-8 -*-


"""Unit tests for SVReferenceMerger."""


import pytest

from sveqtl.genotyping.sv_reference_merger import SVReferenceMerger

_Short = "short_read"
_Long = "long_read"


def _base_variant(**kwargs):
    """Return a baseline variant dict; override fields via kwargs."""
    return {
        "chrom": "chr1",
        "pos": 1000,
        "end": 1200,
        "size": 200,
        "svtype": "DEL",
        "source": _Short,
        "priority": 1,
    } | kwargs


def large_del(**kw):
    return _base_variant(**kw)


@pytest.fixture()
def del_small_a():
    return _base_variant()


@pytest.fixture()
def del_small_b_close(del_small_a):
    return _base_variant(
        pos=1050,
        end=1250,
        source="srcB",
        priority=2,
    )


@pytest.fixture()
def del_small_c_diff_size(del_small_a):
    return _base_variant(
        pos=1080,
        end=1500,
        size=420,  # < 0.8 size ratio w.r.t. A
        source="srcC",
        priority=3,
    )


@pytest.fixture()
def inv_a():
    return _base_variant(
        chrom="chr3",
        pos=10_000,
        end=15_000,
        size=5_000,
        svtype="INV",
        priority=1,
    )


@pytest.fixture()
def inv_b_overlap20(inv_a):
    return inv_a | {"pos": 14_000, "end": 19_000, "source": "srcB", "priority": 2}


@pytest.mark.parametrize(
    "start_a,end_a,start_b,end_b,expected",
    [
        (100, 200, 100, 200, 1.0),  # identical
        (100, 200, 300, 400, 0.0),  # disjoint
        (100, 200, 150, 250, 0.5),  # half overlap
        (100, 400, 200, 300, 1 / 3),  # one contained
        (100, 200, 200, 300, 0.0),  # abutting
    ],
)
def test_reciprocal_overlap(start_a, end_a, start_b, end_b, expected):
    """SVReferenceMerger.reciprocal_overlap returns expected float."""
    assert SVReferenceMerger.reciprocal_overlap(
        start_a, end_a, start_b, end_b
    ) == pytest.approx(expected)


def test_should_merge_small_dels(del_small_a, del_small_b_close):
    """Merge small DELs that are breakpoint-close & size-similar."""
    assert SVReferenceMerger._should_merge_svs(del_small_a, del_small_b_close)


def test_should_merge_large_dels_ro50():
    """Large DELs (≥5 kbp) merge when RO ≥ 0.5."""
    sv_1 = large_del(pos=5_000, end=15_000)
    sv_2 = large_del(pos=7_500, end=17_500, source="srcB")
    assert SVReferenceMerger._should_merge_svs(sv_1, sv_2)


def test_should_not_merge_ins_breakpoint_far():
    """INS pairs with breakpoints >100 bp apart must not merge."""
    ins_a = _base_variant(svtype="INS", end=1_001, size=1)
    ins_b = ins_a | {"pos": ins_a["pos"] + 250, "end": ins_a["end"] + 250}
    assert not SVReferenceMerger._should_merge_svs(ins_a, ins_b)


def test_merge_svs_late_variant_higher_priority(del_small_a, del_small_b_close):
    """If low-priority variant comes first, higher-priority later variant
    replaces it.
    """
    earlier_low_pri = del_small_b_close  # priority 2
    later_high_pri = del_small_a  # priority 1
    merger = SVReferenceMerger(short_read_svs=[], long_read_svs=[])
    result = merger._merge_svs(
        [earlier_low_pri, later_high_pri],
        SVReferenceMerger._should_merge_svs,
    )
    assert result == [later_high_pri]


def test_filter_short_keeps_inv_when_ro_lt20():
    """SR INV is kept when LR INV overlap <20 %."""
    sr_inv = _base_variant(
        chrom="chr4",
        pos=10_000,
        end=15_000,
        size=5_000,
        svtype="INV",
        source=_Short,
    )
    # 500 bp overlap -> RO = 0.10
    lr_inv = sr_inv | {
        "pos": 14_500,
        "end": 19_500,
        "source": _Long,
    }
    merger = SVReferenceMerger(short_read_svs=[], long_read_svs=[lr_inv])
    kept = merger._filter_short_read_svs([sr_inv])
    assert kept == [sr_inv]


def test_should_not_merge_size_mismatch(del_small_a, del_small_c_diff_size):
    """Do not merge when size ratio < 0.8."""
    assert not SVReferenceMerger._should_merge_svs(del_small_a, del_small_c_diff_size)


def test_should_merge_inversions(inv_a, inv_b_overlap20):
    """Inversions merge at ≥ 20 % reciprocal overlap."""
    assert SVReferenceMerger._should_merge_svs(inv_a, inv_b_overlap20)


def test_merge_svs_returns_highest_priority(del_small_a, del_small_b_close):
    """Catalogue merge keeps the lower numeric priority."""
    merger = SVReferenceMerger(short_read_svs=[], long_read_svs=[])
    result = merger._merge_svs(
        [del_small_b_close, del_small_a], SVReferenceMerger._should_merge_svs
    )
    assert result == [del_small_a]


def test_filter_short_removes_overlapping_sr_calls():
    """SR call is removed when ≥10 % RO with LR call (short DEL example)."""
    sr = _base_variant()
    lr = sr | {"source": _Long, "pos": 1040, "end": 1240}
    merger = SVReferenceMerger(short_read_svs=[], long_read_svs=[lr])
    kept = merger._filter_short_read_svs([sr])
    assert kept == []


def test_merge_callsets_end_to_end(del_small_a, del_small_b_close):
    """Public pipeline produces a single merged variant, sorted correctly."""
    merger = SVReferenceMerger(
        short_read_svs=[del_small_b_close],
        long_read_svs=[del_small_a],  # higher priority, identical region
    )
    merger.merge_callsets()
    assert merger.merged_variants == [del_small_a]
