#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""Plot some basic statistics for an SV reference panel."""


from collections import Counter
from collections import defaultdict
from typing import Dict, List, Set, Tuple

from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedFormatter
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import NullLocator
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pysam

from sveqtl.figures import _set_matplotlib_publication_parameters

SHORT_READ_SOURCES = ["1kg30x", "ccdg"]
LONG_READ_SOURCES = ["ont100", "ont1k", "hgsvc3"]


def _normalize_source(source_str: str) -> str:
    """Adjust HGSVC3 string for consistency."""
    return "hgsvc3" if source_str.startswith("hgsvc3") else source_str


def _parse_vcf_size_data(vcf_path: str) -> Dict[str, List[int]]:
    """Parse the given VCF file to extract SV sizes keyed by SVTYPE."""
    variants: Dict[str, List[int]] = defaultdict(list)

    with pysam.VariantFile(vcf_path) as vcf_file:
        for record in vcf_file.fetch():
            svtype = record.info.get("SVTYPE", None)
            svsize = record.info.get("SIZE", None)

            if svtype is None or svsize is None:
                continue

            svtype = str(svtype).upper()
            size = int(svsize)
            variants[svtype].append(size)

    return variants


def _parse_vcf_source_data(vcf_path: str) -> Dict[str, List[str]]:
    """Return a dict of lists mapping each SVTYPE to a list of sources, plus an
    'ALL' entry for the total across types.
    """
    svtype_to_sources = defaultdict(list)

    with pysam.VariantFile(vcf_path) as vcf_file:
        for record in vcf_file.fetch():
            svtype = record.info.get("SVTYPE", None)
            source = record.info.get("SOURCE", None)

            if not svtype or not source:
                continue

            svtype = svtype.upper()
            source = _normalize_source(source)

            svtype_to_sources[svtype].append(source)

    # Ensure 'ALL' is included
    all_sources = []
    for s in svtype_to_sources.values():
        all_sources.extend(s)
    svtype_to_sources["ALL"] = all_sources

    return svtype_to_sources


def create_sv_barplot(variants: Dict[str, List[int]]) -> Figure:
    """Create a barplot of variant counts per SV type."""
    counts_by_type = {t: len(sizes) for t, sizes in variants.items()}

    sorted_counts = sorted(counts_by_type.items(), key=lambda x: x[1], reverse=True)
    sorted_types, counts = zip(*sorted_counts)
    unique_types = sorted_types

    cmap = plt.get_cmap("tab10")
    type_to_color = {t: cmap(i % 10) for i, t in enumerate(unique_types)}

    barplot, ax_bar = plt.subplots(figsize=(0.85, 2))
    bar_positions = np.arange(len(unique_types))
    ax_bar.bar(
        bar_positions,
        counts,
        color=[type_to_color[t] for t in unique_types],
        edgecolor="dimgray",
        linewidth=0.25,
        width=0.5,
    )

    for spine in ["top", "right", "bottom"]:
        ax_bar.spines[spine].set_visible(False)

    plt.ylabel("Number of SVs")

    # Remove x-tick lines and color x-tick labels
    ax_bar.tick_params(axis="x", which="both", bottom=False, top=False, pad=-2.5)
    ax_bar.set_xticks(bar_positions)
    ax_bar.set_xticklabels(unique_types, rotation=90)
    for label, bar_color in zip(
        ax_bar.get_xticklabels(), [type_to_color[t] for t in unique_types]
    ):
        label.set_color(bar_color)

    # Annotate each bar
    for idx, c in enumerate(counts):
        if idx == 0:
            ha = "left"
            x_coord = idx - 0.25
        else:
            ha = "center"
            x_coord = idx

        ax_bar.text(
            x_coord,
            c,
            format(c, ",d"),
            ha=ha,
            va="bottom",
            color=type_to_color[unique_types[idx]],
        )

    # Add total count as title
    total_count = sum(counts)
    ax_bar.set_title(f"n={format(total_count, ',d')}")
    return barplot


def create_sv_size_plots(
    variants: Dict[str, List[int]],
    min_val: int = 50,
    max_val: int = 1_000_000,
) -> Figure:
    """Create a stacked set of SV size plots (hist + KDE) per SV type, all
    sharing the same X-axis range.
    """
    # For consistent ordering
    unique_types = sorted(variants.keys())

    if not unique_types:
        fig_empty, ax_empty = plt.subplots()
        ax_empty.set_title("No SV data available.")
        return fig_empty

    fig, ax = plt.subplots(figsize=(2.25, 1.85))

    cmap = plt.get_cmap("tab10")
    type_to_color = {t: cmap(i % 10) for i, t in enumerate(unique_types)}
    bins = np.logspace(np.log10(min_val), np.log10(max_val), 500)

    window_size = 6
    window = np.ones(window_size) / window_size
    for t in unique_types:
        these_sizes = np.array(variants[t])
        if len(these_sizes) > 1:
            counts, edges = np.histogram(these_sizes, bins=bins)
            midpoints = 0.5 * (edges[:-1] + edges[1:])
            smoothed_counts = np.convolve(counts, window, mode="same")
            ax.plot(
                midpoints,
                smoothed_counts,
                label=t,
                color=type_to_color[t],
                linewidth=0.75,
                alpha=0.85,
            )

    ax.set_xscale("log")
    ax.set_yscale("log")

    x_labels = ["100bp", "1kb", "10kb", "100kb", "1Mb"]
    x_ticks = [100, 1e3, 1e4, 1e5, 1e6]
    y_ticks = [10, 100, 1_000, 10_000]

    ax.set(
        xlim=(min_val, max_val),
        ylim=(0, 10_000),
        xticks=x_ticks,
        xticklabels=x_labels,
        yticks=y_ticks,
        xlabel="SV size",
        ylabel="Number of SVs",
    )

    ax.xaxis.set_major_locator(FixedLocator(x_ticks))
    ax.xaxis.set_major_formatter(FixedFormatter(x_labels))
    ax.yaxis.set_major_formatter(ScalarFormatter())
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_minor_locator(NullLocator())

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.legend(frameon=False, handlelength=0.75, labelspacing=0.175)
    fig.tight_layout()

    return fig


def create_sv_source_stacked_barplot(
    svtype_to_sources: Dict[str, List[str]],
) -> Figure:
    """Creates a stacked horizontal barplot where each bar is an SV type
    (ALL, DEL, DUP, INS, INV) showing the proportion of each source
    contributing to that SV type.
    """
    legend_name_mapping = {
        "hgsvc3": "Logsdon et al.",
        "ont1k": "Schloissnig et al.",
        "ont100": "Gustafson et al.",
        "1kg30x": "Byrska-Bishop et al.",
        "ccdg": "Abel et al.",
    }
    desired_order = ["ALL", "DEL", "DUP", "INS", "INV"]
    short_read_colors = ["#9B59B6", "#BB8FCE"]  # purples
    long_read_colors = ["#27AE60", "#52BE80", "#82E0AA"]  # greens
    fallback_color = "#95A5A6"  # gray

    counts_per_type: Dict[str, Counter] = {}
    for svtype in desired_order:
        if svtype not in svtype_to_sources:
            counts_per_type[svtype] = Counter()
        else:
            counts_per_type[svtype] = Counter(svtype_to_sources[svtype])

    def get_source_color_map(all_sources: Set[str]) -> Dict[str, str]:
        """Defines a color map for sources."""
        source_to_color = {}
        sr_i, lr_i = 0, 0
        for src in all_sources:
            if src in SHORT_READ_SOURCES:
                source_to_color[src] = short_read_colors[sr_i % len(short_read_colors)]
                sr_i += 1
            elif src in LONG_READ_SOURCES:
                source_to_color[src] = long_read_colors[lr_i % len(long_read_colors)]
                lr_i += 1
            else:
                source_to_color[src] = fallback_color

        return source_to_color

    unique_sources: Set[str] = set()
    for counts in counts_per_type.values():
        unique_sources.update(counts.keys())
    ordered_sources: List[str] = sorted(unique_sources)
    source_to_color = get_source_color_map(unique_sources)

    fig, ax = plt.subplots(figsize=(3.25, 1.35))
    fig.subplots_adjust(left=0.26, right=0.82)

    y_positions = np.arange(len(desired_order))[::-1]
    ax.set_yticks(y_positions)
    ax.set_yticklabels(desired_order)
    ax.set_ylim(-0.5, len(desired_order) - 0.5)

    # Make the left spine line up with the top and bottom ticks:
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_bounds(0, len(desired_order) - 1)
    ax.spines["left"].set_position(("outward", 25))

    # Build bars for each SV type
    for y_pos, svtype in zip(y_positions, desired_order):
        counts = counts_per_type[svtype]
        total_count = sum(counts.values())
        if total_count == 0:
            continue

        longread_svs = [s for s in counts if s in LONG_READ_SOURCES]
        shortread_svs = [s for s in counts if s in SHORT_READ_SOURCES]
        other_svs = [
            s
            for s in counts
            if s not in LONG_READ_SOURCES and s not in SHORT_READ_SOURCES
        ]

        longread_svs.sort(key=lambda s: counts[s], reverse=True)
        shortread_svs.sort(key=lambda s: counts[s], reverse=True)
        other_svs.sort(key=lambda s: counts[s], reverse=True)
        grouped_sources = longread_svs + shortread_svs + other_svs

        left_start = 0.0
        for src in grouped_sources:
            proportion = counts[src] / total_count
            ax.barh(
                y=y_pos,
                width=proportion,
                left=left_start,
                color=source_to_color[src],
                edgecolor="white",
                linewidth=0.4,
            )
            left_start += proportion

        long_read_count = sum(counts[s] for s in counts if s in LONG_READ_SOURCES)
        short_read_count = sum(counts[s] for s in counts if s in SHORT_READ_SOURCES)

        ax.text(
            -0.02,
            y_pos,  # type: ignore
            f"{long_read_count:,d}",
            ha="right",
            va="center",
            fontsize=5,
        )

        ax.text(
            1.02,
            y_pos,  # type: ignore
            f"{short_read_count:,d}",
            ha="left",
            va="center",
            fontsize=5,
        )

    ax.set_xlim(0, 1)
    ax.set_xlabel("Proportion of SVs")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Build the legend with long-read first, then short-read, then others
    longread_in_all = [s for s in ordered_sources if s in LONG_READ_SOURCES]
    shortread_in_all = [s for s in ordered_sources if s in SHORT_READ_SOURCES]
    others_in_all = [
        s
        for s in ordered_sources
        if s not in LONG_READ_SOURCES and s not in SHORT_READ_SOURCES
    ]
    legend_order = longread_in_all + shortread_in_all + others_in_all

    # Ensure color mapping is consistent
    legend_handles = []
    seen_colors = {}
    for src in legend_order:
        color = source_to_color[src]
        if color not in seen_colors:
            patch = plt.Rectangle((0, 0), 1, 1, color=color, ec="none")  # type: ignore
            legend_label = legend_name_mapping.get(src, src)
            legend_handles.append((patch, legend_label))
            seen_colors[color] = True

    ax.legend(
        [h for (h, _) in legend_handles],
        [lbl for (_, lbl) in legend_handles],
        frameon=False,
        bbox_to_anchor=(1.18, 1),
        loc="upper left",
        handlelength=0.8,
        labelspacing=0.2,
        fontsize=5,
    )

    return fig


def main(vcf_path: str = "sv_genotype_reference.vcf") -> None:
    """Main function to execute the plotting."""
    _set_matplotlib_publication_parameters()
    variants = _parse_vcf_size_data(vcf_path)
    svtype_to_sources = _parse_vcf_source_data(vcf_path)

    barplot = create_sv_barplot(variants)
    ridgelineplot = create_sv_size_plots(variants)
    stacked_barplot = create_sv_source_stacked_barplot(svtype_to_sources)

    barplot.savefig("sv_barplot.svg", bbox_inches="tight")
    ridgelineplot.savefig("sv_ridgeline.svg", bbox_inches="tight")
    stacked_barplot.savefig("sv_stacked_barplot.svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
