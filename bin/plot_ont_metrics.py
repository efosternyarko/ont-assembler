#!/usr/bin/env python3
"""
plot_ont_metrics.py — visualise ONT read QC metrics across all samples.

Generates an 8-panel figure (2 rows × 4 columns):
  Row 1: Read N50 histogram | Read N50 boxplot | Read depth histogram | Read depth boxplot
  Row 2: Total reads hist   | Total reads box  | Mean read length hist | Mean read length box

A vertical dashed line on the depth histogram marks --min_depth, and the
panel title reports how many samples pass the depth threshold.

Usage:
    plot_ont_metrics.py --stats sample1_read_stats.tsv [sample2 ...] \\
                        --min_depth 20 \\
                        --output ont_read_metrics.png \\
                        --summary ont_read_metrics_summary.tsv
"""

import argparse
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import gaussian_kde

# ── Colour palette (matches reference figure: blue for N50, green for depth) ──
COLOURS = {
    'n50':    '#6a9fca',  # steel blue
    'depth':  '#7bbf7b',  # muted green
    'reads':  '#e08c6a',  # warm orange
    'length': '#9b7cbc',  # purple
}

PANEL_LABELS = 'ABCDEFGH'


# ── Helpers ───────────────────────────────────────────────────────────────────

def _format_axis_label(label: str, data: np.ndarray) -> str:
    """Append SI suffix to axis label if values are large."""
    if len(data) == 0:
        return label
    mx = np.max(data)
    if mx >= 1e9:
        return f'{label} (Gb)'
    if mx >= 1e6:
        return f'{label} (Mb)'
    if mx >= 1e3:
        return f'{label} (kb)'
    return label


def _scale_data(data: np.ndarray) -> tuple[np.ndarray, str]:
    """Return (scaled_data, suffix_label) so axis values are readable."""
    mx = np.max(data) if len(data) > 0 else 0
    if mx >= 1e9:
        return data / 1e9, 'Gb'
    if mx >= 1e6:
        return data / 1e6, 'Mb'
    if mx >= 1e3:
        return data / 1e3, 'kb'
    return data, ''


def _n_bins(data: np.ndarray) -> int:
    """Reasonable bin count for small-to-medium sample sets."""
    n = len(data)
    if n <= 5:   return n
    if n <= 20:  return max(5, n // 2)
    return min(30, max(10, int(np.ceil(np.sqrt(n)))))


def plot_histogram(ax, data: np.ndarray, colour: str, xlabel: str,
                   title: str, panel_label: str, vline: float | None = None,
                   vline_label: str = '', scale: bool = True) -> None:
    """Histogram + KDE overlay, matching the reference figure style."""
    if len(data) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                transform=ax.transAxes, color='#999999')
        ax.set_title(f'{panel_label}   {title}', loc='left', fontweight='bold')
        return

    if scale:
        scaled, suffix = _scale_data(data)
    else:
        scaled, suffix = data, ''

    bins = _n_bins(scaled)

    ax.hist(scaled, bins=bins, color=colour, alpha=0.70,
            edgecolor='white', linewidth=0.5)

    # KDE overlay (needs ≥ 2 distinct values)
    if len(scaled) >= 2 and np.std(scaled) > 0:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            kde = gaussian_kde(scaled)
        xs = np.linspace(scaled.min(), scaled.max(), 300)
        bin_width = (scaled.max() - scaled.min()) / bins
        ax.plot(xs, kde(xs) * len(scaled) * bin_width,
                color=colour, linewidth=2.0)

    # Depth threshold line
    if vline is not None:
        scaled_vline = vline / (1e9 if suffix == 'Gb' else
                                 1e6 if suffix == 'Mb' else
                                 1e3 if suffix == 'kb' else 1)
        ax.axvline(scaled_vline, color='#d62728', linestyle='--',
                   linewidth=1.5, label=vline_label)
        ax.legend(fontsize=8, framealpha=0.8)

    x_label = f'{xlabel} ({suffix})' if suffix else xlabel
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel('Frequency', fontsize=9)
    ax.set_title(f'{panel_label}   {title}', loc='left', fontweight='bold', fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.7)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=8)


def plot_boxplot(ax, data: np.ndarray, colour: str, xlabel: str,
                 title: str, panel_label: str, scale: bool = True) -> None:
    """Horizontal boxplot matching the reference figure style."""
    if len(data) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                transform=ax.transAxes, color='#999999')
        ax.set_title(f'{panel_label}   {title}', loc='left', fontweight='bold')
        return

    if scale:
        scaled, suffix = _scale_data(data)
    else:
        scaled, suffix = data, ''

    bp = ax.boxplot(
        scaled,
        vert=False,
        patch_artist=True,
        widths=0.5,
        flierprops=dict(marker='D', markerfacecolor='#888888',
                        markeredgecolor='none', markersize=5, alpha=0.7),
        medianprops=dict(color='#333333', linewidth=2),
        boxprops=dict(facecolor=colour, alpha=0.70),
        whiskerprops=dict(color='#555555', linewidth=1.2),
        capprops=dict(color='#555555', linewidth=1.2),
    )

    x_label = f'{xlabel} ({suffix})' if suffix else xlabel
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_yticks([])
    ax.set_title(f'{panel_label}   {title}', loc='left', fontweight='bold', fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.7, axis='x')
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.tick_params(labelsize=8)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stats',   nargs='+', required=True,
                        help='Per-sample *_read_stats.tsv files')
    parser.add_argument('--min_depth', type=float, default=20.0,
                        help='Minimum depth threshold to mark on the depth histogram')
    parser.add_argument('--output',  required=True,
                        help='Output PNG file path')
    parser.add_argument('--summary', required=True,
                        help='Output merged summary TSV path')
    args = parser.parse_args()

    # ── Load and merge per-sample stats ──────────────────────────────────────
    frames = []
    for f in args.stats:
        try:
            df = pd.read_csv(f, sep='\t')
            frames.append(df)
        except Exception as e:
            print(f"WARNING: could not read {f}: {e}", file=sys.stderr)

    if not frames:
        sys.exit("ERROR: no readable stats files provided")

    data = pd.concat(frames, ignore_index=True)

    # Save merged summary
    data.to_csv(args.summary, sep='\t', index=False)

    n = len(data)

    # ── Extract metric arrays ─────────────────────────────────────────────────
    read_n50   = data['read_N50'].to_numpy(dtype=float)
    depth      = data['estimated_depth'].to_numpy(dtype=float)
    num_reads  = data['num_reads'].to_numpy(dtype=float)
    mean_len   = data['mean_length'].to_numpy(dtype=float)

    n_pass = int(np.sum(depth >= args.min_depth))

    # ── Build figure ──────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    fig.suptitle(f'ONT Read Metrics  (n = {n} samples)',
                 fontsize=14, fontweight='bold', y=0.98)

    # Row 1 — Read N50 + Depth
    plot_histogram(axes[0, 0], read_n50,  COLOURS['n50'],
                   'Read N50', 'Read N50 Distribution',   'A')
    plot_boxplot  (axes[0, 1], read_n50,  COLOURS['n50'],
                   'Read N50', 'Read N50 Boxplot',         'B')

    depth_title = (f'Avg Read Depth Distribution  '
                   f'({n_pass}/{n} ≥ {args.min_depth:.0f}×)')
    plot_histogram(axes[0, 2], depth, COLOURS['depth'],
                   'Average Read Depth (×)', depth_title, 'C',
                   vline=args.min_depth,
                   vline_label=f'Min depth ({args.min_depth:.0f}×)',
                   scale=False)
    plot_boxplot  (axes[0, 3], depth, COLOURS['depth'],
                   'Average Read Depth (×)', 'Avg Read Depth Boxplot', 'D',
                   scale=False)

    # Row 2 — Total reads + Mean read length
    # Pre-scale reads to millions so axis shows e.g. "1.2 M reads"
    # (avoids _scale_data tagging the unit as "Mb" = megabases, which is wrong)
    reads_m = num_reads / 1e6
    plot_histogram(axes[1, 0], reads_m, COLOURS['reads'],
                   'Total Reads (M)', 'Total Reads Distribution',  'E', scale=False)
    plot_boxplot  (axes[1, 1], reads_m, COLOURS['reads'],
                   'Total Reads (M)', 'Total Reads Boxplot',       'F', scale=False)

    plot_histogram(axes[1, 2], mean_len, COLOURS['length'],
                   'Mean Read Length', 'Mean Read Length Distribution', 'G')
    plot_boxplot  (axes[1, 3], mean_len, COLOURS['length'],
                   'Mean Read Length', 'Mean Read Length Boxplot',      'H')

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(args.output, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {args.output}  ({n} samples, {n_pass} pass depth ≥ {args.min_depth}×)")


if __name__ == '__main__':
    main()
