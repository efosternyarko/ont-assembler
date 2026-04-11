#!/usr/bin/env python3
"""
parse_read_stats.py — convert seqkit stats -a -T output to a clean per-sample
read metrics TSV, computing estimated sequencing depth.

Usage:
    parse_read_stats.py <sample_id> <genome_size> <seqkit_stats.txt>

genome_size accepts:   5m  5M  5000000  4.5m  500k  etc.
Output (stdout): one-row TSV with columns:
    sample_id  num_reads  total_bases  read_N50  mean_length
    mean_quality  gc_pct  estimated_depth
"""

import sys
import re


def parse_genome_size(gsize_str: str) -> int:
    """Parse genome size string (e.g. '5m', '5M', '500k', '5000000') → int bp."""
    m = re.match(r'^([0-9]+(?:\.[0-9]+)?)([mkMK]?)$', gsize_str.strip())
    if not m:
        raise ValueError(f"Cannot parse genome_size: '{gsize_str}'. "
                         "Use e.g. '5m', '5000000', '500k'.")
    num  = float(m.group(1))
    unit = m.group(2).lower()
    if   unit == 'm': return int(num * 1_000_000)
    elif unit == 'k': return int(num * 1_000)
    else:             return int(num)


def safe_float(val: str, default: float = 0.0) -> float:
    try:
        v = val.replace(',', '').strip()
        return float(v) if v not in ('', 'N/A', '-', 'NA') else default
    except (ValueError, AttributeError):
        return default


def safe_int(val: str, default: int = 0) -> int:
    try:
        v = val.replace(',', '').strip()
        return int(float(v)) if v not in ('', 'N/A', '-', 'NA') else default
    except (ValueError, AttributeError):
        return default


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <sample_id> <genome_size> <seqkit_stats.txt>")

    sample_id   = sys.argv[1]
    genome_size = parse_genome_size(sys.argv[2])
    stats_file  = sys.argv[3]

    with open(stats_file) as fh:
        lines = [l.rstrip('\n') for l in fh if l.strip()]

    if len(lines) < 2:
        sys.exit(f"ERROR: unexpected seqkit stats output in {stats_file} "
                 f"(expected ≥2 lines, got {len(lines)})")

    header = lines[0].split('\t')
    vals   = lines[1].split('\t')
    d      = dict(zip(header, vals))

    num_reads    = safe_int(d.get('num_seqs', '0'))
    total_bases  = safe_int(d.get('sum_len',  '0'))
    read_n50     = safe_int(d.get('N50',      '0'))
    mean_length  = safe_float(d.get('avg_len', '0'))
    # AvgQual only present for FASTQ; seqkit outputs 0 or N/A for FASTA
    mean_quality = safe_float(d.get('AvgQual', '0'))
    # seqkit may label GC as 'GC(%)' or 'GC'
    gc_pct       = safe_float(d.get('GC(%)', d.get('GC', '0')))

    estimated_depth = round(total_bases / genome_size, 1) if genome_size > 0 else 0.0

    print('\t'.join([
        'sample_id', 'num_reads', 'total_bases', 'read_N50',
        'mean_length', 'mean_quality', 'gc_pct', 'estimated_depth'
    ]))
    print('\t'.join([
        sample_id,
        str(num_reads),
        str(total_bases),
        str(read_n50),
        f'{mean_length:.1f}',
        f'{mean_quality:.2f}',
        f'{gc_pct:.2f}',
        str(estimated_depth),
    ]))


if __name__ == '__main__':
    main()
