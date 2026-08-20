#!/usr/bin/env python3
"""
select_spaced_blocks.py

Greedily selects non-overlapping blocks of a fixed length from a callable
BED file, such that consecutive selected blocks are separated by at least
`--spacing` bp. This is the step that enforces "blocks are independent data
points": composite-likelihood inference over blocks assumes each block is
(approximately) an independent draw, so blocks that sit close together on
the same chromosome (and are therefore linked) must not both be selected.

This replaces the original pipeline's per-gene intron sampler
(bed_pruner.v1.py), which picked one region per gene - that approach
doesn't generalise to intergenic-only filtering, where there is no "one
feature per block" structure to key off. Here selection is purely
positional: walk each chromosome's callable intervals in order, and take a
block whenever one fits, then skip ahead by (block_length + spacing)
before considering the next.

IMPORTANT: any columns beyond chrom/start/end on an input line are
preserved verbatim on every block carved out of that line. This matters
because gIMble's own BED format (as produced by gimbleprep, and expected
by `gimble parse -b`) is 5 columns - chrom, start, end, num_samples,
samples (a comma-separated per-interval callable-sample list) - not a
plain 3-column BED; `gimble parse` reads columns 0,1,2,4 specifically
(see LohseLab/gIMble's parse_intervals) and errors if column 4 doesn't
exist. A sub-window carved from a larger callable interval was, by
construction of that interval, callable for the exact same set of samples
across its entire span, so inheriting the parent interval's extra columns
unchanged is correct, not just convenient.

Usage:
    select_spaced_blocks.py --bed callable.bed --block-length 200 \
        --spacing 10000 --out spaced_blocks.bed

Input BED must be sorted by (chrom, start); this is enforced by
re-sorting internally so pre-sorting is not required, only that all rows
for a given chromosome are handled together (this is guaranteed once
sorted here).
"""
import argparse
import sys
from collections import defaultdict


def read_bed(path):
    """Returns {chrom: [(start, end, extra_fields_str_or_None), ...]}."""
    intervals = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            extra = "\t".join(fields[3:]) if len(fields) > 3 else None
            if end > start:
                intervals[chrom].append((start, end, extra))
    for chrom in intervals:
        intervals[chrom].sort(key=lambda t: (t[0], t[1]))
    return intervals


def select_blocks(intervals, block_length, spacing):
    """
    Yields (chrom, start, end, extra) for each selected block, where
    `extra` is whatever extra columns (if any) the source interval had,
    inherited unchanged, or None if the input had no columns past 0,1,2.
    """
    for chrom in sorted(intervals):
        cursor = None  # earliest position the next block may start at
        for (istart, iend, extra) in intervals[chrom]:
            pos = istart if cursor is None else max(cursor, istart)
            while pos + block_length <= iend:
                yield (chrom, pos, pos + block_length, extra)
                pos = pos + block_length + spacing
                cursor = pos


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", required=True, help="input callable BED (chrom, start, end, ...)")
    p.add_argument("--block-length", type=int, required=True,
                   help="length in bp of each block")
    p.add_argument("--spacing", type=int, required=True,
                   help="minimum bp gap between the end of one block and the start of the next")
    p.add_argument("--out", required=True, help="output BED of selected blocks")
    args = p.parse_args()

    if args.block_length <= 0 or args.spacing < 0:
        sys.exit("[X] --block-length must be > 0 and --spacing must be >= 0")

    intervals = read_bed(args.bed)
    n = 0
    with open(args.out, "w") as out:
        for chrom, start, end, extra in select_blocks(intervals, args.block_length, args.spacing):
            if extra is not None:
                out.write(f"{chrom}\t{start}\t{end}\t{extra}\n")
            else:
                out.write(f"{chrom}\t{start}\t{end}\n")
            n += 1

    print(f"[=] Selected {n} blocks of length {args.block_length} bp, "
          f">= {args.spacing} bp apart, written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()

