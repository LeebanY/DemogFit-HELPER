#!/usr/bin/env python3
"""
validate_blocks.py

Independent, post-hoc validation that a pair's final selected blocks
actually satisfy the properties appropriate to the region strategy in use
(config.yaml's region_strategy: "intergenic", "intron", or "mixed"):

  - "intergenic": zero overlap with the gene+buffer zone (checking
    against the ALREADY-BUFFERED gene BED in one step confirms both "not
    inside a gene" and "the full configured buffer distance is
    respected"), zero overlap with repeats, block spacing respected.
  - "intron": zero overlap with EXONS (genuinely bug-catching - a
    coordinate/strand-handling bug in extract_introns.py could bleed
    into real coding/UTR sequence), zero overlap with repeats, block
    spacing respected. Overlapping the gene+buffer zone is EXPECTED and
    NOT checked - that is the entire point of this strategy.
  - "mixed": the intron-appropriate checks (exon overlap, repeat overlap,
    spacing) apply to the whole combined block set; the gene+buffer
    overlap check does not apply, since some blocks are legitimately
    drawn from inside gene bodies.

This is deliberately NOT a re-derivation of the same bedtools logic that
built the blocks in the first place - it exists because upstream failures
(an empty gene BED, a stale/incorrect repeat BED, a swapped manifest
column...) can silently produce a technically-valid but scientifically
wrong block set with no error anywhere else in the pipeline. Re-running
the same construction logic wouldn't catch a bug in that logic; checking
the actual output file against the stated properties independently will.

Exits non-zero (and the report says FAIL) if any check fails - this is
meant to be a hard gate before the blocks are ever handed to gIMble.
"""
import argparse
import os
import subprocess
import sys
import tempfile
from collections import defaultdict


def run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def check_no_overlap(blocks_bed, source_cmd, label):
    """source_cmd: a shell command that writes the comparison BED to stdout
    (lets us transparently zcat -f a possibly-gzipped repeat file)."""
    fd, tmp_path = tempfile.mkstemp(suffix=".bed")
    os.close(fd)
    try:
        r = run(f"{source_cmd} > {tmp_path}")
        if r.returncode != 0:
            return False, f"failed to prepare {label} for comparison: {r.stderr.strip()}"
        r = run(f"bedtools intersect -a {blocks_bed} -b {tmp_path} -u")
        if r.returncode != 0:
            return False, f"bedtools intersect against {label} failed: {r.stderr.strip()}"
        lines = [l for l in r.stdout.splitlines() if l.strip()]
        if lines:
            return False, f"{len(lines)} block(s) overlap {label} (expected 0). First offender: {lines[0]}"
        return True, f"0 blocks overlap {label}"
    finally:
        os.unlink(tmp_path)


def check_spacing(blocks_bed, min_spacing):
    by_chrom = defaultdict(list)
    with open(blocks_bed) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, start, end = line.split("\t")[:3]
            by_chrom[chrom].append((int(start), int(end)))

    worst = None
    for chrom, intervals in by_chrom.items():
        intervals.sort()
        for i in range(1, len(intervals)):
            gap = intervals[i][0] - intervals[i - 1][1]
            if gap < 0:
                return False, f"OVERLAPPING blocks on {chrom}: {intervals[i-1]} and {intervals[i]}"
            if worst is None or gap < worst[0]:
                worst = (gap, chrom, intervals[i - 1], intervals[i])

    if worst is None:
        return True, "fewer than 2 blocks on any single chromosome - nothing to compare"
    gap, chrom, a, b = worst
    if gap < min_spacing:
        return False, (f"minimum inter-block gap is {gap}bp on {chrom} "
                        f"(between {a} and {b}), expected >= {min_spacing}bp")
    return True, f"minimum inter-block gap is {gap}bp (>= {min_spacing}bp required)"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--blocks", required=True, help="final selected blocks BED to validate")
    p.add_argument("--strategy", required=True, choices=["intergenic", "intron", "mixed"])
    p.add_argument("--genes-buffered",
                   help="gene BED AFTER buffering/merging - required and checked only for "
                        "strategy=intergenic")
    p.add_argument("--exons",
                   help="exons_used.bed from extract_introns.py - required and checked only "
                        "for strategy=intron/mixed")
    p.add_argument("--repeats", required=True, help="the repeat BED (gzip or plain)")
    p.add_argument("--block-spacing", type=int, required=True)
    p.add_argument("--report", required=True, help="where to write the human-readable report")
    args = p.parse_args()

    if args.strategy == "intergenic" and not args.genes_buffered:
        sys.exit("[X] --genes-buffered is required for --strategy intergenic")
    if args.strategy in ("intron", "mixed") and not args.exons:
        sys.exit(f"[X] --exons is required for --strategy {args.strategy}")

    checks = []
    if args.strategy == "intergenic":
        checks.append((
            "no overlap with gene+buffer zone",
            *check_no_overlap(args.blocks, f"cat {args.genes_buffered}", "the gene+buffer zone")
        ))
    if args.strategy in ("intron", "mixed"):
        checks.append((
            "no overlap with exons",
            *check_no_overlap(args.blocks, f"cat {args.exons}", "exons")
        ))
    checks.append((
        "no overlap with repeats",
        *check_no_overlap(args.blocks, f"zcat -f {args.repeats}", "repeats")
    ))
    checks.append((
        "block spacing respected",
        *check_spacing(args.blocks, args.block_spacing)
    ))

    passed = all(c[1] for c in checks)
    with open(args.report, "w") as out:
        out.write(f"Block validation report (strategy={args.strategy})\n" + "=" * 40 + "\n")
        for name, ok, msg in checks:
            out.write(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}\n")
        out.write(f"\nOverall: {'PASS' if passed else 'FAIL'}\n")

    for name, ok, msg in checks:
        print(f"[{'=' if ok else 'X'}] {name}: {msg}", file=sys.stderr)

    if not passed:
        print("\n[X] Block validation FAILED - see report for details. "
              "The bSFS produced from these blocks should NOT be trusted.", file=sys.stderr)
        sys.exit(1)
    print("\n[=] All block validation checks passed.", file=sys.stderr)


if __name__ == "__main__":
    main()
