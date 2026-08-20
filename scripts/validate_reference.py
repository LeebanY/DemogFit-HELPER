#!/usr/bin/env python3
"""
validate_reference.py

Independent, post-hoc validation of the REFERENCE-LEVEL intergenic
construction (gene buffering + repeat removal). Run immediately after
these BED files are built and BEFORE any per-sample work (read QC,
mapping, variant calling) begins - this step depends only on the
reference genome, GFF, and repeat BED, never on any sample's reads, so
there is no reason to wait until after expensive per-sample processing to
find out the region construction was wrong. That is what motivated this
script: an empty gene BED (from a GFF/FASTA sequence-name mismatch)
previously went undetected until deep into a real run, after mapping and
variant calling had already completed, silently reclassifying the entire
genome as "intergenic". This validates the actual output files, not the
same bedtools logic that built them - a bug in that logic wouldn't be
caught by re-deriving it the same way.

Checks:
  1. genes.bed and genes.buffered.merged.bed are non-empty. (This is also
     independently enforced as a hard failure inside the gene_bed rule
     itself - checking again here is cheap, and catches the case where
     genes.bed exists from a stale prior run but a later step silently
     produced nothing.)
  2. intergenic.repeat_free.bed has ZERO overlap with the gene+buffer
     zone - a structural guarantee of `bedtools complement`, checked
     independently in case of a coordinate or sorting bug in the chain
     that built it.
  3. intergenic.repeat_free.bed has ZERO overlap with repeats.
  4. intergenic.repeat_free.bed's total bp is strictly less than the
     genome's total bp. If they are equal, gene buffering excluded
     nothing at all - the same silent-failure signature as the empty-
     genes.bed incident this script exists to catch.

Exits non-zero (and the report says FAIL) if any check fails.
"""
import argparse
import os
import subprocess
import sys
import tempfile


def run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def count_lines(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.strip())


def total_bp(path):
    total = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.split("\t")
            total += int(fields[2]) - int(fields[1])
    return total


def check_non_empty(path, label):
    n = count_lines(path)
    if n == 0:
        return False, f"{label} is empty (0 intervals) - see gene_bed's own error for likely cause"
    return True, f"{label} has {n} interval(s)"


def check_no_overlap(target_bed, source_cmd, label):
    fd, tmp_path = tempfile.mkstemp(suffix=".bed")
    os.close(fd)
    try:
        r = run(f"{source_cmd} > {tmp_path}")
        if r.returncode != 0:
            return False, f"failed to prepare {label} for comparison: {r.stderr.strip()}"
        r = run(f"bedtools intersect -a {target_bed} -b {tmp_path} -u")
        if r.returncode != 0:
            return False, f"bedtools intersect against {label} failed: {r.stderr.strip()}"
        lines = [l for l in r.stdout.splitlines() if l.strip()]
        if lines:
            return False, (f"{len(lines)} intergenic interval(s) overlap {label} "
                            f"(expected 0). First offender: {lines[0]}")
        return True, f"0 intergenic intervals overlap {label}"
    finally:
        os.unlink(tmp_path)


def check_meaningfully_smaller_than_genome(intergenic_bed, genome_file):
    genome_bp = 0
    with open(genome_file) as fh:
        for line in fh:
            if not line.strip():
                continue
            genome_bp += int(line.split("\t")[1])
    intergenic_bp = total_bp(intergenic_bed)
    if intergenic_bp >= genome_bp:
        return False, (f"intergenic.repeat_free.bed covers {intergenic_bp}bp, which is >= "
                        f"the genome's {genome_bp}bp total - gene buffering/repeat removal "
                        f"excluded nothing at all, which is the same failure signature as "
                        f"an empty genes.bed")
    pct_excluded = 100 * (1 - intergenic_bp / genome_bp) if genome_bp else 0
    return True, (f"intergenic.repeat_free.bed covers {intergenic_bp}bp of {genome_bp}bp "
                   f"genome ({pct_excluded:.2f}% excluded by gene buffer + repeats)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True)
    p.add_argument("--genes-buffered", required=True)
    p.add_argument("--repeats", required=True)
    p.add_argument("--intergenic-repeat-free", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--report", required=True)
    args = p.parse_args()

    checks = [
        ("genes.bed non-empty", *check_non_empty(args.genes, "genes.bed")),
        ("genes.buffered.merged.bed non-empty",
         *check_non_empty(args.genes_buffered, "genes.buffered.merged.bed")),
        ("intergenic is smaller than the genome",
         *check_meaningfully_smaller_than_genome(args.intergenic_repeat_free, args.genome)),
        ("no overlap with gene+buffer zone",
         *check_no_overlap(args.intergenic_repeat_free, f"cat {args.genes_buffered}",
                            "the gene+buffer zone")),
        ("no overlap with repeats",
         *check_no_overlap(args.intergenic_repeat_free, f"zcat -f {args.repeats}", "repeats")),
    ]

    passed = all(c[1] for c in checks)
    with open(args.report, "w") as out:
        out.write("Reference-level intergenic construction validation\n" + "=" * 55 + "\n")
        for name, ok, msg in checks:
            out.write(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}\n")
        out.write(f"\nOverall: {'PASS' if passed else 'FAIL'}\n")

    for name, ok, msg in checks:
        print(f"[{'=' if ok else 'X'}] {name}: {msg}", file=sys.stderr)

    if not passed:
        print("\n[X] Reference validation FAILED - stopping before any sample-level work "
              "(read QC, mapping, variant calling) begins. See report for details.", file=sys.stderr)
        sys.exit(1)
    print("\n[=] All reference validation checks passed.", file=sys.stderr)


if __name__ == "__main__":
    main()
