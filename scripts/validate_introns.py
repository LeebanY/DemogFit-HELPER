#!/usr/bin/env python3
"""
validate_introns.py

Independent, post-hoc validation of the REFERENCE-LEVEL intron
construction (extract_introns.py's output, after repeat removal),
run immediately after it is built and BEFORE any sample-level work
begins - same placement rationale as validate_reference.py for the
intergenic construction.

Unlike the intergenic case, "overlaps the gene+buffer zone" is not a
failure here - it is the entire point of an intron-based strategy, so
that check does not apply. What DOES matter for introns:

  1. introns.bed and the post-repeat-removal version are both non-empty.
  2. Zero overlap with repeats (same reasoning as intergenic).
  3. Zero overlap with EXONS specifically - genuinely bug-catching,
     unlike gene-body overlap: a coordinate or strand-handling bug in
     extract_introns.py could in principle produce an interval that
     bleeds into real exonic (coding/UTR) sequence, which would be a
     serious problem (strong selective constraint) unlike intronic
     sequence itself.

Exits non-zero (and the report says FAIL) if any check fails.
"""
import argparse
import os
import subprocess
import sys


def run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def count_lines(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.strip())


def check_non_empty(path, label):
    n = count_lines(path)
    if n == 0:
        return False, f"{label} is empty (0 intervals)"
    return True, f"{label} has {n} interval(s)"


def check_no_overlap(target_bed, source_cmd, label):
    fd, tmp_path = tempfile_mkstemp()
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
            return False, (f"{len(lines)} interval(s) overlap {label} "
                            f"(expected 0). First offender: {lines[0]}")
        return True, f"0 intervals overlap {label}"
    finally:
        os.unlink(tmp_path)


def tempfile_mkstemp():
    import tempfile
    return tempfile.mkstemp(suffix=".bed")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--introns", required=True, help="raw candidate introns (before repeat removal)")
    p.add_argument("--introns-repeat-free", required=True)
    p.add_argument("--exons", required=True, help="exons_used.bed from extract_introns.py")
    p.add_argument("--repeats", required=True)
    p.add_argument("--report", required=True)
    args = p.parse_args()

    checks = [
        ("introns.bed non-empty", *check_non_empty(args.introns, "introns.bed")),
        ("intron.repeat_free.bed non-empty", *check_non_empty(args.introns_repeat_free, "intron.repeat_free.bed")),
        ("no overlap with exons",
         *check_no_overlap(args.introns_repeat_free, f"cat {args.exons}", "exons")),
        ("no overlap with repeats",
         *check_no_overlap(args.introns_repeat_free, f"zcat -f {args.repeats}", "repeats")),
    ]

    passed = all(c[1] for c in checks)
    with open(args.report, "w") as out:
        out.write("Reference-level intron construction validation\n" + "=" * 50 + "\n")
        for name, ok, msg in checks:
            out.write(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}\n")
        out.write(f"\nOverall: {'PASS' if passed else 'FAIL'}\n")

    for name, ok, msg in checks:
        print(f"[{'=' if ok else 'X'}] {name}: {msg}", file=sys.stderr)

    if not passed:
        print("\n[X] Intron validation FAILED - stopping before any sample-level work "
              "begins. See report for details.", file=sys.stderr)
        sys.exit(1)
    print("\n[=] All intron validation checks passed.", file=sys.stderr)


if __name__ == "__main__":
    main()
