#!/usr/bin/env python3
"""
extract_introns.py

For each gene in a GFF3, selects ONE candidate intron: the LARGEST intron
that is NOT the first intron (in transcription order, i.e. 5'->3' along
the mRNA - this is coordinate-DECREASING for minus-strand genes, not
coordinate order on the chromosome), trims a fixed number of bp off each
end of it, and emits it as a BED interval. Genes with fewer than 2 introns
on their canonical transcript (0 or 1 introns - excluding the first would
leave nothing to choose from) contribute no interval and are skipped.

Rationale for "largest non-first intron": the first intron is
disproportionately likely to carry promoter-proximal or 5' regulatory
elements (and is systematically shorter/more constrained in many
species); a gene's LARGEST intron, among the remainder, is the one most
likely to have room for a full block after trimming, while still being a
real, single, well-defined interval rather than pooling multiple introns
together.

Multi-transcript genes: there is no single unambiguous "the" intron
structure for a gene with multiple annotated isoforms. This script uses
the LONGEST transcript (by genomic span, i.e. max(exon end) - min(exon
start)) as the gene's canonical transcript - a documented, defensible
choice, not the only possible one. If you disagree with it for your
species/annotation, this is the one place in this script that encodes a
judgment call rather than a mechanical rule from the spec.

GFF3 assumptions (standard, but check against your actual file if this
produces unexpectedly few/many introns - see the "How to check" section
of the pipeline README):
  - gene features have Parent-less top-level IDs
  - transcript features (default feature type: "mRNA") have Parent=<gene ID>
  - exon features (default feature type: "exon") have Parent=<transcript ID>
  - GFF3 attributes are semicolon-separated key=value pairs; Parent may be
    a comma-separated list (a shared exon across isoforms) - each parent
    is handled

Output BED (0-based half-open, matching all other BED files in this
pipeline): chrom, start, end, gene_id - the 4th column is for traceability/
debugging only and is stripped by the same zcat-based awk filters used
elsewhere in the pipeline, so it does not break anything downstream that
expects a plain 3-column BED.

Output: TWO BED files (0-based half-open, matching every other BED file in
this pipeline):
  --out          candidate intron intervals: chrom, start, end, gene_id
  --exons-out    ALL exons in the GFF, from every transcript of every gene
                 (not restricted to canonical transcripts): chrom, start,
                 end, gene_id. This is deliberately comprehensive, not
                 just "the exons of the gene this intron came from" - an
                 intron candidate never overlaps its OWN gene's exons by
                 construction, but CAN legitimately overlap a DIFFERENT
                 gene's exon if that gene is nested inside (or overlaps)
                 the host gene's intron, which is real and reasonably
                 common (confirmed directly against a real genome during
                 this pipeline's development, not a hypothetical - see
                 README). This file is used both to subtract exonic
                 sequence from intron candidates genome-wide and for
                 downstream validation that construction did so correctly.

Usage:
    extract_introns.py --gff genes.gff3[.gz] --trim 10 \
        --transcript-feature mRNA --exon-feature exon \
        --out introns.bed --exons-out exons_used.bed
"""
import argparse
import gzip
import sys
from collections import defaultdict


def open_maybe_gzip(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "r")


def parse_attributes(attr_str):
    attrs = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field or "=" not in field:
            continue
        key, val = field.split("=", 1)
        attrs[key] = val
    return attrs


def parse_gff(path, transcript_feature, exon_feature):
    """
    Returns:
      transcript_parent: {transcript_id: gene_id}
      transcript_strand: {transcript_id: '+' or '-'}
      exons_by_transcript: {transcript_id: [(start, end), ...]} (1-based, GFF-inclusive)
    """
    transcript_parent = {}
    transcript_strand = {}
    exons_by_transcript = defaultdict(list)

    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attr_str = fields[:9]
            attrs = parse_attributes(attr_str)

            if feature == transcript_feature:
                tid = attrs.get("ID")
                parent = attrs.get("Parent")
                if tid is None or parent is None:
                    continue
                # a transcript should have exactly one gene parent; if a
                # GFF ever lists more than one, take the first rather than fail
                gene_id = parent.split(",")[0]
                transcript_parent[tid] = gene_id
                transcript_strand[tid] = strand

            elif feature == exon_feature:
                parent = attrs.get("Parent")
                if parent is None:
                    continue
                for tid in parent.split(","):
                    exons_by_transcript[tid].append((int(start), int(end), chrom))

    return transcript_parent, transcript_strand, exons_by_transcript


def canonical_transcript_per_gene(transcript_parent, exons_by_transcript):
    """Longest transcript (by genomic span) per gene."""
    genes_transcripts = defaultdict(list)
    for tid, gene_id in transcript_parent.items():
        genes_transcripts[gene_id].append(tid)

    canonical = {}
    for gene_id, tids in genes_transcripts.items():
        best_tid, best_span = None, -1
        for tid in tids:
            exons = exons_by_transcript.get(tid)
            if not exons:
                continue
            span = max(e[1] for e in exons) - min(e[0] for e in exons)
            if span > best_span:
                best_tid, best_span = tid, span
        if best_tid is not None:
            canonical[gene_id] = best_tid
    return canonical


def introns_from_exons(exons, strand):
    """
    exons: list of (start, end) 1-based GFF-inclusive coordinates, any order.
    Returns introns as a list of (start, end) 1-based GFF-inclusive gaps,
    ORDERED IN TRANSCRIPTION DIRECTION (5'->3' along the mRNA) - increasing
    genomic coordinate for '+' strand, decreasing genomic coordinate for
    '-' strand. This ordering is what "first intron" means biologically,
    and is NOT the same as genomic left-to-right order for minus-strand genes.
    """
    uniq_exons = sorted(set((s, e) for s, e, *_ in exons))
    introns = []
    for i in range(1, len(uniq_exons)):
        prev_end = uniq_exons[i - 1][1]
        this_start = uniq_exons[i][0]
        if this_start > prev_end + 1:
            introns.append((prev_end + 1, this_start - 1))  # 1-based inclusive gap
    if strand == "-":
        introns = list(reversed(introns))
    return introns


def select_intron(introns):
    """Excludes the first (transcription-order) intron; returns the
    largest of the remainder, or None if fewer than 2 introns exist."""
    if len(introns) < 2:
        return None
    remainder = introns[1:]
    return max(remainder, key=lambda iv: iv[1] - iv[0])


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff", required=True)
    p.add_argument("--trim", type=int, default=10,
                   help="bp trimmed from EACH end of the selected intron (default 10)")
    p.add_argument("--transcript-feature", default="mRNA")
    p.add_argument("--exon-feature", default="exon")
    p.add_argument("--out", required=True)
    p.add_argument("--exons-out", required=True,
                   help="ALL exons from ALL transcripts in the GFF (not just canonical "
                        "ones) - used both to subtract exonic sequence from intron "
                        "candidates (nested/overlapping genes are real and common; an "
                        "intron candidate can legitimately fall inside a DIFFERENT "
                        "gene's exon even though it never overlaps its OWN gene's exons "
                        "by construction) and for downstream validation")
    args = p.parse_args()

    transcript_parent, transcript_strand, exons_by_transcript = parse_gff(
        args.gff, args.transcript_feature, args.exon_feature)

    canonical = canonical_transcript_per_gene(transcript_parent, exons_by_transcript)

    n_genes = len(canonical)
    n_no_introns = 0        # 0 or 1 introns on the canonical transcript
    n_too_short_after_trim = 0
    n_written = 0

    # ALL exons, from every transcript of every gene - deliberately not
    # restricted to canonical transcripts, since this is used to subtract
    # exonic sequence genome-wide, and missing even a non-canonical
    # isoform's exon would let a block land on real coding/UTR sequence.
    all_exon_rows = set()
    for tid, exons in exons_by_transcript.items():
        for e_start, e_end, e_chrom in exons:
            gene_id = transcript_parent.get(tid, tid)
            all_exon_rows.add((e_chrom, e_start, e_end, gene_id))

    with open(args.exons_out, "w") as exons_out:
        for e_chrom, e_start, e_end, gene_id in sorted(all_exon_rows):
            exons_out.write(f"{e_chrom}\t{e_start - 1}\t{e_end}\t{gene_id}\n")

    with open(args.out, "w") as out:
        for gene_id, tid in canonical.items():
            exons = exons_by_transcript[tid]
            chrom = exons[0][2]
            strand = transcript_strand.get(tid, "+")
            introns = introns_from_exons(exons, strand)
            chosen = select_intron(introns)
            if chosen is None:
                n_no_introns += 1
                continue
            start_1based, end_1based = chosen
            # trim from each end, then convert 1-based inclusive -> 0-based half-open
            trimmed_start_1based = start_1based + args.trim
            trimmed_end_1based = end_1based - args.trim
            if trimmed_end_1based < trimmed_start_1based:
                n_too_short_after_trim += 1
                continue
            bed_start = trimmed_start_1based - 1
            bed_end = trimmed_end_1based
            out.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_id}\n")
            n_written += 1

    print(f"[=] {n_genes} genes with a canonical transcript", file=sys.stderr)
    print(f"[=] {n_no_introns} skipped (fewer than 2 introns - no non-first intron to choose)", file=sys.stderr)
    print(f"[=] {n_too_short_after_trim} skipped (selected intron shorter than 2x trim after trimming)", file=sys.stderr)
    print(f"[=] {n_written} candidate intron intervals written to {args.out}", file=sys.stderr)
    print(f"[=] {len(all_exon_rows)} exons (all transcripts, all genes) written to {args.exons_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
