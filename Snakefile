# =============================================================================
# Species-pair demographic pipeline
#
# Driven by two manifests (samples_manifest, pairs_manifest in config.yaml)
# rather than one config block per pair - works unchanged for a batch of 2
# pairs or a batch of 200. Reference-derived processing (gene buffer,
# repeat removal, genome chunking, indexing) is deduplicated per distinct
# reference genome; mapping is deduplicated per (sample, reference) pair -
# so if the same individual or the same reference genome is reused across
# many pairs, that work only happens once, however many pairs use it.
#
# See README.md for the full rationale and the CLI audit against gIMble's
# actual (modern) source.
# =============================================================================

import re
from pathlib import Path
import pandas as pd

configfile: "config.yaml"

OUT = config["output_dir"]

REGION_STRATEGY = config.get("region_strategy", "mixed")
if REGION_STRATEGY not in ("intergenic", "intron", "mixed"):
    raise WorkflowError(
        f"config.yaml's region_strategy must be 'intergenic', 'intron', or 'mixed' "
        f"- got {REGION_STRATEGY!r}")
USE_INTERGENIC = REGION_STRATEGY in ("intergenic", "mixed")
USE_INTRON = REGION_STRATEGY in ("intron", "mixed")

# -----------------------------------------------------------------------------
# Manifests
# -----------------------------------------------------------------------------
samples_df = pd.read_csv(config["samples_manifest"], sep="\t", dtype=str).fillna("")
samples_df = samples_df.set_index("sample_id", drop=False)
SAMPLES = list(samples_df.index)

pairs_df = pd.read_csv(config["pairs_manifest"], sep="\t", dtype=str).fillna("")
pairs_df = pairs_df.set_index("pair_id", drop=False)

def _sanitize(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)

def _ref_id(fasta_path):
    return _sanitize(Path(fasta_path).stem)

pairs_df["ref_id"] = pairs_df["reference_fasta"].apply(_ref_id)
PAIRS = list(pairs_df.index)

# One row per DISTINCT reference genome (dedup key = ref_id, i.e. the
# fasta's own filename - pairs sharing a reference file share this row).
REF_TABLE = pairs_df.drop_duplicates(subset="ref_id").set_index("ref_id", drop=False)
REF_IDS = list(REF_TABLE.index)

def ref_fasta(ref_id):    return REF_TABLE.loc[ref_id, "reference_fasta"]
def ref_gff(ref_id):      return REF_TABLE.loc[ref_id, "reference_gff"]
def ref_repeats(ref_id):  return REF_TABLE.loc[ref_id, "repeat_bed"]
def ref_include_sequences(ref_id):
    v = REF_TABLE.loc[ref_id, "include_sequences"]
    return v if v else None

def ref_id_of(pair):      return pairs_df.loc[pair, "ref_id"]
def sample_a_of(pair):    return pairs_df.loc[pair, "sampleA"]
def sample_b_of(pair):    return pairs_df.loc[pair, "sampleB"]

def sample_r1(sample): return samples_df.loc[sample, "r1"]
def sample_r2(sample):
    v = samples_df.loc[sample, "r2"]
    return v if v else None
def is_paired(sample): return sample_r2(sample) is not None

# -----------------------------------------------------------------------------
# Manifest validation - runs at parse time, before any job is scheduled.
#
# Catches the mistakes that would otherwise only surface deep into a run
# (after minutes of mapping/variant calling): missing files, sample IDs
# that don't exist in samples.tsv, and - most importantly - files that
# exist but are the WRONG KIND of file. Column-shift errors in pairs.tsv
# are easy to make (a single tab where two were needed shifts every
# subsequent field left by one) and otherwise fail late and cryptically,
# e.g. a chromosome allowlist landing in the repeat_bed column and being
# handed to `bedtools merge`, which needs >= 3 columns.
# -----------------------------------------------------------------------------
def _open_maybe_gzip(path):
    import gzip
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "r")

def _first_data_line(path, comment_prefixes=("#", "track", "browser")):
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith(comment_prefixes):
                return s
    return None

def _validate_manifests():
    errors = []

    required_sample_cols = {"sample_id", "r1"}
    missing = required_sample_cols - set(samples_df.columns)
    if missing:
        errors.append(f"samples manifest is missing column(s): {', '.join(sorted(missing))}")

    required_pair_cols = {"pair_id", "sampleA", "sampleB",
                          "reference_fasta", "reference_gff", "repeat_bed"}
    missing = required_pair_cols - set(pairs_df.columns)
    if missing:
        errors.append(f"pairs manifest is missing column(s): {', '.join(sorted(missing))}")
    if errors:
        raise WorkflowError("Manifest problems found:\n  - " + "\n  - ".join(errors))

    # --- samples: referenced fastqs exist -----------------------------------
    for sid in samples_df.index:
        for col in ("r1", "r2"):
            val = samples_df.loc[sid, col] if col in samples_df.columns else ""
            if val and not Path(val).exists():
                errors.append(f"sample '{sid}': {col} file not found: {val}")

    # --- pairs --------------------------------------------------------------
    for pair in pairs_df.index:
        row = pairs_df.loc[pair]

        for col in ("sampleA", "sampleB"):
            if row[col] not in samples_df.index:
                errors.append(
                    f"pair '{pair}': {col} '{row[col]}' is not a sample_id in the samples manifest")

        for col in ("reference_fasta", "reference_gff", "repeat_bed"):
            val = row[col]
            if not val:
                errors.append(f"pair '{pair}': {col} is empty (all three are required)")
            elif not Path(val).exists():
                errors.append(f"pair '{pair}': {col} file not found: {val}")

        inc = row.get("include_sequences", "")
        if inc and not Path(inc).exists():
            errors.append(f"pair '{pair}': include_sequences file not found: {inc}")

        # --- file-KIND checks: right path, wrong type of file ---------------
        rb = row["repeat_bed"]
        if rb and Path(rb).exists():
            line = _first_data_line(rb)
            if line is None:
                errors.append(f"pair '{pair}': repeat_bed appears to be empty: {rb}")
            else:
                fields = line.split("\t")
                if len(fields) < 3:
                    errors.append(
                        f"pair '{pair}': repeat_bed does not look like a BED file (needs at least "
                        f"3 tab-separated columns: chrom, start, end). Got {len(fields)} column(s) "
                        f"in first data line of {rb}: {line[:80]!r}. "
                        f"NOTE: if this is your chromosome allowlist, it belongs in the "
                        f"'include_sequences' column, not 'repeat_bed' - check for a missing tab "
                        f"shifting your columns.")
                else:
                    try:
                        int(fields[1]); int(fields[2])
                    except ValueError:
                        errors.append(
                            f"pair '{pair}': repeat_bed columns 2 and 3 should be numeric "
                            f"(start, end). Got {fields[1]!r}, {fields[2]!r} in {rb}")

        if inc and Path(inc).exists():
            line = _first_data_line(inc)
            if line is not None and len(line.split("\t")) > 1:
                errors.append(
                    f"pair '{pair}': include_sequences should be ONE sequence name per line, but "
                    f"the first line of {inc} has multiple tab-separated columns: {line[:80]!r}. "
                    f"NOTE: if this is a BED file, you may have swapped it with 'repeat_bed'.")

        gff = row["reference_gff"]
        if gff and Path(gff).exists():
            line = _first_data_line(gff, comment_prefixes=("#",))
            if line is not None and len(line.split("\t")) < 9:
                errors.append(
                    f"pair '{pair}': reference_gff does not look like a GFF (expects 9 "
                    f"tab-separated columns). First data line of {gff} has "
                    f"{len(line.split(chr(9)))}: {line[:80]!r}")

    if errors:
        raise WorkflowError(
            "Manifest validation failed ({} problem(s)):\n  - ".format(len(errors))
            + "\n  - ".join(errors))

_validate_manifests()

# -----------------------------------------------------------------------------
# 0. Top-level target: every pair's final bSFS + info file
# -----------------------------------------------------------------------------
rule all:
    input:
        expand(f"{OUT}/bsfs/{{pair}}.bsfs.csv", pair=PAIRS),
        expand(f"{OUT}/gimble/{{pair}}/info.txt", pair=PAIRS),
        expand(f"{OUT}/gimble/{{pair}}/block_validation.txt", pair=PAIRS)

# -----------------------------------------------------------------------------
# 1. Read QC/trimming - keyed by {sample} only (no reference dependency,
#    so shared automatically across every pair that uses this sample)
# -----------------------------------------------------------------------------
ruleorder: fastp_pe > fastp_se

rule fastp_pe:
    input:
        r1 = lambda wc: sample_r1(wc.sample),
        r2 = lambda wc: sample_r2(wc.sample),
        ref_validation = lambda wc: [ancient(f) for f in ref_validations_for_sample(wc.sample)],  # hard gate (ancient: gate only, ignore its mtime)
    output:
        r1 = f"{OUT}/cleaned/{{sample}}_1.cleaned.fastq.gz",
        r2 = f"{OUT}/cleaned/{{sample}}_2.cleaned.fastq.gz",
        html = f"{OUT}/cleaned/{{sample}}.fastp.html",
        json = f"{OUT}/cleaned/{{sample}}.fastp.json",
    log: f"{OUT}/logs/fastp.{{sample}}.log"
    conda: "envs/mapping.yaml"
    threads: 4
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "-h {output.html} -j {output.json} -w {threads} > {log} 2>&1"

rule fastp_se:
    input:
        r1 = lambda wc: sample_r1(wc.sample),
        ref_validation = lambda wc: [ancient(f) for f in ref_validations_for_sample(wc.sample)],  # hard gate (ancient: gate only, ignore its mtime)
    output:
        r1 = f"{OUT}/cleaned/{{sample}}_1.cleaned.fastq.gz",
        html = f"{OUT}/cleaned/{{sample}}.fastp.html",
        json = f"{OUT}/cleaned/{{sample}}.fastp.json",
    log: f"{OUT}/logs/fastp.{{sample}}.log"
    conda: "envs/mapping.yaml"
    threads: 4
    shell:
        "fastp -i {input.r1} -o {output.r1} "
        "-h {output.html} -j {output.json} -w {threads} > {log} 2>&1"

# -----------------------------------------------------------------------------
# 2. Reference staging and indexing - keyed by {ref_id}, so a reference
#    used by many pairs is only staged/indexed once.
# -----------------------------------------------------------------------------
rule stage_reference:
    # Symlinks (not copies) the reference into our own results tree, so
    # indexing never needs write access to wherever your original genome
    # file lives (often a shared, read-only directory on HPC systems).
    input: lambda wc: ref_fasta(wc.ref_id)
    output: f"{OUT}/reference/{{ref_id}}/ref.fasta"
    shell: "ln -sf $(realpath {input}) {output}"

rule bwa_index:
    input: rules.stage_reference.output
    output: multiext(f"{OUT}/reference/{{ref_id}}/ref.fasta", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123")
    conda: "envs/mapping.yaml"
    log: f"{OUT}/logs/bwa_index.{{ref_id}}.log"
    shell: "bwa-mem2 index {input} > {log} 2>&1"

rule samtools_faidx:
    input: rules.stage_reference.output
    output: f"{OUT}/reference/{{ref_id}}/ref.fasta.fai"
    conda: "envs/mapping.yaml"
    shell: "samtools faidx {input}"

rule genome_file:
    # 2-column (chrom, length) file bedtools needs for slop/complement.
    # If this reference's include_sequences is set, restricts to only
    # those sequence IDs - excluding e.g. unplaced scaffolds/contigs from
    # block selection. Every BED used downstream (gene_bed, the repeat
    # BED) is filtered against this same list.
    input: rules.samtools_faidx.output
    output: f"{OUT}/reference/{{ref_id}}/genome.txt"
    params:
        include_filter = lambda wc: (f"grep -Fwf {ref_include_sequences(wc.ref_id)}"
                                      if ref_include_sequences(wc.ref_id) else "cat")
    shell: "cut -f1,2 {input} | sort -k1,1 | {params.include_filter} > {output}"

# -----------------------------------------------------------------------------
# 3. Mapping, sorting, dedup, read groups - keyed by {ref_id}/{sample}, so
#    the same individual mapped against the same reference in multiple
#    pairs is only ever mapped once.
# -----------------------------------------------------------------------------
rule bwa_map_pe:
    input:
        r1 = f"{OUT}/cleaned/{{sample}}_1.cleaned.fastq.gz",
        r2 = f"{OUT}/cleaned/{{sample}}_2.cleaned.fastq.gz",
        idx = rules.bwa_index.output,
        ref = rules.stage_reference.output,
        ref_validation = lambda wc: [ancient(f) for f in ref_validations_for_ref(wc.ref_id)],  # hard gate (ancient: gate only, ignore its mtime)
    output: temp(f"{OUT}/mapped/{{ref_id}}/{{sample}}.bam")
    conda: "envs/mapping.yaml"
    threads: config["threads"]["mapping"]
    log: f"{OUT}/logs/bwa_map.{{ref_id}}.{{sample}}.log"
    shell:
        "bwa-mem2 mem -t {threads} {input.ref} {input.r1} {input.r2} 2> {log} "
        "| samtools view -b -o {output} -"

rule bwa_map_se:
    input:
        r1 = f"{OUT}/cleaned/{{sample}}_1.cleaned.fastq.gz",
        idx = rules.bwa_index.output,
        ref = rules.stage_reference.output,
        ref_validation = lambda wc: [ancient(f) for f in ref_validations_for_ref(wc.ref_id)],  # hard gate (ancient: gate only, ignore its mtime)
    output: temp(f"{OUT}/mapped/{{ref_id}}/{{sample}}.se.bam")
    conda: "envs/mapping.yaml"
    threads: config["threads"]["mapping"]
    log: f"{OUT}/logs/bwa_map.{{ref_id}}.{{sample}}.log"
    shell:
        "bwa-mem2 mem -t {threads} {input.ref} {input.r1} 2> {log} "
        "| samtools view -b -o {output} -"

def mapped_bam(wc):
    suffix = "bam" if is_paired(wc.sample) else "se.bam"
    return f"{OUT}/mapped/{wc.ref_id}/{wc.sample}.{suffix}"

rule sort_bam:
    input: mapped_bam
    output: temp(f"{OUT}/sorted/{{ref_id}}/{{sample}}.bam")
    conda: "envs/mapping.yaml"
    threads: 8
    shell: "sambamba sort -t {threads} -o {output} {input}"

rule markdup:
    input: rules.sort_bam.output
    output: temp(f"{OUT}/markdup/{{ref_id}}/{{sample}}.bam")
    conda: "envs/mapping.yaml"
    threads: 8
    shell: "sambamba markdup -t {threads} -r {input} {output}"

rule add_readgroups:
    # Direct picard call (not the Snakemake wrapper repo, which needs
    # internet access at runtime and is a poor fit for offline HPC use).
    input: rules.markdup.output
    output: f"{OUT}/final_bams/{{ref_id}}/{{sample}}.bam"
    conda: "envs/mapping.yaml"
    log: f"{OUT}/logs/readgroups.{{ref_id}}.{{sample}}.log"
    shell:
        "picard AddOrReplaceReadGroups I={input} O={output} "
        "RGID={wildcards.sample} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={wildcards.sample} "
        "> {log} 2>&1"

rule index_bam:
    input: rules.add_readgroups.output
    output: f"{OUT}/final_bams/{{ref_id}}/{{sample}}.bam.bai"
    conda: "envs/mapping.yaml"
    shell: "sambamba index {input} {output}"

# -----------------------------------------------------------------------------
# 4. Intergenic region construction - keyed by {ref_id}
# -----------------------------------------------------------------------------
rule gene_bed:
    input:
        gff = lambda wc: ref_gff(wc.ref_id),
        genome = rules.genome_file.output,
    output: f"{OUT}/reference/{{ref_id}}/genes.bed"
    params: feature = config["intergenic"]["gff_feature"]
    shell:
        # zcat -f transparently handles gzip-compressed OR plain-text input
        # (passes plain files through unchanged) - the GFF may legitimately
        # be either.
        # GFF is 1-based inclusive; BED is 0-based half-open, hence start-1.
        #
        # CRITICAL GUARD: an empty genes.bed is NOT a harmless edge case.
        # `bedtools complement` on empty input reports the ENTIRE genome as
        # intergenic, so the pipeline would silently proceed with zero gene
        # buffering and draw blocks from genic sequence - invalidating the
        # analysis with no error anywhere. We fail hard here instead, and
        # report which of the two possible causes it is.
        r"""tmp=$(mktemp)
        trap 'rm -f "$tmp"' EXIT
        zcat -f {input.gff} \
          | awk -F'\t' -v feat="{params.feature}" '$3==feat {{OFS="\t"; print $1, $4-1, $5}}' \
          | sort -k1,1 -k2,2n > "$tmp"
        n_pre=$(wc -l < "$tmp")
        awk 'NR==FNR{{keep[$1];next}} ($1 in keep)' {input.genome} "$tmp" > {output}
        n_post=$(wc -l < {output})
        if [ "$n_post" -eq 0 ]; then
          echo "" >&2
          echo "[X] ERROR: no '{params.feature}' features survived for reference '{wildcards.ref_id}'." >&2
          echo "    An empty gene set would make bedtools report the WHOLE GENOME as" >&2
          echo "    intergenic, silently removing all gene buffering. Refusing to continue." >&2
          echo "" >&2
          if [ "$n_pre" -eq 0 ]; then
            echo "    CAUSE: the GFF contains no rows with '{params.feature}' in column 3." >&2
            echo "    Check what feature types it actually has, and set intergenic.gff_feature" >&2
            echo "    in config.yaml to match:" >&2
            echo "      zcat -f {input.gff} | awk -F'\t' '!/^#/{{print \$3}}' | sort | uniq -c | sort -rn | head" >&2
          else
            echo "    CAUSE: found $n_pre '{params.feature}' features, but NONE are on sequences" >&2
            echo "    present in the reference. The GFF's sequence names do not match the FASTA's." >&2
            echo "" >&2
            echo "    Sequence names in the GFF:" >&2
            cut -f1 "$tmp" | sort -u | head -5 | sed 's/^/      /' >&2
            echo "    Sequence names in the reference ({input.genome}):" >&2
            cut -f1 {input.genome} | head -5 | sed 's/^/      /' >&2
            echo "" >&2
            echo "    (If include_sequences is set, also confirm those names match both files.)" >&2
          fi
          echo "" >&2
          exit 1
        fi"""

rule gene_buffer_bed:
    input:
        genes = rules.gene_bed.output,
        genome = rules.genome_file.output,
    output: f"{OUT}/reference/{{ref_id}}/genes.buffered.merged.bed"
    params: buffer = config["intergenic"]["gene_buffer"]
    conda: "envs/bedtools.yaml"
    shell:
        "bedtools slop -i {input.genes} -g {input.genome} -b {params.buffer} "
        "| sort -k1,1 -k2,2n | bedtools merge -i - > {output}"

rule intergenic_candidate_bed:
    input:
        buffered_genes = rules.gene_buffer_bed.output,
        genome = rules.genome_file.output,
    output: f"{OUT}/reference/{{ref_id}}/intergenic_candidate.bed"
    conda: "envs/bedtools.yaml"
    shell: "bedtools complement -i {input.buffered_genes} -g {input.genome} > {output}"

rule repeat_free_intergenic_bed:
    input:
        intergenic = rules.intergenic_candidate_bed.output,
        repeats = lambda wc: ref_repeats(wc.ref_id),
        genome = rules.genome_file.output,
    output: f"{OUT}/reference/{{ref_id}}/intergenic.repeat_free.bed"
    conda: "envs/bedtools.yaml"
    shell:
        # zcat -f: same reasoning as gene_bed - repeats.bed may be gzipped or plain.
        "zcat -f {input.repeats} | sort -k1,1 -k2,2n | bedtools merge -i - "
        "| awk 'NR==FNR{{keep[$1];next}} ($1 in keep)' {input.genome} - "
        "| bedtools subtract -a {input.intergenic} -b - > {output}"

rule validate_reference_intergenic:
    # Independent, post-hoc validation of the reference-level intergenic
    # construction (gene buffering + repeat removal) - see
    # scripts/validate_reference.py. Deliberately placed here, immediately
    # after these BED files are built and BEFORE any sample-level work:
    # this construction depends only on the reference/GFF/repeat BED, so
    # there's no reason to wait until after mapping+calling to find out it
    # was wrong. Gates fastp_pe/fastp_se and bwa_map_pe/bwa_map_se below.
    input:
        genes = rules.gene_bed.output,
        genes_buffered = rules.gene_buffer_bed.output,
        repeats = lambda wc: ref_repeats(wc.ref_id),
        intergenic_repeat_free = rules.repeat_free_intergenic_bed.output,
        genome = rules.genome_file.output,
    output: f"{OUT}/reference/{{ref_id}}/reference_validation.txt"
    conda: "envs/bedtools.yaml"
    shell:
        "python3 scripts/validate_reference.py --genes {input.genes} "
        "--genes-buffered {input.genes_buffered} --repeats {input.repeats} "
        "--intergenic-repeat-free {input.intergenic_repeat_free} --genome {input.genome} "
        "--report {output}"

# -----------------------------------------------------------------------------
# 4b. Intron-based candidate region construction (only used if
#     region_strategy is "intron" or "mixed") - mirrors the intergenic
#     chain above, but candidate regions come from inside gene bodies
#     (one intron per gene) instead of outside them. See
#     scripts/extract_introns.py for the exact selection rule (largest
#     non-first intron, strand-aware, trimmed) and why "overlaps the
#     gene+buffer zone" is expected here, not a failure.
# -----------------------------------------------------------------------------
rule extract_introns:
    input:
        gff = lambda wc: ref_gff(wc.ref_id),
        genome = rules.genome_file.output,
    output:
        introns = f"{OUT}/reference/{{ref_id}}/introns.bed",
        exons = f"{OUT}/reference/{{ref_id}}/exons_used.bed",
    params:
        trim = config["intron"]["trim"],
        transcript_feature = config["intron"]["transcript_feature"],
        exon_feature = config["intron"]["exon_feature"],
    shell:
        # zcat -f handled inside extract_introns.py itself (same
        # gzip-or-plain convention as every other GFF/BED read in this
        # pipeline). include_sequences filtering happens in the next rule
        # (repeat_free_introns_bed), same as the intergenic chain does it
        # for gene_bed's output there rather than here.
        "python3 scripts/extract_introns.py --gff {input.gff} --trim {params.trim} "
        "--transcript-feature {params.transcript_feature} --exon-feature {params.exon_feature} "
        "--out {output.introns} --exons-out {output.exons}"

rule repeat_free_introns_bed:
    input:
        introns = rules.extract_introns.output.introns,
        exons = rules.extract_introns.output.exons,
        repeats = lambda wc: ref_repeats(wc.ref_id),
        genome = rules.genome_file.output,
    output: f"{OUT}/reference/{{ref_id}}/introns.repeat_and_exon_free.bed"
    conda: "envs/bedtools.yaml"
    shell:
        # Repeat subtraction: same include_sequences filtering + repeat
        # subtraction as repeat_free_intergenic_bed.
        # Exon subtraction: an intron candidate never overlaps its OWN
        # gene's exons by construction, but real, reasonably common
        # nested/overlapping genes mean it CAN legitimately overlap a
        # DIFFERENT gene's exon - confirmed directly against a real
        # genome during this pipeline's development (see README), not
        # hypothetical. exons_used.bed covers every transcript of every
        # gene genome-wide (not just canonical transcripts), so this
        # subtracts all annotated exonic sequence, not just the host
        # gene's own.
        # Note introns.bed carries a 4th column (gene_id) from
        # extract_introns.py - bedtools subtract preserves -a's extra
        # columns by default (confirmed for intersect earlier in this
        # pipeline's history; subtract behaves the same way), so that
        # traceability column survives into the output here.
        "sort -k1,1 -k2,2n {input.introns} > {output}.sorted.tmp && "
        "zcat -f {input.repeats} | sort -k1,1 -k2,2n | bedtools merge -i - "
        "| awk 'NR==FNR{{keep[$1];next}} ($1 in keep)' {input.genome} - "
        "| bedtools subtract -a {output}.sorted.tmp -b - "
        "| sort -k1,1 -k2,2n "
        "| bedtools subtract -a - -b {input.exons} > {output} && "
        "rm -f {output}.sorted.tmp"

rule validate_reference_introns:
    # Independent, post-hoc validation of the reference-level intron
    # construction - see scripts/validate_introns.py. Same placement
    # rationale as validate_reference_intergenic: depends only on the
    # reference/GFF/repeat BED, so it runs before any sample-level work.
    input:
        introns = rules.extract_introns.output.introns,
        introns_repeat_free = rules.repeat_free_introns_bed.output,
        exons = rules.extract_introns.output.exons,
        repeats = lambda wc: ref_repeats(wc.ref_id),
    output: f"{OUT}/reference/{{ref_id}}/reference_validation_introns.txt"
    conda: "envs/bedtools.yaml"
    shell:
        "python3 scripts/validate_introns.py --introns {input.introns} "
        "--introns-repeat-free {input.introns_repeat_free} --exons {input.exons} "
        "--repeats {input.repeats} --report {output}"

# -----------------------------------------------------------------------------
# 4c. Unify whichever candidate pool(s) region_strategy selects into one
#     region a pair's callable sites get intersected against. Untagged by
#     design (not "which source did this region come from") - see the
#     Snakefile's own commit history / README for why: bedtools intersect
#     only preserves the -a side's extra columns, so a source tag placed
#     on these regions (the -b side, once intersected against gimbleprep's
#     callable sites) would silently vanish. Spacing is enforced ONCE,
#     globally, across this unified pool by select_spaced_blocks.py below
#     - not computed separately per source and concatenated, which would
#     allow an intergenic block and an intron block to end up closer
#     together than block_spacing allows.
# -----------------------------------------------------------------------------
rule region_candidate_bed:
    input:
        intergenic = rules.repeat_free_intergenic_bed.output if USE_INTERGENIC else [],
        introns = rules.repeat_free_introns_bed.output if USE_INTRON else [],
    output: f"{OUT}/reference/{{ref_id}}/region_candidate.bed"
    shell:
        # cut -f1-3: introns.repeat_and_exon_free.bed carries a 4th (gene_id,
        # traceability-only) column that intergenic.repeat_free.bed
        # doesn't; concatenating without normalising column count first
        # produces a BED with inconsistent field counts across rows,
        # which bedtools correctly refuses to parse ("wrong number of
        # fields") - confirmed directly, not assumed, by hitting exactly
        # this error when testing the mixed-strategy path end to end.
        "cat {input.intergenic} {input.introns} | cut -f1-3 | sort -k1,1 -k2,2n > {output}"

def ref_validations_for_sample(sample):
    # Every reference whose pair(s) use this sample - in the common case
    # (no sample reused across pairs against different references) this
    # is a single file; handled generally in case that ever changes.
    # Includes whichever reference-level validation(s) region_strategy
    # actually uses - both, for "mixed".
    used_refs = {ref_id_of(pair) for pair in PAIRS
                 if sample_a_of(pair) == sample or sample_b_of(pair) == sample}
    out = []
    if USE_INTERGENIC:
        out += [f"{OUT}/reference/{r}/reference_validation.txt" for r in used_refs]
    if USE_INTRON:
        out += [f"{OUT}/reference/{r}/reference_validation_introns.txt" for r in used_refs]
    return out

def ref_validations_for_ref(ref_id):
    # Same as above, but for rules keyed directly by {ref_id} (bwa_map_pe/se)
    # rather than by {sample}.
    out = []
    if USE_INTERGENIC:
        out.append(f"{OUT}/reference/{ref_id}/reference_validation.txt")
    if USE_INTRON:
        out.append(f"{OUT}/reference/{ref_id}/reference_validation_introns.txt")
    return out

# -----------------------------------------------------------------------------
# 5. Variant calling (freebayes), chunked per reference genome for
#    scalability to large (non-Drosophila-sized) genomes.
# -----------------------------------------------------------------------------
rule freebayes_regions:
    input: rules.samtools_faidx.output
    output: f"{OUT}/reference/{{ref_id}}/regions.txt"
    params: n_chunks = config["freebayes"]["n_chunks"]
    shell:
        """
        python3 -c "
import sys
fai = [l.split() for l in open('{input}')]
contigs = [(c[0], int(c[1])) for c in fai]
total = sum(l for _, l in contigs)
chunk_size = max(1, total // {params.n_chunks})
with open('{output}', 'w') as out:
    for chrom, length in contigs:
        start = 0
        while start < length:
            end = min(start + chunk_size, length)
            out.write(f'{{chrom}}:{{start}}-{{end}}\\n')
            start = end
"
        """

CHUNKS = list(range(config["freebayes"]["n_chunks"]))

rule freebayes_call_chunk:
    input:
        ref = lambda wc: f"{OUT}/reference/{ref_id_of(wc.pair)}/ref.fasta",
        fai = lambda wc: f"{OUT}/reference/{ref_id_of(wc.pair)}/ref.fasta.fai",
        bamA = lambda wc: f"{OUT}/final_bams/{ref_id_of(wc.pair)}/{sample_a_of(wc.pair)}.bam",
        bamB = lambda wc: f"{OUT}/final_bams/{ref_id_of(wc.pair)}/{sample_b_of(wc.pair)}.bam",
        baiA = lambda wc: f"{OUT}/final_bams/{ref_id_of(wc.pair)}/{sample_a_of(wc.pair)}.bam.bai",
        baiB = lambda wc: f"{OUT}/final_bams/{ref_id_of(wc.pair)}/{sample_b_of(wc.pair)}.bam.bai",
        regions = lambda wc: f"{OUT}/reference/{ref_id_of(wc.pair)}/regions.txt",
    output: temp(f"{OUT}/calls/{{pair}}/chunk_{{chunk}}.vcf")
    params:
        ploidy = config["freebayes"]["ploidy"],
        theta = config["freebayes"]["theta"],
    conda: "envs/variant_calling.yaml"
    threads: config["threads"]["freebayes_chunk"]
    shell:
        """
        region=$(awk -v n={wildcards.chunk} 'NR==n+1{{print}}' {input.regions})
        freebayes -f {input.ref} --haplotype-length -1 --no-population-priors \
          --hwe-priors-off --use-mapping-quality --ploidy {params.ploidy} \
          --theta {params.theta} --region "$region" {input.bamA} {input.bamB} > {output}
        """

rule concat_vcf:
    input: expand(f"{OUT}/calls/{{{{pair}}}}/chunk_{{chunk}}.vcf", chunk=CHUNKS)
    output: f"{OUT}/calls/{{pair}}.combined.vcf"
    conda: "envs/variant_calling.yaml"
    shell: "bcftools concat {input} | bcftools sort -o {output} -"

rule bgzip_index_vcf:
    input: rules.concat_vcf.output
    output:
        gz = f"{OUT}/calls/{{pair}}.combined.vcf.gz",
        tbi = f"{OUT}/calls/{{pair}}.combined.vcf.gz.tbi",
    conda: "envs/variant_calling.yaml"
    shell: "bgzip -c {input} > {output.gz} && tabix -p vcf {output.gz}"

# -----------------------------------------------------------------------------
# 6. gimbleprep, then a pair-specific BAM directory for it - gimbleprep's
#    -b flag takes a whole directory, so if two pairs share a reference
#    but not both samples (e.g. ambigua_obscura and ambigua_persimilis
#    both use the same reference), a shared final_bams/{ref_id}/ directory
#    would incorrectly hand gimbleprep a third, irrelevant sample's BAM.
#    This stages a clean, pair-specific view via symlinks - no extra copies.
# -----------------------------------------------------------------------------
rule stage_pair_bams:
    input:
        bamA = lambda wc: f"{OUT}/final_bams/{ref_id_of(wc.pair)}/{sample_a_of(wc.pair)}.bam",
        baiA = lambda wc: f"{OUT}/final_bams/{ref_id_of(wc.pair)}/{sample_a_of(wc.pair)}.bam.bai",
        bamB = lambda wc: f"{OUT}/final_bams/{ref_id_of(wc.pair)}/{sample_b_of(wc.pair)}.bam",
        baiB = lambda wc: f"{OUT}/final_bams/{ref_id_of(wc.pair)}/{sample_b_of(wc.pair)}.bam.bai",
    output: directory(f"{OUT}/gimble/{{pair}}/bam_dir")
    run:
        import os
        os.makedirs(output[0], exist_ok=True)
        for src, name in [(input.bamA, f"{sample_a_of(wildcards.pair)}.bam"),
                           (input.baiA, f"{sample_a_of(wildcards.pair)}.bam.bai"),
                           (input.bamB, f"{sample_b_of(wildcards.pair)}.bam"),
                           (input.baiB, f"{sample_b_of(wildcards.pair)}.bam.bai")]:
            dst = os.path.join(output[0], name)
            if not os.path.exists(dst):
                os.symlink(os.path.realpath(src), dst)

rule gimbleprep:
    input:
        ref = lambda wc: f"{OUT}/reference/{ref_id_of(wc.pair)}/ref.fasta",
        vcf = rules.bgzip_index_vcf.output.gz,
        bam_dir = rules.stage_pair_bams.output,
    output:
        vcf = f"{OUT}/gimble/{{pair}}/prep.vcf.gz",
        bed = f"{OUT}/gimble/{{pair}}/prep.bed",
        genomefile = f"{OUT}/gimble/{{pair}}/prep.genomefile",
        samples = f"{OUT}/gimble/{{pair}}/prep.samples.csv",
        coverage = f"{OUT}/gimble/{{pair}}/prep.coverage_summary.csv",
    params:
        min_depth = config["gimble"]["min_depth"],
        max_depth_flag = (f"-M {config['gimble']['max_depth']}"
                          if config["gimble"].get("max_depth") is not None else ""),
        min_qual = config["gimble"]["min_qual"],
        snpgap = config["gimble"]["snpgap"],
        outprefix = f"{OUT}/gimble/{{pair}}/prep",
    conda: "envs/gimbleprep.yaml"
    threads: config["threads"]["gimble"]
    shell:
        "gimbleprep -f {input.ref} -v {input.vcf} -b {input.bam_dir} "
        "-m {params.min_depth} {params.max_depth_flag} -q {params.min_qual} "
        "-g {params.snpgap} -t {threads} -o {params.outprefix}"

rule filter_prep_bed_to_allowed_sequences:
    # gimbleprep necessarily runs genome-wide (mapping/calling correctly
    # still cover unplaced scaffolds, to avoid reads that belong there
    # mismapping onto real chromosomes instead), so its own prep.bed
    # output is NOT restricted to include_sequences at all. The final
    # callable_intergenic_bed step was already implicitly correct (a
    # chromosome present only in prep.bed and absent from the
    # already-restricted intergenic.repeat_free.bed simply cannot
    # intersect with anything - confirmed directly, not assumed), but
    # that restriction was invisible: you could not look at any single
    # file and confirm scaffolds were actually excluded, and prep.bed
    # itself stayed bloated with scaffold-derived callable sites that
    # never contributed to the final blocks anyway. This makes the
    # restriction explicit, inspectable, and reduces file size directly.
    input:
        bed = rules.gimbleprep.output.bed,
        genome = lambda wc: f"{OUT}/reference/{ref_id_of(wc.pair)}/genome.txt",
    output: f"{OUT}/gimble/{{pair}}/prep.allowed_sequences.bed"
    shell:
        "awk 'NR==FNR{{keep[$1];next}} ($1 in keep)' {input.genome} {input.bed} > {output}"

rule callable_regions_bed:
    # Renamed from callable_intergenic_bed - with region_strategy possibly
    # including introns, "intergenic" would no longer be an accurate name
    # for what this file can contain.
    input:
        callable = rules.filter_prep_bed_to_allowed_sequences.output,
        region = lambda wc: f"{OUT}/reference/{ref_id_of(wc.pair)}/region_candidate.bed",
    output: f"{OUT}/gimble/{{pair}}/callable_regions.bed"
    conda: "envs/bedtools.yaml"
    shell: "bedtools intersect -a {input.callable} -b {input.region} > {output}"

rule select_spaced_blocks:
    # Enforces block spacing (>= config's block_spacing between any two
    # selected blocks, GLOBALLY across whichever region(s) region_strategy
    # includes) - see scripts/select_spaced_blocks.py. The >=100kb-from-gene
    # criterion is only meaningful/applied for intergenic-derived regions
    # (already baked into region_candidate_bed upstream via gene buffering);
    # intron-derived regions are deliberately exempt - see config.yaml.
    input: bed = rules.callable_regions_bed.output
    output: f"{OUT}/gimble/{{pair}}/spaced_blocks.bed"
    params:
        block_length = config["gimble"]["block_length"],
        spacing = config["intergenic"]["block_spacing"],
    shell:
        "python3 scripts/select_spaced_blocks.py --bed {input.bed} "
        "--block-length {params.block_length} --spacing {params.spacing} --out {output}"

rule validate_blocks:
    # Independent, post-hoc check that the final blocks actually satisfy
    # the properties appropriate to region_strategy (see
    # scripts/validate_blocks.py for exactly what differs between
    # strategies and why). This is a hard gate: gimble_parse depends on
    # this rule's output, so a failure here stops the pipeline before any
    # bad blocks reach gIMble.
    input:
        blocks = rules.select_spaced_blocks.output,
        genes_buffered = (lambda wc: f"{OUT}/reference/{ref_id_of(wc.pair)}/genes.buffered.merged.bed") if REGION_STRATEGY == "intergenic" else [],
        exons = (lambda wc: f"{OUT}/reference/{ref_id_of(wc.pair)}/exons_used.bed") if REGION_STRATEGY in ("intron", "mixed") else [],
        repeats = lambda wc: ref_repeats(ref_id_of(wc.pair)),
    output: f"{OUT}/gimble/{{pair}}/block_validation.txt"
    params:
        block_spacing = config["intergenic"]["block_spacing"],
        strategy = REGION_STRATEGY,
        genes_buffered_flag = (lambda wc, input: f"--genes-buffered {input.genes_buffered}"
                                if REGION_STRATEGY == "intergenic" else ""),
        exons_flag = (lambda wc, input: f"--exons {input.exons}"
                      if REGION_STRATEGY in ("intron", "mixed") else ""),
    conda: "envs/bedtools.yaml"
    shell:
        "python3 scripts/validate_blocks.py --blocks {input.blocks} --strategy {params.strategy} "
        "{params.genes_buffered_flag} {params.exons_flag} --repeats {input.repeats} "
        "--block-spacing {params.block_spacing} --report {output}"

rule samples_csv:
    # gimbleprep writes a single-column sample list; gimble parse requires
    # a second column of population labels (exactly two populations).
    input: rules.gimbleprep.output.samples
    output: f"{OUT}/gimble/{{pair}}/samples.csv"
    run:
        with open(input[0]) as fh:
            sample_ids = [line.strip() for line in fh if line.strip()]
        if len(sample_ids) != 2:
            raise ValueError(
                f"Expected exactly 2 samples in {input[0]}, found {len(sample_ids)}: {sample_ids}. "
                "Check that both BAMs' read groups (RGSM) matched their VCF sample names."
            )
        pop_label = {sample_ids[0]: "A", sample_ids[1]: "B"}
        with open(output[0], "w") as out:
            for sid in sample_ids:
                out.write(f"{sid},{pop_label[sid]}\n")

rule gimble_parse:
    input:
        vcf = rules.gimbleprep.output.vcf,
        bed = rules.select_spaced_blocks.output,
        # FIX: use our own already-filtered genome.txt, not gimbleprep's
        # genomefile - the latter reflects the WHOLE reference (every
        # unplaced scaffold), and gimble parse's variant-reading step
        # issues one separate VCF query PER SEQUENCE in this file,
        # regardless of whether that sequence has any BED intervals at
        # all (confirmed directly in gIMble's source, lib/gimble.py
        # _read_variants). Using the full genomefile meant thousands of
        # wasted per-sequence VCF queries for every unplaced scaffold -
        # the entire cause of multi-hour gimble_parse runtimes.
        genomefile = lambda wc: f"{OUT}/reference/{ref_id_of(wc.pair)}/genome.txt",
        samples = rules.samples_csv.output,
        validation = ancient(rules.validate_blocks.output),  # hard gate (ancient: gate only, ignore its mtime)
    output: directory(f"{OUT}/gimble/{{pair}}/{{pair}}.z")
    params:
        # FIX: gimble parse's -z is a PREFIX, not the final path - it
        # appends its own ".z" itself (confirmed in gIMble's source:
        # Store.path = "%s.z" % prefix). Passing a value that already
        # ends in ".z" creates "{pair}.z.z", which Snakemake can never
        # find under the "{pair}.z" name it expects -> MissingOutputException.
        # No trailing .z here.
        prefix = f"{OUT}/gimble/{{pair}}/{{pair}}",
    shell:
        "gimble parse -g {input.genomefile} -v {input.vcf} -b {input.bed} "
        "-s {input.samples} -z {params.prefix} -f"

rule gimble_blocks:
    input: rules.gimble_parse.output
    output: touch(f"{OUT}/gimble/{{pair}}/.blocks_done")
    params:
        block_length = config["gimble"]["block_length"],
        block_span = config["gimble"]["block_span"],
    shell: "gimble blocks -z {input} -l {params.block_length} -m {params.block_span} -f"

rule gimble_tally:
    # NOTE: gimble tally's own default kmax is 2 (see config.yaml comment) -
    # this MUST be raised above the largest per-category mutation count in
    # your real data, or gimble will silently marginalize (truncate) it.
    input:
        z = rules.gimble_parse.output,
        done = rules.gimble_blocks.output,
    output: touch(f"{OUT}/gimble/{{pair}}/.tally_done")
    params:
        kmax = ",".join([str(config["gimble"]["kmax"])] * 4),
        label = config["gimble"]["tally_label"],
    shell: "gimble tally -z {input.z} -t blocks -l {params.label} -k {params.kmax} -f"

rule gimble_info:
    input:
        z = rules.gimble_parse.output,
        done = rules.gimble_blocks.output,
    output: f"{OUT}/gimble/{{pair}}/info.txt"
    shell: "gimble info -z {input.z} > {output}"

rule gimble_query_bsfs:
    input:
        z = rules.gimble_parse.output,
        done = rules.gimble_tally.output,
    output: f"{OUT}/gimble/{{pair}}/tally_{config['gimble']['tally_label']}.tsv"
    params: label = config["gimble"]["tally_label"]
    shell: "gimble query -z {input.z} -l tally/{params.label}"

# -----------------------------------------------------------------------------
# 7. Final output: ready for demogfit::bsfs_to_s_distribution()
#    (column order from gimble is count,het_b,het_a,het_ab,fixed_diff - the
#    het_a/het_b swap relative to demogfit's hetA/hetB is harmless, since
#    they are summed together before use - see demogfit source)
# -----------------------------------------------------------------------------
rule bsfs_to_demogfit_csv:
    input: rules.gimble_query_bsfs.output
    output: f"{OUT}/bsfs/{{pair}}.bsfs.csv"
    shell:
        r"""echo "count,hetA,hetB,hetAB,fixed" > {output}
            tail -n +2 {input} | tr '\t' ',' >> {output}"""
