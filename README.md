# Species-pair demographic pipeline

`fastq (speciesA, speciesB) + reference genome + GFF + repeat BED` →
mapping → variant calling → intergenic block selection → blockwise SFS,
ready for the [`demogfit`](../demogfit) R package.

## Running many pairs at once

The pipeline is driven by **two manifests**, not a config block per pair -
this is what lets it scale from 2 pairs to 200 without editing anything
per run:

- **`samples.tsv`** - one row per unique sequenced individual:
  ```
  sample_id  r1  r2
  ```
  (leave `r2` blank for single-end data)

- **`pairs.tsv`** - one row per species pair to analyze:
  ```
  pair_id  sampleA  sampleB  reference_fasta  reference_gff  repeat_bed  include_sequences
  ```
  (`include_sequences` blank = no scaffold filtering for that reference)

See `samples_example.tsv` / `pairs_example.tsv` for a worked example with 6
samples and 3 fully independent pairs - each pair has its own reference
genome, GFF and repeat BED, and no sample appears in more than one pair.
This is the typical usage pattern.

**If a reference genome or a sample ever is reused across pairs** (not the
usual case, but the pipeline handles it correctly either way), that work
is automatically deduplicated rather than redone per pair - gene-buffer/
repeat-removal/genome-indexing runs once per distinct reference file, and
mapping runs once per distinct (sample, reference) combination. Verified
via `snakemake -n` job counts on a deliberately overlapping test case
(2 references, 5 samples, 3 pairs with intentional reuse): reference
indexing ran 2 times rather than 3, mapping ran 5 times rather than 6.
On the actual example above (nothing shared), every one of those rules
correctly runs once per sample/reference/pair as expected, with no
deduplication effect at all - the mechanism only ever kicks in when
something genuinely is shared.

`config.yaml` now holds only parameters that apply uniformly across every
pair in a batch (freebayes settings, depth filters, block length/spacing/
buffer, thread counts) - nothing pair- or reference-specific lives there
anymore.

## Installation

Install in this exact order - `gimble` itself is **not** managed by
`--use-conda` the way everything else is, and needs to already be on
`PATH` before you run anything.

**1. Create one environment and install Snakemake + `gimble` into it:**

```bash
conda create -n demography_pipeline -c conda-forge -c bioconda snakemake mamba gimble
conda activate demography_pipeline
```

`gimble`'s own conda-forge package has repeatedly failed to solve when
Snakemake tries to create it as a separate, isolated per-rule environment
via `--use-conda` (see `envs/gimble.yaml`, kept only as a record/fallback -
it is not actually used by any rule). Installing it directly into the
same environment you run `snakemake` from sidesteps that entirely: the
`gimble_parse`/`blocks`/`tally`/`info`/`query` rules deliberately have no
`conda:` directive, so they just use whatever `gimble` is already active
on your `PATH`, rather than asking Snakemake to manage it.

**2. Let Snakemake build every *other* tool's environment from the `envs/*.yaml` files:**

```bash
snakemake --use-conda --conda-frontend mamba --conda-create-envs-only --cores 1
```

This creates the `mapping`, `variant_calling`, `bedtools` and
`gimbleprep` environments automatically - `gimbleprep` is a separate tool
from `gimble` itself and installs cleanly this way; it does not need the
same manual treatment. Run this once, on a login node with internet
access if you're on a cluster where compute nodes don't have it - see
`run_pipeline.slurm`, which builds this in as its own first step before
ever submitting anything to a compute node.

**3. Confirm both parts worked:**

```bash
gimble --version
snakemake -n --cores 1   # dry run - should resolve without error
```

## Requirements

- [Snakemake](https://snakemake.readthedocs.io) >= 8, with `--use-conda` support
- conda/mamba
- `gimble`, installed as in step 1 above - not via `--use-conda`

## Quick start

1. Copy `samples_example.tsv` -> `samples.tsv` and `pairs_example.tsv` ->
   `pairs.tsv`, and fill them in with your actual samples/pairs (see above).
2. Check `config.yaml`'s depth/block/spacing parameters against your
   species (see below for what to actually change vs. leave alone) - these
   apply to every pair in the batch.
3. Run:
   ```
   snakemake --cores 16 --use-conda
   ```
4. Final output: one `results/bsfs/{pair_id}.bsfs.csv` per row in
   `pairs.tsv`, ready for `demogfit::bsfs_to_s_distribution()` /
   `fit_demography()`.

For a cluster, add a [Snakemake executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/)
for your scheduler (e.g. slurm) rather than hardcoding cluster-specific
flags into the Snakefile itself - see "What changed" below for why.

## What this pipeline does NOT do

- **Genome annotation.** You supply the GFF; this pipeline does not run
  BRAKER or any other annotation tool.
- **Model fitting.** That's `demogfit`'s job - this pipeline stops at the
  blockwise SFS.

## Where blocks are drawn from: `region_strategy`

Set in `config.yaml`, one of three values:

- **`"intergenic"`** - blocks only from gene-buffered, repeat-free
  intergenic sequence (see the three criteria below). The cleanest
  signal, but for gene-dense genomes the ceiling on how much sequence is
  even available this way can be small - measured directly on a real
  *D. melanogaster* genome during this pipeline's development: raw,
  unbuffered gene footprint alone covered ~76% of the accessible genome,
  independently cross-checked straight from the raw GTF outside this
  pipeline entirely (102.8Mb vs. 101.6Mb - within 1.2%, not a pipeline
  bug). The absolute best case with zero buffer at all was ~24% of the
  genome; realistic buffer settings leave meaningfully less.
- **`"intron"`** - blocks only from ONE intron per gene: the **largest
  intron that is NOT the first intron**, in transcription order (strand-
  aware - for a minus-strand gene this is coordinate-*decreasing* order,
  not left-to-right on the chromosome), trimmed by `intron.trim` bp on
  each end (default 10). Genes with fewer than 2 introns on their
  canonical transcript (the longest, by genomic span, for multi-isoform
  genes) contribute nothing and are skipped, not an error. This
  deliberately draws sequence from *inside* gene bodies - that is the
  point of the strategy, not a bug, and `validate_blocks`/
  `validate_reference_introns` do not flag it as one. But it is a real
  tradeoff, not a free upgrade: an intron sits inside its gene's own
  transcript, arguably *more* exposed to that gene's linked/background
  selection than intergenic sequence safely buffered away from it. Both
  this pipeline's own predecessor (`bed_pruner.v1.py`) and DISMaL
  default to introns, so this is a precedented, defensible choice in
  this research lineage - not a fringe one - but it should be a decision
  you make deliberately, not accept because it yields more blocks.
- **`"mixed"`** (default) - the union of both candidate pools above,
  with block spacing enforced **globally** across the combined set: an
  intergenic block and an intron block are never placed closer together
  than `intergenic.block_spacing` allows either, since spacing is
  computed once over the whole unified pool rather than independently
  per source and concatenated afterward (which would allow two blocks
  from different sources to end up linked to each other).

### What's actually checked, per strategy

`validate_reference_intergenic`/`validate_reference_introns` (reference-
level, before any sample data is touched) and `validate_blocks` (per-pair,
final blocks) apply different checks depending on `region_strategy`,
since "overlaps a gene" is a failure for intergenic blocks and an
expected, required property for intron blocks:

| Check | intergenic | intron | mixed |
|---|---|---|---|
| No overlap with gene+buffer zone | checked | not applicable | not checked (some blocks legitimately overlap) |
| No overlap with exons | n/a (already excluded via gene buffer) | checked | checked |
| No overlap with repeats | checked | checked | checked |
| Block spacing | checked | checked | checked, globally across both sources |

**Exon overlap is construction, not just validation** - `repeat_free_introns_bed`
subtracts exonic sequence (via `extract_introns.py`'s `exons_out`, which
covers every transcript of every gene genome-wide, not just canonical
ones) the same way it subtracts repeats, before the validation check
above ever runs. This isn't a hypothetical safeguard: the first real run
of this feature immediately hit exactly this case on real data - 865 of
11,752 intron candidates overlapped an exon. Diagnosed directly (not
assumed) by comparing the intron's own gene ID against the overlapping
exon's gene ID for every violation: 100% were a *different* gene's exon,
zero were the same gene's own exon. That rules out a bug in the intron/
gap-computation logic itself (which cannot, by construction, produce an
interval overlapping its own gene's exons) and confirms real, reasonably
common **nested or overlapping genes** - one gene's exon sitting inside
a different gene's intron, which Drosophila (and many other genomes) has
plenty of. Worth understanding why this is more likely here than it
might seem at first: intron selection deliberately picks the *largest*
non-first intron per gene, and large introns are exactly where a nested
gene has room to fit - so encountering this wasn't a fluke of the test
data, it's a direct, expected consequence of the selection criterion.
Excluding this sequence is still correct even though it isn't a bug: a
block sitting inside a *different* gene's exon is under that gene's own
selective constraint, exactly the contamination this whole framework
exists to avoid, regardless of which gene's intron it happened to be
drawn from.

### A note on how `region_candidate_bed` is built

For `"mixed"`, the intergenic and intron candidate pools are simply
concatenated (after both independently pass repeat removal and
`include_sequences` filtering) and re-sorted - not tagged by source. This
was a deliberate choice after testing the alternative: `bedtools
intersect` only preserves the `-a` side's extra columns, not `-b`'s
(confirmed directly), so a per-region source tag placed on these
candidates would silently vanish once intersected against a pair's
callable sites later in the pipeline. Tagging also isn't needed - the
strategy-aware validation above doesn't require knowing which block came
from which source, only what the final combined set as a whole should
and shouldn't overlap. One concrete bug this surfaced during testing:
`extract_introns.py`'s output carries a 4th (`gene_id`, traceability-only)
column that the intergenic candidates don't, so naively concatenating
them produces a BED with inconsistent column counts across rows, which
`bedtools` correctly refuses to parse - `region_candidate_bed` strips
down to 3 columns before combining, specifically because this was hit
and fixed during real end-to-end testing, not anticipated in the abstract.

## The three intergenic-block filtering criteria

(these apply to the `"intergenic"` half of `region_strategy` specifically
- see above for what governs intron-derived blocks instead: the largest
non-first intron per gene, trimmed, with no gene-buffer or spacing
distinction from intergenic blocks once combined under `"mixed"`)

1. **>= 100kb from any annotated gene** (`intergenic.gene_buffer`) - reduces
   the chance that background/linked selection is acting on a block.
2. **Repeats entirely removed** using your supplied repeat BED - repeat
   regions are frequently mismapped and show inflated, spurious
   heterozygosity that would otherwise distort the composite likelihood.
3. **>= 10kb between consecutive blocks** (`intergenic.block_spacing`) -
   blocks are treated as independent data points by the composite-likelihood
   models in `demogfit`, so blocks that sit close together (and are
   therefore linked) must not both be selected. Enforced by
   `scripts/select_spaced_blocks.py` - a new script (see below), tested on
   synthetic data to confirm the spacing constraint holds exactly.

   **Important format note**: `gimble parse -b` does not accept a plain
   3-column BED. `gimbleprep`'s own callable BED (and everything derived
   from it) is 5 columns - `chrom, start, end, num_samples, samples` (a
   comma-separated per-interval callable-sample list) - and `gimble parse`
   reads columns 0, 1, 2 and 4 specifically, erroring
   (`Too many columns specified`) if column 4 doesn't exist.
   `select_spaced_blocks.py` preserves any columns beyond chrom/start/end
   unchanged on every block it carves out - confirmed directly by
   replicating gIMble's own pandas read call against its output, not just
   assumed from reading the source.

Optionally, a fourth: **unplaced scaffolds/contigs can be excluded from
block selection** via `reference.include_sequences` (a text file listing
one sequence ID per line to keep). These are often assembly artifacts,
collapsed repeats, or contamination rather than real chromosomal sequence.
Left `null` by default, since many draft assemblies (especially for less-
studied species) have no chromosome-level scaffolds at all - forcing a
filter on would discard the entire usable genome for those. Reads still
map against the *full* assembly regardless of this setting - only block
selection is restricted - since excluding unplaced scaffolds from mapping
would make reads that truly belong there mismap onto paralogous regions
of real chromosomes instead, contaminating the signal there.

Every BED file involved in block selection is explicitly filtered against
`include_sequences`, including `gimbleprep`'s own `prep.bed` (via
`filter_prep_bed_to_allowed_sequences`) - not just relying on
`bedtools complement`'s implicit `-g` scoping to exclude scaffolds
downstream. That implicit path was already correct (verified directly:
`bedtools intersect` cleanly excludes a chromosome present only in one of
its two inputs), but invisible - you couldn't inspect any single file and
confirm scaffolds were actually excluded, and `prep.bed` itself stayed
bloated with scaffold-derived sites that never contributed to the final
blocks. This also directly reduces file sizes throughout, since scaffold
content is dropped as early as each file can reasonably drop it rather
than carried through several intermediate files before the final intersect.

### Independent validation of the reference-level construction (before any sample data is touched)

`scripts/validate_reference.py` (run automatically as the
`validate_reference_intergenic` rule) checks the *reference-only*
intergenic construction - gene buffering and repeat removal, which depend
only on the reference genome/GFF/repeat BED, never on any sample's reads -
immediately after those BED files are built, and **before any read QC,
mapping, or variant calling begins**. It gates `fastp_pe`/`fastp_se` and
`bwa_map_pe`/`bwa_map_se` directly (confirmed via `snakemake -n -p`, not
just assumed - both rules' actual input lists include the validation
report as a real dependency edge). A failure here stops the pipeline
before a single fastq is trimmed, let alone mapped, rather than
discovering the problem after hours of mapping and variant calling.

This matters because the region-selection criteria exist specifically to
reduce the risk of erroneously inferring gene flow - from linked blocks,
from selection acting on badly-placed blocks, or from repeat-driven
inflated heterozygosity - so a bug here isn't a minor inconvenience, it's
a direct threat to the validity of the whole analysis. Checks:

1. `genes.bed` and `genes.buffered.merged.bed` are non-empty (also
   independently enforced inside `gene_bed` itself - checked again here
   in case a later step silently produced nothing from valid input).
2. `intergenic.repeat_free.bed` has zero overlap with the gene+buffer zone.
3. `intergenic.repeat_free.bed` has zero overlap with repeats.
4. `intergenic.repeat_free.bed`'s total bp is strictly less than the
   genome's - if equal, gene buffering excluded nothing at all, which is
   exactly the failure signature of an empty `genes.bed` reaching this far.

Tested against the actual "whole genome misclassified as intergenic" bug
this pipeline hit for real: it now fails on three checks simultaneously
(size, gene overlap, and repeat overlap), not just one - genuine defense
in depth. `validate_blocks` (below) still runs later too, since it catches
a different class of bug specific to the per-pair depth-filtering step,
which can't be checked before sample data exists.

### Independent validation of the final per-pair blocks

`scripts/validate_blocks.py` (run automatically as the `validate_blocks`
rule, gating `gimble_parse` - a failure here stops the pipeline before any
bad blocks reach gIMble) checks the *actual final block set* against all
three criteria independently, rather than re-running the same construction
logic that built them:

1. Zero overlap with the gene+buffer zone
2. Zero overlap with repeats
3. >= the configured spacing between any two blocks on the same chromosome

This matters because a bug upstream (an empty gene BED, a stale or
mismatched repeat file, a swapped manifest column) can silently produce a
technically-valid but scientifically wrong block set with no error
anywhere else in the pipeline - re-deriving the same logic wouldn't catch
a bug in that logic, but checking the real output against the stated
properties will. Tested against a gene-overlap violation, a repeat-overlap
violation, a spacing violation, and a valid case, each producing the
expected pass/fail with a specific diagnostic. The report lands at
`results/gimble/{pair_id}/block_validation.txt`.

### Validation gates don't cause unnecessary re-mapping/re-calling

`reference_validation.txt` and `block_validation.txt` are real hard gates
(see above) - `fastp_pe`/`bwa_map_pe`/`gimble_parse` genuinely depend on
them - but they aren't *data* those rules read, just pass/fail checkpoints.
Snakemake's default staleness check is timestamp-based, so without
special handling, any time a validation file gets re-generated (even
re-confirming the exact same PASS result) its timestamp becomes newer
than already-finished downstream outputs like `final_bams/*.bam` or a
completed VCF, and Snakemake would conclude the entire mapping/calling
chain needs to be redone - expensive, and pointless, since nothing about
the actual reference or reads changed. All five of these dependencies are
wrapped in `ancient()`, which tells Snakemake "this must exist, but ignore
its timestamp." Verified directly, not just asserted: simulated a
completed mapping chain, touched `reference_validation.txt` to be
genuinely newer than the finished BAM, and confirmed `snakemake -n` still
reports "nothing to be done" - while a real change (touching the raw
fastq) still correctly triggers the full re-mapping chain, confirming
`ancient()` only suppresses the false-positive case.

## Scaling to non-Drosophila genomes

This was designed and validated on Drosophila (~150-200 Mb genomes), but is
meant to work for anything from an insect to a mammal or fish. Things that
matter as genome size grows:

- **Freebayes is not internally parallel.** For anything past a small
  genome, a single freebayes invocation will not finish in a practical
  time. This pipeline splits the reference into `freebayes.n_chunks`
  region chunks and calls each in parallel, then concatenates the results -
  increase `n_chunks` for larger genomes.
- **bwa-mem2 indexing memory scales with genome size** (rule of thumb:
  budget roughly 28x the FASTA size in RAM). Not something this pipeline
  can fix - just something to plan compute for.
- **Depth filters are kept low** (`gimble.min_depth`, default 2) because
  this pipeline is designed around a single wild-caught individual per
  species, not a population sample - a strict depth filter would discard
  most of the genome for anything short of very high coverage data.
  `gimble.max_depth` is left `null` by default, so `gimbleprep`'s own
  default (2x each BAM's mean coverage) applies unmodified - a sensible
  default on its own, not something this pipeline needs to override.
- **`gimble.kmax` must be set above the largest per-category mutation
  count actually observed in your data** (see the warning in `config.yaml`
  and in the Snakefile's `gimble_tally` rule) - this is a real gIMble
  behavior, not something specific to this pipeline, and it will silently
  truncate/corrupt your bSFS if left at gIMble's own default of 2.

## What changed relative to the original per-script Drosophila pipeline

This is a genuine redesign, not just a wrapper around the old scripts - the
following issues were found and fixed:

1. **Not actually one pipeline.** Variant calling, gIMble processing, and
   region filtering previously lived in separate bash scripts invoked by
   hand after the Snakemake portion finished. Everything is now inside one
   Snakefile, runnable end to end with a single `snakemake` command.
2. **Hardcoded absolute paths** to a personal gIMble install
   (`/home/lyusuf/scratch/.../gIMble`) and SLURM partition names specific
   to one HPC cluster. Removed - tools are resolved via conda envs, and
   the Snakefile itself has no cluster-specific assumptions.
3. **Manual `sed` placeholder substitution before every run** (reference
   genome name, min/max feature length values baked into filenames as
   literal `PLACEHOLDER` strings requiring hand-editing). Replaced with
   `config.yaml` parameters - nothing needs text-editing between runs.
4. **An internet-dependent Snakemake wrapper** for `AddOrReplaceReadGroups`
   (fetched from the Snakemake wrapper repo at runtime - fails on offline
   HPC nodes). Replaced with a direct `picard` call in the conda env.
5. **A bwa/bwa-mem2 inconsistency with an Oxford Nanopore preset**
   (`bwa mem -x ont2d`) applied to what was actually short-read Illumina
   single-end data, alongside `bwa-mem2` for paired-end data in the same
   file. Both paths now consistently use `bwa-mem2` with standard settings.
6. **The gIMble CLI has changed since the original scripts were written**
   (`preprocess` → `gimbleprep`, `setup` → `parse`, plus a now-required
   `tally` step with its own `kmax` truncation parameter). Verified
   directly against the current `LohseLab/gimbleprep` and `LohseLab/gIMble`
   source rather than assumed, since guessing at CLI flags for a tool that
   can't be run here would be a bad way to find out something broke.
7. **Duplicated whole-genome vs. intronic processing tracks** (two parallel
   `gIMble setup`/`blocks`/`info` invocations). Collapsed to a single
   intergenic-only track, since intronic filtering was dropped.
8. **The per-gene intron sampler (`bed_pruner.v1.py`) doesn't generalize**
   to intergenic-only filtering - it picked one region per gene, which
   has no equivalent structure once you're not keying off genes. Replaced
   with `scripts/select_spaced_blocks.py`, a genome-wide greedy spacing
   selector purpose-built for the "≥10kb between any two blocks" rule.
9. **gimbleprep already has a built-in max-depth filter**
   (`-M`/`--max_depth`, expressed as a multiple of each BAM's own mean
   coverage, default 2x) - a custom depth-outlier filter was planned before
   this was discovered in gimbleprep's source; using the tool's own
   parameter instead avoids reimplementing something it already does.
10. **A duplicate `samtools faidx` call** across two of the original
    scripts. Now a single Snakemake rule, run once.
