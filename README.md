# Species-pair demographic pipeline

`fastq (speciesA, speciesB) + reference genome + GFF + repeat BED` →
mapping → variant calling → intergenic block selection → blockwise SFS,
ready for the [`demogfit`](../demogfit) R package.

## Running many pairs at once

The pipeline is driven by **two manifests**, with a set pipeline configuration across pairs,
which is what lets it scale from 2 pairs to 200 without editing anything
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
genome (belonging to one of the species in the pair), GFF and repeat BED. 
This is the typical usage pattern.


## Installation

Install in this exact order - `gimble` itself is **not** managed by
`--use-conda` the way everything else is, and needs to already be on
`PATH` before you run anything.

**1. Create one environment and install Snakemake + `gimble` into it:**

```bash
conda create -n demography_pipeline -c conda-forge -c bioconda snakemake mamba gimble
conda activate demography_pipeline
```

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
4. Final output: one `pipeline/{pair_id}.tally_blocks_tally.tsv` per row in
   `pairs.tsv`, ready for `demogfit::bsfs_to_s_distribution()` /
   `fit_demography()`.

For a cluster, add a [Snakemake executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/)
for your scheduler (e.g. slurm) rather than hardcoding cluster-specific
flags into the Snakefile itself - see "What changed" below for why.

## What this pipeline does NOT do

- **Genome annotation.** You supply the GFF and repeat annotation; this pipeline does not run
  BRAKER or any other annotation tool.
- **Model fitting.** That's `demogfit`'s job - this pipeline stops at the
  blockwise SFS.

## Where blocks are drawn from: `region_strategy`

Set in `config.yaml`, one of three values:

- **`"intergenic"`** - blocks only from gene-buffered (100kb from nearest gene),
  repeat-free intergenic sequence (see the three criteria below). The cleanest
  signal, but for gene-dense genomes the ceiling on how much sequence is
  even available this way can be small - measured directly on a
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
  genes) contribute nothing and are skipped. Non-first introns that are
  trimmed in the correct way can show less functional constraint than other
  genic regions and some intergenic sequence, and is therefore ideal for
  creating blocks to analyse.
- **`"mixed"`** (default) - the union of both candidate pools above,
  with block spacing enforced **globally** across the combined set: an
  intergenic block and an intron block are never placed closer together
  than `intergenic.block_spacing` allows either, since spacing is
  computed once over the whole unified pool rather than independently
  per source and concatenated afterward (which would allow two blocks
  from different sources to end up linked to each other).

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
collapsed repeats, or contamination. 

## Scaling to non-Drosophila genomes

This was designed and validated on Drosophila (~150-200 Mb genomes), but is
meant to work for anything for any eukaryotic species/system. Things that
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

## Citation
If you have used this package for your work (thanks!), please cite the following:

Yusuf, L.H., Laetsch, D.R., Lohse, K. and Ritchie, M.G., 2026. Genomic analyses in Drosophila do not support the classic allopatric model of speciation. Evolution Letters, 10(2), pp.186-194.

Please also cite all of the packages that are used within this Snakemake pipeline (see config.yaml). 

## A small (AI) note: 
Hi all, for transparency sake: much of this pipeline was developed from standalone bash scripts with the help of Claude (Sonnet/Opus) models. This dramatically sped up the time it took to finish putting together this pipeline together and I have tried to sanity check the results as much as possible, but LLMs are of course error-prone. Please check the outputs carefully to see whether they make sense. If you do spot a bug or an issue with outputs, please do raise an issue and I will try to fix it as quickly as I can.
