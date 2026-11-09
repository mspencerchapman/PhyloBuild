# PhyloBuild

Somatic phylogenies from whole-genome sequencing of clonal samples (e.g. single-cell derived colonies), as a
[Snakemake](https://snakemake.readthedocs.io) pipeline. It runs the method described in
Coorens, Spencer Chapman *et al.*, *Nature Protocols* 19, 1866–1886 (2024): read counts at every candidate
site in every sample, statistical filtering of germline variants and artefacts, tree building with MPBoot,
assignment of mutations to branches with treemut, and sensitivity correction of branch lengths.

It replaces the copy-and-paste bash workflow (`bash_analysis_for_*.sh`) and the scripts it called. The R code is
in an R package, `phylobuild/`, which you can also use on its own.

## What it does

For each donor in the config:

| Step | Rule(s) | Replaces |
|---|---|---|
| Stage BAMs out of iRODS (Sanger only, optional) | `stage_bam` | `stageBam.pl` by hand |
| PASS SNVs (CLPM=0, ASMD above threshold) and indels per sample | `pass_snvs`, `pass_indels` | `import_new_samples_only.R`, `filter_pass_subs*.pl`, `filter_indels_new.pl` |
| Hairpin (LCM artefact) filter per sample | `hairpin_filter`, `hairpin_pass`, `hairpin_sites` | `runlcmfilter.sh`, `filter_for_bed_file.pl` |
| Sites for cgpVAF | `site_bed` | `cut ... \| sort \| uniq` |
| cgpVAF read counts, in batches × chromosomes, then merged | `cgpvaf_*` | `createVafCmd.pl`, `create_split_config_ini.R`, `sed` on `run_bsub.sh`, the `cut`/`paste` merge, `rerun_cgpvaf_fails.sh` |
| Sex check | `sex_check` | (new) |
| Filter parameters, with/without upfront colony filtering | `filter_params` | `Mutation_filtering_get_parameters_v2.R` (`-m -v`) |
| Keep hairpin-pass SNVs; drop near-certain germline | `reduce` | `Reducing_mutset_from_MSfilters.R` |
| Mean-depth window | `depth_cutoffs` | the commented-out `hist()` check |
| Per-sample sensitivity | `sensitivity` | `Sensitivity_analysis_from_SNPs.R` |
| Tree building | `build_tree` | `filtering_from_table_mix_remove_GRCh38.R` |

Snakemake works out what needs running, submits jobs to LSF as soon as their inputs are ready, retries failures
(with more memory), and on a rerun only does what is missing or out of date. There are no "once X has finished"
steps.

## Quick start (Sanger farm)

[GUIDE.md](GUIDE.md) has the full version: installing, the sample sheet, the config file, and running.

```bash
# 1. Get the pipeline and a Python environment with Snakemake (once)
git clone https://github.com/mspencerchapman/PhyloBuild.git && cd PhyloBuild
module load python/3.12.3         # Snakemake 9 needs Python >= 3.11; the farm's default python3 is 3.10
python3 -m venv ~/phylobuild-venv
~/phylobuild-venv/bin/pip install "snakemake>=9" snakemake-executor-plugin-cluster-generic
# The R package, into your personal R library (its other dependencies are already on the farm's R 4.4)
module load R/4.4.0
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)'
Rscript -e 'if (!requireNamespace("treemut", quietly = TRUE) || packageVersion("treemut") < "1.2") remotes::install_github("nangalialab/treemut")'
R CMD INSTALL phylobuild

# 2. A study directory, with a copy of the config template in it. Edit the copy: study_dir, and the donors with
#    their CanApps projects (projects: [1234]), ages, sex and settings
S=/lustre/scratch125/casm/teams/team273/users/$USER/my_study
mkdir -p $S
cp config/config.yaml $S/config.yaml

# 3. Check what will run (this also builds $S/samples.tsv from nst_links; check it), then run it
#    (from the PhyloBuild folder, as a basement-queue job; see GUIDE.md)
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf --configfile $S/config.yaml -n
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf --configfile $S/config.yaml
```

Results go to `<study_dir>/<donor>/`, one folder per donor. Rerunning the same command picks up where it stopped.
To run part of the pipeline, name an output before `--configfile`, e.g.
`snakemake --profile profiles/sanger-lsf $S/PD12345/qc/sex_check.tsv --configfile $S/config.yaml`.

## Configuration

Each study has its own config file, copied from the template `config/config.yaml`; the comments there explain
each setting. The main ones:

- **`study_dir`**: the study directory (it must exist). Each donor's results go to `<study_dir>/<donor>/`; the
  sample sheet (`samples.tsv`), staged BAMs (`staged_bams/`) and LSF logs (`logs/lsf/`) default to it too.
- **`donors`**: one entry per donor (or transplant pair); anything set here overrides `defaults`.
- **`projects`**, **`pd_ids`**, **`include_samples`**: where each donor's samples are in nst_links, used to build the
  sample sheet when there isn't one. **`exclude_samples`** / **`exclude_pattern`** leave samples out, by name or by
  regular expression (e.g. `"_lo[0-9]+$"`).
- **`sex`**: `male`, `female`, or `unknown` to use the sex inferred from the data. The pipeline checks the
  configured sex against the data and stops if they disagree (`sex_check: warn` to carry on instead).
- **`upfront_colony_filter`**: e.g. `[none, 0.4]`. Each value gives its own set of trees: `none` builds from all
  colonies; a number first removes colonies whose peak VAF is below it (the old `_m40` runs).
- **`tree_filtering`**: `[pval, vaf]`, one tree per value.
- **`cgpvaf_snv_sites`**: `hairpin` (genotype only SNV sites that passed the hairpin filter; same result, fewer
  sites) or `caveman` (the old order: all CaVEMan PASS sites, reduced afterwards).
- **`bam_staging`**: `sanger` (stage BAMs with `stageBam.pl`; uses the sample sheet's `project` column) or `none`
  (use the sample sheet's `bam` column as given; at Sanger this reads live BAMs straight from nst_links).
- **`tree.depth_cutoffs`**: `auto` sets the mean-depth window from the data, as median ± `depth_mad_k` (default 4)
  × MAD of the sites' mean depth (separately for X/Y in males), so it follows the donor's coverage and number of
  colonies. For a 47-colony test donor at ~15× this gives 6.8–23.5, close to the old fixed 6–23. The histogram
  with the cut-offs is in `<study_dir>/<donor>/qc/depth_cutoffs.pdf`.

The **sample sheet** (`<study_dir>/samples.tsv`; `config/samples_example.tsv` shows the format) is
tab-separated with columns `donor, sample, project, bam, caveman_vcf, pindel_vcf`, and optionally `irods_db`
(`test` for BAMs in a test project, which `stageBam.pl` stages with `-l test`). cgpVAF expects `.bai` and `.bas`
files next to each BAM. At Sanger, the pipeline builds it from the donors' `projects` if it doesn't exist. For more
complex cases, or outside Sanger, make it yourself (by hand or with your own script; with `bam_staging: none`
outside Sanger): a sheet you make is used as it is, and the build step is skipped. See
[GUIDE.md](GUIDE.md#2-the-study-directory-and-sample-sheet).

A config given with `--configfile` replaces the template `config/config.yaml` (it isn't merged with it), so each
study's config is complete on its own.

Farm-specific settings (LSF group, queues) are in `profiles/sanger-lsf/config.yaml`; jobs are submitted by
`profiles/sanger-lsf/lsf_submit.sh`.

## Outputs

Under `<study_dir>/<donor>/`:

| Path | Contents |
|---|---|
| `qc/sex_check.tsv`, `qc/sex_check_per_colony.pdf` | sex measures; X dosage per colony (flags X loss or colonies of the other sex) |
| `qc/depth_cutoffs.pdf` | mean-depth histograms with the window used |
| `qc/sensitivity_<donor>.txt` | per-sample SNV/indel sensitivity |
| `cgpvaf/<donor>_merged_{snv,indel}.tsv` | read count matrices (same layout as the old `merged_SNVs_*.tsv`) |
| `mats_and_params/mats_and_params_<donor>_<cf>[_postMS[_reduced]].rds` | matrices + filter parameters |
| `mats_and_params/colony_vaf_plots_*.pdf`, `colony_peak_vaf_*.tsv` | upfront colony filtering |
| `trees/<donor>_<cf>/` | tree-building outputs, with the same folders and file names as before |
| `trees/<donor>_<cf>/plots/*_poor_fit_mutations.pdf` | heatmap of mutations that fit the tree poorly |
| `logs/` | a log per step |

`<cf>` is `none` or e.g. `cf0.4`. Objects are saved as `.rds`: read them with
`x = readRDS(path)` (e.g. `x$COMB_mats`, `x$filter_params`), or `phylobuild::read_mats_and_params(path)`,
which also reads old `save()` files. Set `legacy_rdata: true` to also write the old format.

## Changes from the bash workflow

Results match the old scripts (filter parameters for a test donor reproduce its existing `mats_and_params`
to floating-point precision; the PASS-call filtering reproduces the old bed files), except where these changes
were deliberate:

- **Sex** comes from the config and is checked against the data. The old rule (Y calls / X calls > 1%)
  labelled 7 of 9 female donors tested as male, which removed nearly all chrX mutations from their trees.
- **Upfront colony filtering** no longer removes every sample when no colony fails (the old `[, -which()]`).
- **cgpVAF merge** finds columns by name, so it works for any batch size.
- **`-y` (filter_y)** uses the correct X/Y names and pseudo-autosomal regions for GRCh38.
- **`-j` (SNVs only) and the copy-number table** are now also applied when the tree is rebuilt after removing
  non-clonal samples (previously only before).
- **MPBoot** runs with a fixed seed, so trees are reproducible.
- **Male X PARs** are treated as diploid (expected VAF 0.5, autosomal depth window); previously all of X/Y was
  treated as haploid, which lost real PAR mutations.
- **Duplicate samples**: the threshold is configurable (`duplicate_min_private`, `duplicate_prop`) and the mean
  burden no longer counts the Ancestral tip (same for the default age scaling).
- **Robustness**: CN-table rows for removed samples are skipped and chromosome names match with or without `chr`;
  missing sensitivities stop the run instead of giving infinite branches; the ultrametric tree handles
  zero-length branches; donor/recipient IDs match from the start of sample names; mutations the user excluded
  are never added back by the germline add-back.
- **New QC**: `plots/<run>_poor_fit_mutations.pdf` (+ `.tsv`) shows every mutation that fits the tree poorly, as a
  VAF heatmap under the tree with the samples its assigned branch implies marked.
- **New options, off by default** (so results stay comparable): `remove_poor_fit`, `addback_require_filters`,
  `clonality_exclude_cn_xy`.
- **Depth window** is set from the data (median ± 4 MAD) instead of fixed 6–23 / 3–10, which only suited ~15× data.
- **Speed**: the binomial tests are vectorised (same values), so the filter-parameter step takes minutes and a
  few GB instead of hours and 64 GB.

## Development

```bash
Rscript -e 'testthat::test_local("phylobuild")'    # R unit tests
snakemake -n                                         # check the workflow builds
```

## Licence

MIT; see [LICENSE](LICENSE).
