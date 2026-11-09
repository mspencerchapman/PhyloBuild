# Installing and running PhyloBuild

This guide covers installing PhyloBuild, setting up a study (the study directory, sample sheet and config file),
and running it, mainly on the Sanger farm. For what the pipeline does and what it produces, see the [README](README.md).

## 1. Installing

PhyloBuild has three parts:

- **the workflow** (`workflow/`, `profiles/`, `config/`): what you run. You get it by cloning the repository;
- **Snakemake**, which runs the workflow. This is a Python program, installed once into its own environment;
- **the `phylobuild` R package** (`phylobuild/`), which does the filtering and tree building. Every R step loads it.

### On the Sanger farm

Everything else (cgpVAF, the hairpin filter container, MPBoot, VAGrENT, R 4.4) is already on the farm and gets
loaded as environment modules by the pipeline.

```bash
# 1. Clone the repository (results don't go in it, so anywhere you like; a lustre folder is fine)
cd /lustre/scratch125/casm/teams/team273/users/$USER
git clone https://github.com/mspencerchapman/PhyloBuild.git
cd PhyloBuild

# 2. A Python environment for Snakemake (once; it can be shared by all your projects)
module load python/3.12.3      # or any Python >= 3.11 (`module avail python` lists them)
python3 -m venv ~/phylobuild-venv
~/phylobuild-venv/bin/pip install "snakemake>=9" snakemake-executor-plugin-cluster-generic

# 3. The R package, into your own R library
module load R/4.4.0
mkdir -p ~/R/phylobuild-lib
export R_LIBS=~/R/phylobuild-lib       # add this line to ~/.bashrc so every job sees the package
Rscript -e 'remotes::install_github("nangalialab/treemut", lib = "~/R/phylobuild-lib")'   # the R module adds its own paths to R_LIBS, so name the folder
R CMD INSTALL --library=$R_LIBS phylobuild
```

Instead of cloning first, you can install just the R package straight from GitHub (useful if you only want the
functions in your own R sessions):

```r
remotes::install_github("mspencerchapman/PhyloBuild", subdir = "phylobuild")
```

To check the installation:

```bash
Rscript -e 'library(phylobuild); packageVersion("treemut")'   # treemut should be 1.1 or later
~/phylobuild-venv/bin/snakemake --version
```

### Elsewhere

Use the conda environments in `workflow/envs/`: `snakemake.yaml` for Snakemake, and `phylobuild-r.yaml` for R and
MPBoot (then install treemut and the R package as in step 3). You also need cgpVAF (`cgpVaf.pl`), Singularity with
the LCM filtering container, and optionally VAGrENT. Point the `tools:` section of the config at them, and empty the
`modules:` lists. Set `bam_staging: none` and give BAM paths in the sample sheet (section 2). You'll also need a
profile for your cluster: `profiles/sanger-lsf/` is the one to copy (it's written for LSF).

### Updating

```bash
cd PhyloBuild
git pull
R CMD INSTALL --library=$R_LIBS phylobuild    # only needed if anything in phylobuild/ changed
```

Snakemake notices when a step's code or settings have changed. A dry run (section 4) shows what would rerun.

## 2. The study directory and sample sheet

Each study (a set of donors analysed together) has its own directory, which you make. It holds the study's config
and sample sheet, and the pipeline puts each donor's results in a subdirectory of it:

```
my_study/
  config.yaml          your copy of the config template (section 3)
  samples.tsv          the sample sheet
  PD12345/             results for each donor: calls, hairpin, cgpvaf, qc, mats_and_params, trees, logs
  PD23456/
  staged_bams/         BAMs staged from iRODS (deleted as they're finished with if keep_staged_bams: false)
  logs/lsf/            LSF output of every job
```

The sample sheet is a tab-separated file with one row per sample (`config/samples_example.tsv` shows the format):

| Column | What it is |
|---|---|
| `donor` | the donor key used in the config (usually the PD ID; for a transplant pair, an experiment ID) |
| `sample` | sample name, e.g. `PD12345a_lo0001` |
| `project` | CanApps project number (needed when `bam_staging: sanger`) |
| `bam` | path to the BAM (needed when `bam_staging: none`; `.bai` and `.bas` must sit next to it) |
| `caveman_vcf`, `pindel_vcf` | the `caveman_c.annot.vcf.gz` and `pindel.annot.vcf.gz` files |
| `irods_db` | optional: `live` (default) or `test`, the iRODS database the BAM is staged from |

There are two ways to get one.

### At Sanger: let the pipeline build it

Give each donor in the config the CanApps projects its samples are in (section 3), and leave out `samples.tsv`:

```yaml
donors:
  PD12345:
    projects: [1234]
  PD23456:
    projects: [1234, 1235]
    exclude_samples: [PD23456a_lo0010]
```

The first time you run the pipeline (a dry run is enough), it looks in each project's nst_links folder and writes
`samples.tsv` into the study directory, with every sample whose name starts with the donor ID. It prints how many
samples it found for each donor, and lists any whose VCFs aren't in nst_links yet (these are left out). Check the
sheet before the real run. These donor settings choose the samples:

| Setting | What it does |
|---|---|
| `projects` | CanApps project numbers to look in, e.g. `[1234, 1235]` |
| `pd_ids` | take samples whose names start with these (default: the donor key), e.g. `[PD34567, PD45678]` for a transplant pair under an experiment ID |
| `include_samples` | extra samples whose names don't start with `pd_ids`, e.g. samples that were labelled with another individual's PD ID |
| `exclude_samples` | samples to leave out, by name |
| `exclude_pattern` | leave out samples whose names match this regular expression, e.g. `"_lo[0-9]+$"` |
| `nst_links` | default `/software/SGP/nst_links/live`; `/software/SGP/nst_links/test` for a test project (section 5) |

`exclude_samples` and `exclude_pattern` are applied every time the sheet is read, so you can change them at any time.
The sheet's first line records the `projects`, `pd_ids` and `include_samples` it was built from. If you change
those, the pipeline stops and asks you to delete `samples.tsv`, and then builds it again on the next run. Delete it
too to pick up samples that have arrived in nst_links since it was built.

`bam_staging` (section 3) decides how the BAMs are reached: `none` reads them straight from nst_links (fine for
live projects, and it uses no disk space), `sanger` stages them out of iRODS first (needed for test projects). The
sheet has both the BAM path and the project, so it works either way.

### Or make it yourself

For anything the settings above can't express, or outside Sanger, make `samples.tsv` yourself: by hand, in R, or
with your own script. **If a `samples.tsv` exists without the pipeline's first line, it is used exactly as it is,
and the pipeline doesn't build one.** `projects`, `pd_ids` and `include_samples` are then ignored, but
`exclude_samples` and `exclude_pattern` still apply. If you edit a sheet the pipeline built, delete its first line
(the one starting `# Built by PhyloBuild`) so that your edits are kept.

At Sanger, `workflow/scripts/sanger_sample_sheet.py` (which is what the pipeline uses) is a good starting point. Run
it from the PhyloBuild folder:

```bash
S=/lustre/scratch125/casm/teams/team273/users/$USER/my_study
python3 workflow/scripts/sanger_sample_sheet.py --donor PD12345 --projects 1234 > $S/samples.tsv
python3 workflow/scripts/sanger_sample_sheet.py --donor PD23456 --projects 1234,1235 --no-header >> $S/samples.tsv
```

It takes `--pd-ids`, `--include`, `--exclude`, `--exclude-pattern` and `--nst-links`, which work like the settings of
the same names above.

## 3. The config file

All settings for a study are in one YAML file. `config/config.yaml` in the repository is a template: copy it into
the study directory and edit the copy,

```bash
cp config/config.yaml $S/config.yaml
```

and give it to every run with `--configfile $S/config.yaml` (section 4). YAML is plain text: `key: value` pairs,
indentation (spaces, not tabs) to nest settings, `#` for comments, and lists written as `[a, b]`. The file has five
parts.

### Study-wide settings (top of the file)

```yaml
project: my_study
study_dir: /lustre/scratch125/casm/teams/team273/users/me/my_study   # must exist; results go to <study_dir>/<donor>/
# samples: ...                       # default <study_dir>/samples.tsv (built from `projects` if missing; section 2)
bam_staging: sanger                  # sanger: stage BAMs from iRODS; none: use the bam column as given
# staging_dir: ...                   # default <study_dir>/staged_bams
keep_staged_bams: true               # false: delete each staged BAM once cgpVAF and hairpin have used it
```

`study_dir` is the one setting you must change: the run stops if it isn't an existing directory.

`keep_staged_bams: false` saves a lot of space (about 25 GB per colony): each BAM is deleted once the last job that
reads it has finished. If you rerun from scratch later, the BAMs are staged again. Even so, all of a donor's BAMs
are on disk at once while cgpVAF runs, so a donor with hundreds of colonies can need more space than the team has.
For live projects, `bam_staging: none` avoids this by reading the BAMs straight from nst_links.

### `reference`, `tools` and `modules`

These are the reference genome files, the in silico normal and the software. You rarely need to change them within
the same genome build. The in silico normal is staged from `normal_project` (with `normal_irods_db`) like the
samples, or read from `normal_bam` if `bam_staging: none`.

### `defaults`: settings for every donor

These are the settings that used to be command-line flags in the bash scripts. The ones you're most likely to
change:

| Setting | What it does |
|---|---|
| `sex` | `male`, `female` or `unknown`. With `unknown`, the sex inferred from the data is used. Otherwise it's checked against the data, and the run stops if they disagree (`sex_check: warn` carries on instead). |
| `projects`, `pd_ids`, `include_samples`, `nst_links` | where to find the donor's samples when the pipeline builds the sample sheet (section 2); usually set per donor |
| `exclude_samples`, `exclude_pattern` | samples to leave out of every step, by name or by regular expression |
| `asmd_threshold` | CaVEMan ASMD cut-off for PASS SNVs (the bash workflow used 120) |
| `upfront_colony_filter` | `[none, 0.4]`: build trees from all colonies, and also after removing colonies with peak VAF below 0.4 (the old `_m40` runs). `[none]` for just the first. |
| `tree_filtering` | `[pval, vaf]`: one tree per filtering strategy |
| `cgpvaf_batch_size` | samples per cgpVAF job |
| `age` | age at sampling, for the ultrametric tree (`null` scales to the mean burden) |
| `tree:` | tree-building options (`filtering_from_table` flags). The comments in the file explain each one. |

Each combination of `upfront_colony_filter` × `tree_filtering` gives one set of tree outputs, so the defaults give
four trees per donor: `none`/`pval`, `none`/`vaf`, `cf0.4`/`pval` and `cf0.4`/`vaf`.

### `donors`: one entry per donor

Anything set under a donor replaces the default for that donor only. Nested settings (under `tree:`) are merged, so
you only write the ones that differ:

```yaml
donors:
  PD12345:
    age: 20
    projects: [1234]
  PD23456:
    age: 20
    projects: [1234]
    exclude_samples: [PD23456a_lo0010]
    tree:
      filter_y: true            # only this tree option changes; the others keep their defaults
  TX_001:                     # a transplant pair under an experiment ID
    projects: [1234, 1235]
    pd_ids: [PD34567, PD45678]
    age: 42
    upfront_colony_filter: [none]
    tree:
      group_pattern: PD34567    # donor samples (the default is the donor key)
      recipient_pattern: PD45678
```

Every donor in the config must have samples in the sample sheet. To run only some donors, comment the others out
(put `#` at the start of their lines) or name the outputs you want (section 4).

`group_pattern` defaults to the donor key, and it's only used to colour the tree plots. If the key isn't the start
of the sample names (an experiment ID, or a key like `PD23456_remapped`), set it to the PD ID.

### Several studies

Each study has its own directory and config, so one copy of PhyloBuild serves them all. The config given with
`--configfile` replaces the template; it isn't merged with it, so each study's config is complete on its own.

You can also change a single setting for one run from the command line, e.g.
`--config keep_staged_bams=False` (after any output names, like `--configfile`).

## 4. Running

Always run from the top of the PhyloBuild folder, because the profile refers to its scripts by relative path. The
results still go to the study directory.

```bash
cd PhyloBuild
export R_LIBS=~/R/phylobuild-lib                    # if it isn't in ~/.bashrc
S=/lustre/scratch125/casm/teams/team273/users/$USER/my_study

# What would run? (a dry run: lists the jobs, runs nothing)
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf --configfile $S/config.yaml -n

# Run everything for every donor in the config
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf --configfile $S/config.yaml

# Or just some outputs, e.g. one donor's sex check, or one tree. Output names go before --configfile (anything
# after it is read as another config file)
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf $S/PD12345/qc/sex_check.tsv --configfile $S/config.yaml
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf \
    $S/PD12345/trees/PD12345_none/tree_files/tree_PD12345_none_pval_post_mix_post_dup.tree --configfile $S/config.yaml
```

The Snakemake process itself (the "controller") does little work; it submits every step as its own LSF job. But it
must keep running until the end, which takes hours (about 3.5 hours for 48 colonies from scratch). So submit it as
a job on the basement queue:

```bash
cd $S
cat > run_phylobuild.sh <<'EOF'
#!/bin/bash
set -euo pipefail
S=/lustre/scratch125/casm/teams/team273/users/$USER/my_study
cd /lustre/scratch125/casm/teams/team273/users/$USER/PhyloBuild
export R_LIBS=~/R/phylobuild-lib
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf --configfile $S/config.yaml
EOF
chmod +x run_phylobuild.sh
bsub -q basement -G team273 -M 4000 -R "select[mem>=4000] rusage[mem=4000]" \
    -o run_phylobuild.log -e run_phylobuild.log -J phylobuild ./run_phylobuild.sh
```

Then follow it with `tail -f $S/run_phylobuild.log`, where it prints `N of M steps (x%) done`.

**If something fails**:
- **Finding the error:** Snakemake retries a failed step once, and keeps going with everything that doesn't depend
  on it. At the end it lists the failed steps. Each step's log is under `<study_dir>/<donor>/logs/`, and the LSF
  output is under `<study_dir>/logs/lsf/`.
- **Carrying on:** fix the problem and run the same command again. Only what's missing or out of date runs.
- **Lock errors:** if the controller was killed, the next run may say the directory is locked. Run the command once
  with `--unlock`, then again without it. (Snakemake keeps its own small bookkeeping folder, `.snakemake/`, in the
  PhyloBuild folder; that's the only thing a run writes there.)
- **Rerunning a step:** to rerun a step that worked, add `--forcerun <rule>` (e.g. `--forcerun build_tree`) or
  delete its outputs.

Runs for different studies can go at the same time from the same PhyloBuild folder. Two runs of the same study can
only if they write different files (e.g. different donors).

## 5. Example: a second analysis of a donor from a test project

Say PD23456 has been remapped to a modified reference and the new results are in a CanApps test project, 9999
(`/software/SGP/nst_links/test/9999/`). Two things differ from a normal run:

1. **The BAMs are in the iRODS test database.** `stageBam.pl` only finds them with `-l test`. The BAMs in nst_links
   are links into iRODS that can't be read directly, so they must be staged.
2. **It's a separate analysis of the same donor,** so it needs its own outputs: treat it as its own study.

- **Study directory:** `R=/lustre/scratch125/casm/teams/team273/users/$USER/my_study_remapped; mkdir -p $R`
- **Config:** copy `config/config.yaml` (or the standard study's config) to `$R/config.yaml`, set `study_dir` to
  `$R` and `bam_staging: sanger`, and give only PD23456 under `donors`, with the test project:
  ```yaml
  donors:
    PD23456:
      projects: [9999]
      nst_links: /software/SGP/nst_links/test
  ```
  The pipeline builds `$R/samples.tsv` from that folder (section 2). Because it's a test folder, `irods_db` is filled
  in as `test`, so the BAMs are staged from the test database.
- **Reference genome:** set `reference.genome` to the FASTA the BAMs were mapped to. If the modified reference has
  the same contigs and coordinates (e.g. a region replaced with Ns), the other reference files can stay as they are.
- **In silico normal:** if there's no normal mapped to the modified reference, keep the standard one. In the
  pipeline the normal is only a column in the cgpVAF matrices; the filtering and trees use the colonies. If the
  normal is in a test project too, set `normal_project` and `normal_irods_db: test`.

Run it with:

```bash
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf --configfile $R/config.yaml -n
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf --configfile $R/config.yaml
```

The results go to `my_study_remapped/PD23456/`, separate from the standard run in `my_study/PD23456/`.

## 6. Using the R package on its own

Every step's output is an `.rds` file that you can load in R:

```r
library(phylobuild)
d = "/path/to/my_study/PD12345"
x = read_mats_and_params(file.path(d, "mats_and_params/mats_and_params_PD12345_none.rds"))
names(x)                 # COMB_mats, filter_params
tree = ape::read.tree(file.path(d, "trees/PD12345_none/tree_files/tree_PD12345_none_pval_post_mix_post_dup.tree"))
muts = readRDS(file.path(d, "trees/PD12345_none/annotated_muts/annotated_mut_set_PD12345_none_pval_post_mix_post_dup.rds"))
```

`read_mats_and_params()` also reads the old `save()` files from the bash workflow. The tree-plotting functions you
know (`plot_tree`, `add_annotation`, `plot_sharing_info`, and so on) are in the package too.
