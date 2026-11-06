# Installing and running PhyloBuild

This guide covers installing PhyloBuild, setting up a project (the sample sheet and the config file), and running
it, mainly on the Sanger farm. For what the pipeline does and what it produces, see the [README](README.md).

## 1. Installing

PhyloBuild has three parts:

- **the workflow** (`workflow/`, `profiles/`, `config/`): what you run. You get it by cloning the repository;
- **Snakemake**, which runs the workflow. This is a Python program, installed once into its own environment;
- **the `phylobuild` R package** (`phylobuild/`), which does the filtering and tree building. Every R step loads it.

### On the Sanger farm

Everything else (cgpVAF, the hairpin filter container, MPBoot, VAGrENT, R 4.4) is already on the farm and gets
loaded as environment modules by the pipeline.

```bash
# 1. Clone the repository somewhere on lustre (the results go inside it by default, and they're big)
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
Rscript -e 'remotes::install_github("nangalialab/treemut", lib = Sys.getenv("R_LIBS"))'
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

## 2. The sample sheet

The sample sheet is a tab-separated file with one row per sample:

| Column | What it is |
|---|---|
| `donor` | the donor key used in the config (usually the PD ID; for a transplant pair, an experiment ID) |
| `sample` | sample name, e.g. `PD12345a_lo0001` |
| `project` | CanApps project number (needed when `bam_staging: sanger`) |
| `bam` | path to the BAM (needed when `bam_staging: none`; `.bai` and `.bas` must sit next to it) |
| `caveman_vcf`, `pindel_vcf` | the `caveman_c.annot.vcf.gz` and `pindel.annot.vcf.gz` files |
| `irods_db` | optional: `live` (default) or `test`, the iRODS database the BAM is staged from |

At Sanger, build it from nst_links rather than by hand:

```bash
python3 workflow/scripts/sanger_sample_sheet.py --donor PD12345 --projects 1234 > config/samples.tsv
python3 workflow/scripts/sanger_sample_sheet.py --donor PD23456 --projects 1234 --no-header >> config/samples.tsv
```

- `--projects 1234,1235` collects samples spread over several projects.
- `--include` adds samples whose names don't start with the donor ID.
- `--exclude` leaves samples out.
- For a transplant pair, `--donor TX_001 --pd-ids PD34567,PD45678` puts both individuals' samples under one key.

Samples whose VCFs aren't in nst_links yet are listed and left out, so rerun it once they're ready. For a test
project, point it at the test folder (`--nst-links /software/SGP/nst_links/test`); it then fills in `irods_db` as
`test` (section 5 has an example).

## 3. The config file

All settings for a project are in one YAML file, `config/config.yaml` by default. YAML is plain text: `key: value`
pairs, indentation (spaces, not tabs) to nest settings, `#` for comments, and lists written as `[a, b]`. The file
has five parts.

### Project-wide settings (top of the file)

```yaml
project: my_project
outdir: results                      # everything goes to results/<donor>/
samples: config/samples.tsv
bam_staging: sanger                  # sanger: stage BAMs from iRODS; none: use the bam column as given
staging_dir: results/staged_bams
keep_staged_bams: true               # false: delete each staged BAM once cgpVAF and hairpin have used it
```

`keep_staged_bams: false` saves a lot of space (about 25 GB per colony): each BAM is deleted once the last job that
reads it has finished. If you rerun from scratch later, the BAMs are staged again.

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
| `exclude_samples` | samples to leave out of every step |
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
  PD23456:
    age: 20
    exclude_samples: [PD23456a_lo0010]
    tree:
      filter_y: true            # only this tree option changes; the others keep their defaults
  TX_001:                     # a transplant pair under an experiment ID
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

### Several projects

Use one config file per project (each with its own sample sheet, `outdir` and, if runs overlap in time,
`staging_dir`), and choose it with `--configfile`:

```bash
snakemake --profile profiles/sanger-lsf --configfile config/config_PD23456_remapped.yaml
```

A config given this way replaces `config/config.yaml`; it isn't merged with it. So the file must be complete:
copy `config/config.yaml` and edit the copy.

You can also change a single setting for one run from the command line, e.g.
`--config keep_staged_bams=False` (put it at the end of the command, after any output names).

## 4. Running

Always run from the top of the repository, because the profile refers to its scripts by relative path.

```bash
cd PhyloBuild
export R_LIBS=~/R/phylobuild-lib                    # if it isn't in ~/.bashrc

# What would run? (a dry run: lists the jobs, runs nothing)
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf -n

# Run everything for every donor in the config
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf

# Or just some outputs, e.g. one donor's sex check, or one tree
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf results/PD12345/qc/sex_check.tsv
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf \
    results/PD12345/trees/PD12345_none/tree_files/tree_PD12345_none_pval_post_mix_post_dup.tree
```

The Snakemake process itself (the "controller") does little work; it submits every step as its own LSF job. But it
must keep running until the end, which takes hours (about 3.5 hours for 48 colonies from scratch). So submit it as
a job on the basement queue:

```bash
cat > run_project.sh <<'EOF'
#!/bin/bash
set -euo pipefail
cd /lustre/scratch125/casm/teams/team273/users/$USER/PhyloBuild
export R_LIBS=~/R/phylobuild-lib
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf
EOF
chmod +x run_project.sh
bsub -q basement -G team273 -M 4000 -R "select[mem>=4000] rusage[mem=4000]" \
    -o run_project.log -e run_project.log -J phylobuild ./run_project.sh
```

Then follow it with `tail -f run_project.log`, where it prints `N of M steps (x%) done`.

**If something fails**:
- **Finding the error:** Snakemake retries a failed step once, and keeps going with everything that doesn't depend
  on it. At the end it lists the failed steps. Each step's log is under `results/<donor>/logs/`, and the LSF output
  is under `logs/lsf/`.
- **Carrying on:** fix the problem and run the same command again. Only what's missing or out of date runs.
- **Lock errors:** if the controller was killed, the next run may say the directory is locked. Run the command once
  with `--unlock`, then again without it.
- **Rerunning a step:** to rerun a step that worked, add `--forcerun <rule>` (e.g. `--forcerun build_tree`) or
  delete its outputs.

Two runs can go at the same time from the same folder only if they write different files, e.g. the same config
with different donors, or two configs with different `outdir`.

## 5. Example: a second analysis of a donor from a test project

Say PD23456 has been remapped to a modified reference and the new results are in a CanApps test project, 9999
(`/software/SGP/nst_links/test/9999/`). Two things differ from a normal run:

1. **The BAMs are in the iRODS test database.** `stageBam.pl` only finds them with `-l test`. The BAMs in nst_links
   are links into iRODS that can't be read directly, so they must be staged.
2. **It's a separate analysis of the same donor,** so it needs its own outputs.

Set it up with its own sample sheet and config (`config/local/` is ignored by git, so it's a good place for
project files):

- **Sample sheet:**
  ```bash
  mkdir -p config/local
  python3 workflow/scripts/sanger_sample_sheet.py --donor PD23456 --projects 9999 \
      --nst-links /software/SGP/nst_links/test > config/local/samples_PD23456_remapped.tsv
  ```
  Because the path is a test folder, `irods_db` is filled in as `test`, so the BAMs are staged from the test
  database.
- **Config:** copy `config/config.yaml` to `config/local/config_PD23456_remapped.yaml` and set `samples` to the new
  sheet, `outdir: results_remapped`, `staging_dir: results_remapped/staged_bams` (so it doesn't share staged files
  with a standard run going at the same time), and only PD23456 under `donors`.
- **Reference genome:** set `reference.genome` to the FASTA the BAMs were mapped to. If the modified reference has
  the same contigs and coordinates (e.g. a region replaced with Ns), the other reference files can stay as they are.
- **In silico normal:** if there's no normal mapped to the modified reference, keep the standard one. In the
  pipeline the normal is only a column in the cgpVAF matrices; the filtering and trees use the colonies. If the
  normal is in a test project too, set `normal_project` and `normal_irods_db: test`.

Run it with:

```bash
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf --configfile config/local/config_PD23456_remapped.yaml -n
~/phylobuild-venv/bin/snakemake --profile profiles/sanger-lsf --configfile config/local/config_PD23456_remapped.yaml
```

The results go to `results_remapped/PD23456/`, next to the standard run in `results/PD23456/`, so the two are easy
to compare.

## 6. Using the R package on its own

Every step's output is an `.rds` file that you can load in R:

```r
library(phylobuild)
x = read_mats_and_params("results/PD12345/mats_and_params/mats_and_params_PD12345_none.rds")
names(x)                 # COMB_mats, filter_params
tree = ape::read.tree("results/PD12345/trees/PD12345_none/tree_files/tree_PD12345_none_pval_post_mix_post_dup.tree")
muts = readRDS("results/PD12345/trees/PD12345_none/annotated_muts/annotated_mut_set_PD12345_none_pval_post_mix_post_dup.rds")
```

`read_mats_and_params()` also reads the old `save()` files from the bash workflow. The tree-plotting functions you
know (`plot_tree`, `add_annotation`, `plot_sharing_info`, and so on) are in the package too.
