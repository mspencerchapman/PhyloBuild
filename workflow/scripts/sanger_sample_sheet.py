#!/usr/bin/env python3
"""Build PhyloBuild sample sheet rows from the Sanger nst_links folders.

Example (one donor, samples spread over two projects, one extra sample, one left out):
    python3 workflow/scripts/sanger_sample_sheet.py --donor PD12345 --projects 1234,1235 \\
        --include PD99999a_lo0001 --exclude PD12345a_lo0010 >> config/samples.tsv

Prints a header line too unless --no-header is given (use it when adding a second donor to an existing sheet).
For a test project, point --nst-links at the test folder; the irods_db column is then set to test, so that
stageBam.pl looks in the test database:
    python3 workflow/scripts/sanger_sample_sheet.py --donor PD23456 --projects 9999 \
        --nst-links /software/SGP/nst_links/test > config/samples_PD23456_remapped.tsv
Samples whose CaVEMan/Pindel VCFs aren't there yet are listed on stderr and left out.
"""
import argparse
import os
import sys

NST = "/software/SGP/nst_links/live"
COLUMNS = ["donor", "sample", "project", "bam", "caveman_vcf", "pindel_vcf", "irods_db"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donor", required=True, help="donor ID used in the config (e.g. PD12345)")
    ap.add_argument("--projects", required=True, help="comma-separated CanApps project numbers")
    ap.add_argument("--pd-ids", help="comma-separated PD IDs whose samples to include (default: the donor ID)")
    ap.add_argument("--include", default="", help="comma-separated extra sample names, even if they don't match")
    ap.add_argument("--exclude", default="", help="comma-separated sample names to leave out")
    ap.add_argument("--nst-links", default=NST, help=f"nst_links folder (default {NST})")
    ap.add_argument("--irods-db", choices=["live", "test"],
                    help="iRODS database for stageBam.pl (default: test if the nst_links folder is a test one, else live)")
    ap.add_argument("--no-header", action="store_true")
    args = ap.parse_args()

    pd_ids = (args.pd_ids or args.donor).split(",")
    irods_db = args.irods_db or ("test" if "/test" in os.path.normpath(args.nst_links) else "live")
    include = set(filter(None, args.include.split(",")))
    exclude = set(filter(None, args.exclude.split(",")))

    rows, not_ready = [], []
    for project in args.projects.split(","):
        pdir = os.path.join(args.nst_links, project)
        if not os.path.isdir(pdir):
            sys.exit(f"Project folder {pdir} not found")
        for sample in sorted(os.listdir(pdir)):
            matches = any(sample.startswith(p) for p in pd_ids) or sample in include
            if not matches or sample in exclude:
                continue
            files = {
                "bam": f"{sample}.sample.dupmarked.bam",
                "caveman_vcf": f"{sample}.caveman_c.annot.vcf.gz",
                "pindel_vcf": f"{sample}.pindel.annot.vcf.gz",
            }
            paths = {k: os.path.join(pdir, sample, v) for k, v in files.items()}
            missing = [files[k] for k in ("caveman_vcf", "pindel_vcf") if not os.path.exists(paths[k])]
            if missing:
                not_ready.append(f"{project}/{sample} (missing {', '.join(missing)})")
                continue
            rows.append({"donor": args.donor, "sample": sample, "project": project, **paths, "irods_db": irods_db})

    found = {r["sample"] for r in rows}
    for s in sorted(include - found):
        print(f"WARNING: --include sample {s} not found or not ready", file=sys.stderr)
    if not_ready:
        print(f"{len(not_ready)} matching samples are not ready yet and were left out:", file=sys.stderr)
        for s in not_ready:
            print(f"  {s}", file=sys.stderr)
    print(f"{len(rows)} samples written for donor {args.donor}", file=sys.stderr)

    if not args.no_header:
        print("\t".join(COLUMNS))
    for r in rows:
        print("\t".join(r[c] for c in COLUMNS))


if __name__ == "__main__":
    main()
