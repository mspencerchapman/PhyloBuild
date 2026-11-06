# PASS calls per sample, and the union of sites per donor (the bed files cgpVAF genotypes every sample at).
# Replaces import_new_samples_only.R + filter_pass_subs_set_threshold.pl / filter_indels_new.pl.

rule pass_snvs:
    """CaVEMan SNVs with FILTER=PASS, CLPM=0.00 and ASMD above the threshold."""
    input:
        vcf = lambda w: sample_row(w.donor, w.sample)["caveman_vcf"],
    output:
        f"{OUT}/{{donor}}/calls/snv/{{sample}}.pass.vcf",
    params:
        asmd = lambda w: donor_cfg(w.donor)["asmd_threshold"],
    resources:
        mem_mb = 1000,
    shell:
        r"""
        zcat {input.vcf} | awk -F'\t' -v t={params.asmd} '
            !/^#/ && $7 == "PASS" && $8 ~ /(^|;)CLPM=0\.00(;|$)/ {{
                if (match($8, /(^|;)ASMD=[0-9.]+/)) {{
                    split(substr($8, RSTART, RLENGTH), a, "="); if (a[2] + 0 > t) print
                }}
            }}' > {output}
        """

rule pass_indels:
    """Pindel indels with FILTER=PASS."""
    input:
        vcf = lambda w: sample_row(w.donor, w.sample)["pindel_vcf"],
    output:
        f"{OUT}/{{donor}}/calls/indel/{{sample}}.pass.vcf",
    resources:
        mem_mb = 1000,
    shell:
        r"""zcat {input.vcf} | awk -F'\t' '!/^#/ && $7 == "PASS"' > {output}"""

def site_bed_inputs(w):
    calls = expand(f"{OUT}/{{donor}}/calls/{{vartype}}/{{sample}}.pass.vcf", donor=w.donor, vartype=w.vartype, sample=samples_of(w.donor))
    if w.vartype == "snv" and donor_cfg(w.donor)["cgpvaf_snv_sites"] == "hairpin":
        return {"calls": calls, "hairpin": f"{OUT}/{w.donor}/hairpin/{w.donor}_hairpin_pass_sites.bed"}
    return {"calls": calls}

rule site_bed:
    """Sites for cgpVAF: everything called PASS in at least one sample (Chrom, Pos, Ref, Alt).
    With cgpvaf_snv_sites: hairpin, SNV sites must also have passed the hairpin filter, so cgpVAF genotypes
    only the sites that would be kept anyway."""
    input:
        unpack(site_bed_inputs),
    output:
        f"{OUT}/{{donor}}/calls/{{donor}}_{{vartype}}_sites.bed",
    resources:
        mem_mb = 4000,
    run:
        shell("cat {input.calls} | cut -f 1,2,4,5 | LC_ALL=C sort -u > {output}.tmp")
        if "hairpin" in input.keys():
            shell("LC_ALL=C sort -u {input.hairpin} | LC_ALL=C comm -12 {output}.tmp - > {output} && rm {output}.tmp")
        else:
            shell("mv {output}.tmp {output}")
