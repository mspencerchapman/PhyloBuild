#Simulate read counts for a set of clonal colonies: shared mutations, private mutations and germline SNPs.
#Sex sets X/Y depth and genotypes; chromosome names use the "chr" prefix (GRCh38 style).
simulate_COMB_mats = function(sex = "female", n_samples = 12, depth = 20, seed = 1) {
  set.seed(seed)
  chroms = c(rep(paste0("chr", 1:22), each = 40), rep("chrX", 60), if (sex == "male") rep("chrY", 10))
  n = length(chroms)
  pos = sample(3e6:1.5e8, n)
  haploid = (chroms == "chrX" | chroms == "chrY") & sex == "male"
  type = sample(c("germline", "shared", "private"), n, replace = TRUE, prob = c(0.3, 0.2, 0.5))
  present = t(sapply(seq_len(n), function(i) switch(type[i],
    germline = rep(TRUE, n_samples),
    shared = seq_len(n_samples) <= n_samples / 2,
    private = seq_len(n_samples) == sample.int(n_samples, 1))))
  vaf = ifelse(haploid, 0.95, 0.5)
  NR = matrix(rpois(n * n_samples, ifelse(haploid, depth / 2, depth)), nrow = n)
  NV = matrix(rbinom(n * n_samples, NR, ifelse(present, vaf, 0.001)), nrow = n)
  samples = sprintf("PD00001a_lo%04d", seq_len(n_samples))
  mat = data.frame(Chrom = chroms, Pos = pos, Ref = "A", Alt = "C", stringsAsFactors = FALSE)
  mat$mut_ref = paste(mat$Chrom, mat$Pos, mat$Ref, mat$Alt, sep = "-")
  mat$Mut_type = "SNV"
  NV = as.data.frame(NV); NR = as.data.frame(NR)
  colnames(NV) = colnames(NR) = samples
  rownames(NV) = rownames(NR) = mat$mut_ref
  list(mat = mat, NV = NV, NR = NR)
}

#Write a fake cgpVAF *_vaf.tsv (with "##" header lines and the in silico normal first)
write_fake_cgpvaf = function(file, sites, samples, normal = "PDv38is_wgs", seed = 1) {
  set.seed(seed)
  cols = c("#Normal", "VariantID", "Chrom", "Pos", "Ref", "Alt")
  body = data.frame(normal, paste0("id", seq_len(nrow(sites))), sites$Chrom, sites$Pos, sites$Ref, sites$Alt, check.names = FALSE)
  colnames(body) = cols
  for (s in c(normal, samples)) {
    body[[paste0(s, "_FAZ")]] = 0
    body[[paste0(s, "_MTR")]] = rpois(nrow(sites), 5)
    body[[paste0(s, "_WTR")]] = 0
    body[[paste0(s, "_DEP")]] = body[[paste0(s, "_MTR")]] + rpois(nrow(sites), 10)
    body[[paste0(s, "_VAF")]] = 0
  }
  writeLines(c("##fileformat=VCFv4.1", "##source=test"), file)
  suppressWarnings(write.table(body, file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
  body
}
