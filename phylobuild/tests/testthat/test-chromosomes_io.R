test_that("chromosome helpers work with and without the chr prefix", {
  expect_equal(is_xy(c("X", "chrX", "chrY", "1", "chr1")), c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_true(in_par("chrX", 2000000, "hg38"))
  expect_false(in_par("chrX", 2781480, "hg38"))
  expect_true(in_par("X", 2699520, "hg19"))
  expect_false(in_par("chr1", 2000000, "hg38"))
  expect_equal(is.snv(c("chr1-10-A-C", "chr1-10-AT-A")), c(TRUE, FALSE))
})

test_that("objects round-trip as .rds, and legacy save() files can be written and read", {
  x = list(COMB_mats = list(a = 1), filter_params = data.frame(b = 2))
  p = file.path(tempdir(), "mats_and_params_TEST.rds")
  save_object(x, p, legacy_name = TRUE)
  expect_equal(read_mats_and_params(p), x)
  expect_equal(read_mats_and_params(sub("\\.rds$", "", p)), x)
})

test_that("mutations are assigned to the branches they were simulated on (treemut)", {
  tree = ape::read.tree(text = "((s1:1,s2:1):1,(s3:1,(s4:1,s5:1):1):1);")
  tree = add_ancestral_outgroup(tree)
  tree$edge.length = rep(1, nrow(tree$edge))
  df = treemut::reconstruct_genotype_summary(tree)
  #one 0/1 genotype profile per branch, over df$samples; simulate two mutations on each branch carried by >=1 sample
  profiles = lapply(df$df$profile, function(p) as.integer(strsplit(p, "")[[1]]))
  branches = which(sapply(profiles, sum) > 0 & sapply(profiles, function(p) p[df$samples == "Ancestral"] == 0))
  truth = rep(branches, each = 2)
  mtr = do.call(rbind, lapply(truth, function(b) profiles[[b]] * 15))
  depth = matrix(30, nrow = length(truth), ncol = length(df$samples))
  colnames(mtr) = colnames(depth) = df$samples
  depth[, "Ancestral"] = 10
  res = assign_mutations(mtr, depth, df, error_rate = ifelse(df$samples == "Ancestral", 1e-6, 0.01))
  expect_equal(res$summary$edge_ml, truth)
  expect_true(all(c("edge_ml", "pval") %in% names(res$summary)))
})

test_that("X PAR sites are diploid in males, non-PAR X and Y haploid", {
  cm = list(gender = "male", mat = data.frame(Chrom = c("chrX", "chrX", "chrY", "chr1"), Pos = c(2000000, 50000000, 3000000, 100)))
  expect_equal(haploid_rows(cm), c(FALSE, TRUE, TRUE, FALSE))
  cm$gender = "female"
  expect_equal(haploid_rows(cm), rep(FALSE, 4))
})

test_that("mean burden for duplicates/age scaling leaves out the Ancestral tip", {
  t = ape::read.tree(text = "((A:10,B:20):5,Ancestral:0);")
  expect_equal(mean_sample_burden(t), 20)
})
