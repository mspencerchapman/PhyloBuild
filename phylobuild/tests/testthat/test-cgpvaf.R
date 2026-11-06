test_that("cgpVAF batches merge by column name, whatever the row order", {
  sites = data.frame(Chrom = c("chr1", "chr2", "chrX"), Pos = c(100, 200, 300), Ref = c("A", "G", "T"), Alt = c("C", "T", "A"))
  f1 = tempfile(fileext = ".tsv"); f2 = tempfile(fileext = ".tsv")
  b1 = write_fake_cgpvaf(f1, sites, c("S1", "S2"), seed = 1)
  b2 = write_fake_cgpvaf(f2, sites[c(3, 1, 2), ], c("S3"), seed = 2)
  merged = merge_cgpvaf_files(c(f1, f2))
  expect_equal(colnames(merged), c("Chrom", "Pos", "Ref", "Alt", "PDv38is_wgs_MTR", "PDv38is_wgs_DEP", "S1_MTR", "S1_DEP", "S2_MTR", "S2_DEP", "S3_MTR", "S3_DEP"))
  expect_equal(merged$S3_MTR, b2$S3_MTR[c(2, 3, 1)]) #rows matched by site, not position in the file
  expect_equal(merged$S1_DEP, b1$S1_DEP)
})

test_that("merging fails loudly if batches were run on different sites", {
  sites = data.frame(Chrom = c("chr1", "chr2"), Pos = c(100, 200), Ref = "A", Alt = "C")
  f1 = tempfile(); f2 = tempfile()
  write_fake_cgpvaf(f1, sites, "S1")
  write_fake_cgpvaf(f2, sites[1, ], "S2")
  expect_error(merge_cgpvaf_files(c(f1, f2)), "Variant sites differ")
})

test_that("import drops the in silico normal from SNV and indel matrices", {
  sites = data.frame(Chrom = "chr1", Pos = 1:5, Ref = "A", Alt = "C")
  snv = tempfile(); indel = tempfile()
  write_merged_cgpvaf(merge_cgpvaf_files({f = tempfile(); write_fake_cgpvaf(f, sites, c("S1", "S2")); f}), snv)
  write_merged_cgpvaf(merge_cgpvaf_files({f = tempfile(); write_fake_cgpvaf(f, sites, c("S1", "S2")); f}), indel)
  cm = import_cgpvaf_SNV_and_INDEL(snv, NULL, ref_ID = "PDv38is_wgs")
  expect_equal(colnames(cm$NV), c("S1", "S2"))
  cm2 = import_cgpvaf_SNV_and_INDEL(snv, indel, ref_ID = "PDv38is_wgs")
  expect_equal(nrow(cm2$mat), 10)
  expect_equal(table(cm2$mat$Mut_type)[["INDEL"]], 5)
})
