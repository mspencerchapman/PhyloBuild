# Chromosome naming and sex-chromosome helpers.
# Everything works with or without the "chr" prefix, so GRCh37 ("X") and GRCh38 ("chrX") data
# go through the same code.

#Pseudo-autosomal regions on X, by genome build
PAR_REGIONS = list(
  hg38 = data.frame(start = c(10001, 155701383), end = c(2781479, 156030895)),
  hg19 = data.frame(start = c(60001, 154931044), end = c(2699520, 155260560))
)

strip_chr = function(chrom) sub("^chr", "", chrom)

is_x = function(chrom) strip_chr(chrom) == "X"
is_y = function(chrom) strip_chr(chrom) == "Y"
is_xy = function(chrom) strip_chr(chrom) %in% c("X", "Y")
is_autosome = function(chrom) strip_chr(chrom) %in% as.character(1:22)

#Guess the build from the chromosome names: GRCh38 data at Sanger uses the "chr" prefix
guess_build = function(chrom) ifelse(any(grepl("^chr", chrom)), "hg38", "hg19")

in_par = function(chrom, pos, build) {
  if (!build %in% names(PAR_REGIONS)) stop("Unknown genome build '", build, "' - use one of ", paste(names(PAR_REGIONS), collapse = ", "))
  r = PAR_REGIONS[[build]]
  is_x(chrom) & ((pos >= r$start[1] & pos <= r$end[1]) | (pos >= r$start[2] & pos <= r$end[2]))
}

#Loci that behave as haploid in a male: Y, and X outside the PARs
is_haploid_in_male = function(chrom, pos, build) {
  is_y(chrom) | (is_x(chrom) & !in_par(chrom, pos, build))
}

#Rows treated as haploid in this donor: in males, Y and non-PAR X (PAR sites behave like autosomes)
haploid_rows = function(COMB_mats) {
  COMB_mats$gender == "male" & is_haploid_in_male(COMB_mats$mat$Chrom, COMB_mats$mat$Pos, guess_build(COMB_mats$mat$Chrom))
}

#Split "Chrom-Pos-Ref-Alt" mut_refs into a data frame
split_mut_ref = function(mut_ref) {
  parts = stringr::str_split(mut_ref, pattern = "-", simplify = TRUE)
  data.frame(mut_ref = mut_ref, Chrom = parts[, 1], Pos = as.integer(parts[, 2]), Ref = parts[, 3], Alt = parts[, 4], stringsAsFactors = FALSE)
}

is.snv = function(mut_ref) {
  parts = split_mut_ref(mut_ref)
  nchar(parts$Ref) == 1 & nchar(parts$Alt) == 1
}
