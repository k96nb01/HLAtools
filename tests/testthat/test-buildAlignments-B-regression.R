# Regression test for the HLA-B cDNA build failure.
#
# A few IMGT alleles (e.g. B*44:568Q, B*51:197) appear in trailing alignment
# blocks that the reference allele does not span, so their split sequences are
# longer or shorter than the reference. do.call(rbind, ...) silently recycled
# these ragged rows and the build later died with
#   "number of items to replace is not a multiple of replacement length"
# which crashed the default alignmentFull(loci = "all"). buildAlignments() now
# normalizes every cDNA/gDNA sequence to the reference length (trim beyond the
# reference, pad short partial sequences with "."), leaving all other loci
# unchanged.
#
# These tests require network access to the ANHIG/IMGTHLA GitHub repo, so they
# are skipped on CRAN and when offline. The version is pinned for reproducibility
# (IMGT release branches are immutable), which makes the exact dimensions stable.

test_that("HLA-B cDNA builds despite ragged trailing-block alleles", {
  skip_on_cran()
  skip_if_offline()

  b <- suppressWarnings(suppressMessages(
    buildAlignments("B", "cDNA", version = "3.64.0")))

  # The cDNA build returns, under the locus, list(codon, cDNA, Version).
  cdna <- b[["B"]][["cDNA"]]
  expect_s3_class(cdna, "data.frame")

  # Rectangular table: 4 metadata columns + the reference (row 1) position
  # columns. The bug produced a width mismatch and aborted before this point.
  expect_equal(ncol(cdna), 1477L)   # 4 + 1473 reference positions (v3.64.0)
  expect_equal(nrow(cdna), 11110L)

  # The two formerly-ragged alleles are present and conform to the table width.
  expect_true(all(c("B*44:568Q", "B*51:197") %in% cdna$allele_name))
  expect_equal(sum(cdna$allele_name == "B*44:568Q"), 1L)
  expect_equal(sum(cdna$allele_name == "B*51:197"), 1L)

  # The reference allele (row 1) is unchanged by the normalization.
  expect_equal(cdna$allele_name[1], "B*07:02:01:01")
})

test_that("buildAlignments reports when it normalizes ragged HLA-B cDNA sequences", {
  skip_on_cran()
  skip_if_offline()
  expect_message(
    suppressWarnings(buildAlignments("B", "cDNA", version = "3.64.0")),
    "did not match the reference length"
  )
})
