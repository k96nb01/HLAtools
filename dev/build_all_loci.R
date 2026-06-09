# Build alignments for ALL loci, all alignment types, and time it.
# Used to measure the optimization at full scale and to prove the optimized
# output is identical() to the original across every locus.
#
# Usage:  Rscript build_all_loci.R <label>
#   <label> is "original" or "optimized" -- selects the output filename so the
#   two runs can be compared afterward. The R source files are swapped between
#   runs via git (see the surrounding workflow), so the SAME script builds both.

suppressMessages(pkgload::load_all("C:/GitHub/HLAtools_fast", quiet = TRUE))

label <- commandArgs(trailingOnly = TRUE)
label <- if (length(label)) label[1] else "run"
VERSION <- "3.64.0"
out <- sprintf("C:/GitHub/HLAtools_fast/dev/all_loci_%s.rds", label)

cat(sprintf("[%s] building ALL loci, alignType=all, version=%s ...\n", label, VERSION))
t <- system.time({
  res <- alignmentFull(loci = "all", alignType = "all", version = VERSION)
})
saveRDS(res, out)

cat(sprintf("\n[%s] DONE\n", label))
cat("elapsed (s):", round(unname(t["elapsed"]), 1),
    " user/CPU (s):", round(unname(t["user.self"]), 1), "\n")
cat("loci per type: ",
    paste(sprintf("%s=%d", names(res)[1:4],
                  vapply(res[1:4], length, integer(1))), collapse="  "), "\n")
cat("object size (MB):", round(as.numeric(object.size(res)) / 1024^2, 1), "\n")
cat("saved ->", out, "\n")
