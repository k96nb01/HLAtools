# Compare the ALL-loci builds from the original vs optimized code and assert
# they are identical() across every locus and alignment type. This is the
# full-scale proof that the optimization preserves behavior.
orig <- readRDS("C:/GitHub/HLAtools_fast/dev/all_loci_original.rds")
opt  <- readRDS("C:/GitHub/HLAtools_fast/dev/all_loci_optimized.rds")

cat("Comparing all-loci builds (original vs optimized)...\n")
if (identical(orig, opt)) {
  cat("\nRESULT: IDENTICAL across all loci and alignment types. Behavior preserved at full scale.\n")
  # Summarize coverage
  for (typ in c("prot","codon","nuc","gen")) {
    cat(sprintf("  %-6s : %d loci\n", typ, length(opt[[typ]])))
  }
  quit(status = 0)
} else {
  cat("\nRESULT: *** MISMATCH ***\n")
  for (typ in names(orig)) {
    if (!identical(orig[[typ]], opt[[typ]])) {
      cat(sprintf("  [%s] differs\n", typ))
      if (is.list(orig[[typ]])) {
        for (lc in names(orig[[typ]])) {
          if (!identical(orig[[typ]][[lc]], opt[[typ]][[lc]]))
            cat(sprintf("      -> locus %s differs\n", lc))
        }
      }
    }
  }
  quit(status = 1)
}
