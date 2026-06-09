# Check whether the position labels in corr_table[1,] (used as the match key in
# the INDEL/EXONB labeling loops) are unique. If unique, a match()-based
# vectorisation is exactly equivalent to the per-indel loop.
suppressMessages(pkgload::load_all("C:/GitHub/HLAtools_fast", quiet = TRUE))
VERSION <- "3.64.0"

for (lc_src in list(c("A","cDNA"), c("A","gDNA"), c("DRB1","gDNA"), c("DPB1","cDNA"))) {
  ct <- suppressWarnings(buildAlignments(lc_src[1], lc_src[2], version = VERSION,
                                         return_corr_table = TRUE))
  ctdf <- ct[[1]][[1]]                         # the corr_table data.frame
  row1 <- as.character(unlist(ctdf[1, ]))      # position labels
  n    <- length(row1)
  ndup <- n - length(unique(row1))
  cat(sprintf("%-5s/%-4s : ncol=%d, duplicate labels=%d %s\n",
              lc_src[1], lc_src[2], n, ndup,
              if (ndup == 0) "(UNIQUE -> match() safe)" else "(*** has dups ***)"))
}
