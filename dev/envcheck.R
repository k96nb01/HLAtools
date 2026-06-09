# Environment check for the performance-optimization project.
# Verifies required packages are installed and that R can reach the
# ANHIG/IMGTHLA GitHub repo (alignmentFull depends on this network access).
pkgs <- c("DescTools","dplyr","fmsb","rvest","stringr","tibble","xfun",
          "testthat","devtools","pkgload","profvis","bench","curl","memoise")
inst <- rownames(installed.packages())
for (p in pkgs) cat(sprintf("%-12s %s\n", p, if (p %in% inst) "OK" else "MISSING"))
cat("---NETWORK TEST---\n")
invisible(tryCatch({
  u <- url("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/release_version.txt")
  on.exit(close(u))
  x <- readLines(u, n = 5, warn = FALSE)
  cat("network OK, first lines:\n"); cat(x, sep="\n"); cat("\n")
}, error = function(e) { cat("NETWORK FAILED:", conditionMessage(e), "\n") }))
