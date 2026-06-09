# Test whether a faster splitter reproduces strsplit(x, split="*") EXACTLY,
# and benchmark candidates on a realistic locus sequence vector.
suppressMessages(pkgload::load_all("C:/GitHub/HLAtools_fast", quiet = TRUE))
suppressWarnings(suppressMessages(library(stringi)))

# Build a realistic vector of aligned sequences the way buildAlignments does
# right before line 413: pull from a real cDNA build's internal stage is hard,
# so synthesize representative sequences (mixed letters, -, ., *, |) of equal
# and unequal lengths, plus pull real ones if available.
set.seed(1)  # NB: deterministic, Math.random not used
alphabet <- c("A","C","G","T","-",".","*","|","R","N")
mkseq <- function(L) paste(sample(alphabet, L, replace=TRUE), collapse="")
seqs <- vapply(rep(c(300, 301, 1579), c(2000,5,2000)), mkseq, character(1))

ref  <- lapply(seqs, function(s) strsplit(s, split="*")[[1]])   # base, the truth
cand_base   <- strsplit(seqs, split="*")
cand_empty  <- strsplit(seqs, split="")
cand_stri   <- stri_split_regex(seqs, "")
cand_stri_b <- stri_split_boundaries(seqs, type="character")

cat("strsplit('*')   == truth:", identical(cand_base, ref), "\n")
cat("strsplit('')    == truth:", identical(cand_empty, ref), "\n")
cat("stri_split_regex('') == truth:", identical(cand_stri, ref), "\n")
cat("stri_boundaries  == truth:", identical(cand_stri_b, ref), "\n")

# Show a mismatch example if any
if (!identical(cand_stri, ref)) {
  for (k in seq_along(ref)) if (!identical(cand_stri[[k]], ref[[k]])) {
    cat("first stri mismatch at", k, "\n")
    cat(" base head:", paste(head(ref[[k]]),collapse="|"),
        " tail:", paste(tail(ref[[k]]),collapse="|"), "len", length(ref[[k]]), "\n")
    cat(" stri head:", paste(head(cand_stri[[k]]),collapse="|"),
        " tail:", paste(tail(cand_stri[[k]]),collapse="|"), "len", length(cand_stri[[k]]), "\n")
    break
  }
}

cat("\n--- benchmark (median) ---\n")
b <- bench::mark(
  base_star  = strsplit(seqs, split="*"),
  base_empty = strsplit(seqs, split=""),
  stri_regex = stri_split_regex(seqs, ""),
  check = FALSE, iterations = 5
)
print(b[, c("expression","min","median","itr/sec","mem_alloc")])
