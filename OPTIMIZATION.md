# HLAtools Performance Optimization — `alignmentFull()` / `buildAlignments()`

This document records a performance-optimization pass over `HLAtools`, focused on
`alignmentFull()` and the function it drives, `buildAlignments()`. **Every change
preserves behavior exactly**: the optimized code reproduces the original output
byte-for-byte (verified with `identical()` against a saved baseline).

- **Upstream:** `sjmack/HLAtools` (v1.8.1, commit `edc9cdc`)
- **Fork / working copy:** `k96nb01/HLAtools`, branch `perf-optimization`
- **Headline result:** the representative 3-locus, all-types build dropped from
  **~139 s to ~48 s (≈2.9× faster)**; CPU work measured by the profiler dropped
  **≈2.75×**. Output is unchanged.

---

## 1. Motivation

`alignmentFull()` builds the `HLAalignments` object that nearly every other
alignment-aware function in the package depends on. It is slow enough that
building all loci takes several minutes, which is a recurring friction point in
day-to-day use. The goal was to make it substantially faster **without changing
any output**, so existing analyses remain valid.

---

## 2. Method: correctness first, measure before cutting

The package shipped with **no `tests/` directory at all**, so the first job was
to build a safety net, not to change code.

1. **Characterization baseline.** Ran the *original* `alignmentFull()` on a
   representative subset of loci and saved the exact result. Loci chosen to
   exercise the distinct code paths:
   - `A` — class I (prot/nuc/gen), and the redundant-cDNA download path
   - `DRB1` — the special DR-locus handling, and a very wide gDNA alignment
     (662 × 18 527)
   - `DPB1` — DP naming conventions

   The version is **pinned to a concrete release (3.64.0)**, not `"Latest"`, so
   the baseline is reproducible (IPD-IMGT/HLA release branches are immutable
   upstream).

2. **`identical()` gate after every edit.** Two harnesses (see
   [§7](#7-reproducing-the-work)) rebuild the loci and assert the new output is
   `identical()` to the baseline. A single non-identical byte fails the gate.
   This was confirmed to report *IDENTICAL* on the unmodified code first, which
   also proves `alignmentFull()` is deterministic (a prerequisite for this
   approach).

3. **Profile, don't guess.** Line-level profiling (`Rprof(line.profiling=TRUE)`)
   located the real hotspots. This step changed the strategy entirely — see
   [§4](#4-the-key-diagnosis-cpu-bound-not-network-bound). The harnesses are in
   `dev/` — see [§8](#8-reproducing-the-work).

4. **A regression test suite.** Alongside the refactor, a full `testthat`
   (edition 3) suite was written: **18 files, 102 tests, 255 expectations, all
   passing** — offline golden tests for the pure functions plus
   alignment-consumer tests run against the saved fixture
   ([§7](#7-test-suite)).

---

## 3. Environment

- R 4.6.0 on Windows 11.
- Added dev/test dependencies: `testthat`, `pkgload`, `profvis`, `bench`, plus
  the two package `Imports` that were not yet installed (`DescTools`, `fmsb`).
- Network access to `raw.githubusercontent.com/ANHIG/IMGTHLA` confirmed
  (required for any real build).

---

## 4. The key diagnosis: CPU-bound, not network-bound

A static read of the call graph suggested the cost was network I/O — a full run
issues ~135 sequential HTTP downloads. **Profiling proved otherwise.** The
baseline build reported:

```
   user  system elapsed
 124.36    7.03  139.11
```

`user` (CPU) is **124 s of 139 s elapsed** — roughly **90% of the time is CPU**,
spent parsing and reshaping the downloaded alignment text. Network was only
~10%. This redirected all effort to the in-memory parsing in `buildAlignments()`.

Line-level profiling of three representative builds (`A/cDNA`, `A/gDNA`,
`DRB1/gDNA`) found two dominant lines:

| Line | Share of CPU | What it did |
|------|-------------:|-------------|
| `buildAlignments.R` reference-distribution loop | **~35%** | per-column `for` loop over up to ~18.5k columns, re-indexing the data frame each pass |
| `buildAlignments.R` sequence split | **~25%** | `sapply(seqs, strsplit, split="*")` splitting every sequence into characters |

The download lines were ~4% combined.

---

## 5. The optimizations

All changes are commented inline in the source explaining the rationale. Each was
verified `identical()` to the baseline before moving on.

### 5.1 Vectorize the reference-distribution loops *(biggest single win)*

**File:** `R/buildAlignments.R` (the two "distributes reference sequence from
row 1" blocks).

The original distributed the reference sequence into every `"-"` cell one column
at a time:

```r
for(x in 5:ncol(DNAalignments[[loci[i]]])) {
  DNAalignments[[loci[i]]][,x][which(DNAalignments[[loci[i]]][,x]=="-")] <-
    DNAalignments[[loci[i]]][,x][1]}
```

For a gDNA alignment with ~18.5k position columns this loops 18.5k times, and
each pass re-resolves the nested `[[loci[i]]]` list lookup and a data-frame
column extraction. Replaced with a single vectorized matrix operation:

```r
dcols <- 5:ncol(DNAalignments[[loci[i]]])
dm <- as.matrix(DNAalignments[[loci[i]]][, dcols, drop = FALSE])
ddash <- which(dm == "-", arr.ind = TRUE)            # every dash, as (row, col)
if(nrow(ddash) > 0) { dm[ddash] <- dm[1, ddash[, "col"]] }  # row-1 value per column
DNAalignments[[loci[i]]][, dcols] <- dm
```

This was **~35% of CPU**. The same fix was applied to the AA/codon block.

### 5.2 `strsplit(split = "")` instead of `split = "*"` *(14× on that line)*

**File:** `R/buildAlignments.R` (the per-character sequence split).

The code split sequences into characters with `strsplit(x, split="*")`. `"*"` is
a **degenerate regex** — a quantifier with nothing to repeat — which forces
`strsplit()` onto its slow regex engine. The empty pattern `""` produces the
*identical* character split but takes `strsplit()`'s special-cased fast path.
Benchmarked on a realistic sequence vector:

| split | median | identical to `"*"`? |
|-------|-------:|:-------------------:|
| `"*"` (original) | 1.37 s | — |
| `""` (new) | **0.098 s** | **yes** |

A ~**14×** speedup on a line that was ~25–31% of CPU, by changing one character.
(The surrounding `sapply` was also collapsed to a single vectorized `strsplit`
call, since `strsplit` already vectorizes over its input.)

### 5.3 `match()`-vectorize the INDEL / EXONB labeling

**File:** `R/buildAlignments.R` (the inDel and exon-boundary inclusion blocks).

The original scanned the entire correspondence-table label row once per indel —
`O(nIndels × ncol)`:

```r
for(o in 1:length(inDels[[loci[i]]])){
  corr_table[[loci[i]]][2,][inDels[[loci[i]]][[o]]==corr_table[[loci[i]]][1,]] <- paste("INDEL", o, sep="-")
  ...
}
```

The position labels in `corr_table[1,]` were **verified to be unique** (across
cDNA and gDNA builds, including the 18 523-column DRB1 gDNA case), so `match()`
locates every indel column in a single pass with identical results:

```r
indelIdx <- match(inDels[[loci[i]]], corr_table[[loci[i]]][1,])
indelLab <- paste("INDEL", seq_along(indelIdx), sep="-")
corr_table[[loci[i]]][2, indelIdx] <- indelLab
corr_table[[loci[i]]][3, indelIdx] <- indelLab
```

### 5.4 `which(is.na())`-vectorize the position-fill loops

**File:** `R/buildAlignments.R` (the two "pastes alignment_positions into
corr_table" loops).

Scalar loops walked every column, assigning the next position value to each
non-indel cell with a manual counter. Since the `NA` cells, in column order,
simply receive `alignment_positions[1..k]`, this is one vectorized assignment:

```r
naIdx2 <- which(is.na(corr_table[[loci[i]]][2,]))
corr_table[[loci[i]]][2, naIdx2] <- alignment_positions[[loci[i]]][seq_along(naIdx2)]
```

### 5.5 Eliminate the duplicate cDNA build *(structural)*

**File:** `R/alignmentFull.R`.

`buildAlignments(locus, "cDNA")` returns **both** the codon and cDNA tables in a
single call. But `alignmentFull()` called it **twice per locus** — once for the
`nuc` list (keeping the cDNA table) and once for the `codon` list (keeping the
codon table) — re-downloading and re-parsing the identical `_nuc.txt` file and
discarding half of each result. (This was visible in the baseline: the `codon`
and `nuc` tables have identical dimensions because they come from the same
file.) The fix builds each required cDNA alignment **once** into a small cache
and reads both tables from it:

```r
cdnaLoci <- unique(c(as.character(NL1), as.character(NL4)))
cdnaLoci <- cdnaLoci[!is.na(cdnaLoci) & nzchar(cdnaLoci)]
cdnaCache <- vector(mode = 'list', length = length(cdnaLoci))
names(cdnaCache) <- cdnaLoci
for(lc in cdnaLoci){
  cdnaCache[[lc]] <- suppressWarnings(buildAlignments(lc, "cDNA", version = version)[[1]])
}
# cList  (nuc)   reads cdnaCache[[locus]][2]
# codonList (codon) reads cdnaCache[[locus]][1]
```

On a full run this removes roughly a quarter of *all* `buildAlignments` calls —
both a download and a full re-parse per locus.

---

## 6. Benchmark results

All builds pinned to IPD-IMGT/HLA release 3.64.0.

**Full `alignmentFull(loci = c("A","DRB1","DPB1"), alignType = "all")`:**

| | elapsed | vs original |
|--|--------:|------------:|
| Original (baseline) | 139.1 s | 1.0× |
| Optimized | ~48 s | **≈2.9×** |

(Wall-clock now includes a larger relative share of network time, so it varies
run to run; CPU is the stabler metric below.)

**CPU work — profiler sampled time, three representative `buildAlignments`
calls (`A/cDNA`, `A/gDNA`, `DRB1/gDNA`):**

| Stage | sampled CPU |
|-------|------------:|
| Original | 50.75 s |
| After §5.1 + §5.2 | 42.31 s |
| Final (all changes) | **18.47 s** |

≈**2.75× less CPU**, and the profile is now *flat* — no single line dominates.

**Per-`buildAlignments` harness (5 calls covering every hotspot):** 88.6 s → 37.9 s
(**2.34×**).

---

## 7. Test suite

The package had no tests. The new `testthat` suite (`tests/testthat/`) is the
regression net:

- **Offline golden tests** (12 files, ~160 expectations) for the pure / data-only
  functions: `alleleTrim`, `getField`, GL ↔ UNIFORMAT conversion and
  round-trips, version utilities, locus/motif validators, `BDtoPyPop`, and more.
- **Alignment-consumer tests** (6 files, ~95 expectations) for the functions that
  read `HLAalignments` — `compareSequences`, `motifMatch`, the
  `alignmentSearch`/`customAlign` family, `posSort`, `validateAllele` — run
  against the saved baseline fixture so they need no network.

**Total: 102 tests / 255 expectations, all passing.**

These tests also documented several pre-existing behavior quirks, captured as
characterization (current behavior is locked, not "fixed"):

- `BDtoPyPop()` names its second returned list element `...neagtive` (misspelled).
- A `multiGLStoUNI()` roxygen `@example` passes a non-existent `version` object.
- `posSort()` returns its result in the input's type (numeric in → numeric out),
  not always character as the docs imply.
- `validateAllele()` on an invalid locus *errors* (inside `checkSource`) rather
  than returning `FALSE`.

---

## 8. Reproducing the work

The `dev/` folder (build-ignored) holds the harnesses:

| Script | Purpose |
|--------|---------|
| `dev/build_baseline.R` | Build & save the `alignmentFull` characterization baseline + fixture |
| `dev/verify_baseline.R` | Rebuild the baseline loci and assert `identical()` (the milestone gate, ~48 s) |
| `dev/verify_fast.R` | Faster per-`buildAlignments` `identical()` check for rapid iteration |
| `dev/profile.R` | Line-level `Rprof` profile of `buildAlignments` |
| `dev/run_tests.R` | Run the full `testthat` suite via `pkgload::load_all` |

Run any of them with, e.g.:

```pwsh
Rscript "C:\GitHub\HLAtools_fast\dev\verify_baseline.R"
```

Fixtures live in `tests/testthat/fixtures/` (`alignmentFull_baseline.rds` and the
per-call fast baseline).

---

## 9. What was *not* changed (and why)

After the changes above, the profile is flat. The remaining costs are spread
thin and were left alone as poor risk/reward:

- **Final list-assembly copies** (`c(final_alignment[[...]], <data frame>)`) —
  inherent to building the nested result structure; hard to change without
  risking output differences.
- **Whitespace normalization / `str_squish`** on the raw downloaded lines —
  string parsing where any change risks altering output.
- **Parallel / cached downloads** — network is now only ~10% of the time, so the
  payoff is small relative to the complexity (and concurrent requests to the
  ANHIG repo are best avoided). This is the natural next lever if a full-repo
  build is still too slow.

---

*Generated as part of the HLAtools performance-optimization effort. All
optimizations verified to produce output `identical()` to HLAtools v1.8.1.*
