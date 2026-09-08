# hitsr
High-performance C++ OpenMP implementation of discrete tabu search for maximum score estimator as an R package

# `hitsr`: High-Performance Discrete Tabu Search in R

[![R-CMD-check](https://img.shields.io/badge/R--build-passing-brightgreen.svg)](#)
[![License: GPL v2+](https://img.shields.io/badge/License-GPL%20v2%2B-blue.svg)](#)

`hitsr` is a high-performance R package providing OpenMP-accelerated C++ implementations of discrete Tabu Search for feature selection and high-dimensional optimization problems. It seamlessly bridges fast low-level computational routines with an intuitive R interface.

---

## Key Features

- **OpenMP Parallelization:** Fast parallel neighborhood evaluations powered by multi-core computational backends.
- **Strict Matrix Memory Alignment:** Direct zero-copy data passing between R and C++ runtime buffers for minimal overhead.
- **Exact Numerical Equivalence:** Fully aligned indexing, initialization schemes, and linear algebra routines guarantee identical outputs across implementations.

---

## Installation

### Prerequisites

Ensure you have a C++ compiler supporting C++11 and OpenMP installed:

- **Windows:** Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version.
- **macOS:** Install Xcode Command Line Tools (`xcode-select --install`) and an OpenMP-enabled `g++`/`clang`.
- **Linux:** Install `r-base-dev` and `g++`.

### Installing `hitsr`

Install the package directly from your local directory using `devtools`:

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install hitsr from local source
devtools::install("path/to/hitsr")
```

---

## Main Functionality

```R
run_tabu_search(X, y, b0 = -1.0, d = 10, iSeed = 123456)
```

### Arguments

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `X` | `matrix` | *Required* | Feature matrix of size $T \times p_{\text{in}}$. |
| `y` | `vector` | *Required* | Integer response/target vector of length $T$. |
| `b0` | `numeric` | `-1.0` | Intercept parameter / threshold value. |
| `d` | `integer` | `10` | Neighborhood depth / search parameter. |
| `iSeed` | `integer` | `123456` | Pseudo-random number generator seed. |

### Return Value

A named `list` containing:
- `attributes`: Vector of selected feature indices (1-based R indexing).
- `coeffs`: Estimated model parameters/coefficients for the selected subset.
- `score`: Evaluation score associated with the optimal feature set.

---

## Complete Example: Multi-Seed Benchmark (20 Random Starts)

This benchmark demonstrates how to execute 20 distinct random restarts using `hitsr`, summarize overall solution quality, and compute selection frequency statistics for each attribute.

```R
library(hitsr)

# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------
# Expecting standard tabular input files X.txt and y.txt
X_raw <- as.matrix(read.table("X.txt"))
y_raw <- as.matrix(read.table("y.txt"))

# Extract design matrix X and target y
X <- X_raw[, 2:ncol(X_raw)]
y <- as.integer(y_raw[, 2])

# -----------------------------------------------------------------------------
# 2. Execution Setup
# -----------------------------------------------------------------------------
b0 <- -1.0
d  <- 10
num_runs <- 20

# Generate 20 reproducible seeds
set.seed(42)
seeds <- sample.int(1e6, num_runs)

# Preallocate container
results_list <- vector("list", num_runs)

cat("Running 20 Tabu Search iterations across different random seeds...

")

# -----------------------------------------------------------------------------
# 3. Multi-Seed Loop
# -----------------------------------------------------------------------------
for (i in seq_len(num_runs)) {
  current_seed <- seeds[i]
  
  # Run discrete tabu search
  res <- run_tabu_search(X = X, y = y, b0 = b0, d = d, iSeed = current_seed)
  
  # Store per-run statistics
  results_list[[i]] <- data.frame(
    Run            = i,
    Seed           = current_seed,
    Score          = res$score,
    Num_Attributes = length(res$attributes),
    Attributes     = paste(res$attributes, collapse = ", "),
    stringsAsFactors = FALSE
  )
}

# Combine into a single structured summary table
results_df <- do.call(rbind, results_list)

# Display tabular output
print(results_df, row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Global Performance Summary
# -----------------------------------------------------------------------------
cat("
================ Summary Across 20 Runs ================
")
cat(sprintf("Best Score Found:  %10.4f
", max(results_df$Score)))
cat(sprintf("Mean Score:        %10.4f
", mean(results_df$Score)))
cat(sprintf("Score Std. Dev.:   %10.4f
", sd(results_df$Score)))
cat("========================================================
")

# -----------------------------------------------------------------------------
# 5. Attribute Selection Frequency Analysis
# -----------------------------------------------------------------------------
# Unlist all selected attributes across all 20 runs
all_selected <- unlist(lapply(results_list, function(r) {
  as.numeric(strsplit(r$Attributes, ", ")[[1]])
}))

# Calculate selection counts and percentages
attr_counts <- table(all_selected)
attr_summary <- data.frame(
  Attribute  = as.numeric(names(attr_counts)),
  Frequency  = as.numeric(attr_counts),
  Percentage = round(as.numeric(attr_counts) / num_runs * 100, 1)
)

# Sort by frequency (descending)
attr_summary <- attr_summary[order(-attr_summary$Frequency), ]

cat("
Top Selected Attributes Across Runs:
")
print(attr_summary, row.names = FALSE)
```

---

## License

Distributed under the GPL-2 / GPL-3 License. See `LICENSE` for details.
