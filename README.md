# hitsr
High-performance C++ OpenMP implementation of discrete tabu search for maximum score estimator as an R package

# `hitsr`: High-Performance Tabu Search in R - max score estimator

[![R-CMD-check](https://img.shields.io/badge/R--build-passing-brightgreen.svg)](#)
[![License: GPL v2+](https://img.shields.io/badge/License-GPL%20v2%2B-blue.svg)](#)

`hitsr` is a high-performance R package providing OpenMP-accelerated C++ implementations of discrete Tabu Search for maximum score estimation and high-dimensional optimization problems. It seamlessly bridges fast low-level computational routines with an intuitive R interface.

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

### Installing `hitsr` : Way 1

Install the package directly from your local directory using `devtools`:

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install hitsr from local source
devtools::install("path/to/hitsr")
```

### Installing `hitsr` : Way 2 (preferred)

Install the package directly from the Internet - GitHub (preferred):
Open R/Rstudio as an Administrator.

```R
# Install devtools and remotes if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install hitsr GitHub
devtools::install_github("kflorios/hitsr")
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
| `b0` | `numeric` | `-1.0` | Coefficient for first column of feature matrix. |
| `d` | `numeric` | `10` | Neighborhood depth / search parameter. |
| `iSeed` | `integer` | `123456` | Pseudo-random number generator seed. |

### Return Value

A named `list` containing:
- `attributes`: Vector of selected observations indices (1-based R indexing).
- `coeffs`: Estimated model parameters/coefficients for the selected subset.
- `score`: Evaluation score associated with the optimal observations set.

---

## Complete Example: Multi-Seed Benchmark (20 Random Starts)

This benchmark demonstrates how to execute 20 distinct random restarts using `hitsr`, summarize overall solution quality, and compute selection frequency statistics for each attribute.

```R
library(hitsr)

# 1. Load and preprocess data
X_raw <- as.matrix(read.table("X.txt"))
y_raw <- as.matrix(read.table("y.txt"))

X <- X_raw[, 2:ncol(X_raw)]
y <- as.integer(y_raw[, 2])

# 2. Parameters
b0 <- -1.0
d <- 1e4
num_runs <- 20

# Generate 20 distinct random seeds
set.seed(42) # For reproducible seed generation
seeds <- sample.int(1e6, num_runs)

# 3. Storage for results
results_list <- vector("list", num_runs)

cat("Starting 20 runs...\n")

# 4. Run loop
for (i in seq_len(num_runs)) {
  current_seed <- seeds[i]
  
  # Execute Tabu Search
  res <- run_tabu_search(X = X, y = y, b0 = b0, d = d, iSeed = current_seed)
  
  # Format selected attributes as a string (e.g., "1, 4, 7")
  attr_str <- paste(res$attributes, collapse = ", ")
  
  # Save metrics
  results_list[[i]] <- data.frame(
    Run = i,
    Seed = current_seed,
    Score = res$score,
    Coeffs = res$coeffs,
    Num_Attributes = length(res$attributes),
    Attributes = attr_str,
    stringsAsFactors = FALSE
  )
}

# 5. Combine into a single data frame
results_df <- do.call(rbind, results_list)

# 6. Display results table
print(results_df, row.names = FALSE)

# 7. Summary statistics
cat("\n=== Summary Across 20 Runs ===\n")
cat("Best Score Found: ", max(results_df$Score), "\n")
cat("Mean Score:       ", mean(results_df$Score), "\n")
cat("Score Std Dev:    ", sd(results_df$Score), "\n")
```

---

## License

---
Distributed under the GPL-2 / GPL-3 License. See `LICENSE` for details.

## Citation

If you use `hitsr` in your research, please cite the underlying methodology papers:

```bibtex
@article{florios2025hits,
  title={HITS: Hyperplanes intersection tabu search for maximum score estimation},
  author={Florios, Kostas and Louka, Alexandros and Bilias, Yannis},
  journal={SoftwareX},
  volume={30},
  pages={102164},
  year={2025},
  publisher={Elsevier}
}

@article{manski1975maximum,
  title={Maximum score estimation of the stochastic utility model of choice},
  author={Manski, Charles F},
  journal={Journal of econometrics},
  volume={3},
  number={3},
  pages={205--228},
  year={1975},
  publisher={Elsevier}
}
  



