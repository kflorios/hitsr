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
