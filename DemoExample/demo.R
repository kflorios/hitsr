library(hitsr)

# 1. Load Data
# Adjust these paths if X.txt and y.txt are stored in a different directory
X_raw <- as.matrix(read.table("X.txt"))
y_raw <- as.matrix(read.table("y.txt"))

# 2. Preprocess Data (Matching Python demo logic)
# Drop the first column of X and extract the target column from y
X <- X_raw[, 2:ncol(X_raw)]
y <- as.integer(y_raw[, 2])

# 3. Execution Parameters
b0 <- -1.0
d <- 10
iSeed <- 123456

# 4. Run Tabu Search
result <- run_tabu_search(X = X, y = y, b0 = b0, d = d, iSeed = iSeed)

# 5. Display Outputs
cat("\n=== Tabu Search Results ===\n")
cat("Score:", result$score, "\n\n")

cat("Selected Attributes (1-based R indices):\n")
print(result$attributes)

cat("\nCoefficients:\n")
print(result$coeffs)
