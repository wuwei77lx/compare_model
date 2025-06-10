# Predict peak-gene results and truth peak-gene results: A data frame with n rows and 3 variables.
# Column1: {gene}; column2: {peak}; column3: {score}.

# AUPRC
library(PRROC) # Load the required libraries
#' Compute AUPRC given predicted and ground truth peak-gene scores
#' @param predicted Data frame with columns: gene, peak, score (model prediction score)
#' @param ground_truth Data frame with columns: gene, peak, score (binary 0/1 true labels)
#' @return A list containing PR curve and AUPRC
AUP <- function(predicted, ground_truth) {
  # Merge on gene-peak pairs
  data <- merge(ground_truth, predicted, by = c("gene", "peak"), all.x = TRUE)
  
  # Replace missing prediction scores with 0
  data$score.y[is.na(data$score.y)] <- 0
  
  # PRROC expects: scores.class0 = predicted scores; weights.class0 = 1 for positives, 0 for negatives
  pr_curve_data <- PRROC::pr.curve(
    scores.class0 = abs(as.numeric(data$score.y)),
    weights.class0 = as.numeric(data$score.x),
    curve = TRUE
  )
  
  return(pr_curve_data)
}

# AUPRC ratio
# pr <- AUP(results, truth)
# ratio <- pr[["auc.integral"]] / mean(truth$score)


# Early precision
#' @param predicted A data frame containing predicted peak-gene interactions with columns: gene, peak, score
#' @param ground_truth A data frame containing ground truth interactions with binary score column: 1 (positive), 0 (negative)
#' @return Early precision score
EPR <- function(predicted, ground_truth) {
  # k = number of true peak-gene links
  k <- sum(ground_truth$score == 1)
  
  # Merge on gene and peak
  merged <- merge(ground_truth, predicted, by = c("gene", "peak"), all.x = TRUE)
  merged$score.y <- as.numeric(merged$score.y)
  merged$score.y[is.na(merged$score.y)] <- 0
  
  # Order by predicted score descending
  merged <- merged[order(abs(merged$score.y), decreasing = TRUE), ]
  
  # Select top-k predictions
  if (sum(merged$score.y != 0) >= k) {
    top_k <- merged[1:k, ]
  } else {
    message("Warning: Fewer than k predicted interactions available, using all non-zero predictions.")
    top_k <- merged[merged$score.y != 0, ]
  }
  
  # Compute early precision
  early_precision <- sum(top_k$score.x == 1) / nrow(top_k)
  return(early_precision)
}

# Early precision ratio
# ep <- EPR(results, truth)
# ratio <- ep / mean(truth$score)
