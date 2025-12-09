# Calculate fraction overlap
fraction_overlap <- function(copies_query, copies_target) {
  # Normalize the target OG copies by median copies in the outgroup
  shared_species <- colnames(copies_query)[which(copies_query > 0 & copies_target > 0)]
  if (length(shared_species) == 0) {
    return(0)
  }
  copies_target_median <- median(as.numeric(copies_target[,copies_target > 0]))
  if (copies_target_median > 0) {
    copies_target <- copies_target - copies_target_median
    copies_target[copies_target == 0] <- 1
    copies_target[copies_target < 0] <- 0
  }
  
  # Calculate overlap coefficient
  intersect_copies_query_and_target <- sum(copies_target[,copies_query > 0] > 0)
  overlap_coef <- intersect_copies_query_and_target / min(sum(copies_query > 0), sum(copies_target > 0))
  
  return(overlap_coef)
}
