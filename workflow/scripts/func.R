pick_cells_MB <- function(log_unique_reads, frip,min_threshold = 100,threshold_reads = 1000,threshold_frip = 0.1) {
  library(mclust)
  library(dplyr)

  d <- as.data.frame(cbind(log_unique_reads,frip))
  d_small <- subset(d, log_unique_reads > log10(min_threshold))

  d.model <- Mclust(d_small$log_unique_reads, 1:9)
  d_small$class <- d.model$classification

  # Now trnasfer class labels from d_small to d
  d$class <- NA
  d$class[d$log_unique_reads > log10(min_threshold)] <- d_small$class

  cluster_stats <- d %>%
    group_by(class) %>%
    summarise(
      mean_log_unique_reads = mean(log_unique_reads),
      mean_frip = mean(frip),
      .groups = "drop"
    )
  
  clusters_pass <- cluster_stats$class[cluster_stats$mean_log_unique_reads > log10(threshold_reads) & cluster_stats$mean_frip > threshold_frip]

  d$pass_model <- FALSE
  d$pass_model[d$class %in% clusters_pass] <- TRUE
  return(d)
}
