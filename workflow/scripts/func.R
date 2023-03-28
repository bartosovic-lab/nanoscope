pick_cells_MB <- function(log_unique_reads, frip){
  library(mclust)
  library(dplyr)

  d <- as.data.frame(cbind(log_unique_reads,frip))
  d.model <- Mclust(d,3)
  d$class <- d.model$classification
  top.cluster <- d %>% group_by(class) %>% summarise(mean(log_unique_reads)) %>% top_n(1) %>% dplyr::select(class)
  top.cluster <- as.character(top.cluster)

  d$pass_model <- FALSE
  d$pass_model[d$class == top.cluster] <- TRUE
  return(d)
}
