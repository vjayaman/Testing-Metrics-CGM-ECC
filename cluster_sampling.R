

tp2ClusterSizes <- function(td) {
  t1clusters <- td[[1]] %>% rownames_to_column("Strain") %>% as_tibble() %>% set_colnames(c("Strain", "Threshold"))
  t2clusters <- td[[2]] %>% rownames_to_column("Strain") %>% as_tibble() %>% set_colnames(c("Strain", "Threshold"))
  
  t2clusters$Threshold %>% table() %>% 
    as.data.frame() %>% as_tibble() %>% 
    set_colnames(c("Threshold", "Freq")) %>% 
    mutate(across(Threshold, as.character)) %>% 
    mutate(across(Threshold, as.integer)) %>% 
    return()
}

datasetSizeBelow <- function(td, maxsize) {
  currsize <- 0
  i <- 0
  y <- c()
  
  t2frequencies <- tp2ClusterSizes(td) %>% filter(Freq < maxsize)
  # t2frequencies %>% slice_sample(.data = ., n = 10)
  clusters <- sample(t2frequencies$Threshold, nrow(t2frequencies), replace = FALSE)

  while (currsize < maxsize & i < length(clusters)) {
    i <- i + 1
    newsize <- t2frequencies %>% filter(Threshold == clusters[i]) %>% pull(Freq)
    newsum <- currsize + newsize
    
    if (newsum < maxsize) {
      y <- append(y, clusters[i])
      currsize <- newsum
    }
  }
  
  return(y)
}
