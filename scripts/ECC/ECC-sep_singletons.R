# Vasena's edits - for speed improvement on Windows, where we have less control over cores 

epi_cohesion_sep <- function(g_cuts, epi_matrix, cpus){
  
  genome_names <- g_cuts %>% select(genome) %>% pull()

  epi_melt <- melt(as.matrix(1-epi_matrix)) # create similarity values from epi distance matrix
  
  epi_melt_joined <- 
    expand_grid(genome_names, genome_names, .name_repair = function(x) {c("Var1", "Var2")}) %>% 
    left_join(epi_melt, by = c("Var1", "Var2")) %>% filter(!is.na(value))
  
  calculate_s1 <- function(k) {
    epi_melt_joined %>% 
      filter(Var1 %in% k & Var2 %in% k) %>%
      select(value) %>% 
      pull() %>% 
      sum()
  }

  cut_cluster_members <-
    g_cuts %>%
    pivot_longer(-genome, names_to = "cut", values_to = "cluster") %>%  
    group_by(cut, cluster) %>% 
    summarise(members = list(cur_data()$genome), .groups = "drop_last")
  
  sums <- cut_cluster_members %>% 
    mutate(
      s1 = map_dbl(members, calculate_s1),
      cluster_size = map_int(members, length),
      ECC = (s1 - cluster_size) / (cluster_size * (cluster_size - 1))
    ) %>% 
    ungroup()
  
  sums %>% 
    select(-cut, -members, -s1) %>% 
    set_colnames(c(names(g_cuts)[2], paste0(names(g_cuts)[2], "_Size"), 
                   paste0(names(g_cuts)[2], "_ECC")))
}
