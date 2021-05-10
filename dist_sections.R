### Incorporating the allele data with the epidemiological data 
actualResults <- function(source_pw, strain_data, sigma, tau, gamma, typing_datum) {
  
  # look into a method that runs the distances when given two separate vectors (then we can 
  # avoid the extra work that comes with overlapping)
  temp_pw <- generateDistances(strain_data, "temp", "Date", "Temp.Dist")
  geog_pw <- generateDistances(strain_data, "geo", c("Latitude", "Longitude"), "Geog.Dist")
  
  geog_temp_tr <- left_join(geog_pw$transformed, temp_pw$transformed, by = c("Strain.1", "Strain.2"))
  
  ## This generates a table of your comparisons based on your epidemiological data (source, time, geographical) 
  ## with the assigned weights of s, t and g and then computes the similarity/distance and generates a matrix
  epi.table <- EpiTable(strain_data, source_pw, sigma, tau, gamma, geog_temp_tr)
  # saveRDS(epi.table, "epitable.Rds")
  epi.matrix <- EpiMatrix(epi.table)
  
  # typing_datum <- typing_data[[2]]
  
  # ### Section 3: Incorporating the allele data with the epidemiological data - typing_data
  # # Calculate ECC in parallel; this may not work on Windows, but should work out of the box on Linux and OSX
  # eccs <- lapply(typing_data, function(typing_datum) {
  g_cuts <- typing_datum %>% rownames_to_column("genome") %>% as_tibble()
  # Method with singletons merged separately, for speed/efficiency:
  eccs <- epi_cohesion_sep(g_cuts, epi.matrix, cpus) %>% as.data.table()
  # beep(3)
  # Original method: epi_cohesion_calc(g_cuts, epi.matrix, cpus = cpus)
  # })
  
  return(eccs)
}

# ---------------------------------------------------------------------------------------------------
# Note: dr stands for data representative
testResults <- function(strain_data, tau, gamma, typing_datum) {

  # in example: strain_data has 35,627 rows (strains), assignments has 5,504 rows (> 6-fold smaller)
  assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
    unique() %>% rownames_to_column("dr")
  
  # Temporal distances - all possible date pair distances
  dm_temp <- assignments %>% select(dr, Date) %>% pairwiseDists(., "temp", "Date", c("dr1", "dr2", "Temp.Dist"))
  
  # Geographical distances - all possible lat-long pair distances
  dm_geo <- assignments %>% select(dr, Latitude, Longitude) %>% 
    pairwiseDists(., "geo", c("Latitude", "Longitude"), c("dr1", "dr2", "Geog.Dist"))
  
  epi_table <- merge.data.table(dm_temp, dm_geo) %>% 
    mutate(Total.Dist = sqrt( ((Temp.Dist^2)*tau) + ((Geog.Dist^2)*gamma) )) %>% 
    select(dr1, dr2, Total.Dist) %>% as_tibble()
  
  epi_matrix <- dcast(epi_table, formula = dr1 ~ dr2, value.var = "Total.Dist")
  epi_matrix <- as.matrix(epi_matrix[,2:ncol(epi_matrix)]) 
  rownames(epi_matrix) <- colnames(epi_matrix)

  # create similarity values from epi distance matrix:
  epi_melt <- melt(as.matrix(1-epi_matrix)) %>%
    mutate(across(c(Var1, Var2), as.character)) %>% as.data.table()
  
  # Identifying which strains match with which non-redundant data representatives
  dr_matches <- strain_data %>% 
    left_join(., assignments, by = c("Latitude", "Longitude", "Date", "Location")) %>% 
    select(Strain, dr)
  
  dr_td1 <- typing_datum %>% rownames_to_column("Strain") %>% as_tibble() %>% 
    left_join(., dr_matches, by = "Strain") %>% 
    mutate(across(dr, as.character)) %>% select(-Strain)
  
  # Counting data representatives (so we know how much to multiply each ECC value by to represent all strains)
  tallied_reps <- dr_td1 %>% group_by(T0) %>% count(dr) %>% ungroup()
  g_cuts <- left_join(dr_td1, tallied_reps, by = intersect(colnames(tallied_reps), colnames(dr_td1))) %>% 
    unique() %>% mutate(across(dr, as.character))
  
  ecc_results <- epi_cohesion_new(g_cuts, epi_melt)
  return(ecc_results)
}
