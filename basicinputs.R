libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr","reshape2","fossil","tidyr","purrr")
library(data.table)
y <- lapply(libs, require, character.only = TRUE)
assert("All packages loaded correctly", all(unlist(y)))

source("scripts/ECC/classes_ecc.R")
source("cluster_sampling.R")

checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

option_list <- list(
  make_option(c("-a", "--source"), metavar = "file", default = "inputs/processed/source_data.tsv", help = "Source data"),
  make_option(c("-b", "--strains"), metavar = "file", default = "inputs/Strain_data.txt", help = "Strain data"),
  make_option(c("-c", "--tp1"), metavar = "file", default = "inputs/processed/tp1_clusters.txt", help = "TP1 cluster assignments"), 
  make_option(c("-d", "--tp2"), metavar = "file", default = "inputs/processed/tp2_clusters.txt", help = "TP2 cluster assignments"), 
  make_option(c("-x", "--heights"), metavar = "character", default = "0", 
              help = "Comma-delimited string of heights to collect ECCs for"), 
  make_option(c("-p", "--cpus"), metavar = "numeric", default = 1, help = "CPUs"),
  make_option(c("-t", "--trio"), metavar = "character", default = "010", 
              help = "source, temporal, geographic coefficients"))

params <- parse_args(OptionParser(option_list=option_list))

combos <- params$trio %>% strsplit(., "-") %>% unlist()
z <- vector("list", length = length(combos)) %>% set_names(combos)

hx <- params$heights %>% strsplit(split = ",") %>% unlist() %>% tibble(h = ., th = paste0("T", .))

tp1 <- Timepoint$new(params$tp1, "tp1")$Process(hx)$listHeights(hx)
tp2 <- Timepoint$new(params$tp2, "tp2")$Process(hx)$listHeights(hx)

c1 <- strsplit(combos[1], split = "") %>% unlist() %>% as.numeric() %>% 
  as.list() %>% set_names(c("sigma", "tau", "gamma"))

tau <- c1$tau
gamma <- c1$gamma
sigma <- c1$sigma
source_file <- params$source