pacman::p_load(readxl, tidyverse, plotly, heatmaply, grid)

files <- list.files("../data/metabolic/taxonomy", full.names = T)

read_files <- function(file){
  file %>% read_delim(delim = "\t")
}

taxonomy <- lapply(files, read_files) %>%
  reduce(bind_rows) %>%
  select(user_genome, classification) %>%
  separate(user_genome, c("site", "genome"), sep = "_") %>%
  separate(classification, c("domain","phylum", "class", "order", "family", "genus", "species"), sep = ";p__|;c__|;o__|;f__|;g__|;s__") %>%
  mutate(domain = str_remove(domain, "d__"),
         site = str_remove(site, "eMMO"))

selected_archaea <- read_csv("../data/Archaea_to_plot.csv")

selected_cpr <- read_csv("../data/Patescibacteria_to_plot.csv")

genome_abundance <- taxonomy %>%
  select( -genome) %>%
  group_by_all() %>%
  summarise(abundance = n())

write_csv(genome_abundance, "../data/genome_abundance.csv")

taxonomy %>%
  select(phylum, class) %>%
  distinct() %>%
  write_csv("../data/distinct_taxonomies.csv")
