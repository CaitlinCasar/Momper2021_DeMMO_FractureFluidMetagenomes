pacman::p_load(readxl, tidyverse, plotly, heatmaply)

#load gene counts produced from .faa files in metabolic output
gene_counts <- read_delim("data/geneCounts.txt", delim = " ", col_names = F) %>%
  rename(gene_count = "X2") %>%
  separate(X1, c("site", "genome"), "_") %>%
  mutate(site = str_remove(site, "eMMO"),
         gene_count = as.numeric(gene_count))

#load metadata
path <- "../submission/data/DeMMO_genome_master.xlsx"
metadata  <- path %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = path) %>%
  reduce(bind_rows) %>%
  select(user_genome, Completeness, Contamination) %>%
  separate(user_genome, c("site", "genome"), sep = "_") %>%
  mutate(site = str_remove(site, "eMMO")) %>%
  left_join(gene_counts) %>%
  left_join(gene_counts_old) %>%
  mutate(gene_diff = `gene_count`/`gene_count_old`*100)

#load GTDBtk taxonomy from metabolic output
files <- list.files("../data/metabolic/taxonomy/", full.names = T)

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


# load metabolic annotations 
files <- list.files("../submission/data/metabolic/annotations", full.names = T)

read_files <- function(file){
  site <- str_extract(file, "(?<=data/metabolic/)(.*)(?=_METABOLIC_result[.]xlsx)")
  map <- file %>%
    read_excel(sheet = "map") %>%
    rename(genome_id = "...2") %>%
    separate(`Genome ID map`, c("site", "genome"), "_") %>%
    mutate(site = str_remove(site, "eMMO"))
  
  presence <- file %>%
    read_excel(sheet = "HMMHitNum") %>%
    select(-contains("Hits"), -contains("Hit numbers")) %>%
    gather(genome_id, presence, contains("Hmm presence")) %>%
    mutate(genome_id = str_remove(genome_id, " Hmm presence")) 
  
  file %>%
    read_excel(sheet = "HMMHitNum") %>%
    select(-contains("Hits"), -contains("Hmm presence")) %>%
    gather(genome_id, hits, contains("Hit numbers")) %>%
    mutate(genome_id = str_remove(genome_id, " Hit numbers")) %>%
    left_join(presence) %>%
    full_join(map) 
}

file_list = lapply(files, read_files)
data <- reduce(file_list, bind_rows) 

carbon_fixation <- data %>%
  filter(Category == "Carbon fixation") %>%
  mutate(presence = if_else(presence == "Present", 1, 0)) %>%
  group_by(site, Function) %>%
  summarise(presence = sum(presence)) %>%
  left_join(n_genomes) %>%
  mutate(presence = (presence/n_genomes)*100) %>%
  ggplot() +
  geom_tile(aes(x=site, y=Function, fill = presence)) +
  scale_fill_viridis()


plotly::ggplotly(carbon_fixation)
