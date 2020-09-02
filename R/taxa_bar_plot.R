pacman::p_load(tidyverse, cowplot)

taxa_color_dict <- read_csv("../data/distinct_taxonomies.csv") 

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

genome_abundance <- taxonomy %>%
  select(site, phylum, class) %>%
  group_by_all() %>%
  summarise(`Number of Genomes` = n()) %>%
  mutate(site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC")))) %>%
  left_join(taxa_color_dict) %>%
  mutate(phylum = if_else(!is.na(new_taxonomy), new_taxonomy, phylum),
         class = if_else(!is.na(new_taxonomy), new_taxonomy, class),
         color_subplot = if_else(is.na(color_subplot), color, color_subplot),
         color = if_else(is.na(color), "#000000", color))

#calculate genome abundance at phylum level
phylum_abundance <- genome_abundance %>%
  group_by(site, phylum, color) %>%
  summarize(`Number of Genomes` = sum(`Number of Genomes`))

#create color palette for phyla
phylum_colors <- phylum_abundance %>%
  ungroup() %>%
  select(phylum, color) %>%
  distinct() %>%
  select(color) %>%
  pull()
names(phylum_colors) <- phylum_abundance %>%
  ungroup() %>%
  select(phylum, color) %>%
  distinct() %>%
  select(phylum) %>%
  pull()

phylum_abundance_plot <- phylum_abundance %>%
  ggplot(aes(site, `Number of Genomes`, fill=phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) + 
  coord_flip() +
  ggtitle("Phylum-level Taxonomy")

class_colors <- genome_abundance %>%
  ungroup() %>%
  select(class, color) %>%
  distinct() %>%
  select(color) %>%
  pull()

names(class_colors) <- genome_abundance %>%
  ungroup() %>%
  select(class, color) %>%
  distinct() %>%
  select(class) %>%
  pull()

subplot_a <- genome_abundance %>%
  filter(subplot == "A") %>%
  ggplot(aes(site, `Number of Genomes`, fill = phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) + 
  coord_flip() +
  theme(legend.position = "none") +
  ggtitle("Archaea")

subplot_p <- genome_abundance %>%
  filter(subplot == "P") %>%
  ggplot(aes(site, `Number of Genomes`, fill = phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) + 
  coord_flip() +
  theme(legend.position = "none") +
  ggtitle("Proteobacteria")

subplot_c <- genome_abundance %>%
  filter(subplot == "C") %>%
  ggplot(aes(site, `Number of Genomes`, fill = phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) + 
  coord_flip() + 
  theme(legend.position = "none") +
  ggtitle("Candidate Phyla Radiation")

bottom_row <- plot_grid(subplot_a, subplot_p, subplot_c, labels = c('B.', 'C.' ,'D.'), label_size = 12, nrow=1)

plot_grid(phylum_abundance_plot, bottom_row, labels = c('A.', ''), label_size = 12, ncol = 1, rel_heights = c(2, 1))
  
  