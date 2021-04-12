pacman::p_load(tidyverse)

data <- read_csv("data/metabolic_MAG_data.csv")
# metadata ----------------------------------------------------------------

taxonomy <- read_csv("../data/taxonomy.csv")


metabolic_MAGS <- data %>% select(site, genome) %>% distinct()

#these gene counts come from Caitlin's METABOLIC run from the amino acid files
gene_counts <- read_delim("data/geneCounts_19Feb2021.txt", delim = "\t", col_names = F) %>%
  separate(X1, c("id", "gene_count"), " ") %>%
  separate(id, c("site", "genome"), "_") %>%
  mutate(site = str_remove(site, "eMMO"),
         gene_count = as.numeric(gene_count)) %>%
  inner_join(metabolic_MAGS)


metadata  <- read_csv("data/checkm.csv") %>%
  mutate(genome = as.character(genome)) %>%
  right_join(gene_counts) %>%
  inner_join(taxonomy)

write_csv(metadata, "data/metadata.csv")

site_color_dict <- c("#1520d2", "#39d0ec", "#48d016", "#eec314", "#ec7f39", "#f71027", "#000000", "#f5f2f2" )
names(site_color_dict) <- c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC")

#completeness bar plot
completeness_plot <- metadata %>%
  mutate(bins = cut_width(Completeness, 10, boundary = 0)) %>%
  group_by(site, bins) %>% 
  summarize(n=n()) %>%
  mutate(bins = str_replace(bins, ",", "TO"),
         bins = str_replace_all(bins, "[[:punct:]]", " "),
         bins = str_replace(bins, "TO", "-")) %>%
  ggplot(aes(bins, n, fill = site)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = site_color_dict) +
  coord_flip() +
  guides(fill=guide_legend(title="Site")) +
  labs(x="% Completeness") +
  labs(y="# of MAGs") 

