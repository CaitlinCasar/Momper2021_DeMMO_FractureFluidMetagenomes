pacman::p_load(tidyverse)

metabolic_MAGS <- data %>% select(site, genome) %>% distinct()

#these gene counts come from Caitlin's METABOLIC run from the amino acid files
gene_counts <- read_delim("data/geneCounts_19Feb2021.txt", delim = "\t", col_names = F) %>%
  separate(X1, c("id", "gene_count"), " ") %>%
  separate(id, c("site", "genome"), "_") %>%
  mutate(site = str_remove(site, "eMMO"),
         gene_count = as.numeric(gene_count)) %>%
  inner_join(metabolic_MAGS)

MAG_stats_plot <- gene_counts %>% 
  mutate(site = factor(site, levels = rev(c("D1", "D2", "D3", "D4", "D5", "D6", "SW", "WC")))) %>%
  ggplot(aes(site, gene_count)) + 
  geom_violin(fill = "gray", draw_quantiles = c(0.25, 0.5, 0.75)) +
  stat_summary(fun.y=mean, geom="point", size=1, color="white") +
  coord_flip() +
  xlab("Gene Count") +
  theme(axis.title.y = element_blank())

MAG_stats <- gene_counts %>% 
  group_by(site) %>%
  summarize(mean = mean(gene_count), median =  median(gene_count), max = max(gene_count),
            min = min(gene_count),
            quantile_25 = quantile(gene_count, 0.25), quantile_50 = quantile(gene_count, 0.5),
            quantile_75 = quantile(gene_count, 0.75))

write_csv(MAG_stats, "data/mag_stats.csv")
