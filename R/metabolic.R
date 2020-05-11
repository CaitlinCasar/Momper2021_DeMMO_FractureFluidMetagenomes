pacman::p_load(readxl, tidyverse, plotly)


metadata <- read_csv("../data/metadata.csv") %>%
  mutate(site = paste0("D", str_extract(`Genome Name / Sample Name`, "(?<=DeMMO)(.*)(?=_)")),
         genome = str_extract(`Genome Name / Sample Name`, "(?<=DeMMO\\d_)(.*)")) %>%
  filter(!is.na(genome))


# load metabolic files
files <- list.files("../data/metabolic", full.names = T)

read_files <- function(file){
  site <- str_extract(file, "(?<=data/metabolic/)(.*)(?=_METABOLIC_result[.]xlsx)")
  map <- file %>%
    read_excel(sheet = "map") %>%
    rename(genome_id = "...2") %>%
    separate(`Genome ID map`, c("site", "genome"), "_") %>%
    mutate(site = str_remove(site, "eMMO"))
    
  file %>%
    read_excel(sheet = "HMMHitNum") %>%
    select(-contains("Hmm presence"), -contains("Hits")) %>%
    rename_at(vars(matches(".*Hit numbers")), ~ str_remove(., " Hit numbers")) %>%
    gather(genome_id, hits, -c(Category,Function,`Gene abbreviation`,`Gene name`,`Hmm file`,`Corresponding KO`, Reaction,Substrate,Product)) %>%
    full_join(map) 
}

file_list = lapply(files, read_files)
data <- reduce(file_list, bind_rows) 


metabolism <- c("Nitrogen cycling", "Sulfur cycling", "As cycling", "Hydrogenases", "Methane metabolism", "Carbon fixation", "Chlorite reduction", "Selenate reduction", "Metal reduction", "Perchlorate reduction")

bubble_plot <- data %>%
  left_join(metadata %>% dplyr::select(site, genome, `Gene Count   * assembled`)) %>%
  filter(!is.na(`Gene Count   * assembled`)) %>%
  filter(hits > 0) %>%
  group_by(site) %>%
  mutate(gene_count = sum(`Gene Count   * assembled`),
         hits = (hits/gene_count)*100) %>%
  group_by(site, Function, Category, `Gene abbreviation`) %>%
  summarise(hits = sum(hits)) %>%
  filter(Category %in% metabolism) %>%
  ggplot(aes(site, `Gene abbreviation`, color = Function, label=Category)) +
  geom_point(ggplot2::aes(size = hits)) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(axis.title.x=ggplot2::element_blank(), 
        axis.title.y=ggplot2::element_blank(),
        legend.position = "none",
        #strip.background = element_blank(), 
        panel.spacing = unit(0,"line"), 
        panel.border = element_rect(size = 0.25, color = "black"))  +
  facet_grid(rows=vars(Category), scales = "free")
  #scale_color_manual(values=metabolism_color_dict) +
  
plotly::ggplotly(bubble_plot)

