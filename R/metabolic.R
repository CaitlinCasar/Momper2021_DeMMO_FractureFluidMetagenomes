pacman::p_load(readxl, tidyverse, plotly, heatmaply, grid)

#load old metadata from JGI annotations for gene counts per genome - need to update this
gene_counts_old <- read_csv("../submission/data/genome_metadata.csv") %>%
  select(`Genome Name / Sample Name`, `Gene Count   * assembled`) %>%
  rename(genome_sample = `Genome Name / Sample Name`) %>%
  rename(gene_count_old = `Gene Count   * assembled`) %>%
  mutate(site = if_else(str_detect(genome_sample, regex("DeMMO", ignore.case = TRUE)),
                        paste0("D", str_extract(genome_sample, "(?<=DeMMO)(.*)(?=_)")),
                        if_else(str_detect(genome_sample, regex("White_creek", ignore.case = TRUE)), "White Creek", "Service Water")),
         genome = str_extract(genome_sample, "(?=[^_]+$)(\\d.*)")) %>%
  select(-genome_sample)

gene_counts <- read_delim("data/geneCounts.txt", delim = " ", col_names = F) %>%
  rename(gene_count = "X2") %>%
  separate(X1, c("site", "genome"), "_") %>%
  mutate(site = str_remove(site, "eMMO"),
         gene_count = as.numeric(gene_count))

#load old metadata - this is not from metabolic output, need total # genes per genome. 
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
files <- list.files("../data/metabolic/annotations/", full.names = T)

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



# bubble plot -------------------------------------------------------------

#metabolism <- c("Nitrogen cycling", "Sulfur cycling", "As cycling", "Hydrogenases", "Methane metabolism", "Carbon fixation", "Chlorite reduction", "Selenate reduction", "Metal reduction", "Perchlorate reduction")
metabolism_color_dict <- read_csv("../data/metabolism_color_dict.csv")
element_cycling <- metabolism_color_dict %>%
  filter(group == "Element Cycling") %>%
  distinct()
element_cycling_colors <- element_cycling$color
names(element_cycling_colors) <- element_cycling$Lump

bubble_plot <- data %>%
  left_join(metadata %>% dplyr::select(site, genome, gene_count)) %>%
  filter(!is.na(gene_count)) %>%
  filter(hits > 0) %>%
  group_by(site) %>%
  mutate(gene_count = sum(na.omit(gene_count)),
         hits = (hits/gene_count)*100) %>%
  group_by(site, Function, Category, `Gene abbreviation`) %>%
  summarise(hits = sum(hits)) %>%
  inner_join(element_cycling) %>%
  ungroup() %>%
  mutate(Category = factor(Category, levels = c("Nitrogen cycling", "Urea utilization","Sulfur cycling", "Hydrogenases","Oxidative phosphorylation",                               
                                                   "Oxygen metabolism (Oxidative phosphorylation Complex IV)",                                         
                                                   "Halogenated compound utilization","Perchlorate reduction", "Chlorite reduction","As cycling",                                           
                                                   "Selenate reduction","Metal reduction"))) %>%
  ggplot(aes(site, `Gene abbreviation`, color = Lump, label=Category)) +
  geom_point(ggplot2::aes(size = hits)) +
  scale_x_discrete(position = "top") +
  scale_color_manual(values = element_cycling_colors, name = "Metabolism Category") +
  theme_bw() +
  guides(col = guide_legend(ncol = 1)) +
  #guides(fill=guide_legend(title="Category")) +
  theme(axis.title.x=ggplot2::element_blank(), 
        axis.title.y=ggplot2::element_blank(),
        #legend.position = "none",
        #strip.background = element_blank(), 
        panel.spacing = unit(0,"line"), 
        panel.border = element_rect(size = 0.25, color = "black"),
        strip.text.y = element_text(angle = 180, size=8, lineheight=1))  +
  facet_grid(rows=vars(Category), switch = "y", scales = "free", space = "free_y",
             labeller = labeller(Category = label_wrap_gen(10)))
  
plotly::ggplotly(bubble_plot)



gt = ggplot_gtable(ggplot_build(bubble_plot))

#show layout of plot to figure out which rows to manually resize
gtable::gtable_show_layout(gt)

for(i in c(9, 19, 21, 23, 25, 27, 29)){
  gt$heights[i] = 1.5*gt$heights[i]
}
grid.draw(gt)


n_genomes <- data %>%
  select(site, genome) %>%
  distinct() %>%
  group_by(site) %>%
  summarise(n_genomes = n())

n_genomes_meta <- metadata %>%
  select(site, genome) %>%
  distinct() %>%
  group_by(site) %>%
  mutate(n_genomes_meta = n()) %>%
  full_join(n_genomes)


