pacman::p_load(readxl, tidyverse, plotly, heatmaply, grid)

#load old metadata from JGI annotations for gene counts per genome - need to update this
gene_counts_old <- read_csv("../../submission/data/genome_metadata.csv") %>%
  select(`Genome Name / Sample Name`, `Gene Count   * assembled`, Phylum) %>%
  rename(genome_sample = `Genome Name / Sample Name`) %>%
  rename(gene_count_old = `Gene Count   * assembled`) %>%
  mutate(site = if_else(str_detect(genome_sample, regex("DeMMO", ignore.case = TRUE)),
                        paste0("D", str_extract(genome_sample, "(?<=DeMMO)(.*)(?=_)")),
                        if_else(str_detect(genome_sample, regex("White_creek", ignore.case = TRUE)), "WC", "SW")),
         genome = str_extract(genome_sample, "(?=[^_]+$)(\\d.*)")) %>%
  select(-genome_sample)

#missing ORFs for white creek and service water 
gene_counts <- read_delim("data/geneCounts.txt", delim = " ", col_names = F) %>%
  rename(gene_count = "X2") %>%
  separate(X1, c("site", "genome"), "_") %>%
  mutate(site = str_remove(site, "eMMO"),
         gene_count = as.numeric(gene_count))

#load old metadata - this is not from metabolic output, need total # genes per genome. 
path <- "../data/DeMMO_genome_master.xlsx"

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


# load metabolic annotations - control files formatted differently, no map sheet
files <- list.files("../data/metabolic/annotations", full.names = T)

read_files <- function(file){
  site <- str_extract(file, "(?<=data/metabolic/annotations/)(.*)(?=_METABOLIC_result[.]xlsx)")
  message(site)
  if(site %in% c("D6", "SW", "WC")){
    file %>%
      read_excel(sheet = "HMMHitNum") %>%
      select(-contains("Hits"), -contains("Hmm.presence")) %>%
      gather(genome_id, hits, contains("Hit.numbers")) %>%
      mutate(site = site, 
             genome = str_extract(genome_id, "(?<=_)(.*)(?=[.]Hit)"),
             presence = if_else(hits > 0, "Present", "Absent")) %>%
      rename_at(vars(contains(".")),funs(gsub("\\.", " ", .)))
  }else{
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
  
}

file_list = lapply(files, read_files)
data <- reduce(file_list, bind_rows) 



# bubble plot -------------------------------------------------------------

#metabolism <- c("Nitrogen cycling", "Sulfur cycling", "As cycling", "Hydrogenases", "Methane metabolism", "Carbon fixation", "Chlorite reduction", "Selenate reduction", "Metal reduction", "Perchlorate reduction")
metabolism_color_dict <- read_csv("../data/metabolism_color_dict.csv")
element_cycling <- metabolism_color_dict %>%
  filter(group == "Element Cycling" | Category %in% c("Fermentation", "C1 metabolism")) %>%
  distinct() %>%
  filter(!Category %in% c("As cycling", "Oxidative phosphorylation")) 
element_cycling_colors <- element_cycling$color
names(element_cycling_colors) <- element_cycling$Lump

bubble_plot <- data %>%
  filter(!`Gene abbreviation` %in% c("E3.8.1.2", "ccoN", "ccoO", "ccoP", "nifK", "nifH", "nifD", "octR", "cydA", "cydB", "mauB"))%>%
  mutate(`Gene abbreviation` = str_remove(`Gene abbreviation`, 'group-')) %>%
  left_join(metadata %>% dplyr::select(site, genome, gene_count_old, gene_count)) %>%
  mutate(gene_count_old = if_else(is.na(gene_count_old), gene_count, gene_count_old)) %>%
  select(-gene_count) %>%
  rename(gene_count = "gene_count_old") %>%
  filter(hits > 0) %>%
  filter(!is.na(gene_count)) %>%
  group_by(site) %>%
  mutate(gene_count = sum(na.omit(gene_count)),
         hits = (hits/gene_count)*100) %>%
  group_by(site, Function, Category, `Gene abbreviation`) %>%
  summarise(hits = sum(hits)) %>%
  inner_join(element_cycling) %>%
  ungroup() %>%
  mutate(Category = if_else(Category == "Metal reduction", "Fe/Mn reduction", Category),
         Category = if_else(Category == "Urea utilization", "Nitrogen cycling", Category),
         Category = if_else(Category == "Halogenated compound utilization", "Halogen cycling", Category),
         Category = if_else(Category == "Perchlorate reduction", "Halogen cycling", Category),
         Category = if_else(Category == "Chlorite reduction", "Halogen cycling", Category),
         Category = factor(Category, levels = c("Nitrogen cycling", "Sulfur cycling", "Hydrogenases",                               
                                                   "Oxygen metabolism (Oxidative phosphorylation Complex IV)",                                         
                                                  "Halogen cycling",                                           
                                                   "Selenate reduction","Fe/Mn reduction", "C1 metabolism", "Fermentation")),
         Lump = as.factor(Lump),
         `Gene abbreviation` = factor(`Gene abbreviation`, levels = unique(`Gene abbreviation`[order(Lump)]))) %>%
  ggplot(aes(site, `Gene abbreviation`, color = Lump, label=Category)) +
  geom_point(ggplot2::aes(size = hits)) +
  scale_size_continuous(breaks = c(2e-05, 5e-05, 1e-04, 5e-04, 1e-03), name = "% metagenome") +
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
#gtable::gtable_show_layout(gt)

for(i in c(28)){
  gt$heights[i] = 2*gt$heights[i]
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




# c-fixation heat map --------------------------------------------------

c_fixation <- metabolism_color_dict %>%
  filter(Category == "Carbon fixation") %>%
  distinct()

c_fixation_colors <- c_fixation$color
names(c_fixation_colors) <- c_fixation$Lump

c_fix_pathways <- read_csv("../data/c_fixation_pathways.csv") 

c_fix_n_genes <- c_fix_pathways %>%
  select(Category, Function) %>%
  group_by(Category, Function) %>%
  summarise(n_pathway_genes = n())

c_fix_plot <- data %>%
  filter(hits > 0 & Category == "Carbon fixation") %>%
  group_by(site, Category, Function, `Gene abbreviation`) %>%
  summarize(hits = sum(hits)) %>%
  left_join(metadata %>% group_by(site) %>% summarise(n_genomes = n())) %>%
  left_join(c_fixation) %>%
  mutate(hits_per_genome = hits/n_genomes) %>%
  ggplot(aes(site, `Gene abbreviation`, color = Lump, label=Category)) +
  geom_point(ggplot2::aes(size = hits_per_genome)) +
  #scale_size_continuous(breaks = c(2e-05, 5e-05, 1e-04, 5e-04, 1e-03), name = "% metagenome") +
  scale_x_discrete(position = "top") +
  scale_color_manual(values = c_fixation_colors, name = "Pathway") +
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
  facet_grid(rows=vars(Function), switch = "y", scales = "free", space = "free_y",
             labeller = labeller(Function = label_wrap_gen(10)))

  

  metabolism_data <- data %>%
    group_by(Category, Function, site, genome) %>%
    summarise(hits = sum(hits)) %>%
    filter(hits > 0)
