pacman::p_load(readxl, tidyverse, plotly)

#load native rock otu data
otu_metadata <- read_csv("../../../DeMMO_NativeRock/DeMMO_NativeRock/data/metadata.csv") 
otu_table <- read_delim("../../../DeMMO_NativeRock/DeMMO_NativeRock/data/DeMMO_NativeRock_noChimera_otuTable_withTaxa_46822.txt", delim = '\t', comment = "# ")

otu_taxonomy <- otu_table %>%
  select(`#OTU ID`, taxonomy) %>%
  mutate(tax = gsub("Gammaproteobacteria; D_3__Betaproteobacteriales", "Betaproteobacteria; D_3__Betaproteobacteriales", taxonomy), #fix taxonomy for Beta's,
         taxonomy = str_remove_all(tax, "D_0__| D_1__| D_2__| D_3__| D_4__| D_5__| D_6__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family", "genus", "species"))

# abundance_table <- otu_table %>%
#   select(-taxonomy) %>%
#   mutate_at(vars(-`#OTU ID`), funs(./sum(.)*100)) %>% #normalize to relative abundance 
#   gather(sample_id, abundance, `9.DeMMO3.steri.11June2019`:`17.DeMMO1.tube8.12Dec2019`) %>%
#   left_join(otu_taxonomy) %>%
#   group_by(sample_id, tax) %>%
#   summarise(abundance = sum(abundance)) %>%
#   mutate(taxonomy = tax,
#          taxonomy = str_remove_all(tax, "D_0__| D_1__| D_2__| D_3__| D_4__| D_5__| D_6__")) %>%
#   separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family", "genus", "species")) %>%
#   filter(abundance >= 5 & !domain == "Unassigned") %>%
#   mutate(taxa = family,
#          #taxa = if_else(is.na(taxa) | str_detect(taxa, "uncultured"), family, taxa),
#          taxa = if_else(is.na(taxa) | str_detect(taxa, "uncultured"), order, taxa),
#          taxa = if_else(is.na(taxa) | str_detect(taxa, "uncultured"), class, taxa),
#          taxa = if_else(is.na(taxa) | str_detect(taxa, "uncultured"), phylum, taxa)) %>%
#   left_join(otu_metadata %>% select(sample_id, site, substrate)) %>%
#   mutate(substrate = if_else(substrate == "fluid", "fluid", "biofilm")) %>%
#   group_by(site, taxa, substrate) %>%
#   summarise(abundance = sum(abundance)) 
# 
# substrate <- abundance_table %>%
#   group_by(site, taxa) %>%
#   summarise(n = n()) %>%
#   left_join(abundance_table) %>%
#   mutate(substrate = if_else(n == 2, "generalist", substrate),
#          taxa = if_else(taxa == "Omnitrophicaeota", "Omnitrophota", taxa),
#          #taxa = if_else(taxa == "Sideroxydans", "Gallionellaceae", taxa),
#          taxa = tolower(taxa)) %>%
#   distinct()



abundance_table <- otu_table %>%
  select(-taxonomy) %>%
  mutate_at(vars(-`#OTU ID`), funs(./sum(.)*100)) %>% #normalize to relative abundance
  gather(sample_id, abundance, `9.DeMMO3.steri.11June2019`:`17.DeMMO1.tube8.12Dec2019`)


taxon_abundance <- function(level, name){
  otu_table %>%
    select(`#OTU ID`, taxonomy) %>%
    mutate(taxonomy = gsub("Gammaproteobacteria; D_3__Betaproteobacteriales", "Betaproteobacteria; D_3__Betaproteobacteriales", taxonomy),
           taxa = str_extract(taxonomy, level),
           taxa = if_else(is.na(taxa), taxonomy, taxa)) %>%
    right_join(abundance_table) %>%
    ungroup() %>%
    group_by(sample_id, taxa) %>%
    summarise(abundance = sum(abundance)) %>%
    inner_join(otu_metadata) %>% #add metadata
    group_by(sample_id, substrate, taxa,`Date Sampled`, experiment_type, site) %>% 
    summarise(abundance = sum(abundance))
}

family_level <- "(.*)(?=; D_5__)"
class_level <- "(.*)(?=; D_3__)"
phylum_level <- "(.*)(?=; D_2__)"

taxon_abundance_table <- taxon_abundance(family_level, "family") %>%
  group_by(taxa) %>%
  filter(max(abundance) >= 5) %>%
  ungroup() %>%
  mutate(taxonomy = str_remove_all(taxa, "D_0__| D_1__| D_2__| D_3__| D_4__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family")) %>%
  mutate(taxa = family, 
         taxa = if_else(is.na(taxa) | str_detect(taxa, "uncultured"), order, taxa),
         taxa = if_else(is.na(taxa) | str_detect(taxa, "uncultured"), class, taxa),
         taxa = if_else(is.na(taxa) | str_detect(taxa, "uncultured"), phylum, taxa),
         taxa = if_else(taxa == "GWF2-40-263", order, taxa),
         taxa = if_else(taxa == "mle1-8", phylum, taxa),
         substrate = if_else(substrate == "fluid", "fluid", "biofilm")) %>%
  select(site, substrate, taxa, abundance) 

substrate_n <- taxon_abundance_table %>%
  filter(abundance >= 5) %>%
  select(-abundance) %>%
  distinct() %>%
  group_by(site, taxa) %>%
  summarise(n = n()) 

substrate <- taxon_abundance_table %>%
  anti_join(substrate_n) %>%
  select(-abundance) %>%
  distinct() %>%
  group_by(site, taxa) %>%
  summarise(n = n()) %>%
  bind_rows(substrate_n) %>%
  distinct() %>%
  left_join(taxon_abundance_table %>%
              filter(abundance >= 5) %>%
              select(-abundance) %>%
              distinct()) %>%
  mutate(substrate = if_else(n == 2, "generalist", substrate),
         taxa = if_else(taxa == "Omnitrophicaeota", "Omnitrophota", taxa),
         #taxa = if_else(taxa == "Sideroxydans", "Gallionellaceae", taxa),
         taxa = tolower(taxa)) %>%
  distinct() %>%
  filter(!is.na(taxa))

#load old metadata from JGI annotations for gene counts per genome - need to update this
gene_counts <- read_csv("../data/genome_metadata.csv") %>%
  select(`Genome Name / Sample Name`, `Gene Count   * assembled`) %>%
  rename(genome_sample = `Genome Name / Sample Name`) %>%
  rename(gene_count = `Gene Count   * assembled`) %>%
  mutate(site = if_else(str_detect(genome_sample, regex("DeMMO", ignore.case = TRUE)),
                        paste0("D", str_extract(genome_sample, "(?<=DeMMO)(.*)(?=_)")),
                        if_else(str_detect(genome_sample, regex("White_creek", ignore.case = TRUE)), "White Creek", "Service Water")),
         genome = str_extract(genome_sample, "(?=[^_]+$)(\\d.*)")) %>%
  select(-genome_sample)


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
  left_join(gene_counts)

#load GTDBtk taxonomy from metabolic output
files <- list.files("../data/metabolic/taxonomy/", full.names = T)

read_files <- function(file){
  file %>% read_delim(delim = "\t")
}

taxonomy <- lapply(files, read_files) %>%
  reduce(bind_rows) %>%
  select(user_genome, classification) %>%
  mutate(taxonomy = classification) %>%
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

#taxa_of_interest <- paste(c("thermodesulfovibrionia", "sideroxydans", "desulfobulbia", "desulfobacteraceae", "thiobacillus", "pseudomonadaceae", "omnitrophota", "woesarchaea", 'gallionella'), collapse = "|")
taxa_of_interest <- paste(substrate$taxa, collapse = "|")

taxa_genes <- data %>%
  left_join(metadata) %>% 
  #filter(Completeness >= 50) %>%
  left_join(taxonomy) %>%
  mutate(taxonomy = tolower(taxonomy)) %>%
  filter(str_detect(taxonomy, taxa_of_interest) & site %in% c("D1", "D3") & hits > 0 ) %>%
  group_by(site, taxonomy) %>%
  #mutate(gene_count = sum(na.omit(gene_count)),
         #hits = (hits/gene_count)*100) %>%
  group_by(site, Function, Category, `Gene abbreviation`, taxonomy) %>%
  summarise(hits = sum(hits)) %>%
  mutate(classification = taxonomy) %>%
  #filter(Category %in% metabolism) %>%
  mutate(id = paste(site,str_extract(taxonomy, taxa_of_interest), sep = "_"),
         taxa = str_extract(taxonomy, taxa_of_interest)) %>%
  left_join(substrate) %>%
  filter(!is.na(substrate)) 

bubble_plot <- taxa_genes %>%
  ggplot(aes(`Gene abbreviation`, id, color = Function, label=Category)) +
  geom_point() +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(axis.title.x=ggplot2::element_blank(), 
        axis.title.y=ggplot2::element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        #legend.position = "none",
        #strip.background = element_blank(), 
        panel.spacing = unit(0,"line"), 
        panel.border = element_rect(size = 0.25, color = "black"))  +
  facet_grid(rows = vars(substrate), cols=vars(Category), scales = "free", labeller = label_wrap_gen(width = 3, multi_line = TRUE))
             #scale_color_manual(values=metabolism_color_dict) +

taxa_bubble_plot <- plotly::ggplotly(bubble_plot)


taxa_genes_hits <- taxa_genes %>%
  ungroup() %>%
  select(id, hits, `Gene abbreviation`) %>%
  group_by(id, `Gene abbreviation`) %>%
  summarise(hits = sum(hits)) %>%
  #mutate(hits = if_else(na.omit(hits) >0, 1, 0)) %>%
  pivot_wider(names_from = `Gene abbreviation`, values_from = hits, values_fill = list(hits=0)) %>%
  column_to_rownames("id")


taxa_genes_NMDS <- taxa_genes_hits %>%
  vegan::metaMDS(k=2)

gene_vectors <-vegan::envfit(taxa_genes_NMDS, taxa_genes_hits, perm=1000)

NMDS_coords <- taxa_genes_NMDS[["points"]] %>%
  as_tibble(rownames = "id") %>%
  left_join(taxa_genes %>% ungroup() %>% select(id, substrate, site) %>% distinct())

vector_coords <- data.frame(gene_vectors[["vectors"]][["arrows"]]*sqrt(gene_vectors[["vectors"]][["r"]])) %>%
  as_tibble(rownames = "Gene abbreviation") %>%
  bind_cols(as_tibble(gene_vectors[["vectors"]][["pvals"]])) %>%
  rename(pval = value) %>%
  left_join(taxa_genes %>% ungroup() %>% select(`Gene abbreviation`, Category) %>% distinct()) %>%
  filter(pval <= 0.05) 



#NMDS plot with controls 
NMDS_plot <- NMDS_coords %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(size=2, alpha=0.8, aes(shape = site, color=substrate, label = id)) +
  geom_segment(data=vector_coords,inherit.aes = FALSE, aes(x=0,xend=NMDS1,y=0,yend=NMDS2, color=Category, label = `Gene abbreviation`), alpha=0.3)+
  #geom_text(data=vector_coords,inherit.aes = FALSE, aes(x = NMDS1, y = NMDS2, label=Category,hjust = 1, vjust = 1),  size=2, color="black") +  
  theme(legend.key.size = unit(.5, "cm"))

#visualize interactive plot
plotly::ggplotly(NMDS_plot)



metal_reducers <- data %>%
  filter(Category == "Metal reduction" & site %in% c("D1", "D3", "D6") & hits > 0) %>%
  left_join(metadata) %>%
  left_join(taxonomy)


# barplot -----------------------------------------------------------------

abundance_table <- otu_table %>%
  select(-taxonomy) %>%
  mutate_at(vars(-`#OTU ID`), funs(./sum(.)*100)) %>% #normalize to relative abundance 
  gather(sample_id, abundance, `9.DeMMO3.steri.11June2019`:`17.DeMMO1.tube8.12Dec2019`) 


taxon_abundance <- function(level, name){
  otu_table %>%
    select(`#OTU ID`, taxonomy) %>%
    mutate(taxonomy = gsub("Gammaproteobacteria; D_3__Betaproteobacteriales", "Betaproteobacteria; D_3__Betaproteobacteriales", taxonomy),
           taxa = str_extract(taxonomy, level),
           taxa = if_else(is.na(taxa), taxonomy, taxa)) %>%
    right_join(abundance_table) %>%
    ungroup() %>%
    group_by(sample_id, taxa) %>%
    summarise(abundance = sum(abundance)) %>%
    left_join(otu_metadata) %>% #add metadata
    group_by(site, `Date Sampled`, substrate, experiment_type, taxa) %>% 
    summarise(abundance = sum(abundance))
}

family_level <- "(.*)(?=; D_5__)"
class_level <- "(.*)(?=; D_3__)"
phylum_level <- "(.*)(?=; D_2__)"

less_abundant_taxa <- taxon_abundance(family_level, "family") %>%
  group_by(taxa) %>%
  filter(max(abundance) < 5) %>% #filter for families that represent less than 5% 
  group_by(site, `Date Sampled`, substrate, experiment_type) %>%
  summarise(abundance = sum(abundance)) %>%
  mutate(taxa = 'Less abundant taxa')


taxon_abundance_table <- taxon_abundance(family_level, "family") %>%
  group_by(taxa) %>%
  filter(max(abundance) >= 5) %>%
  bind_rows(less_abundant_taxa)

#palette
palette <- distinctColorPalette(k = length(unique(bar_plot$family)), altCol = FALSE, runTsne = FALSE)
names(palette) <- unique(unique(bar_plot$family)) 


#bar plot figure 
bar_plot <- taxon_abundance_table %>%
  ungroup() %>%
  mutate(taxonomy = str_remove_all(taxa, "D_0__| D_1__| D_2__| D_3__| D_4__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family")) %>%
  mutate(family = if_else(is.na(family) | str_detect(family, "uncultured"), order, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), class, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), phylum, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), domain, family),
         family = if_else(family == "GWF2-40-263", order, family),
         family = if_else(family == "mle1-8", phylum, family),
         experiment_type = factor(experiment_type, levels = c("rep", "exp")),
         substrate = factor(substrate, levels = rev(c("fluid", "sand", "Yates", "Homestake", "Poorman", "Ellison")))) %>%
  #family = factor(family, levels = family_color$family)) %>% 
  group_by(site, `Date Sampled`, substrate, experiment_type, family) %>%
  summarise(abundance = sum(abundance) *100) %>%
  ggplot(aes(fill=family, y=abundance, x=experiment_type, label = `Date Sampled`)) +
  geom_bar(stat='identity', position='fill') +
  scale_fill_manual(values=palette) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() + 
  theme(axis.title = element_blank(),
        legend.title = ggplot2::element_blank(), 
        legend.text = ggplot2::element_text(size = 8),
        legend.key.size = unit(0.25, "cm")) +
  facet_grid(cols = vars(site), rows = vars(substrate)) + 
  guides(fill = guide_legend(ncol = 1))

taxa_bar_plot <- plotly::ggplotly(bar_plot)


browsable(
  tagList(list(
    tags$div(
      style = 'width:100%;display:block;float:left;',
      taxa_bubble_plot
    ),
    tags$div(
      style = 'width:100%;display:block;float:left;',
      taxa_bar_plot
    )
  ))
)

