pacman::p_load(ggraph, tidyverse, igraph)

files <- list.files("../data/metabolic/energy_flow", full.names = T, pattern = "Metabolic_network")

read_files <- function(file){
  site <- str_extract(file, "(?<=data/metabolic/energy_flow/)(.*)(?=_Metabolic_network_input[.]txt)")
  file %>%
    read_delim(delim="\t") %>%
    mutate(site = site)
}


file_list = lapply(files, read_files)
data <- reduce(file_list, bind_rows)

phylum_colors <- read_csv("../data/phylum_color_data.csv")
  
  
phylum_color_dict <- phylum_colors$hex.color
names(phylum_color_dict) <- phylum_colors$full.name

#find missing phyla in color dict
missing_phyla <- data.frame(unique(data$`Taxonomic Group`)[which(!unique(data$`Taxonomic Group`) %in% names(phylum_color_dict))])

metabolism <- data.frame(Step1 = c("C-S-01:Organic carbon oxidation",        
                "C-S-04:Acetate oxidation",               
                 "C-S-06:Fermentation",
                "C-S-02:Carbon fixation",                 
                "C-S-08:Methanotrophy",
                "C-S-07:Methanogenesis",
                "N-S-01:Nitrogen fixation",
                 "N-S-04:Nitrate reduction",                
                 "N-S-05:Nitrite reduction",               
                 "N-S-06:Nitric oxide reduction",           
                 "N-S-08:Nitrite ammonification", 
                "N-S-03:Nitrite oxidation",               
                "N-S-07:Nitrous oxide reduction",
                 "O-S-02:Arsenate reduction",              
                 "O-S-03:Arsenite oxidation", 
                "O-S-04:Selenate reduction",
                "O-S-01:Metal reduction",  
                 "C-S-09:Hydrogen oxidation",
                "C-S-05:Hydrogen generation",
                "S-S-01:Sulfide oxidation",
                 "S-S-03:Sulfur oxidation",              
                 "S-S-04:Sulfite oxidation",              
                 "S-S-05:Sulfate reduction",              
                 "S-S-07:Thiosulfate oxidation",         
                 "S-S-08:Thiosulfate disproportionation 1",
                "S-S-09:Thiosulfate disproportionation 2"), order = c(1:26),
                color = c(rep("#000", 6), rep("#34ccc0", 7), rep("#B2671B", 4),
                          rep("#787a78", 2), rep("#f4d98c", 7)))
metabolism_colors <- metabolism$color
names(metabolism_colors) <- metabolism$Step1

network_plot <- data %>%
  filter(site == "D1") %>%
  select(-`#Genome`) %>%
  left_join(metabolism) %>%
  arrange(order) %>%
  mutate(Category = str_extract(Step1, "[^-]+"),
         Category = recode(Category, "C" = "Carbon", "S" = "Sulfur", "N" = "Nitrogen", "O" = "Oxygen")) %>%
ggraph(layout = 'linear', circular = TRUE) + 
  geom_edge_arc(aes(colour = as.factor(`Taxonomic Group`), edge_width = `Coverage value(average)`), alpha = 0.2) +
  scale_edge_color_manual(values = phylum_color_dict) +
  geom_node_point(aes(color = as.factor(name)), size = 3)+ 
  scale_color_manual(values = metabolism_colors, guide = F) +
  geom_node_text(aes(label = name)) +
  theme_bw()
  


