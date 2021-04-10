pacman::p_load(tidyverse,  ggalluvial, plotly)

files <- list.files("data/metabolic/energy_flow", full.names = T, pattern = "energy_flow")

read_files <- function(file){
  site <- str_extract(file, "(?<=data/metabolic/energy_flow/)(.*)(?=_Metabolic_energy_flow_input[.]txt)")
  file %>%
    read_delim(delim="\t", col_names = c("Phylum", "Function", "Abundance")) %>%
    mutate(site = site)
}


file_list = lapply(files, read_files)
data <- reduce(file_list, bind_rows)

phylum_colors <- read_csv("data/phylum_color_data.csv") 
  
  
phylum_color_dict <- phylum_colors$hex.color
names(phylum_color_dict) <- phylum_colors$full.name


energy_plot <- data %>%
  filter(site == "D1") %>%
  mutate(Category = str_extract(Function, "[^-]+"),
         Category = recode(Category, "C" = "Carbon", "S" = "Sulfur", "N" = "Nitrogen", "O" = "Oxygen")) %>%
  ggplot(aes(y = Abundance, axis1 = Phylum, axis2 = reorder(Function, Category), axis3 = Category)) +
  geom_alluvium(aes(fill = Phylum), width = 1/12) +
  scale_fill_manual(values = phylum_color_dict) +
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_label(stat = "stratum", infer.label = TRUE) +
  scale_x_discrete(limits = c("Phylum", "Function"), expand = c(.05, .05)) 
  #scale_fill_brewer(type = "qual", palette = "Set1")