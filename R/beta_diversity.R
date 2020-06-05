otu_table <- read_delim("../../submission/data/16S/metagenome_samples_noChimera_otuTable_withTaxa74000.txt", delim = '\t', comment = "# ")

unweighted_unifrac <- read_delim("../../submission/data/16S/beta_diversity/unweighted_unifrac_coords.txt",delim = '\t') %>%
  ggplot(aes(PC1, PC2, color = site)) + 
  geom_point() +
  ggtitle("Unweighted Unifrac")

weighted_unifrac <- read_delim("../../submission/data/16S/beta_diversity/weighted_unifrac_coords.txt",delim = '\t') %>%
  ggplot(aes(PC1, PC2, color = site)) + 
  geom_point() +
  ggtitle("Weighted Unifrac") 


browsable(
  tagList(list(
    tags$div(
      style = 'width:50%;display:block;float:left;',
      ggplotly(unweighted_unifrac)
    ),
    tags$div(
      style = 'width:50%;display:block;float:left;',
      ggplotly(weighted_unifrac)
    )
  ))
)