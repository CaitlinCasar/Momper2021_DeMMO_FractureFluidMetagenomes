pacman::p_load(tidyverse, readxl, stringi)

ncbi_organism_names <- read_xlsx("data/PRJNA563685orgnames.xlsx") %>%
  rename(organism_name = `proposed organism name`,
         organism = `NCBI recommended organism name`)

gene_counts_Feb2021 <- read_delim("data/geneCounts_19Feb2021.txt", delim = "\t", col_names = F) %>%
  separate(X1, c("id", "gene_counts"), " ") %>%
  separate(id, c("site", "genome"), "_") %>%
  mutate(site = str_remove(site, "eMMO"),
         genome = as.numeric(genome),
         gene_counts = as.numeric(gene_counts))

nucleotide_counts <- read_delim("data/nucleotideCounts.txt", delim = "\t", col_names = F) %>%
  separate(X1, c("id", "nucleotide_counts"), " ") %>%
  separate(id, c("site", "genome"), "_") %>%
  mutate(genome = as.numeric(genome),
         sample_name = paste(site, genome, sep = "_"),
         site = str_remove(site, "eMMO"),
         nucleotide_counts = as.numeric(nucleotide_counts)) %>%
  filter(nucleotide_counts >= 100000) %>%
  left_join(gene_counts_Feb2021) %>%
  left_join(taxonomy) %>%
  mutate(NCBI_phylum = recode(phylum, Gemmatimonadota = "Gemmatimonadetes",
                                         Omnitrophota = "Omnitrophica",
                                         Actinobacteriota = "Actinobacteria",
                                         Hadesarchaeota = "Hadesarchaea",
                                         Armatimonadota = "Armatimonadetes",
                                         Planctomycetota = "Planctomycetes",
                                         Chloroflexota = "Chloroflexi",
                                         Latescibacterota = "Latescibacteria",
                                         Bacteroidota = "Bacteroidetes",
                                         Nitrospirota = "Nitrospirae",
                                         Fusobacteriota = "Fusobacteria"),
         type = if_else(domain == "Bacteria", " bacterium", " archaeon")) %>%
  mutate(organism_name =if_else(!str_detect(species, "(?=.*[0-9]).*") & !is.na(species), species,
         if_else(!str_detect(genus, "(?=.*[0-9]).*") & !is.na(genus), paste0(genus, " sp."),
                 if_else(!str_detect(order, "(?=.*[0-9]).*") & !is.na(order), paste0(order, type),
                         if_else(!str_detect(class, "(?=.*[0-9]).*") & !is.na(class), paste0(class, type),
                                 if_else(!is.na(NCBI_phylum), paste0(NCBI_phylum, type), type)))))) %>%
  select(sample_name, organism_name) %>%
  left_join(ncbi_organism_names) %>%
  mutate(organism = if_else(is.na(organism), "Proteobacteria bacterium",organism)) %>%
  select(-organism_name)

#create batch submission for organism biosamples - note these can be submitted with max 400 MAGs
ncbi_mag_batch <- read_delim("data/MIMAG.water.5.0.tsv", delim = "\t", comment = "#") %>%
  select(contains("*")) %>%
  rename_all(.funs = funs(sub("\\*", "", .))) %>%
  right_join(nucleotide_counts) %>%
  mutate(isolate = sample_name,
         collection_date = "Apr-18",
         env_broad_scale = "continental subsurface",
         env_local_scale = "Deep Mine Microbial Observatory",
         env_medium = if_else(str_detect(sample_name, "DeMMO"), "fracture fluid",
                              if_else(str_detect(sample_name, "SW_"), "mine service water", "surficial creek")),
         geo_loc_name = "USA: Sanford Underground Research Facility, Lead, SD",
         lat_lon = "44.35 N 103.76 W",
         isolation_source = if_else(str_detect(sample_name, "DeMMO|SW_"), "DeMMO","surficial creek"),
         depth = if_else(str_detect(sample_name, "DeMMO1|DeMMO2"), "800 feet",
                         if_else(str_detect(sample_name, "DeMMO3"), "2000 feet",
                                 if_else(str_detect(sample_name, "DeMMO4|SW_"), "4100 feet",
                                         if_else(str_detect(sample_name, "DeMMO5|DeMMO6"), "4850 feet", "0 feet")))))

write_delim(ncbi_mag_batch[1:400,], "data/Momper2021_organism_biosample_batch.tsv", delim = "\t")
write_delim(ncbi_mag_batch[401:nrow(ncbi_mag_batch),], "data/Momper2021_organism_biosample_batch2.tsv", delim = "\t")




# MAG batch submission ----------------------------------------------------

MAG_attributes <- read_delim("MAG_attributes.tsv", delim = "\t") %>%
  bind_rows(read_delim("MAG_attributes2.tsv", delim = "\t")) %>%
  select(accession, sample_name) %>%
  rename(biosample_accession = "accession")

ncbi_mag_batch_files <- read_delim("data/sample_batch_genome_accs_fsa.91965d47e287.tsv", delim = "\t") %>%
  select(assembly_method, assembly_method_version, genome_coverage, sequencing_technology, filename) %>%
  right_join(ncbi_mag_batch %>% select(sample_name) %>% mutate(filename = sample_name)) %>%
  left_join(MAG_attributes) %>%
  mutate(assembly_method = "MEGAHIT", 
         assembly_method_version = "v1.0.6.1",
         genome_coverage = "1x",
         sequencing_technology = "Illumina HiSeq 2500",
         filename = paste0(filename, ".fa"))
         
         
mag_attribute_supp_file <- ncbi_mag_batch_files %>%
  mutate(biosample_accession = if_else(str_detect(sample_name, "DeMMO1"), "SAMN18064095",
                              if_else(str_detect(sample_name, "DeMMO2"), "SAMN18064236",
                                      if_else(str_detect(sample_name, "DeMMO3"), "SAMN18064310",
                                              if_else(str_detect(sample_name, "DeMMO4"), "SAMN18064413",
                                                      if_else(str_detect(sample_name, "DeMMO5"), "SAMN18064496",
                                                              if_else(str_detect(sample_name, "DeMMO6"), "SAMN18064575",
                                                                      if_else(str_detect(sample_name, "WC"), "SAMN18004502",
                                                                              if_else(str_detect(sample_name, "SW"), "SAMN18005272", "missing")))))))),
         bioproject_accession = "PRJNA563685") %>%
  rename(MAG_accession = "accession", ncbi_taxonomy = "organism") %>%
  write_delim("momper2021_MAG_attributes.tsv", delim = "\t")


#had to submit in two chunks due to 400 MAG limit per submission
write_delim(ncbi_mag_batch_files %>%
              filter(!str_detect(sample_name, "SW|DeMMO3")), "momper2021_batch_genome_accs_fsa.91965d47e287.tsv", delim = "\t")
write_delim(ncbi_mag_batch_files %>%
              filter(str_detect(sample_name, "SW|DeMMO3")), "momper2021_batch2_genome_accs_fsa.91965d47e287.tsv", delim = "\t")



biosample_dict <- c("SAMN18064095", "SAMN18064236", "SAMN18064310", "SAMN18064413", "SAMN18064496", "SAMN18064575", "SAMN18004502", "SAMN18005272")
names(biosample_dict) <-  c("D1", "D2", "D3", "D4", "D5", "D6", "WC", "SW")


#write files for submission to taxonomist for organism name approval
# unique_taxa <- nucleotide_counts %>%
#   select(organism_name) %>%
#   distinct()
# 
# write_csv(nucleotide_counts, "DeMMO_MAG_organism_names.csv")
# write_csv(unique_taxa, "DeMMO_MAG_organism_names_distinct.csv")
# 
# ncbi_batch_genome_accs <- read_delim("data/sample_batch_genome_accs_fsa.91965d47e287.tsv", delim = "\t")

# correct contaminant seqs ------------------------------------------------

#D5_50 
D5_50_contam_seq <- read_lines("data/NCBI_submission/genbank_errors/D5_50_contaminated.txt", skip = 1) %>%
  str_c(collapse = "") 

D5_50_clean_seq1 <- str_sub(D5_50_contam_seq, 1, 622) %>%
  stri_extract_all_regex('.{1,70}') %>%
  as.data.frame(col.names = c(">lcl|c_000000022732 Hadesarchaea archaeon isolate DeMMO5_50"), check.names =F) %>%
  write_delim("data/NCBI_submission/genbank_errors/DeMMO5_50.fsa", delim = "\n", append = T, col_names = T)

D5_50_clean_seq2 <- str_sub(D5_50_contam_seq, 4798, nchar(D5_50_contam_seq)) %>%
  stri_extract_all_regex('.{1,70}') %>%
  as.data.frame(col.names = c(">lcl|c_000000022733 Hadesarchaea archaeon isolate DeMMO5_50"), check.names =F) %>%
  write_delim("data/NCBI_submission/genbank_errors/DeMMO5_50.fsa", delim = "\n", append = T, col_names = T)

#D3_61
D3_61_contam_seq <- read_lines("data/NCBI_submission/genbank_errors/D3_61_contaminated.txt", skip = 1) %>%
  str_c(collapse = "") 

D3_61_clean_seq1 <- str_sub(D3_61_contam_seq, 1, 1193) %>%
  stri_extract_all_regex('.{1,70}') %>%
  as.data.frame(col.names = c(">lcl|c_000000013178 Candidatus Dojkabacteria bacterium isolate DeMMO3_61"), check.names =F) %>%
  write_delim("data/NCBI_submission/genbank_errors/DeMMO3_61.fsa", delim = "\n", append = T, col_names = T)

D3_61_clean_seq2 <- str_sub(D3_61_contam_seq, 1225, nchar(D3_61_contam_seq)) %>%
  stri_extract_all_regex('.{1,70}') %>%
  as.data.frame(col.names = c(">lcl|c_000000013179 Candidatus Dojkabacteria bacterium isolate DeMMO3_61"), check.names =F) %>%
  write_delim("data/NCBI_submission/genbank_errors/DeMMO3_61.fsa", delim = "\n", append = T, col_names = T)


#D6_55
D6_55_contam_seq <- read_lines("data/NCBI_submission/genbank_errors/D6_55_contaminated.txt", skip = 1) %>%
  str_c(collapse = "") 

D6_55_clean_seq1 <- str_sub(D6_55_contam_seq, 1, 3354) %>%
  stri_extract_all_regex('.{1,70}') %>%
  as.data.frame(col.names = c(">lcl|c_000000003906 Bacteroidales bacterium isolate DeMMO6_55"), check.names =F) %>%
  write_delim("data/NCBI_submission/genbank_errors/DeMMO6_55.fsa", delim = "\n", append = T, col_names = T)

D6_55_clean_seq2 <- str_sub(D6_55_contam_seq, 3414, nchar(D6_55_contam_seq)) %>%
  stri_extract_all_regex('.{1,70}') %>%
  as.data.frame(col.names = c(">lcl|c_000000003907 Bacteroidales bacterium isolate DeMMO6_55"), check.names =F) %>%
  write_delim("data/NCBI_submission/genbank_errors/DeMMO6_55.fsa", delim = "\n", append = T, col_names = T)

