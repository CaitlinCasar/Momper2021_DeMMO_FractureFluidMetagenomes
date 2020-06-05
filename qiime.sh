#interactive job
srun --account=p30777 --time=04:00:00 --partition=short --mem=18G --pty bash

#to activate local qiime environment
module load anaconda3
source activate qiime1

#navigate to projects folder
cd /projects/p30777/DeMMO_16s/DeMMO_Mar2020

###join paired end reads
~/pear-0.9.11-linux-x86_64/bin/pear -f Undetermined_S0_L001_R1_001.fastq -r Undetermined_S0_L001_R2_001.fastq -j 8 -o 08032020_OsburnLab_joined.PEAR|tee pear_assembly_log.txt

-----output-----
Forward reads file.................: Undetermined_S0_L001_R1_001.fastq
Reverse reads file.................: Undetermined_S0_L001_R2_001.fastq
PHRED..............................: 33
Using empirical frequencies........: YES
Statistical method.................: OES
Maximum assembly length............: 999999
Minimum assembly length............: 50
p-value............................: 0.010000
Quality score threshold (trimming).: 0
Minimum read size after trimming...: 1
Maximal ratio of uncalled bases....: 1.000000
Minimum overlap....................: 10
Scoring method.....................: Scaled score
Threads............................: 8

Allocating memory..................: 200,000,000 bytes
Computing empirical frequencies....: DONE
  A: 0.224617
  C: 0.263137
  G: 0.251739
  T: 0.260507
  7 uncalled bases
Assemblying reads: 100%

Assembled reads ...................: 17,174,963 / 20,450,706 (83.982%)
Discarded reads ...................: 0 / 20,450,706 (0.000%)
Not assembled reads ...............: 3,275,743 / 20,450,706 (16.018%)
Assembled reads file...............: 08032020_OsburnLab_joined.PEAR.assembled.fastq
Discarded reads file...............: 08032020_OsburnLab_joined.PEAR.discarded.fastq
Unassembled forward reads file.....: 08032020_OsburnLab_joined.PEAR.unassembled.forward.fastq
Unassembled reverse reads file.....: 08032020_OsburnLab_joined.PEAR.unassembled.reverse.fastq
-----output-----

###remove the discarded reads from the barcodes file
~/fastq-barcode.pl Undetermined_S0_L001_I1_001.fastq 08032020_OsburnLab_joined.PEAR.assembled.fastq > barcodes.join.PEAR.fastq

###check for errors in mapping file###
validate_mapping_file.py -m 200302_Osburn_16S_KEH_200227_edited.txt

###demultiplex the reads 
split_libraries_fastq.py -i 08032020_OsburnLab_joined.PEAR.assembled.fastq -o split_libraries_fastq_output -b barcodes.join.PEAR.fastq --barcode_type 12 -m 200302_Osburn_16S_KEH_200227_edited.txt

###copy files from RDSS to Quest
smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/DeMMO_16s;  cd DeMMO/DNA_Sequence_Data/Environmental_16s/Orig_data/Sep2018_June2019/; mget *' -U "ads\cpc7770"
smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/Osburn2020;  cd /DeMMO/Publications/Osburn2020/16s_data/; mget *' -U "ads\cpc7770"


###split sequences up into smaller groups###
filter_fasta.py -f split_libraries_fastq_output/seqs.fna --sample_id_fp metagenome_samples_map.txt -o ../metagenome_samples.fasta


###check that sample names are formatted properly 
~/usearch-v11/usearch -fastx_get_sample_names metagenome_samples.fasta -sample_delim _ -output sample_names.txt

-----output-----
00:04 37Mb    100.0% 8 samples found
-----output-----

###dereplicate sequences
~/usearch-v11/usearch -fastx_uniques metagenome_samples.fasta -fastaout metagenome_samples_derep.fasta -sizeout

-----output-----
00:02 344Mb   100.0% Reading metagenome_samples.fasta
00:02 310Mb  CPU has 24 cores, defaulting to 10 threads

WARNING: Max OMP threads 1

00:03 348Mb   100.0% DF
00:03 359Mb  770702 seqs, 303474 uniques, 252425 singletons (83.2%)
00:03 359Mb  Min size 1, median 1, max 12717, avg 2.54
00:05 346Mb   100.0% Writing metagenome_samples_derep.fasta
-----output-----

###sort sequences
~/usearch-v11/usearch -sortbysize metagenome_samples_derep.fasta -fastaout metagenome_samples_derep_sorted.fasta

-----output-----
00:00 163Mb   100.0% Reading metagenome_samples_derep.fasta
00:00 129Mb  Getting sizes
00:01 132Mb  Sorting 303474 sequences
00:02 133Mb   100.0% Writing output
-----output-----

###pick otus
~/usearch-v11/usearch -cluster_otus metagenome_samples_derep_sorted.fasta -otus metagenome_samples_derep_sorted_otus.fasta -relabel OTU_ -uparseout metagenome_samples_derep_sorted_otus_results.txt -minsize 2

-----output-----
00:19 64Mb    100.0% 6258 OTUs, 1921 chimeras
-----output-----

###generate otu tables
~/usearch-v11/usearch -usearch_global metagenome_samples.fasta -sample_delim _ -db metagenome_samples_derep_sorted_otus.fasta -sample_delim _ -strand plus -id 0.97 -otutabout metagenome_samples_noChimera.otutable.txt -biomout metagenome_samples_otuTable.biom

-----output-----
00:00 42Mb    100.0% Reading metagenome_samples_derep_sorted_otus.fasta
00:00 8.8Mb   100.0% Masking (fastnucleo)
00:00 9.6Mb   100.0% Word stats
00:00 9.6Mb   100.0% Alloc rows
00:00 16Mb    100.0% Build index
00:00 49Mb   CPU has 24 cores, defaulting to 10 threads

WARNING: Max OMP threads 1

00:53 51Mb    100.0% Searching, 91.9% matched
708279 / 770702 mapped to OTUs (91.9%)
00:53 51Mb   Writing metagenome_samples_noChimera.otutable.txt
00:53 51Mb   Writing metagenome_samples_noChimera.otutable.txt ...done.
00:53 51Mb   Writing metagenome_samples_otuTable.biom
00:53 51Mb   Writing metagenome_samples_otuTable.biom ...done.
-----output-----

###convert otu table to qiime compatible format 
biom convert -i metagenome_samples_noChimera.otutable.txt -o  metagenome_samples_hdf5.biom --table-type="OTU table" --to-hdf5

###assign taxonomy using Silava132 database
assign_taxonomy.py -i metagenome_samples_derep_sorted_otus.fasta -t ~/SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt -r ~/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna --similarity 0.9 --uclust_max_accepts 10 --min_consensus_fraction 0.90 -o uclust_taxa_0.9_10_9.0_chimerasRemoved/

###add metadata to otu table 
biom add-metadata -i metagenome_samples_hdf5.biom -o metagenome_samples_noChimera_otuTable_withTaxa.biom --observation-metadata-fp uclust_taxa_0.9_10_9.0_chimerasRemoved/metagenome_samples_derep_sorted_otus_tax_assignments.txt --sc-separated taxonomy --observation-header OTUID,taxonomy

###check add-metadata output
biom convert -i metagenome_samples_noChimera_otuTable_withTaxa.biom -o metagenome_samples_noChimera_otuTable_withTaxa.txt --to-tsv --header-key taxonomy

###summarize otu data 
biom summarize-table -i metagenome_samples_noChimera_otuTable_withTaxa.biom -o metagenome_samples_noChimera_otuTable_withTaxa_table_summary.txt

-----output-----
Num samples: 8
Num observations: 6,258
Total count: 708,279
Table density (fraction of non-zero values): 0.175

Counts/sample summary:
 Min: 74,080.000
 Max: 137,620.000
 Median: 82,686.500
 Mean: 88,534.875
 Std. dev.: 19,467.361
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
D4: 74,080.000
D6: 74,109.000
WC: 77,521.000
D3: 80,485.000
D5: 84,888.000
D1: 88,588.000
D2: 90,988.000
SW: 137,620.000
-----output-----

###align representative set of sequences (otus)
align_seqs.py -i metagenome_samples_derep_sorted_otus.fasta -e 250

###remove gaps from alignment
filter_alignment.py -i pynast_aligned/metagenome_samples_derep_sorted_otus_aligned.fasta -o pynast_aligned/

###make phylogenetic tree
make_phylogeny.py -i pynast_aligned/metagenome_samples_derep_sorted_otus_aligned_pfiltered.fasta -o metagenome_samples.tre

###perform multiple rarefactions to generate data for rarefaction curves, cutoff 62,294 is median read depth 
multiple_rarefactions.py -i metagenome_samples_noChimera_otuTable_withTaxa.biom -m 500 -x 74000 -s 1000 -n 10 -o rarefied_otu_tables/

###calculate alpha diversity on rarefied data
alpha_diversity.py -i rarefied_otu_tables/ -m chao1,chao1_ci,PD_whole_tree,observed_otus,shannon,simpson,simpson_e -o alpha_div/ -t  metagenome_samples.tre

###concatenate alpha diversity results for generating rarefaction curves 
collate_alpha.py -i alpha_div/ -o collated_alpha/

###move files from Quest to RDSS for processing in R 
smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd /projects/p30777/DeMMO_16s/NativeRock/collated_alpha/; cd /DeMMO/; mput *' -U "ads\cpc7770"

###beta diversity
beta_diversity_through_plots.py -i metagenome_samples_noChimera_otuTable_withTaxa.biom -m ../DeMMO_Mar2020/metagenome_samples_map.txt -o beta_diversity/ -t metagenome_samples.tre --seqs_per_sample 74000 --force

###rarefy otu table to normalize data to a uniform read depth based on lowest read count from dataset 
single_rarefaction.py -i metagenome_samples_noChimera_otuTable_withTaxa.biom -d 74000 -o metagenome_samples_noChimera_otuTable_withTaxa_74000.biom

###summarize otu table
summarize_taxa.py -i metagenome_samples_noChimera_otuTable_withTaxa_74000.biom  -o tax_summaries -a --suppress_biom_table_output

#convert rarefied biom table to tsv
biom convert -i metagenome_samples_noChimera_otuTable_withTaxa_74000.biom -o tax_summaries/metagenome_samples_noChimera_otuTable_withTaxa74000.txt --to-tsv --header-key taxonomy








