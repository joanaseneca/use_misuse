###### SCRIPT SUMMARY  #######

# Author: Joana Séneca, November 2020
# Title: Micron_amplicon_analysis.R
#
# This is the raw R script. A summary of all these analyses can be found on the html Rmarkdown file.
# 
# This script imports, filters and exports OTU tables for the 16S and ITS1 amplicon sequencing dataset originally published by:
# Zheng, Q. et al. Soil multifunctionality is affected by the soil environment and by microbial community composition and diversity. Soil Biol. Biochem. 136, 107521 (2019).  
# 
# The raw fastq files are deposited on NCBI under Bioproject number PRJNA551019.
# 
# Shortly, both 16S and ITS datasets were 
# 1. Taxonomically filtered: OTUs not affiliated with the groups of interests were excluded;
# 2. Abundance filtered: OTUs with an overall abundance < 0.01%  were excluded.

###### LOAD PACKAGES AND FUNCTIONS #####

library(dplyr)
library(tidyverse)
library(tibble)
library(phyloseq)


#https://rdrr.io/github/vmikk/metagMisc/src/R/phyloseq_filter.R

#' @title Remove taxa with small mean relative abundance.
#'
#' @param physeq A phyloseq-class object
#' @param frac The minimum cutoff for the relative OTU abundance
#' @details This function searches for taxa with small mean relative abundance and removes them. Result will be returned with original counts in the abundance table.
#' @return Phyloseq object with a subset of taxa.
#' @export
#'
#' @examples
#' data("esophagus")
#' phyloseq_filter_taxa_rel_abund(esophagus, frac = 0.01)                  
#' phyloseq_filter_taxa_rel_abund(esophagus, frac = 0.1)
#'
phyloseq_filter_taxa_rel_abund <- function(physeq, frac = 1e-4){
  
  # require(phyloseq)
  
  ## Transform OTU counts to relative abundance
  rel <- phyloseq::transform_sample_counts(physeq, function(x) x / sum(x) )
  
  ## Filter OTUs
  rel.subs <- phyloseq::filter_taxa(rel, function(x){ mean(x) > frac }, prune = FALSE)
  
  ## if prune = TRUE
  # tn <- taxa_names(rel.subs)              # OTUs to preserve
  # tr <- setdiff(taxa_names(physeq), tn)   # OTUs to remove
  
  ## Taxa to remove
  tr <- names(rel.subs)[ which(rel.subs == FALSE) ]
  
  ## If all taxa should be removed
  if(length(tr) == phyloseq::ntaxa(physeq)){
    stop("Error: all taxa will be removed with the specified 'frac' cutoff.\n")
  }
  
  ## If there is nothing to remove
  if(length(tr) == 0){
    res <- physeq
    cat("Warning: no taxa removed.\n")
  }
  
  ## Keep taxa which satisfies the truncation threshold
  if(length(tr) > 0){
    res <- phyloseq::prune_taxa(taxa = rel.subs, physeq)
  }
  
  return(res)
}



###### 16S IMPORT ########

#we will load the metadata file and the OTU table. The last columns of the OTU table document contain the taxonomy
metadata_16S <- read.delim("~/Desktop/MicroN/01_input/metadata_16S.txt")
otu_table_16S <- read.delim("~/Desktop/MicroN/01_input/otuCountTable_16S.txt")

rownames(metadata_16S) <- metadata_16S$Sequencing_ID
dim(metadata_16S)
dim(otu_table_16S) 

# The samples and fastq files deposited on NCBI were more than the ones used for this experiment
# therefore we need to subset the samples we are interested in

my_samples <- metadata_16S$Sequencing_ID
my_samples <- as.character(my_samples)

sub_otu_table_16S <- select_(otu_table_16S, .dots = my_samples)
sub_otu_table_16S <- cbind(otu_table_16S$OTUnumber, sub_otu_table_16S)
rownames(sub_otu_table_16S) <- sub_otu_table_16S[,1]
sub_otu_table_16S$`otu_table_16S$OTUnumber` <- NULL

taxonomy_16S <- otu_table_16S %>% select(OTUnumber, domain_id, phylum_id, class_id, order_id, family_id, genus_id ) 
rownames(taxonomy_16S) <- taxonomy_16S$OTUnumber
taxonomy_16S$OTUnumber <- NULL
taxonomy_16S <- as.matrix(taxonomy_16S)

#final check. Do the rownames match?
dim(taxonomy_16S)
dim(sub_otu_table_16S)
dim(metadata_16S)

colnames(sub_otu_table_16S) == rownames(metadata_16S) # should be TRUE for all samples

#create phyloseq object for 16S
phy_16S <-phyloseq(otu_table(sub_otu_table_16S, taxa_are_rows=T), tax_table(taxonomy_16S), sample_data(metadata_16S))



###### ITS IMPORT ########

#we will load the metadata file and the OTU table. The last columns of the OTU table document contain the taxonomy
metadata_ITS <- read.delim("~/Desktop/MicroN/01_input/metadata_ITS.txt")
otu_table_ITS <- read.delim("~/Desktop/MicroN/01_input/otuCountTable_ITS.txt")
rownames(metadata_ITS) <- metadata_ITS$Sequencing_ID
dim(metadata_ITS)
dim(otu_table_ITS) 

# The samples and fastq files deposited on NCBI were more than the ones used for this experiment
# therefore we need to subset the samples we are interested in

my_samples <- metadata_ITS$Sequencing_ID
my_samples <- as.character(my_samples)

sub_otu_table_ITS <- select_(otu_table_ITS, .dots = my_samples)
sub_otu_table_ITS <- cbind(otu_table_ITS$OTUnumber, sub_otu_table_ITS)
rownames(sub_otu_table_ITS) <- sub_otu_table_ITS[,1]
sub_otu_table_ITS$`otu_table_ITS$OTUnumber` <- NULL

taxonomy_ITS <- otu_table_ITS %>% select(OTUnumber, domain_id, phylum_id, class_id, order_id, family_id, genus_id ) 
rownames(taxonomy_ITS) <- taxonomy_ITS$OTUnumber
taxonomy_ITS$OTUnumber <- NULL
taxonomy_ITS <- as.matrix(taxonomy_ITS)

#final check. Do the rownames match?
dim(taxonomy_ITS)
dim(sub_otu_table_ITS)
dim(metadata_ITS)

colnames(sub_otu_table_ITS) == rownames(metadata_ITS) # should be TRUE for all samples

#create phyloseq object for ITS
phy_ITS <-phyloseq(otu_table(sub_otu_table_ITS, taxa_are_rows=T), tax_table(taxonomy_ITS), sample_data(metadata_ITS))


#clean up to keep only the phyloseq objects and the function. Any subset of data can be extracted from the phyloseq object
rm(list=setdiff(ls(), c("phy_16S", "phy_ITS","phyloseq_filter_taxa_rel_abund")))


#the filtering parameters used are summarized in an RMarkdown
###### FILTERING ITS#####


# initial phy objects
phy_ITS # 3060 OTUs, 24 samples

rank_names(phy_ITS)
#"domain_id" "phylum_id" "class_id"  "order_id"  "family_id" "genus_id" 

plot_bar(phy_ITS, fill='phylum_id') #there's a sample that failed and must be removed
sample_names(phy_ITS) # get sample names
sample_sums(phy_ITS) # get reads per sample

phy_ITS_pruned <- subset_samples(phy_ITS, Sequencing_ID !="SAMPLE.767.GF.II.8.60")
phy_ITS_pruned # 3609 OTUs, 23 samples

get_taxa_unique(phy_ITS_pruned, taxonomic.rank = 'domain_id') # good, all OTUs are Eukaryotes
get_taxa_unique(phy_ITS_pruned, taxonomic.rank = 'phylum_id') # good, all OTUs are Fungi, no need to remove non-specific OTUs

# removing OTUs with an overall abundance of < 0.01% in the whole dataset
phy_ITS_f1 <- phyloseq_filter_taxa_rel_abund(phy_ITS_pruned, frac=0.0001)
phy_ITS_f1 # 1072 taxa, 23 samples

plot_bar(phy_ITS_f1, fill='phylum_id') 

#obtain minimum, average and maximum numer of reads per sample (just to check)
min(sample_sums(phy_ITS_f1)) # the minimum number of reads per sample is 717
mean(sample_sums(phy_ITS_f1)) # the average number of reads per sample is 2646.304
max(sample_sums(phy_ITS_f1))# the maximum number of reads per sample is 5321

rm(phy_ITS_pruned)



###### FILTERING 16S#####

# initial phy objects
phy_16S # 7550 OTUs, 24 samples

rank_names(phy_16S)
#"domain_id" "phylum_id" "class_id"  "order_id"  "family_id" "genus_id" 

plot_bar(phy_16S, fill='phylum_id') #there's a sample that failed, and must be removed
sample_names(phy_16S) # get sample names
sample_sums(phy_16S) # get reads per sample

phy_16S_pruned <- subset_samples(phy_16S, Sequencing_ID !="SAMPLE.325.GF.II.5.60")
phy_16S_pruned # 7550 OTUs, 23 samples


get_taxa_unique(phy_16S_pruned, taxonomic.rank = 'domain_id')
#"Bacteria"  "Archaea"   "unknown"   "Eukaryota"

# We will want to remove Eukaryota but also maybe the unknowns. Let's see how abundant they are.
sub <- subset_taxa(phy_16S_pruned, domain_id %in% c("unknown", "Eukaryota"))
ntaxa(sub) # 39 OTUs total
plot_bar(sub, fill='domain_id') # how are these OTUs distributed and how abundant are they?

# ok, they're very low abundant (max 20 reads in one sample)
# All OTUs with this taxonomic classification will be removed.

phy_16S_f1 <- subset_taxa(phy_16S_pruned, !domain_id %in% c("unknown", "Eukaryota"))

phy_16S_f1 # 7511 OTUs, 23 samples


# removing OTUs with an overall abundance of < 0.01% in the whole dataset
phy_16S_f2 <- phyloseq_filter_taxa_rel_abund(phy_16S_f1, frac=0.0001)
phy_16S_f2 # 1486 taxa, 23 samples

plot_bar(phy_16S_f2, fill='phylum_id') 

#obtain minimum, average and maximum numer of reads per sample (just to check)
min(sample_sums(phy_16S_f2)) # the minimum number of reads per sample is 1302
mean(sample_sums(phy_16S_f2)) # the average number of reads per sample is 8823.522
max(sample_sums(phy_16S_f2))# the maximum number of reads per sample is 16795

rm(phy_16S_f1, phy_16S_pruned,sub)


###### EXPORTING FILTERED OTU TABLES FOR POWER ANALYSIS #####

# saving original phyloseq objects. This makes things easier afterwards, as we don't need to create them again.
save(phy_16S, phy_16S_f2, phy_ITS, phy_ITS_f1, file="~/Desktop/MicroN/02_data/MicroN_phyloseq_objects.RData")

#exporting the ITS filtered OTU table
filt_otu_table_ITS <- otu_table(phy_ITS_f1)
filt_otu_table_ITS <- as.data.frame(filt_otu_table_ITS)

str(filt_otu_table_ITS) 

#convert integers to numbers
filt_otu_table_ITS[1:ncol(filt_otu_table_ITS)] <- lapply(filt_otu_table_ITS[1:ncol(filt_otu_table_ITS)], as.numeric)
str(filt_otu_table_ITS)

write.table(filt_otu_table_ITS, file="~/Desktop/MicroN/03_output/ITS_OTU_table_filt.txt", quote=F)


#exporting the 16S filtered OTU table
filt_otu_table_16S <- otu_table(phy_16S_f2)
filt_otu_table_16S <- as.data.frame(filt_otu_table_16S)

str(filt_otu_table_16S) 

#convert integers to numbers
filt_otu_table_16S[1:ncol(filt_otu_table_16S)] <- lapply(filt_otu_table_16S[1:ncol(filt_otu_table_16S)], as.numeric)
str(filt_otu_table_16S)

write.table(filt_otu_table_16S, file="~/Desktop/MicroN/03_output/16S_OTU_table_filt.txt", quote=F)

