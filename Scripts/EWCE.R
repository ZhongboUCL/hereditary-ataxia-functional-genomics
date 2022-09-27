# Script for Expression Weighted Cell-type Enrichment(EWCE)-------------------------------------------------------------------------------------
##- Skene, et al. Identification of Vulnerable Cell Types in Major Brain Disorders Using Single Cell Transcriptomes and Expression Weighted Cell Type Enrichment. Front. Neurosci, 2016.
##- Zeisel, et al. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science, 2015.
##- Level 3 data 

library(tidyverse)
library(readxl)
install_github("nathanskene/ewce") # EWCE from Skene Github
library(EWCE)

# Load data -------------------------------------------------------------------------------------------

# Load ataxia gene lists
ataxia_green <- read_delim("/home/zchen/ML/panelapp/new_gene_list/reclassified_ataxia_genes_25062020csv.csv", delim = ",")

# Functions -------------------------------------------------------------------------------------------

# Main ------------------------------------------------------------------------------------------------

# 1. Generate gene lists for testing in EWCE

ataxia_green <- ataxia_green %>% mutate(Panel = reclassification_category)

# Three groups:
adult_only <- ataxia_green %>% 
  filter(Panel == "Adult-onset")

child_only <- ataxia_green %>% 
  filter(Panel == "Childhood-onset")

adult_child_overlap <- ataxia_green %>% 
  filter(Panel == "Overlap")

# Background gene list:
# All human genes with mouse orthologs: as accounting for transcript length and GC content, use human gene IDs and not mouse

data("mouse_to_human_homologs")
m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])

human.hits.adult = unique(m2h[m2h$HGNC.symbol %in% adult_only$gene_symbol,"HGNC.symbol"])
human.bg.adult = unique(c(human.hits.adult,m2h$HGNC.symbol))

human.hits.child <- unique(m2h[m2h$HGNC.symbol %in% child_only$gene_symbol,"HGNC.symbol"])
human.bg.child = unique(c(human.hits.child,m2h$HGNC.symbol))

human.hits.overlap <- unique(m2h[m2h$HGNC.symbol %in% adult_child_overlap$gene_symbol,"HGNC.symbol"])
human.bg.overlap = unique(c(human.hits.overlap,m2h$HGNC.symbol))

#mouse.bg  = unique(setdiff(m2h$MGI.symbol,mouse.hits))
#mouse.bg = unique(m2h$MGI.symbol)

# 2. Setting analysis parameters
# Running EWCE analysis on genetic data
# ctd_Zeisel2018_Skene
load("/home/zchen/ML/predictors/gene_feature_predictors/ctd_Zeisel2018_Skene.rda")

# Zeisel Data - all cell types Level three

level=3 #  (i.e. Interneurons)
reps=10000

full_results_adults_level3 = bootstrap.enrichment.test(sct_data=ctd_Zeisel2018,
                                                       hits=human.hits.adult,
                                                       bg=human.bg.adult,
                                                       reps=reps,       
                                                       annotLevel=3,
                                                       geneSizeControl=TRUE,
                                                       genelistSpecies="human",
                                                       sctSpecies="mouse")

full_results_child_level3 = bootstrap.enrichment.test(sct_data=ctd_Zeisel2018,
                                                      hits=human.hits.child,
                                                      bg=human.bg.child,
                                                      reps=reps,       
                                                      annotLevel=3,
                                                      geneSizeControl=TRUE,
                                                      genelistSpecies="human",
                                                      sctSpecies="mouse")

full_results_overlap_level3 = bootstrap.enrichment.test(sct_data=ctd_Zeisel2018,
                                                        hits=human.hits.overlap,
                                                        bg=human.bg.overlap,
                                                        reps=reps,       
                                                        annotLevel=3,
                                                        geneSizeControl=TRUE,
                                                        genelistSpecies="human",
                                                        sctSpecies="mouse")


# Most significant results
# Adult-onset
full_results_adults_results_level3 <- as.data.frame(full_results_adults_level3$results)

# Most significant results
# Childhiid-onset
full_results_child_results_level3  <- as.data.frame(full_results_child_level3$results)

# Most significant results
# Overlap-onset
full_results_overlap_results_level3  <- as.data.frame(full_results_overlap_level3$results)

# Merge results
full_results_adult_merge_level3 = data.frame(full_results_adults_results_level3,list="Adult-onset")
full_results_child_merge_level3 = data.frame(full_results_child_results_level3,list="Childhood-onset")
full_results_overlap_merge_level3 = data.frame(full_results_overlap_results_level3,list="Overlap")

merged_results_level3_all_Zeisel = rbind(full_results_adult_merge_level3, full_results_child_merge_level3,full_results_overlap_merge_level3)


# Change SD <0 to 0 and apply significance to those significant
# Apply BH correction by number of conditions (i.e. celltypes x 3)
merged_results_level3_all_Zeisel_noneg <- merged_results_level3_all_Zeisel %>%
  mutate(sd_from_mean = ifelse(sd_from_mean <=0, 0, sd_from_mean)) %>%
  mutate(BH_corrected_p = p.adjust(p, method = "BH", n = length(p)))%>%
  mutate(p_value = ifelse(BH_corrected_p < 0.05, "significant", "non-significant")) %>%
  mutate(class = ifelse(str_detect(CellType, "neuron"), "Neuron",
                        ifelse(str_detect(CellType, "immune"), "Immune cells",
                               ifelse(str_detect(CellType, "vascular"), "Vascular cells", 
                                      ifelse(CellType %in% c("Oligodendrocytes", "Astroependymal cells", "Neural crest-like glia"), "Glia", "Other")))))%>%
  mutate(CellType = fct_reorder(CellType, class))


## Which genes are contributing to cerebellar neuronal expression in overlap?

cerebellar_neurons_bootstrap <- full_results_overlap_level3$bootstrap_data %>% as.data.frame() %>%
  dplyr::select(`Cerebellum neurons`)%>%  
  mutate(percent_expression = (`Cerebellum neurons`/sum(`Cerebellum neurons`))*100)

# mean specificity
mean_specificity_across_overlap_level3 <- ctd_Zeisel2018[[3]]$specificity %>% 
  as.data.frame() %>%
  dplyr::select(`Cerebellum neurons`) 

mean_specificity_across_overlap_level3$gene <- rownames(mean_specificity_across_overlap_level3)

mouse.hits.overlap<- unique(m2h[m2h$HGNC.symbol %in% adult_child_overlap$gene_symbol,"MGI.symbol"])

mean_specificity_across_overlap1_level3<- mean_specificity_across_overlap_level3%>%
  dplyr::filter(gene %in% mouse.hits.overlap) %>%
  mutate(expression_cell_type = (`Cerebellum neurons`/sum(`Cerebellum neurons`))*100) %>%
  left_join(m2h, by = c("gene" = "MGI.symbol")) 



# Save data -------------------------------------------------------------------------------------------
