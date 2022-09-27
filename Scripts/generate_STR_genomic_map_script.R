# - Script to annotate HipSTR STR map with genic locations and form predictor -------------------------

library(tidyverse)
library(readxl)
library(ggpubr)
library(rtracklayer)
library(data.table)
library(biomaRt)
library(GenomicRanges)

# Load data -------------------------------------------------------------------------------------------
## Load GRCh38 option of HipSTR Reference Genome

STR_map <- fread("/home/zchen/ML/predictors/STR/hipstr/resources/GRCh38.hipstr_reference.bed.gz")

# Functions -------------------------------------------------------------------------------------------

# Main ------------------------------------------------------------------------------------------------

## Load GRCh38 option

STR_map_granges <- GRanges(seqnames = STR_map$V1,
                             ranges = IRanges(start = STR_map$V2,
                                              end = STR_map$V3),
                             number_nucleotides = STR_map$V4,
                             repeat_number = STR_map$V5,
                             HipSTR_number = STR_map$V6,
                             STR_sequence = STR_map$V7)

## Annotate STR map /w/ Ensembl ID

Homo_sapiens.GRCh38.97 <- import("/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf")

gene_list <- Homo_sapiens.GRCh38.97[Homo_sapiens.GRCh38.97$type == "gene"]

STR_overlap_gtf_1 <- findOverlaps(query = STR_map_granges, subject = gene_list, 
                                  maxgap = -1, minoverlap = 1)

STR_overlap_gtf_1_w_annot <- 
  STR_overlap_gtf_1 %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(subject_annot = (gene_list$gene_id)[c(subjectHits(STR_overlap_gtf_1))])

STR_overlap_gtf_1_w_annot_no_dup <- 
  STR_overlap_gtf_1_w_annot %>%
  filter(!duplicated(queryHits, fromLast = F)) 

overlap_STR_map <- 
  STR_overlap_gtf_1 %>% 
  as.data.frame() %>%
  as_tibble() %>%
  mutate(overlapped_genes = (STR_map_granges$HipSTR_number %>% as.character())[queryHits(STR_overlap_gtf_1)])

annotated_STR_map <- left_join(STR_map, overlap_STR_map, by = c("V6" = "overlapped_genes"))

annotate_STR_map_Ensembl <- annotated_STR_map %>% 
  left_join(STR_overlap_gtf_1_w_annot_no_dup, by = "queryHits") %>%
  dplyr::select(-(c(subjectHits.x, subjectHits.y, queryHits))) %>%
  mutate(region = ifelse(is.na(subject_annot), "intergenic", NA)) %>% unique()

# How many STRs have no Ensembl ID? 631902/1638945 STRs are not annotated with Ensembl ID, i.e. "intergenic"
annotate_STR_map_Ensembl$subject_annot %>% is.na() %>% sum()

colnames(annotate_STR_map_Ensembl) <- c("chr", "start", "end", "nucleotides_per_repeat", "number_of_repeats", 
                                        "HipSTR_ID", "repeat_sequence", "gene_id", "region")

## Double check NA in regions = 1133733
annotate_STR_map_Ensembl$region %>% is.na() %>% sum()

# Annotate by priority order of intragenic location
## 1. Annotate with exon

exon_list <- Homo_sapiens.GRCh38.97[Homo_sapiens.GRCh38.97$type == "exon"]

STR_overlap_gtf_exon <- findOverlaps(query = STR_map_granges, subject = exon_list, 
                                     maxgap = -1, minoverlap = 1)

STR_overlap_gtf_1_w_annot_exon <- 
  STR_overlap_gtf_exon %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(subject_annot = (exon_list$gene_id %>% 
                            as.character())[subjectHits(STR_overlap_gtf_exon)])

STR_overlap_gtf_1_w_annot_no_dup_exon <- 
  STR_overlap_gtf_1_w_annot_exon %>%
  filter(!duplicated(queryHits, fromLast = F)) 

overlap_STR_map_exon <- 
  STR_overlap_gtf_exon %>% 
  as.data.frame() %>%
  as_tibble() %>%
  mutate(overlapped_genes = (STR_map_granges$HipSTR_number %>% as.character())[queryHits(STR_overlap_gtf_exon)])

annotated_STR_map_exon <- left_join(STR_map, overlap_STR_map_exon, by = c("V6" = "overlapped_genes"))

annotate_STR_map_Ensembl_exon <- annotated_STR_map_exon %>% 
  left_join(STR_overlap_gtf_1_w_annot_no_dup_exon, by = "queryHits") %>%
  dplyr::select(-(c(subjectHits.x, subjectHits.y, queryHits))) %>%
  mutate(subject_annot = ifelse(!is.na(subject_annot), "exon", "NA")) %>%
  dplyr::select(c(V6, subject_annot))

annotate_STR_map_Ensembl_exon_joined <- annotate_STR_map_Ensembl %>% 
  left_join(annotate_STR_map_Ensembl_exon, by = c("HipSTR_ID" = "V6")) %>%
  mutate(region = ifelse(subject_annot == "exon", "exon", region)) %>%
  dplyr::select(-subject_annot) %>% unique()

# exon: 48995, intergenic: 631902 
table(annotate_STR_map_Ensembl_exon_joined$region)

## 2. Annotate with 5'UTR

five_utr_list <- Homo_sapiens.GRCh38.97[Homo_sapiens.GRCh38.97$type == "five_prime_utr"]

STR_overlap_gtf_five_utr <- findOverlaps(query = STR_map_granges, subject = five_utr_list, 
                                         maxgap = -1, minoverlap = 1)

STR_overlap_gtf_1_w_annot_5utr <- 
  STR_overlap_gtf_five_utr %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(subject_annot = (five_utr_list$gene_id %>% 
                            as.character())[subjectHits(STR_overlap_gtf_five_utr)])

STR_overlap_gtf_1_w_annot_no_dup_5utr <- 
  STR_overlap_gtf_1_w_annot_5utr %>%
  filter(!duplicated(queryHits, fromLast = F)) 

overlap_STR_map_5utr <- 
  STR_overlap_gtf_five_utr %>% 
  as.data.frame() %>%
  as_tibble() %>%
  mutate(overlapped_genes = (STR_map_granges$HipSTR_number %>% as.character())[queryHits(STR_overlap_gtf_five_utr)])

annotated_STR_map_5utr <- left_join(STR_map, overlap_STR_map_5utr, by = c("V6" = "overlapped_genes"))

annotate_STR_map_Ensembl_5utr <- annotated_STR_map_5utr %>% 
  left_join(STR_overlap_gtf_1_w_annot_no_dup_5utr, by = "queryHits") %>%
  dplyr::select(-(c(subjectHits.x, subjectHits.y, queryHits))) %>%
  mutate(subject_annot = ifelse(!is.na(subject_annot), "5UTR", "NA")) %>%
  dplyr::select(c(V6, subject_annot))

# Not 5UTR overlaps with exons but in this case, let 5UTR take priority
## annotate_STR_map_Ensembl_5utr_joined <- annotate_STR_map_Ensembl_exon_joined%>% 
  # left_join(annotate_STR_map_Ensembl_5utr, by = c("HipSTR_ID" = "V6")) %>%
  # mutate(region = ifelse(!is.na(region), region, 
                         # ifelse(subject_annot == "5UTR", "5UTR", "NA"))) %>%
  # dplyr::select(-subject_annot) %>% unique()

annotate_STR_map_Ensembl_5utr_joined <- annotate_STR_map_Ensembl_exon_joined%>% 
  left_join(annotate_STR_map_Ensembl_5utr, by = c("HipSTR_ID" = "V6")) %>%
  mutate(region = ifelse(subject_annot == "5UTR", "5UTR", region)) %>%
  dplyr::select(-subject_annot) %>% unique()
  
table(annotate_STR_map_Ensembl_5utr_joined$region)

## 3. Annotate with 3'UTR

three_utr_list <- Homo_sapiens.GRCh38.97[Homo_sapiens.GRCh38.97$type == "three_prime_utr"]

STR_overlap_gtf_3utr <- findOverlaps(query = STR_map_granges, subject = three_utr_list, 
                                     maxgap = -1, minoverlap = 1)

STR_overlap_gtf_1_w_annot_3utr <- 
  STR_overlap_gtf_3utr %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(subject_annot = (three_utr_list$gene_id %>% 
                            as.character())[subjectHits(STR_overlap_gtf_3utr)])

STR_overlap_gtf_1_w_annot_no_dup_3utr <- 
  STR_overlap_gtf_1_w_annot_3utr %>%
  filter(!duplicated(queryHits, fromLast = F)) 

overlap_STR_map_3utr <- 
  STR_overlap_gtf_3utr %>% 
  as.data.frame() %>%
  as_tibble() %>%
  mutate(overlapped_genes = (STR_map_granges$HipSTR_number %>% as.character())[queryHits(STR_overlap_gtf_3utr)])

annotated_STR_map_3utr <- left_join(STR_map, overlap_STR_map_3utr, by = c("V6" = "overlapped_genes"))

annotate_STR_map_Ensembl_3utr <- annotated_STR_map_3utr %>% 
  left_join(STR_overlap_gtf_1_w_annot_no_dup_3utr, by = "queryHits") %>%
  dplyr::select(-(c(subjectHits.x, subjectHits.y, queryHits))) %>%
  mutate(subject_annot = ifelse(!is.na(subject_annot), "3UTR", "NA")) %>%
  dplyr::select(c(V6, subject_annot))

annotate_STR_map_Ensembl_3utr_joined <- annotate_STR_map_Ensembl_5utr_joined%>% 
  left_join(annotate_STR_map_Ensembl_3utr, by = c("HipSTR_ID" = "V6")) %>%
  mutate(region = ifelse(subject_annot == "3UTR", "3UTR", region)) %>%
  dplyr::select(-subject_annot) %>% unique()

# Again, 3'UTR takes over from exon
table(annotate_STR_map_Ensembl_3utr_joined$region)

## 4. Annotation with intron
### Estimating introns via splice table: can get intron annotation by working out junction coordinates

gtf <- "/data/references/ensembl/gtf_gff3/v97/Homo_sapiens.GRCh38.97.gtf"
ref <- refGenome::ensemblGenome()
refGenome::basedir(ref) <- dirname(gtf)
refGenome::read.gtf(ref, gtf %>% stringr::str_replace(".*/", "")) 
ref_junc <- refGenome::getSpliceTable(ref) 
ref_junc <- ref_junc@ev$gtf

## intron start/end = lend + 1/rstart - 1

ref_junc_intron <- ref_junc %>% 
  mutate(intron_start = lend + 1, intron_end = rstart - 1, 
         intron = "intron", intron_width = intron_end - intron_start)


## convert ref_junc to granges object for finding overlap
ref_junc_intron_granges <- GRanges(seqnames = ref_junc_intron$seqid,
                                   ranges = IRanges(start = ref_junc_intron$intron_start,
                                                    end = ref_junc_intron$intron_end),
                                   intron = ref_junc_intron$intron,
                                   intron_width = ref_junc_intron$intron_width,
                                   gene_id = ref_junc_intron$gene_id,
                                   gene_name = ref_junc_intron$gene_name,
                                   transcript_id = ref_junc_intron$transcript_id)

STR_overlap_gtf_intron <- findOverlaps(query = STR_map_granges, subject = ref_junc_intron_granges, 
                                       maxgap = -1, minoverlap = 1)

STR_overlap_gtf_1_w_annot_intron <- 
  STR_overlap_gtf_intron %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  mutate(subject_annot = (ref_junc_intron_granges$gene_id %>% 
                            as.character())[subjectHits(STR_overlap_gtf_intron)])

STR_overlap_gtf_1_w_annot_no_dup_intron <- 
  STR_overlap_gtf_1_w_annot_intron %>%
  filter(!duplicated(queryHits, fromLast = F)) 

overlap_STR_map_intron <- 
  STR_overlap_gtf_intron %>% 
  as.data.frame() %>%
  as_tibble() %>%
  mutate(overlapped_genes = (STR_map_granges$HipSTR_number %>% as.character())[queryHits(STR_overlap_gtf_intron)])

annotated_STR_map_intron <- left_join(STR_map, overlap_STR_map_intron, by = c("V6" = "overlapped_genes"))

annotate_STR_map_Ensembl_intron <- annotated_STR_map_intron %>% 
  left_join(STR_overlap_gtf_1_w_annot_no_dup_intron, by = "queryHits") %>%
  dplyr::select(-(c(subjectHits.x, subjectHits.y, queryHits))) %>%
  mutate(subject_annot = ifelse(!is.na(subject_annot), "intron", "NA")) %>%
  dplyr::select(c(V6, subject_annot)) %>% unique()

annotate_STR_map_Ensembl_intron_joined <- annotate_STR_map_Ensembl_3utr_joined%>% 
  left_join(annotate_STR_map_Ensembl_intron, by = c("HipSTR_ID" = "V6")) %>%
  mutate(region = ifelse(subject_annot == "intron" & is.na(region), "intron", region)) %>% 
  dplyr::select(-subject_annot) %>% unique()

# 3UTR: 19334, 5UTR:4855, exon: 24806, intron: 957046, intergenic: 631902
table(annotate_STR_map_Ensembl_intron_joined$region)

## How many regions are unannotated? 1002
annotate_STR_map_Ensembl_intron_joined$region %>% is.na() %>% sum()

## Finalising STR table
genomic_STR_map <- annotate_STR_map_Ensembl_intron_joined %>%
  mutate(STR_present = ifelse(nucleotides_per_repeat>= 1, 1, 0)) %>%
  mutate(STR_width = nucleotides_per_repeat*number_of_repeats)

## Adding gene symbols (1072 unique genes in total annotated)
Human_genes <- Homo_sapiens.GRCh38.97 %>% as.data.frame()%>% dplyr::select(gene_id, gene_name) %>% unique()

genomic_STR_map_symbol <- left_join(genomic_STR_map, Human_genes, 
                                    by = "gene_id")

genomic_STR_map_symbol <- genomic_STR_map_symbol %>% 
                                mutate(region = ifelse(is.na(region), "unannotated", region)) %>% unique()

genomic_STR_map_symbol$region %>% is.na() %>% sum()

table(genomic_STR_map_symbol$region)

# Forming STR Predictor ----------------------------------------------------------------------------------------------------------------

# Working out the number of STRs per gene

number_STRs <- aggregate(genomic_STR_map_symbol$STR_present, 
                         by=list(gene_id = genomic_STR_map_symbol$gene_id), FUN=sum)
colnames(number_STRs) <- c("gene_id", "number_STRs")

# Working out the STR width per gene

STR_width <- aggregate(genomic_STR_map_symbol$STR_width, 
                       by=list(gene_id = genomic_STR_map_symbol$gene_id), FUN=sum)
colnames(STR_width) <- c("gene_id", "STR_width")

# Working out nucleotides per repeat per gene

study <- genomic_STR_map_symbol$nucleotides_per_repeat %>% unique()

for(i in seq_along(study)){
  
  Number_nucleotides <- subset(genomic_STR_map_symbol, nucleotides_per_repeat == study[i])
  
  print(i)
  
  Number_nucleotides_count <- aggregate(Number_nucleotides$nucleotides_per_repeat, 
                                        by=list(gene_id = Number_nucleotides$gene_id), FUN=sum)
  
  Number_nucleotides_count_final <- Number_nucleotides_count %>% mutate(count = x/study[i]) %>% dplyr::select(-x)
  
  colnames(Number_nucleotides_count_final) = c("gene_id", study[i])
  
  if(i > 1){
        all_number_nucleotides_count_final <- full_join(all_number_nucleotides_count_final, Number_nucleotides_count_final,
                                                    by = "gene_id") 
    
  }else{
        all_number_nucleotides_count_final <- Number_nucleotides_count_final
    }
}

colnames(all_number_nucleotides_count_final) <- c("gene_id", "6nt", "3nt", "2nt", "1nt", "4nt", "5nt")

  all_number_nucleotides_count_final[is.na(all_number_nucleotides_count_final)] = 0


# Working out number of STRs per region per gene

region_study <- c("exon", "intron", "3UTR", "5UTR", "unannotated")

# loop over each region
for(i in seq_along(region_study)){
  
  total_number_STRs <- subset(genomic_STR_map_symbol, region == region_study[i])
  
  print(i)
  
  total_number_STRs_count <- aggregate(total_number_STRs$STR_present, 
                                       by=list(gene_id = total_number_STRs$gene_id), FUN=sum)
  

  colnames(total_number_STRs_count) = c("gene_id", region_study[i])
  
  if(i > 1){
    total_number_STRs_count_final <- full_join(total_number_STRs_count_final, total_number_STRs_count,
                                                    by = "gene_id") 
    
  }else{
    total_number_STRs_count_final <-  total_number_STRs_count
  }
}

colnames(total_number_STRs_count_final) <- c("gene_id", "number_STRs_in_exon", 
                                             "number_STRs_in_intron", "number_STRs_5UTR", 
                                             "number_STRs_3UTR", "number_STRs_unannotated")

total_number_STRs_count_final1 <- total_number_STRs_count_final %>% 
                                    mutate(number_STRs_5UTR = ifelse(is.na(number_STRs_5UTR), 0, number_STRs_5UTR),
                                           number_STRs_3UTR = ifelse(is.na(number_STRs_3UTR), 0, number_STRs_3UTR),
                                           number_STRs_in_intron = ifelse(is.na(number_STRs_in_intron), 0, number_STRs_in_intron),
                                           number_STRs_unannotated = ifelse(is.na(number_STRs_unannotated), 0, number_STRs_unannotated))

##generating df with width of each gene

gene_list_width <- gene_list %>% as.data.frame() %>%
                    dplyr::select(c(gene_id, gene_name, width))

## generating df with exonic length per gene

exon_list_width <- exon_list %>% as.data.frame() %>% 
                      dplyr::select(c(gene_id, gene_name, width))

exon_list_width_per_gene <- aggregate(exon_list_width$width, by=list(gene_id = exon_list_width$gene_id), FUN=sum)

colnames(exon_list_width_per_gene) <- c("gene_id", "exon_width")


## generating df with width of intron per gene

ref_junc_intron_width <- ref_junc_intron %>% 
                        dplyr::select(c(gene_id, gene_name, intron_width))

intron_list_width_per_gene <- aggregate(ref_junc_intron_width$intron_width, 
                                        by=list(gene_id = ref_junc_intron_width$gene_id), FUN=sum)

colnames(intron_list_width_per_gene) <- c("gene_id", "intron_width")


## Dataframe with gene, exon and intron widths

gene_exon_width <- as.data.frame(left_join(gene_list_width, 
                                    exon_list_width_per_gene, 
                                    by = "gene_id"))

gene_exon_intron_width <- gene_exon_width %>% 
                        left_join(intron_list_width_per_gene, by = "gene_id") %>%
                        mutate(intron_width = ifelse(is.na(intron_width), 0, intron_width))


# Linking STR properties in a single predictor

STR_predictor <- number_STRs %>% left_join(gene_exon_intron_width, by = "gene_id") %>%
                                  left_join(STR_width, by = "gene_id") %>%
                                    left_join(all_number_nucleotides_count_final, by = "gene_id") %>%
                                      left_join(total_number_STRs_count_final1, by = "gene_id")
 
STR_predictor[is.na(STR_predictor)] = 0

table(STR_predictor$number_STRs_unannotated)


# Save data -------------------------------------------------------------------------------------------

## Save predictors as STR_w_annotation
genomic_STR_map_symbol

## save_generic gene GTF DF containing gene/exon/intron widths

gene_exon_intron_width

## save STR_predictor

STR_predictor


