# Script for characterising PanelApp Ataxia genes - work out statistical differences between different gene lists -----------

# Load libraries
library(forcats)
library(ggpubr)
library(readxl)
library(rtracklayer)
library(stringr)
library(tidyverse)

# Load data -------------------------------------------------------------------------------------------

# Load gene lists and genic features
# Load all predictors
clean_gene_features <- read.csv("/home/zchen/ataxia-functional-genomics/Results/input_files/clean_gene_features_matrix.csv")

# Generate a table that labels variables as either binary or continuous

binary_classification <- clean_gene_features %>% dplyr::select(-c(gene, Panel)) %>%
  purrr::map_lgl(~all(.x %in% c(0,1, NA))) %>% 
  .[-1] %>% 
  as.data.frame() %>%  
  setNames("Binary") %>% 
  add_rownames(var = "Variable")

non_binary <- binary_classification %>% dplyr::filter(Binary == "FALSE") %>% 
  dplyr::filter(Variable != "number_STRs_unannotated") %>%
  dplyr::filter(Variable != "number_rpt_unannotated") %>%
  dplyr::filter(Variable != "number_rpt_for_intron_width")

binary <- binary_classification %>% dplyr::filter(Binary == "TRUE")

# Comparisons required

Comparison_list <- c("Adult-onset vs. controls", "Childhood-onset vs. controls", "Overlap-onset vs. controls", 
                     "Adult-onset vs. childhood-onset", "Adult-onset vs. overlap-onset", "Childhood-onset vs. overlap-onset")

Comparison_group1 <- c("Adult-onset", "Childhood-onset", "Overlap", "Adult-onset", "Adult-onset", "Childhood-onset")

Comparison_group2 <- c("Not ataxia", "Not ataxia", "Not ataxia", "Childhood-onset", "Overlap", "Overlap")

comparison_df <- data.frame(Comparison_list, Comparison_group1, Comparison_group2)


# 1. Recursive Wilcoxon rank sum test comparing between all non-binary features ----------------------------------

for(i in seq_along(non_binary$Variable)){
  for(j in seq_along(comparison_df$Comparison_list)){
    
    print(paste("i=", i))
    print(paste("j=", j))
    
    test_variable <- non_binary$Variable[i]
    Annotation <- print(test_variable)
    
    new_df_to_test <- clean_gene_features %>% dplyr::select(c(test_variable), Panel)
    
    Comparison <- comparison_df$Comparison_list[j]
    print(Comparison)
    
    compare_group1 <- new_df_to_test %>% dplyr::filter(Panel == print(comparison_df$Comparison_group1[j]))
    compare_group2 <- new_df_to_test %>% dplyr::filter(Panel == print(comparison_df$Comparison_group2[j]))
    
    test <- wilcox.test(as.numeric(unlist(compare_group1[,1])), 
                        as.numeric(unlist(compare_group2[,1])), 
                        conf.int = TRUE, conf.level = 0.95)
    
    P_value <-  test$p.value
    Lower_CI <- test$conf.int[1]
    Upper_CI <- test$conf.int[2]
    
    output <- data.frame(Comparison, Annotation, P_value, Lower_CI, Upper_CI) 
    
    if(j>1){
      output_final <- bind_rows(output_final, output)
      
    }else{
      output_final <- output
      
    }
    
      output_final$FDR_P <- p.adjust(p = output_final$P_value, method = "fdr") # correcting for all six comparisons per annotation
    
  }
  
  
  if(i>1){
    output_final_final <- bind_rows(output_final_final, output_final)
    
  }else{
    output_final_final <- output_final # output_final_final: comparison table with FDR adjustment for non-binary features
    
  }
  
}

# 2. Binary features - chi-sq comparisons -------------------------------------------------------------------------------

panel_numbers <- as.data.frame(table(clean_gene_features$Panel))
adult_total <- panel_numbers %>% dplyr::filter(Var1 == "Adult-onset")
child_total <- panel_numbers %>% dplyr::filter(Var1 == "Childhood-onset") 
overlap_total <- panel_numbers %>% dplyr::filter(Var1 == "Overlap") 
controls_total <- panel_numbers %>% dplyr::filter(Var1 == "Not ataxia") 


for(i in seq_along(binary$Variable)){
  
  test_variable <- binary$Variable[i]
  name_variable <- print(test_variable)
  
  new_df_to_test <- clean_gene_features %>% dplyr::select(c(test_variable, Panel)) 
  colnames(new_df_to_test) <- c("Variable", "Panel")
  
  # Set different matrices based on age-of-onset or control
  
  adult_onset <- new_df_to_test %>% dplyr::filter(Panel == "Adult-onset") %>% dplyr::filter(Variable == 1)
  childhood_onset <- new_df_to_test %>% dplyr::filter(Panel == "Childhood-onset") %>% dplyr::filter(Variable == 1)
  overlap_onset <- new_df_to_test %>% dplyr::filter(Panel == "Overlap") %>% dplyr::filter(Variable == 1)
  controls <- new_df_to_test %>% dplyr::filter(Panel == "Not ataxia") %>% dplyr::filter(Variable == 1)
  
  # Adult
  Adult <- c(nrow(adult_onset), adult_total$Freq - nrow(adult_onset))
  # Childhood
  Child <- c(nrow(childhood_onset), child_total$Freq - nrow(childhood_onset))
  # Overlap
  Overlap <- c(nrow(overlap_onset), overlap_total$Freq - nrow(overlap_onset))
  # Controls
  Controls <- c(nrow(controls), controls_total$Freq - nrow(controls))
  
  adult_v_control <- data.frame(Adult, Controls)
  child_v_control <- data.frame(Child, Controls)
  overlap_v_control <- data.frame(Overlap, Controls)
  adult_v_child <- data.frame(Adult, Child)
  adult_v_overlap <- data.frame(Adult, Overlap)
  child_v_overlap <- data.frame(Child, Overlap)
  
  # (i) Adult vs controls
  adult_v_control_test <- chisq.test(adult_v_control)
  adult_v_control_pval <-  adult_v_control_test$p.value
  
  # (ii) Child vs controls
  child_v_control_test <- chisq.test(child_v_control)
  child_v_control_test_pval <- child_v_control_test$p.value
  
  # (iii) Overlap vs controls
  overlap_onset_v_control_test <- chisq.test(overlap_v_control)
  overlap_onset_v_control_test_pval <- overlap_onset_v_control_test$p.value
  
  # (iv) Adult vs child
  adult_v_child_test <- chisq.test(adult_v_child)
  adult_v_child_test_pval <- adult_v_child_test$p.value
  
  # (v) Adult vs overlap
  adult_v_overlap_test <- chisq.test(adult_v_overlap)
  adult_v_overlap_test_pval <- adult_v_overlap_test$p.value
  
  # (vi) Child vs overlap
  child_v_overlap_test <- chisq.test(child_v_overlap)
  child_v_overlap_test_pval <- child_v_overlap_test$p.value
  
  Comparison_list <- c("Adult-onset vs. controls", "Childhood-onset vs. controls", "Overlap-onset vs. controls", 
                       "Adult-onset vs. childhood-onset", "Adult-onset vs. overlap-onset", "Childhood-onset vs. overlap-onset")
  
  P_value <- c(adult_v_control_pval,child_v_control_test_pval, overlap_onset_v_control_test_pval,
               adult_v_child_test_pval, adult_v_overlap_test_pval, child_v_overlap_test_pval)
  
  output_chisq <- data.frame(Comparison_list, P_value, name_variable)
  
  output_chisq$FDR_p <- p.adjust(p = output_chisq$P_value, method = "fdr")
  
  
  if(i > 1){
    output_final_chisq <- bind_rows(output_final_chisq, output_chisq)
    
  }else{
    output_final_chisq <- output_chisq
  }
}

