install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("rhdf5", "S4Vectors", "S4Arrays", "IRanges", "XVector", 
                       "GenomeInfoDb", "GenomicRanges", "DelayedArray", 
                       "SingleCellExperiment", "SparseArray", "phyloseq", 
                       "SummarizedExperiment"))
remotes::install_github("jbisanz/qiime2R")


library(qiime2R)
library(tidyverse)
library(ggsignif)
library(ggforce)
library(ggalt)
library(ggpubr)
library(rstatix)
library(vegan) #For permanova 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
metadata <- read.table(file = "sample-metadata.tsv", sep = "\t", header = T, row.names = 1,check.names = FALSE) %>%
  tibble::rownames_to_column(var = "Sample Name") %>%
  dplyr::mutate(`Day of collection` = trimws(`Day of collection`)) %>%
  dplyr::mutate(`Variant Group` = case_when(`Variant Group` == "Hongkong" ~ "Hong Kong",
                                            TRUE ~ `Variant Group`))
  

########################################## FECAL MATTER ##############################################

########################################## OBSERVED FEATURES

raw_observed_features_vector_dataframe = read_qza("diversity-core-metrics-phylogenetic/observed_features_vector.qza")

observed_features_dataframe = raw_observed_features_vector_dataframe$data %>% 
  tibble::rownames_to_column(var = "Sample Name") %>%
  left_join(. , metadata , by = "Sample Name")

#------- Significance compared to Delta --------#

# Fecal Only Day 4
#600:650

stats_observed_features_data <- observed_features_dataframe %>%
  filter(`Sample Type` == "Fecal", `Day of collection` == "DAY 4") %>%
  mutate(`Variant Group` = factor(`Variant Group`, levels = c("Uninfected", "Hong Kong", "Delta", "Omicron"))) %>%
  dplyr::mutate(Variant_Name = `Variant Group`) %>%
  t_test(observed_features ~ Variant_Name) %>%
  add_significance() %>% add_xy_position(x = "Variant_Name") %>%
  dplyr::filter(group1 == "Uninfected" | group1 == "Delta" | group2 == "Delta") %>%
  dplyr::mutate(y.position = c(198 , 201 , 203 , 198 , 198)) #%>%
  #dplyr::filter(group1 == "Delta" | group2 == "Delta")

ggplot(observed_features_dataframe %>%
         dplyr::filter(`Sample Type` == "Fecal") %>%
         dplyr::filter(`Day of collection` == "DAY 4") %>%
         dplyr::mutate(`Variant Group` = factor(`Variant Group` , levels = c("Uninfected" , "Hong Kong" , "Delta" , "Omicron")))) +
  geom_boxplot(aes(x = `Variant Group`, y = observed_features , fill = `Variant Group`)) +
  theme_bw() +
  labs(x = "Variant", y = "Observed Features") +
  scale_fill_manual(values = c("Uninfected" = "black","Hong Kong" = "red","Delta" = "purple","Omicron" = "green")) +
  theme(axis.title = element_text(size = 20 , face = "bold") , 
        strip.text = element_text(size = 20 , face = "bold") , 
        axis.text.x = element_text(size = 20 , angle = 45, hjust = 1, face = "bold" , colour = "black"),
        axis.text.y = element_text(size = 15 , face = "bold" ,colour = "black"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 15)) +
#stat_pvalue_manual(stats_observed_features_data, label = "p.adj.signif", tip.length = 0.01)
stat_pvalue_manual(stats_observed_features_data, label = "p", tip.length = 0.01,fontface = "bold")

########################################## FAITH'S PHYLOGENETIC DIVERSITY

raw_faith_pd_vector_dataframe = read_qza("diversity-core-metrics-phylogenetic/faith_pd_vector.qza")

faith_pd_vector_dataframe = raw_faith_pd_vector_dataframe$data %>% 
  dplyr::rename(`Sample Name` = V1, 
                faith_pd = V2) %>%
  left_join(. , metadata , by = "Sample Name")

stats_faith_pd_features_data <- faith_pd_vector_dataframe %>%
  filter(`Sample Type` == "Fecal", `Day of collection` == "DAY 4") %>%
  mutate(`Variant Group` = factor(`Variant Group`, levels = c("Uninfected", "Hong Kong", "Delta", "Omicron"))) %>%
  dplyr::mutate(Variant_Name = `Variant Group`) %>%
  t_test(faith_pd ~ Variant_Name) %>%
  add_significance() %>% add_xy_position(x = "Variant_Name") %>%
  dplyr::filter(group1 == "Uninfected" | group1 == "Delta" | group2 == "Delta") %>%
  dplyr::mutate(y.position = c(39.5 , 40.3 , 41.2 , 39.5 , 39.5)) #%>%
#dplyr::filter(group1 == "Delta" | group2 == "Delta")

# Fecal Only Day 4
#600:650
ggplot(faith_pd_vector_dataframe %>% dplyr::filter(`Sample Type` == "Fecal") %>%
         dplyr::filter(`Day of collection` == "DAY 4") %>%
         dplyr::mutate(`Variant Group` = factor(`Variant Group` , levels = c("Uninfected" , "Hong Kong" , "Delta" , "Omicron")))) +
  geom_boxplot(aes(x = `Variant Group`, y = faith_pd, fill = `Variant Group`)) +
  #geom_point(aes(x = `Variant Group`, y = faith_pd)) +
  theme_bw() +
  labs(x = "Variant", y = "Faith's Phylogenetic Diversity") +
  scale_fill_manual(values = c("Uninfected" = "black","Hong Kong" = "red","Delta" = "purple","Omicron" = "green")) +
  theme(axis.title = element_text(size = 20 , face = "bold") , 
        strip.text = element_text(size = 20 , face = "bold") , 
        axis.text.x = element_text(size = 20 , angle = 45, hjust = 1, face = "bold" , colour = "black"),
        axis.text.y = element_text(size = 15 , face = "bold" , colour = "black"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 15)) +
  #stat_pvalue_manual(stats_faith_pd_features_data, label = "p.adj.signif", tip.length = 0.01)
  stat_pvalue_manual(stats_faith_pd_features_data, label = "p", tip.length = 0.01,fontface = "bold")


########################################## BETA DIVERSITY ##################################################

########################################## JACCARD

raw_jaccard_pcoa_results_dataframe = read_qza("diversity-core-metrics-phylogenetic/jaccard_pcoa_results.qza")

jaccard_pcoa_results_dataframe = raw_jaccard_pcoa_results_dataframe$data$Vectors %>% 
  dplyr::rename(`Sample Name` = "SampleID") %>%
  left_join(. , metadata , by = "Sample Name") 

jaccard_pcoa_fecal_day4_dataframe = jaccard_pcoa_results_dataframe %>% 
  dplyr::filter(`Sample Type` == "Fecal") %>%
  dplyr::filter(`Day of collection` == "DAY 4") %>%
  dplyr::mutate(`Variant Group` = factor(`Variant Group` , levels = c("Uninfected" , "Hong Kong" , "Delta" , "Omicron")))

jaccard_pcoa_coords <- jaccard_pcoa_fecal_day4_dataframe %>% 
  dplyr::select(PC1, PC2)

jaccard_grouping <- jaccard_pcoa_fecal_day4_dataframe$`Variant Group`

adonis2(jaccard_pcoa_coords ~ jaccard_grouping, data = jaccard_pcoa_fecal_day4_dataframe, method = "euclidean", permutations = 999)

#1500:600
ggplot(jaccard_pcoa_results_dataframe %>% dplyr::filter(`Sample Type` == "Fecal") %>%
         dplyr::filter(`Day of collection` == "DAY 4") %>%
         dplyr::mutate(`Variant Group` = factor(`Variant Group` , levels = c("Uninfected" , "Hong Kong" , "Delta" , "Omicron"))) , 
       aes(x = PC1, y = PC2,color = `Variant Group`)) +
  geom_point(size = 2.5) +
  geom_mark_ellipse(aes(fill = `Variant Group`),expand = 0.05, alpha = 0.2) +
  theme_bw() +
  labs(x = "PCoA1" , y = "PCoA2") +
  theme(axis.title = element_text(size = 25 , face = "bold") , 
        strip.text = element_text(size = 15 , face = "bold") , 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20 , face = "bold"),
        panel.grid = element_blank(),
        legend.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 20 , face = "bold"),
        legend.key.size = unit(1, 'cm')) +
  expand_limits(x = c(min(0, min(jaccard_pcoa_results_dataframe$PC1) * 0.5),
                      max(jaccard_pcoa_results_dataframe$PC1) * 1.5)) +
  expand_limits(y = c(min(0, min(jaccard_pcoa_results_dataframe$PC2) * 0.5),
                      max(jaccard_pcoa_results_dataframe$PC2) * 1.5))
                
############################################### TABLES FOR PUBLICATION #######################################

observed_features_dataframe %>%
  dplyr::filter(`Sample Type` == "Fecal") %>%
  dplyr::filter(`Day of collection` == "DAY 4") %>%
  dplyr::select(`Variant Group` , `Sample Type` , `Day of collection` , observed_features) %>%
  write.table(. , "/media/user/4TB/Rishad/Microbiome_Analysis_Santhosh/Scripts/Observed_Features_Table_For_Publication.tsv" , sep = "\t" ,row.names = FALSE , col.names = TRUE)

faith_pd_vector_dataframe %>% 
  dplyr::filter(`Sample Type` == "Fecal") %>%
  dplyr::filter(`Day of collection` == "DAY 4") %>%
  dplyr::select(`Variant Group` , `Sample Type` , `Day of collection` , faith_pd) %>%
  write.table(. , "/media/user/4TB/Rishad/Microbiome_Analysis_Santhosh/Scripts/Faiths_pd_Table_For_Publication.tsv" , sep = "\t" ,row.names = FALSE , col.names = TRUE)


jaccard_pcoa_results_dataframe %>% 
  dplyr::filter(`Sample Type` == "Fecal") %>%
  dplyr::filter(`Day of collection` == "DAY 4") %>%
  dplyr::select(`Variant Group` , `Sample Type` , `Day of collection` , PC1 , PC2) %>%
  write.table(. , "/media/user/4TB/Rishad/Microbiome_Analysis_Santhosh/Scripts/Jaccard_PCOA_Table_For_Publication.tsv" , sep = "\t" ,row.names = FALSE , col.names = TRUE)


