library(tidyverse)
library(XML)
library(purrr)
library(patchwork)
library(grid)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################################################## EXPRESSION DATA #################################################
# Colon is renamed as large intestine
EXPRESSION_DATAFRAME = read.csv("rna_tissue_consensus.tsv" , header = TRUE , sep = "\t") %>%
  dplyr::filter(Gene.name %in% c("ACE2" , "TMPRSS2")) %>%
  dplyr::filter(Tissue %in% c("lung" , "small intestine" , "stomach" , "colon")) %>%
  dplyr::mutate(Tissue = str_to_title(Tissue))
#dplyr::mutate(Tissue = gsub("^Colon$" , "Large Intestine" , Tissue)) %>%

####################################################### PROTEIN ABUNDANCE DATA #########################################################

string_id_conversion_guide <- read.delim("9606.protein.info.v12.0.txt",
                                         quote = "",
                                         row.names = NULL,
                                         stringsAsFactors = FALSE) %>%
  dplyr::select(-c(annotation)) %>%
  dplyr::rename(String_ID = "X.string_protein_id")   

organs = c("Nasopharynx" , "Lung" ,"Stomach", "Small_Intestine" , "Colon")

protein_abundance_dataframe_list = lapply(c(1:length(organs)), function(i) { 
  
  required_file_names = list.files("9606/" , pattern = organs[i], ignore.case = TRUE)
  
  lapply(required_file_names, function(x){
    dataframes = read.csv(paste0("9606/" , x) , header = FALSE , comment.char = "#", sep =  "\t") 
    if (ncol(dataframes) == 3) {
      names(dataframes) = c("String_ID" , "Abundance" , "Raw_Spectral_Count")
    } else if (ncol(dataframes) == 2) {
      names(dataframes) = c("String_ID" , "Abundance")
    }
    return(dataframes %>%
             dplyr::mutate(Study = x) %>%
             dplyr::mutate(Organ = organs[i]) %>%
             left_join(. , string_id_conversion_guide) %>%
             dplyr::filter(preferred_name %in% c("ACE2" , "TMPRSS2"))) 
  })
}) %>% 
  unlist(., recursive = FALSE) 

protein_abundance_dataframe <- protein_abundance_dataframe_list %>%
  bind_rows()

# 700:400

lapply(c("ACE2" , "TMPRSS2"), function(i) {
  
  EXPRESSION_PLOT = ggplot(EXPRESSION_DATAFRAME  %>% dplyr::filter(Gene.name == i) %>%  
                             dplyr::mutate(Tissue = factor(Tissue , c("Colon" , "Small Intestine" , "Stomach","Lung"))), aes(x = nTPM , y = Tissue)) +
    geom_bar(stat = "identity", color="black", fill = "blue",width = 0.5) +
    labs(x = "Normalised\nTranscipt per Million\n(nTPM)" , y = "Tissue") +
    theme_minimal() +
    scale_x_reverse() +
    theme(legend.key.size = unit(0.75, 'cm'),
          axis.text.y=element_blank(),
          axis.text.x = element_text(size=16,face = "bold", colour = "black") ,
          legend.text = element_text(size=20),
          axis.title.y=element_blank(),
          axis.title=element_text(size=15,face="bold"),
          axis.line.x = element_line(color = "black", size = 0.5), 
          axis.ticks.x = element_line(color = "black", size = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title = element_text(size=25,face = "bold"),
          strip.text.y = element_text(size=20, face="bold")) +
    geom_bar(data = . %>% dplyr::filter(Tissue == "Lung") , stat = "identity", fill = "red",width = 0.5)
  
  PROTEIN_ABUNDANCE_PLOT = ggplot(protein_abundance_dataframe %>% dplyr::filter(preferred_name == i) %>%
                                    dplyr::mutate(Organ = gsub("Small_Intestine" , "Small Intestine" , Organ)) %>% 
                                    group_by(Organ) %>%
                                    summarise(Abundance_mean = mean(Abundance) , 
                                              Abundance_SD = sd(Abundance),
                                              Abundance_SEM = sd(Abundance)/sqrt(n())) %>%
                                    dplyr::mutate(Organ = factor(Organ , c("Colon" , "Small Intestine" , "Stomach","Lung"))) %>%
                                    dplyr::mutate(Organ_padded = str_pad(Organ, width = max(nchar(protein_abundance_dataframe$Organ)), side = "both")) %>%
                                    dplyr::arrange(Organ) %>%
                                    mutate(Organ_padded = factor(Organ_padded, levels = .$Organ_padded)), aes(x = Abundance_mean , y = Organ_padded)) +
    geom_bar(stat = "identity", color="black", color="black",fill = "blue",width = 0.5) +
    labs(x = "Abundance\nParts per Million\n(ppm)" , y = "Tissue") +
    scale_x_continuous(breaks = seq(0, 45, 15)) +
    expand_limits(x = c(0, 45)) +
    theme_minimal() +
    theme(legend.key.size = unit(0.75, 'cm'),
          axis.text.x = element_text(size=16,face = "bold" , colour = "black") ,
          axis.text.y = element_text(size=16,face="bold",colour = "black"),
          legend.text = element_text(size=20),
          axis.title.y=element_blank(),
          axis.title=element_text(size=16,face="bold"),
          axis.line.x = element_line(color = "black", size = 0.5), 
          axis.ticks.x = element_line(color = "black", size = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title = element_text(size=25,face = "bold"),
          strip.text.y = element_text(size=20, face="bold")) +
    geom_bar(data = . %>% dplyr::filter(Organ == "Lung") , stat = "identity", fill = "red",width = 0.5) 
  
  EXPRESSION_PLOT + PROTEIN_ABUNDANCE_PLOT 
} )


DATAFRAME_1 = EXPRESSION_DATAFRAME %>% dplyr::filter(Gene.name %in% c("ACE2" , "TMPRSS2")) %>%
  dplyr::rename(Gene_Name = Gene.name , 
                Normalised_Transcripts_Per_Million = nTPM , 
                Organ = Tissue) %>%
  dplyr::select(Gene_Name , Normalised_Transcripts_Per_Million , Organ)

DATAFRAME_2 = protein_abundance_dataframe %>% dplyr::filter(preferred_name %in% c("ACE2" , "TMPRSS2")) %>%
  dplyr::mutate(Organ = gsub("Small_Intestine" , "Small Intestine" , Organ)) %>% 
  group_by(Organ ,preferred_name) %>%
  summarise(Mean_Protein_Abundance = mean(Abundance) , 
            Abundance_SD = sd(Abundance),
            Abundance_SEM = sd(Abundance)/sqrt(n())) %>%
  dplyr::rename(Gene_Name = preferred_name) %>%
  dplyr::select(Gene_Name , Organ, Mean_Protein_Abundance)

full_join(DATAFRAME_1 , DATAFRAME_2) %>% write.table(. , "ACE2_TMPRSS2_datatable_export.tsv" , sep = "\t" , quote = FALSE , row.names = FALSE , col.names = TRUE)
