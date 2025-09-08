# In the beginning itself a small number of 0.1 is added to the number

BiocManager::install("ComplexHeatmap")
BiocManager::install("lefser")
BiocManager::install("SummarizedExperiment")

library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(ComplexHeatmap)
library(purrr)
library(rstatix)
library(SummarizedExperiment)
library(lefser)
library(cowplot)

# For Colour Purpose
qual_col_pals = brewer.pal.info[brewer.pal.info$colorblind ==TRUE , ] %>% .[row.names(.) != "Greys" , ]
col_vector = unlist(mapply(brewer.pal , qual_col_pals$maxcolors , rownames(qual_col_pals)))


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

level_7_dataframe = read.csv(file = "level-7.csv" , check.names = FALSE) %>%
  dplyr::mutate(`Variant Group` = case_when(`Variant Group` == "Hongkong" ~ "Hong Kong",
                                            TRUE ~ `Variant Group`))

Full_species_classification = level_7_dataframe %>% 
  dplyr::select(-c("index" , "Variant Group" , "Gender" , "No. of animals" , "Sample Type" , "Day of collection" , "DNA conc (ng/uL)")) %>%
  names()

classification_details_dataframe = Full_species_classification %>%
  strsplit(. , ";") %>% # All elements are of the same number
  do.call(`rbind`,.) %>%
  as.data.frame()

classification_details_dataframe = setNames(classification_details_dataframe , c("Division" , "Phylum" , "Class" , "Order" , "Family","Genus" , "Species"))

gather_group_dataframe = level_7_dataframe %>% 
  filter(`Sample Type` == "Fecal",
         `Day of collection` == "DAY 4") %>%
  select(-c(Gender, `No. of animals`, `DNA conc (ng/uL)`, 
            `Sample Type`, `Day of collection`)) %>%
  dplyr::rename(Variant_Group = `Variant Group`) %>% 
  gather(key = "Group" , value = "Number", all_of(Full_species_classification)) %>% 
  separate(Group , into = c("Division" , "Phylum" , "Class" , "Order" , "Family","Genus" , "Species"), sep = ";",remove = FALSE) %>%
  dplyr::mutate(Division = gsub("k__" , "" , Division) , 
                Division = gsub("_" , "" , Division), 
                Division = gsub("^$" , NA , Division)) %>%
  dplyr::mutate(Phylum = gsub("p__" , "" , Phylum) , 
                Phylum = gsub("_" , "" , Phylum), 
                Phylum = gsub("^$" , NA , Phylum)) %>%
  dplyr::mutate(Class = gsub("c__" , "" , Class) , 
                Class = gsub("_" , "" , Class), 
                Class = gsub("^$" , NA , Class)) %>%
  dplyr::mutate(Order = gsub("o__" , "" , Order) , 
                Order = gsub("_" , "" , Order), 
                Order = gsub("^$" , NA , Order)) %>%
  dplyr::mutate(Family = gsub("f__" , "" , Family) , 
                Family = gsub("_" , "" , Family), 
                Family = gsub("^$" , NA , Family)) %>%
  dplyr::mutate(Genus = gsub("g__" , "" , Genus) , 
                Genus = gsub("_" , "" , Genus), 
                Genus = gsub("^$" , NA , Genus)) %>%
  dplyr::mutate(Species = gsub("s__" , "" , Species) , 
                Species = gsub("_" , "" , Species), 
                Species = gsub("^$" , NA , Species)) %>%
  rowwise() %>%
  mutate(New_group_Name = paste(c_across(c(Division, Phylum, Class, Order, Family, Genus, Species))[!is.na(c_across(c(Division, Phylum, Class, Order, Family, Genus, Species)))], collapse = "|")) %>%
  ungroup() %>%
  dplyr::select(Variant_Group , index , New_group_Name , Number) %>%
  group_by(Variant_Group, index, New_group_Name) %>%
  summarise(Number = mean(Number), .groups = "drop")



# a helper function for pairwise LEfSe comparison

run_lefse_pairwise <- function( group1, group2) {
  df_clean <- gather_group_dataframe %>% 
    pivot_wider(names_from = New_group_Name , values_from = Number) %>%
    dplyr::filter(Variant_Group %in% c(group1 , group2))
  
  # assume df_clean is after your filtering/renaming
  features <- df_clean %>% 
    select(-Variant_Group, -index) %>% 
    t() %>%  # transpose: now rows = features, cols = samples
    as.matrix()
  
  metadata <- df_clean %>% select(Variant_Group, index)
  # build SummarizedExperiment
  se <- SummarizedExperiment(assays = list(counts = features),
                             colData = metadata)
  
  # now lefser works
  res <- lefser(se, classCol = "Variant_Group",  "index")
  
  p <- lefserPlot(res,trim.names = FALSE) + theme_bw() +
    labs(fill = "Variants") +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = "right",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.title = element_text(size = 15,face = "bold"),
          legend.text = element_text(size = 15))
  
  return(p)
}

# Run pairwise comparisons
p1 <- run_lefse_pairwise("Uninfected", "Delta")
p2 <- run_lefse_pairwise("Hong Kong", "Delta")
p3 <- run_lefse_pairwise("Delta", "Omicron")
p4 <- run_lefse_pairwise("Uninfected", "Hong Kong")
p5 <- run_lefse_pairwise("Uninfected", "Omicron")

# Combine all plots vertically
plot_grid(p1, p2, p3, ncol = 1, align = "h")

plot_grid(p1, p2, p3, p4, p5, ncol = 1, align = "h")


#------------- DATAFRAMES -------------------#

# a helper function for pairwise LEfSe comparison

run_dataframe_lefse_pairwise <- function( group1, group2) {
  df_clean <- gather_group_dataframe %>% 
    pivot_wider(names_from = New_group_Name , values_from = Number) %>%
    dplyr::filter(Variant_Group %in% c(group1 , group2))
  
  # assume df_clean is after your filtering/renaming
  features <- df_clean %>% 
    select(-Variant_Group, -index) %>% 
    t() %>%  # transpose: now rows = features, cols = samples
    as.matrix()
  
  metadata <- df_clean %>% select(Variant_Group, index)
  # build SummarizedExperiment
  se <- SummarizedExperiment(assays = list(counts = features),
                             colData = metadata)
  
  # now lefser works
  res <- lefser(se, classCol = "Variant_Group",  "index") %>%
    dplyr::mutate(Group1 = group1 , 
                  Group2 = group2)
  
  return(res)
}


list(run_dataframe_lefse_pairwise("Uninfected", "Delta"),
     run_dataframe_lefse_pairwise("Hong Kong", "Delta"),
     run_dataframe_lefse_pairwise("Delta", "Omicron")) %>%
  do.call('rbind' , .) %>%
  write.table(. , "lefse_analysis_result_dataframe.tsv" , sep = "\t" , quote = FALSE , row.names = FALSE , col.names = TRUE)
