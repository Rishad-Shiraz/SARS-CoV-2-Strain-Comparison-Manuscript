# In the beginning itself a small number of 0.1 is added to the number

BiocManager::install("ComplexHeatmap")
BiocManager::install("microbiome")

library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(ComplexHeatmap)
library(purrr)
library(microbiome)
library(UpSetR)

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

gather_group_dataframe = gather(level_7_dataframe , key = "Group" , value = "Number", all_of(Full_species_classification)) %>% 
  dplyr::mutate(Number = Number+0.1) %>% # A small number of 0.1 is added to the number to make it divisible at all times
  dplyr::filter(`Day of collection` == "DAY 4") %>%
  dplyr::mutate(`Day of collection` = case_when(`Variant Group` == "Uninfected" ~ "DAY 0 (preinfection)" , 
                                                TRUE ~ `Day of collection`)) %>%
  dplyr::mutate(`Variant Group` = case_when(`Day of collection` == "DAY 0 (preinfection)" ~ "Uninfected" , 
                                            TRUE ~ `Variant Group`)) %>%
  group_by(`Variant Group` , `Sample Type` , `Day of collection`, Group) %>%
  summarise(Mean_Number = mean(Number)) %>% #When we have multiple samples for the same group, we take average value in those case
  separate(Group , into = c("Division" , "Phylum" , "Class" , "Order" , "Family","Genus" , "Species"), sep = ";",remove = FALSE) %>%
  dplyr::mutate(Division = gsub("k__" , "" , Division) , 
                Division = gsub("_" , "" , Division), 
                Division = gsub("^$" , "Other" , Division)) %>%
  dplyr::mutate(Phylum = gsub("p__" , "" , Phylum) , 
                Phylum = gsub("_" , "" , Phylum), 
                Phylum = gsub("^$" , "Other" , Phylum)) %>%
  dplyr::mutate(Class = gsub("c__" , "" , Class) , 
                Class = gsub("_" , "" , Class), 
                Class = gsub("^$" , "Other" , Class)) %>%
  dplyr::mutate(Order = gsub("o__" , "" , Order) , 
                Order = gsub("_" , "" , Order), 
                Order = gsub("^$" , "Other" , Order)) %>%
  dplyr::mutate(Family = gsub("f__" , "" , Family) , 
                Family = gsub("_" , "" , Family), 
                Family = gsub("^$" , "Other" , Family)) %>%
  dplyr::mutate(Genus = gsub("g__" , "" , Genus) , 
                Genus = gsub("_" , "" , Genus), 
                Genus = gsub("^$" , "Other" , Genus)) %>%
  dplyr::mutate(Species = gsub("s__" , "" , Species) , 
                Species = gsub("_" , "" , Species), 
                Species = gsub("^$" , "Other" , Species)) %>%
  dplyr::mutate(Species = if_else(Species == "Other", Genus, Species)) %>%
  dplyr::mutate(Species = if_else(Species == "Other", Family, Species)) %>%
  dplyr::mutate(Species = if_else(Species == "Other", Order, Species)) %>%
  dplyr::mutate(Species = if_else(Species == "Other", Class, Species)) %>%
  group_by(`Variant Group` , `Sample Type` , `Day of collection`) %>%
  dplyr::mutate(Sum_Mean_Group_Number = sum(Mean_Number))


################################################# Species Based Plots #################################################3

group_infected_phylum_species_dataframe = gather_group_dataframe %>%
  dplyr::filter(`Variant Group` != "Uninfected") %>%
  group_by(`Variant Group` , `Sample Type` , `Day of collection`, Division , Phylum, Class , Order , Family , Genus, Species) %>%
  dplyr::mutate(Sum_Mean_Species_infected_Number = sum(Mean_Number)) %>% 
  dplyr::mutate(Percentage_infected_Species = Sum_Mean_Species_infected_Number/Sum_Mean_Group_Number*100) %>%
  as.data.frame(.) %>%
  dplyr::select(-c("Mean_Number" , "Sum_Mean_Group_Number"))

group_uninfected_phylum_species_dataframe = gather_group_dataframe %>%
  dplyr::filter(`Variant Group` == "Uninfected") %>%
  group_by(`Variant Group` , `Sample Type` , `Day of collection`, Division , Phylum, Class , Order , Family , Genus, Species) %>%
  dplyr::mutate(Sum_Mean_Species_uninfected_Number = sum(Mean_Number)) %>% 
  dplyr::mutate(Percentage_uninfected_Species = Sum_Mean_Species_uninfected_Number/Sum_Mean_Group_Number*100) %>%
  as.data.frame(.) %>%
  dplyr::select(-c("Mean_Number" , "Sum_Mean_Group_Number")) %>%
  dplyr::select(-c("Variant Group" , "Day of collection"))

comparison_infected_uninfected_species_dataframe = full_join(group_infected_phylum_species_dataframe , group_uninfected_phylum_species_dataframe) %>%
  dplyr::mutate(Percentage_Compared_Uninfected = Percentage_infected_Species/Percentage_uninfected_Species) %>%
  dplyr::mutate(log2_percentage_compared_uninfected = log2(Percentage_Compared_Uninfected)) %>%
  dplyr::mutate(log10_percentage_compared_uninfected = log10(Percentage_Compared_Uninfected)) %>%
  dplyr::distinct() 


##################################### Generating significance data ##################################

gather_individual_group_dataframe = gather(level_7_dataframe , key = "Group" , value = "Number", all_of(Full_species_classification)) %>% 
  dplyr::mutate(Number = Number+0.1) %>% # A small number of 0.1 is added to the number to make it divisible at all times
  dplyr::filter(`Day of collection` == "DAY 4") %>%
  dplyr::mutate(`Day of collection` = case_when(`Variant Group` == "Uninfected" ~ "DAY 0 (preinfection)" , 
                                                TRUE ~ `Day of collection`)) %>%
  dplyr::mutate(`Variant Group` = case_when(`Day of collection` == "DAY 0 (preinfection)" ~ "Uninfected" , 
                                            TRUE ~ `Variant Group`)) %>%
  separate(Group , into = c("Division" , "Phylum" , "Class" , "Order" , "Family","Genus" , "Species"), sep = ";",remove = FALSE) %>%
  dplyr::mutate(Division = gsub("k__" , "" , Division) , 
                Division = gsub("_" , "" , Division), 
                Division = gsub("^$" , "Other" , Division)) %>%
  dplyr::mutate(Phylum = gsub("p__" , "" , Phylum) , 
                Phylum = gsub("_" , "" , Phylum), 
                Phylum = gsub("^$" , "Other" , Phylum)) %>%
  dplyr::mutate(Class = gsub("c__" , "" , Class) , 
                Class = gsub("_" , "" , Class), 
                Class = gsub("^$" , "Other" , Class)) %>%
  dplyr::mutate(Order = gsub("o__" , "" , Order) , 
                Order = gsub("_" , "" , Order), 
                Order = gsub("^$" , "Other" , Order)) %>%
  dplyr::mutate(Family = gsub("f__" , "" , Family) , 
                Family = gsub("_" , "" , Family), 
                Family = gsub("^$" , "Other" , Family)) %>%
  dplyr::mutate(Genus = gsub("g__" , "" , Genus) , 
                Genus = gsub("_" , "" , Genus), 
                Genus = gsub("^$" , "Other" , Genus)) %>%
  dplyr::mutate(Species = gsub("s__" , "" , Species) , 
                Species = gsub("_" , "" , Species), 
                Species = gsub("^$" , "Other" , Species)) %>%
  dplyr::mutate(Species = if_else(Species == "Other", Genus, Species)) %>%
  dplyr::mutate(Species = if_else(Species == "Other", Family, Species)) %>%
  dplyr::mutate(Species = if_else(Species == "Other", Order, Species)) %>%
  dplyr::mutate(Species = if_else(Species == "Other", Class, Species)) %>%
  group_by(index , `Variant Group` , `Sample Type` , `Day of collection`, Gender) %>%
  dplyr::mutate(Sum_Group_Number = sum(Number))


################################################# Species Based Plots #################################################3

lung_viral_load_dataframe = read.csv("Lung_Viral_copy_number.csv" , header = TRUE , sep = ",",check.names = FALSE) %>%
  dplyr::select("Variant Group" , "Sample ID" , "Viral_Load_Lung" , "Viral_Load_Small_Intestine" , "Organ") %>%
  dplyr::mutate(`Sample ID` = gsub(" " , "" , `Sample ID`)) %>%
  dplyr::rename(index = `Sample ID`) %>%
  dplyr::mutate(`Variant Group` = case_when(`Variant Group` == "Hongkong" ~ "Hong Kong",
                                            TRUE ~ `Variant Group`))

############################### FECAL #######################################

group_bacteria_viral_load_dataframe = gather_individual_group_dataframe %>%
  dplyr::mutate(group_new_name = paste0(Division ,";", Phylum ,";", Class ,";", Order ,";", Family ,";", Genus ,";", Species)) %>%
  group_by(index , group_new_name) %>%
  dplyr::mutate(Sum_Species_Number = sum(Number)) %>% 
  dplyr::mutate(Percentage_Species = Sum_Species_Number/Sum_Group_Number*100) %>%
  as.data.frame(.) %>%
  dplyr::select(-Number) %>%
  dplyr::distinct() %>%
  left_join(. , lung_viral_load_dataframe)

small_intestine_group_bacteria_viral_load_dataframe = group_bacteria_viral_load_dataframe %>%
  dplyr::filter(Organ == "Fecal")

small_intestine_group_bacteria_viral_load_dataframe_list = split(small_intestine_group_bacteria_viral_load_dataframe, interaction(small_intestine_group_bacteria_viral_load_dataframe$`Variant Group`, small_intestine_group_bacteria_viral_load_dataframe$group_new_name))

lung_spearman_dataframe = lapply(small_intestine_group_bacteria_viral_load_dataframe_list , function(x) {
  spearman_statistic_lung = cor.test(
    x$Percentage_Species, 
    x$Viral_Load_Lung, 
    method = "spearman",
    exact = TRUE
  )
  spearman_statistic_small_intestine = cor.test(
    x$Percentage_Species, 
    x$Viral_Load_Small_Intestine, 
    method = "spearman",
    exact = TRUE
  )
  x %>%
    dplyr::mutate(Lung_spearman_statistic_rho = spearman_statistic_lung$estimate,
                  Lung_spearman_statistic_p_value = spearman_statistic_lung$p.value)
}) %>%
  do.call('rbind' , .)

small_intestine_spearman_dataframe = lapply(small_intestine_group_bacteria_viral_load_dataframe_list , function(x) {
  spearman_statistic_small_intestine = cor.test(
    x$Percentage_Species, 
    x$Viral_Load_Small_Intestine, 
    method = "spearman",
    exact = TRUE
  )
  x %>%
    dplyr::mutate(Small_Intestine_spearman_statistic_rho = spearman_statistic_small_intestine$estimate,
                  Small_Intestine_spearman_statistic_p_value = spearman_statistic_small_intestine$p.value)
}) %>%
  do.call('rbind' , .)


spearman_dataframe = full_join(lung_spearman_dataframe %>%
                                 dplyr::filter(Lung_spearman_statistic_p_value <= 0.05), 
                               small_intestine_spearman_dataframe %>%
                                 dplyr::filter(Small_Intestine_spearman_statistic_p_value <= 0.05))

spearman_filtered_dataframe = spearman_dataframe %>% 
  dplyr::filter(`Day of collection` == "DAY 4") %>%
  dplyr::mutate(Lung_spearman_statistic_rho = case_when(Lung_spearman_statistic_p_value > 0.05 ~ 0 , 
                                                        TRUE ~ Lung_spearman_statistic_rho)) %>%
  dplyr::mutate(Small_Intestine_spearman_statistic_rho = case_when(Small_Intestine_spearman_statistic_p_value > 0.05 ~ 0 , 
                                                                   TRUE ~ Small_Intestine_spearman_statistic_rho)) %>%
  
  dplyr::rename(`Lung vRNA` = Lung_spearman_statistic_rho , 
                `Small Intestine vRNA` = Small_Intestine_spearman_statistic_rho) %>%
  
  pivot_longer(
    cols = c(`Lung vRNA`, `Small Intestine vRNA`),  # Specify columns to gather
    names_to = "Names",                # Name for the new column that indicates the original column name
    values_to = "Spearman Correlation" # Name for the new column containing the values
  )

# 1200:600
ggplot(spearman_filtered_dataframe %>%
         dplyr::filter(!is.na(`Spearman Correlation`)) %>%
         dplyr::mutate(`Variant Group` = factor(`Variant Group` , levels = c("Hong Kong" , "Delta","Omicron")))) +
  geom_tile(aes(x = `Variant Group`, y = Species , fill = `Spearman Correlation`) , colour = "black" , size = 0.1) +
  theme_set(theme_bw()) +
  theme(axis.text=element_text(size=8,colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title=element_text(size=14,face="bold") ,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1, 'cm') ,
        legend.title = element_text(size=12 , face = "bold"),
        legend.text = element_text(size=12) , 
        strip.text.x = element_text(size = 12 , face = "bold"),
        strip.text.y = element_text(size = 10 , angle = 0 , face = "bold"),
        strip.background.y = element_blank(),
        legend.spacing.y = unit(0, 'cm'),legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm')) +
  labs(y = "Smallest Taxonomic Group\n(Species/Genus/Family/Order/Class)\n",fill = "Spearman\nCorrelation") +
  #facet_wrap(~`Day of collection` , scales = "free_x") +
  #facet_wrap(~Phylum , scales = "free_y")
  facet_grid(Phylum ~ Names , scales = "free_y",space = "free_y") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", 
                       midpoint = 0,
                       labels = scales::label_number(scilimits = c(0, 0)))
