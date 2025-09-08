# In the beginning itself a small number of 0.1 is added to the number

BiocManager::install("ComplexHeatmap")

library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(ComplexHeatmap)
library(purrr)
library(rstatix)

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

group_phylum_species_dataframe = gather_individual_group_dataframe %>%
  group_by(index , `Variant Group` , `Sample Type` , `Day of collection`, Division , Phylum, Class , Order , Family , Genus, Species) %>%
  dplyr::mutate(Sum_Species_Number = sum(Number)) %>% 
  dplyr::mutate(Percentage_Species = Sum_Species_Number/Sum_Group_Number*100) %>%
  as.data.frame(.) %>%
  dplyr::mutate(group_new_name = paste0(Division ,";", Phylum ,";", Class ,";", Order ,";", Family ,";", Genus ,";", Species))

anova_dataframe = lapply(unique(group_phylum_species_dataframe$group_new_name), function(x){
  group_phylum_species_dataframe %>% 
    dplyr::filter(`Sample Type` == "Fecal") %>%
    dplyr::filter(group_new_name == x) %>%
    dplyr::mutate(Variant_Group = `Variant Group`) %>%
    #anova_test(Percentage_Species ~ group_new_name) %>%
    tukey_hsd(Percentage_Species ~ Variant_Group) %>%
    dplyr::mutate(Group_1_variable = case_when(group1 == "Uninfected" ~ 1 ,
                                               TRUE ~ 0)) %>%
    dplyr::mutate(Group_2_variable = case_when(group2 == "Uninfected" ~ 1 ,
                                               TRUE ~ 0)) %>%
    dplyr::mutate(Group_Variable = Group_1_variable+Group_2_variable) %>%
    dplyr::filter(Group_Variable == 1) %>%
    dplyr::mutate(group1 = case_when(group1 == "Uninfected" ~ "" , 
                                     TRUE ~ group1)) %>%
    dplyr::mutate(group2 = case_when(group2 == "Uninfected" ~ "" , 
                                     TRUE ~ group2)) %>%
    dplyr::mutate(`Variant Group` = paste0(group1, group2)) %>%
    dplyr::mutate(group_new_name = x) %>%
    dplyr::select(-c(term,group1,group2,Group_1_variable,Group_2_variable,Group_Variable)) 
}) %>% do.call('rbind' , .) %>%
  full_join(. , 
            comparison_infected_uninfected_species_dataframe %>% 
              dplyr::filter(`Sample Type` == "Fecal") %>%
              dplyr::mutate(group_new_name = paste0(Division ,";", Phylum ,";", Class ,";", Order ,";", Family ,";", Genus ,";", Species)))

filtered_group_names = anova_dataframe %>% dplyr::filter(Percentage_Compared_Uninfected >= 5 | Percentage_Compared_Uninfected <= 0.2) %>% .$group_new_name

t_test_dataframe = lapply(unique(group_phylum_species_dataframe$group_new_name), function(x){
  group_phylum_species_dataframe %>% 
    dplyr::filter(`Sample Type` == "Fecal") %>%
    dplyr::filter(group_new_name == x) %>%
    dplyr::mutate(Variant_Group = `Variant Group`) %>%
    #anova_test(Percentage_Species ~ group_new_name) %>%
    #tukey_hsd(Percentage_Species ~ Variant_Group) %>%
    rstatix::pairwise_t_test(
      Percentage_Species ~ Variant_Group,
      p.adjust.method = "BH"   # options: "bonferroni", "holm", "fdr", etc.
    ) %>%
    dplyr::mutate(Group_1_variable = case_when(group1 == "Uninfected" ~ 1 ,
                                               TRUE ~ 0)) %>%
    dplyr::mutate(Group_2_variable = case_when(group2 == "Uninfected" ~ 1 ,
                                               TRUE ~ 0)) %>%
    dplyr::mutate(Group_Variable = Group_1_variable+Group_2_variable) %>%
    dplyr::filter(Group_Variable == 1) %>%
    dplyr::mutate(group1 = case_when(group1 == "Uninfected" ~ "" , 
                                     TRUE ~ group1)) %>%
    dplyr::mutate(group2 = case_when(group2 == "Uninfected" ~ "" , 
                                     TRUE ~ group2)) %>%
    dplyr::mutate(`Variant Group` = paste0(group1, group2)) %>%
    dplyr::mutate(group_new_name = x) %>%
    dplyr::select(-c(".y.",group1,group2,Group_1_variable,Group_2_variable,Group_Variable)) 
}) %>% do.call('rbind' , .) %>%
  full_join(. , 
            comparison_infected_uninfected_species_dataframe %>% 
              dplyr::filter(`Sample Type` == "Fecal") %>%
              dplyr::mutate(group_new_name = paste0(Division ,";", Phylum ,";", Class ,";", Order ,";", Family ,";", Genus ,";", Species)))

# 1000:1000
# Filtration of groups which have the percentage values less than 5 and  
ggplot(t_test_dataframe %>% dplyr::filter(group_new_name %in% unique(filtered_group_names)) %>% dplyr::filter(`Sample Type` == "Fecal") %>% dplyr::mutate(`Variant Group` = factor(`Variant Group` , levels = c("Hong Kong" , "Delta","Omicron")))) +
  geom_tile(aes(x = `Variant Group`, y = Species , fill = Percentage_Compared_Uninfected) , colour = "black" , size = 0.1) +
  geom_text(aes(x = `Variant Group`, y = Species,
                label = ifelse(p.adj.signif != "ns", p.adj.signif, "")),
            size = 4, fontface = "bold", vjust = 0.7, hjust = 0.5) +
  theme_set(theme_bw()) +
  theme(axis.text=element_text(size=7,colour = "black"),
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.title=element_text(size=20,face="bold") ,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1, 'cm') ,
        legend.title = element_text(size=15 , face = "bold"),
        legend.text = element_text(size=15) , 
        strip.text.x = element_blank(), 
        strip.text.y = element_text(size = 12 , angle = 0 , face = "bold"),
        strip.background.y = element_blank(),
        legend.spacing.y = unit(0, 'cm')) +
  labs(y = "Smallest Taxonomic Group\n(Species/Genus/Family/Order/Class)") +
  #facet_wrap(~`Day of collection` , scales = "free_x") +
  #facet_wrap(~Phylum , scales = "free_y")
  facet_grid(Phylum ~ `Day of collection` , scales = "free_y",space = "free_y") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, 
                       trans = "log2", name = "Log2 Fold\nPercentage Change",
                       labels = scales::label_number(scilimits = c(0, 0)))
