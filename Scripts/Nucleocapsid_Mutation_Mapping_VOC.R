library(tidyverse)
library(RColorBrewer)
library(treeio)
library(Biostrings)
library(msa)
library(parallel)
library(readxl)
library(scales)
library(patchwork)

# Mutation Mapping
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# INDEXING
# Indexing would be done as a whole for the whole genome 
# if both sequence nucleotide as well as Wuhan_Hu_1 sequence nucleotide are "-", then those positions are deleted out before giving the final pre indexed numbering

sequences_file = treeio::read.fasta("Sequences_with_wuhan_aligned.fasta")

sequences_before_indexing_dataframe_list = lapply(c(1:length(sequences_file)), function(x) 
{ dummy = data.frame(dummy_column_name = as.character(sequences_file[[x]]) %>% toupper())
colnames(dummy) <- names(sequences_file[x])
dummy}
)  

sequences_before_indexing_dataframe = do.call('cbind',sequences_before_indexing_dataframe_list) %>%
  dplyr::mutate(Position_Before_Indexing = row_number() , .before = Wuhan_Hu_1)

indexing_of_dataframe = dplyr::select(sequences_before_indexing_dataframe , Position_Before_Indexing , Wuhan_Hu_1) %>%
  dplyr::filter(!Wuhan_Hu_1 %in% c("-")) %>%
  dplyr::mutate(Indexed_Position = row_number())  

sequences_after_indexing_dataframe = full_join(indexing_of_dataframe,sequences_before_indexing_dataframe) %>%
  dplyr::arrange(Position_Before_Indexing)  

# Make temporary Position_Before_Indexing and Indexed_Position columns for later purpose
sequences_after_indexing_dataframe = sequences_after_indexing_dataframe %>%
  dplyr::mutate(Temp_Position_Before_Indexing = Position_Before_Indexing , .after = Position_Before_Indexing) %>%
  dplyr::mutate(Temp_Indexed_Position = Indexed_Position , .after = Indexed_Position)

sequences_gathered_after_indexing_dataframe = sequences_after_indexing_dataframe %>% 
  tidyr::gather(key = "Sequence_Name",value = "Nucleotide",6:ncol(sequences_after_indexing_dataframe))

sequences_after_indexing_dataframe_list = sequences_gathered_after_indexing_dataframe %>%
  split(. , .$Sequence_Name)  

sequences_gathered_after_indexing_processed_dataframe_list = lapply(sequences_after_indexing_dataframe_list, function(x) {
  x %>%
    dplyr::filter(!(Wuhan_Hu_1 == "-" & Nucleotide == "-")) %>% # If both Wuhan-Hu-1 and nucleotide are "-" then those positions are avoided from that sequence
    dplyr::mutate(Position_After_Indexing = row_number()) # Position after indexing is generated
}) 

sequences_gathered_after_indexing_processed_dataframe = do.call('rbind',sequences_gathered_after_indexing_processed_dataframe_list)

nucleotide_variations_dataframe = sequences_gathered_after_indexing_processed_dataframe %>%
  dplyr::filter(Wuhan_Hu_1 != Nucleotide)  


# Single_Amino Acid_Polymorphism (Pre-Canonical Stop)

nucleotide_polymorphism_dataframe = nucleotide_variations_dataframe %>% 
  dplyr::filter(!Wuhan_Hu_1 %in% c("-")) %>%
  dplyr::filter(!Nucleotide %in% c("-","N")) %>%
  dplyr::mutate(Mutation = paste0(Wuhan_Hu_1,Indexed_Position,Nucleotide))

if (!is.null(nucleotide_polymorphism_dataframe)) {
  nucleotide_polymorphism_frequency_dataframe = nucleotide_polymorphism_dataframe %>% 
    dplyr::group_by(Mutation) %>%
    dplyr::summarise(Frequency_of_mutation = length(Sequence_Name),
                     Position = dplyr::first(Indexed_Position),
                     Mutation_Type = "Nucleotide_Polymorphism")
  
  
} else if (is.null(nucleotide_polymorphism_dataframe)) {
  nucleotide_polymorphism_frequency_dataframe = data.frame(Mutation = character(),Frequency_of_mutation = numeric(),Position = numeric(), Mutation_Type = character()) 
}


# Deletion_Till_Wuhan_End

deletion_dataframe_list = nucleotide_variations_dataframe %>% 
  dplyr::filter(Nucleotide == "-" | Wuhan_Hu_1 == "-") %>%
  drop_na(Indexed_Position) %>%
  split(. , .$Sequence_Name)

deletion_mutation_dataframe_list = lapply(deletion_dataframe_list , function(x) {
  x %>%
    dplyr::arrange(Position_After_Indexing) %>%
    dplyr::mutate(group = cumsum(c(1, diff(Position_After_Indexing)) != 1)) %>%
    dplyr::group_by(group) %>% 
    dplyr::summarise( Sequence_Name = dplyr::first(Sequence_Name),
                      Mutation = paste0("Del", dplyr::first(Indexed_Position), "-", dplyr::last(Indexed_Position)),
                      Deleted_Index_Positions = paste0(Indexed_Position,collapse = ","),
                      Deleted_unindex_Positions = paste0(Position_After_Indexing,collapse = ","),
                      Indexed_Position = dplyr::first(Indexed_Position),
                      Position_Before_Indexing = dplyr::first(Position_Before_Indexing),
                      Position_After_Indexing = dplyr::first(Position_After_Indexing),
                      Deleted_Wuhan_Nucleotide = paste0(Wuhan_Hu_1,collapse = ","),
                      Deleted_evidence = paste0(Nucleotide,collapse = ","),
                      Number_of_Nucleotide_deletion = dplyr::n()) %>%
    dplyr::mutate(Mutation = case_when(Number_of_Nucleotide_deletion == 1 ~ paste0("Del",Indexed_Position),
                                       Number_of_Nucleotide_deletion != 1 ~ Mutation))
})

deletion_mutation_dataframe = do.call('rbind',deletion_mutation_dataframe_list)


if (!is.null(deletion_mutation_dataframe)) {
  deletion_frequency_dataframe = deletion_mutation_dataframe %>%
    dplyr::group_by(Mutation) %>%
    dplyr::summarise(Frequency_of_mutation = length(Sequence_Name),
                     Position = dplyr::first(Indexed_Position),
                     Mutation_Type = "Deletion")
} else if (is.null(deletion_mutation_dataframe)) {
  deletion_frequency_dataframe = data.frame(Mutation = character(),Frequency_of_mutation = numeric(),Position = numeric(), Mutation_Type = character()) 
}


# Insertion_Till_Wuhan_End

insertion_dataframe_list = sequences_gathered_after_indexing_processed_dataframe %>% 
  dplyr::filter(Wuhan_Hu_1 == "-") %>% 
  split(. , .$Sequence_Name)

insertion_mutation_dataframe_list = lapply(insertion_dataframe_list , function(x) {
  grouped_internal_list = x %>%
    dplyr::arrange(Position_After_Indexing) %>%
    dplyr::mutate(group = cumsum(c(1, diff(Position_After_Indexing)) != 1)) %>%
    split(.,.$group)
  temp_sequence_name = do.call('rbind',grouped_internal_list)$Sequence_Name %>%
    unique()
  required_processed_sequence_dataframe = sequences_gathered_after_indexing_processed_dataframe %>%
    dplyr::filter(Sequence_Name == temp_sequence_name)
  mclapply(grouped_internal_list, function(y) {
    y %>%
      dplyr::arrange(Position_After_Indexing) %>%
      dplyr::mutate(Insertion_Before_Index = required_processed_sequence_dataframe[required_processed_sequence_dataframe$Position_After_Indexing == (as.numeric(dplyr::first(Position_After_Indexing)) - 1) , ]$Indexed_Position) %>%
      dplyr::mutate(Insertion_After_Index =  required_processed_sequence_dataframe[required_processed_sequence_dataframe$Position_After_Indexing == (as.numeric(dplyr::last(Position_After_Indexing)) + 1) , ]$Indexed_Position) %>%
      dplyr::mutate(Inserted_Nucleotide = paste0(Nucleotide,collapse = ",")) %>%
      dplyr::mutate(Insertion = paste0(Nucleotide,collapse = "")) %>%
      dplyr::mutate(Number_of_Nucleotide_inserted = n()) %>%
      dplyr::mutate( Mutation = paste0(Insertion_Before_Index, "_", Insertion_After_Index ,"ins",Insertion)) %>% 
      dplyr::mutate(N_Status = if_else(all(Nucleotide == "N"), "Only_N",
                                       if_else(any(Nucleotide == "N"), "Some_N",
                                               "No_N"))) %>%
      dplyr::select(Sequence_Name , Insertion_Before_Index , Insertion_After_Index , Inserted_Nucleotide , Number_of_Nucleotide_inserted , Mutation , N_Status) %>%
      head(.,1)
  }, mc.cores = 8) %>% do.call('rbind',.) %>% dplyr::filter(!N_Status %in% "Only_N")
})
insertion_mutation_dataframe = do.call('rbind',insertion_mutation_dataframe_list) 


if (!is.null(insertion_mutation_dataframe)) {
  insertion_frequency_dataframe = insertion_mutation_dataframe %>%
    dplyr::filter(!N_Status %in% "Only_N") %>% # We got rid of Only X configurations here before getting frequency of the mutation
    dplyr::group_by(Mutation) %>%
    dplyr::summarise(Frequency_of_mutation = length(Sequence_Name),
                     Position = dplyr::first(Insertion_Before_Index),
                     Mutation_Type = "Insertion") 
} else if (is.null(insertion_mutation_dataframe)) {
  insertion_frequency_dataframe = data.frame(Mutation = character(),Frequency_of_mutation = numeric(),Position = numeric(), Mutation_Type = character()) 
}



if (!is.null(nucleotide_variations_dataframe)) {
  nucleotide_variations_dataframe= nucleotide_variations_dataframe
}

if (!is.null(deletion_mutation_dataframe)) {
  deletion_mutation_dataframe = deletion_mutation_dataframe %>%
    dplyr::mutate(Mutation_Type = "Deletion")
}

if (!is.null(nucleotide_polymorphism_dataframe)) {
  nucleotide_polymorphism_dataframe = nucleotide_polymorphism_dataframe %>%
    dplyr::mutate(Mutation_Type = "Nucleotide Polymorphism")
}
# Only N mutations would be there in this below dataframe list
if (!is.null(insertion_mutation_dataframe)) {
  insertion_mutation_dataframe = insertion_mutation_dataframe %>%
    dplyr::mutate(Mutation_Type = "Insertion")
}

Combined_Mutation_Frequency_Dataframe = rbind(nucleotide_polymorphism_frequency_dataframe,
                                              deletion_frequency_dataframe, 
                                              insertion_frequency_dataframe)


mutations_for_mapping_list = list()



if (!is.null(deletion_mutation_dataframe)) {
  mutations_for_mapping_list[["deletion_mutation"]] = deletion_mutation_dataframe %>%
    tidyr::separate_rows(. , Deleted_Index_Positions , sep = ",") %>%
    dplyr::select(Deleted_Index_Positions , Mutation,Sequence_Name , Mutation_Type) %>%
    dplyr::rename(Position = Deleted_Index_Positions) %>%
    dplyr::mutate(Position = as.numeric(Position)) %>%
    dplyr::mutate(Nucleotide = "-")
}

if (!is.null(nucleotide_polymorphism_dataframe)) {
  mutations_for_mapping_list[["Nucleotide_Polymorphism"]] = nucleotide_polymorphism_dataframe %>%
    dplyr::select(Indexed_Position , Mutation,Sequence_Name , Mutation_Type , Nucleotide) %>%
    dplyr::rename(Position = Indexed_Position) %>%
    dplyr::mutate(Position = as.numeric(Position)) 
  
}

if (!is.null(insertion_mutation_dataframe)) {
  mutations_for_mapping_list[["Insertion_mutation"]] = insertion_mutation_dataframe %>%
    dplyr::filter(!N_Status %in% "Only_N") %>%
    dplyr::select(Insertion_Before_Index ,Mutation,Sequence_Name,Mutation_Type) %>%
    dplyr::rename(Position = Insertion_Before_Index) %>%
    dplyr::mutate(Position = as.numeric(Position)) 
}

mutations_for_mapping_list[["Nucleotide_Polymorphism"]]$Sequence_Name = factor(mutations_for_mapping_list[["Nucleotide_Polymorphism"]]$Sequence_Name, levels = c("HongKong" , "Delta" , "Omicron"))

colors_used = c(rainbow(4))

# 1400 : 2000
NUCLEOTIDE_MUTATION_MAP_PLOT = ggplot(mutations_for_mapping_list[["Nucleotide_Polymorphism"]] %>% dplyr::filter(Nucleotide != "Y") %>% dplyr::mutate(Nucleotide = case_when(Nucleotide == "T" ~ "U" ,
                                                                                                                                                                            TRUE ~ Nucleotide)) , y = Sequence_Name) +
  geom_point(aes(x = Position, y = Sequence_Name , colour = Nucleotide), shape = "|", size = 10) +
  scale_color_manual(values = colors_used) +
  #scale_color_hue(direction = 1) +
  labs(color='Nucleoside',x = "Position") +
  ylab("Sequences") +
  theme_bw() +
  theme(legend.key.size = unit(0.75, 'cm') ,
        axis.text =element_text(size=20),
        legend.text = element_text(size=20),
        axis.title=element_text(size=30,face="bold"),
        legend.title = element_text(size=25,face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  geom_point(data = mutations_for_mapping_list[["deletion_mutation"]] , aes(x = Position, y = Sequence_Name) , shape = "|",size = 10,color = "black") +
  geom_point(data = mutations_for_mapping_list[["Insertion_mutation"]] , aes(x = Position, y = Sequence_Name) , shape = "v",size = 10,color = "dimgray") +
  scale_shape_manual(name = "Shapes",
                     values = c("|" = 6)) +
  guides(shape = guide_legend(title = "Shapes",override.aes = list(size = 5)))

NUCLEOTIDE_MUTATION_MAP_PLOT

################################################### MAKING THE REPRESENTATIVE WITH DIVISIONS IN SPIKE ###############################################
colours_here = c("darkgrey","darkgrey",hue_pal()(12))


positional_dataframe = read.csv("ORF_locations.tsv" , sep = "\t" , header = F) %>%
  dplyr::rename(ORF = V1) %>%
  dplyr::rename(Start = V2) %>%
  dplyr::rename(Stop = V3) %>%
  dplyr::arrange(Start) %>%
  dplyr::mutate(ORF = gsub("Nucleocapsid","N",ORF)) %>% 
  dplyr::mutate(ORF = gsub("Envelope","E",ORF)) %>% 
  dplyr::mutate(ORF = gsub("Membrane","M",ORF)) %>%
  dplyr::mutate(ORF = gsub("_UTR","'UTR",ORF))


ORF_locations = positional_dataframe %>% 
  dplyr::mutate(median_x = Start + floor((Stop-Start)/2)) %>% 
  dplyr::mutate(text_height = c(7 , 7, 2 , 7,7,19,7,19,0,27,7,19,0 , 7)) %>% 
  dplyr::mutate(y_min_value = text_height - 7) %>% 
  dplyr::mutate(y_max_value = text_height + 7)

# 1400 : 200

p1 <-ggplot(ORF_locations) + 
  geom_rect(aes(xmin=Start,xmax=Stop,fill=ORF,ymin=y_min_value,ymax=y_max_value),
            colour="white", 
            size=0.5, 
            alpha=5) + 
  
  ylim(-10,35)+
  scale_fill_manual(values = colours_here) +
  geom_text(data=ORF_locations,aes(x=median_x,y=text_height,label=ORF), size=4.2,color = "black",fontface = "bold") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank())

(p1/NUCLEOTIDE_MUTATION_MAP_PLOT) + plot_layout(heights = c(0.5, 0.5))

