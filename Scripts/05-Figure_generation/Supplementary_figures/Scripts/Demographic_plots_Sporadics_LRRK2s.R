
# Setup
# Load libraries
suppressPackageStartupMessages({
library(magrittr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggdist)
})

conflicted::conflict_prefer_all(winner = "dplyr", quiet = T)

# 0. Setup ----
figure_out_path<-"/home/jbrenton/TDP43_PD_Paper_Figures/Drafting_figures_code/Supplementary_figures/Plots"

main_path<-"/home/jbrenton/LRRK2_PM_HS/Results/Outlier_case_removed_17.10.24"

analysis_colours<-c("#0c5394", "#e06666")

colours<-c("#C3B3A1", "#91D1C2", "#4B8D86", "#125C5D")

lrrk2_group_order<-c("Control", "LRRK2\nG2019S")

# mid_brain_area_order<-c("Substantia nigra", "Caudate", "Putamen", "Ant.\nCingulate",
#                         "Parahippocampal", 'Temporal',
#                         "Frontal", "Parietal", "Combined")
# 
# late_brain_area_order<-c("Ant.\nCingulate", 'Temporal ',
#                          "Frontal", "Parietal", "Combined")

extrafont::font_import(paths = "/home/jbrenton/Wood_ASAP_post_align_files/Final_dataset/Wood_limma_analysis/data/Files_for_Rome_meeting", 
                       pattern = 'Roboto', prompt = F)

# New theme based off Regina theme:

theme_jb <- function(text_size=12) {
  theme_bw(base_family = "Roboto", base_size = text_size) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_text(vjust = 0.6),
      panel.spacing = unit(0.1, "lines")
    )
}


# set jb theme for all plots
text_size2<-12
get.theme.jb <- theme_set(
  theme_bw(base_family = "Roboto",
           base_size = text_size2)+
    theme( #panel.grid.major.x = element_blank(),
      legend.position = "top",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_text(vjust = 0.6),
      panel.spacing = unit(0.1, "lines"))
)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load Braak stage 3-4 metadata
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# outliers removed already
Braak_3_4_path<-"/home/jbrenton/TDP43_PD_Paper_Figures/Drafting_figures_code/Supplementary_figures/Data/metadata_mid_stage_sporadic"

Braak_3_4_RIN <- read_csv(file.path(Braak_3_4_path, "ASSAY_RNAseq.csv")) %>% select(sample_id,RIN)

Braak_3_4_CLINPATH <- read_csv(file.path(Braak_3_4_path, "CLINPATH.csv")) %>% select(subject_id, duration_pmi, age_at_death, path_braak_nft,  path_braak_asyn, path_thal)

Braak_3_4_SAMPLE <- read_csv(file.path(Braak_3_4_path, "SAMPLE.csv")) %>% select(sample_id, subject_id)
  
Braak_3_4_SUBJECT <- read_csv(file.path(Braak_3_4_path, "SUBJECT.csv")) %>% select(subject_id, sex, age_at_collection, primary_diagnosis) %>% 
  mutate(Group=case_when(primary_diagnosis=="No PD nor other neurological disorder" ~ "Control",
                         primary_diagnosis=="Idiopathic PD" ~ "PD",
                         primary_diagnosis=="Hemiparkinson/hemiatrophy syndrome" ~ "PD",
                           )) %>% select(-primary_diagnosis)


Braak_3_4_metadata<-Braak_3_4_SAMPLE %>% left_join(., Braak_3_4_RIN, by="sample_id") %>%
  left_join(., Braak_3_4_CLINPATH, by="subject_id") %>%
  left_join(., Braak_3_4_SUBJECT, by="subject_id") %>% distinct(sample_id,.keep_all = T)



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Braak stage 6
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load Braak stage 3-4 metadata
# 
Braak_6_path<-"/home/jbrenton/TDP43_PD_Paper_Figures/Drafting_figures_code/Supplementary_figures/Data/metadata_late_stage_sporadic"

Braak_6_RIN <- read_csv(file.path(Braak_6_path, "ASSAY_bulkRNAseq.csv")) %>% select(sample_id,RIN)

Braak_6_CLINPATH <- read_csv(file.path(Braak_6_path, "CLINPATH.csv")) %>% select(subject_id, duration_pmi, age_at_death, path_braak_nft,  path_braak_asyn, path_thal)

Braak_6_SAMPLE <- read_csv(file.path(Braak_6_path, "SAMPLE.csv")) %>% select(sample_id, subject_id, alternate_id)
# find outliers
sample_outliers16 <- readRDS("/home/jbrenton/ASAP_Hardy_post_align_files/Final_Dataset/Hardy_Limma_Analysis_Final_version_only/results/sample_outliers16.rds") %>%
  filter(outlier==T) %>% pull(bxp_id_full)
# remove outliers
Braak_6_SAMPLE<-Braak_6_SAMPLE %>% filter(!alternate_id %in% sample_outliers16) %>% select(sample_id, subject_id)

Braak_6_SUBJECT <- read_csv(file.path(Braak_6_path, "SUBJECT.csv")) %>% select(subject_id, sex, age_at_collection, primary_diagnosis) %>% 
  mutate(Group=case_when(primary_diagnosis=="No PD nor other neurological disorder" ~ "Control",
                         primary_diagnosis=="Idiopathic PD" ~ "PD",
                         primary_diagnosis=="Hemiparkinson/hemiatrophy syndrome" ~ "PD",
                         primary_diagnosis=="Other neurological disorder" ~ "Control"
  )) %>% select(-primary_diagnosis)


Braak_6_metadata<-Braak_6_SAMPLE %>% left_join(., Braak_6_RIN, by="sample_id") %>%
  left_join(., Braak_6_CLINPATH, by=c("sample_id"="subject_id")) %>%
  left_join(., Braak_6_SUBJECT, by="subject_id") %>% distinct(sample_id,.keep_all = T)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LRRK2 Samples
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# RIN, Brain areas, Sex and age
LRRK2_metadata <- read_csv("/home/jbrenton/LRRK2_PM_HS/Results/Outlier_case_removed_17.10.24/metadata_without_removal_of_correlated_terms.csv")

LRRK2_metadata %<>% select(Case, Group, `Brain Area`, Sex, Age, Group, RIN)






# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Graphing ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Numeric Variables: ----


# RIN ----
p_RIN_data <- rbind(
  Braak_3_4_metadata %>% select(RIN, Group) %>% mutate(Dataset = "Mid-stage"),
  Braak_6_metadata %>% select(RIN, Group) %>% mutate(Dataset = "Late-stage"),
  LRRK2_metadata %>% select(RIN, Group) %>% mutate(Dataset = "LRRK2 G2019S")
) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("Mid-stage", "Late-stage", "LRRK2 G2019S")),
    
    # Create new variable to differentiate PD across datasets
    Group_Colour = case_when(
      Group == "Control" ~ "Control",
      Group == "PD" & Dataset == "Mid-stage" ~ "PD_Mid",
      Group == "PD" & Dataset == "Late-stage" ~ "PD_Late",
      Group == "LRRK2-G2019S" ~ "LRRK2 G2019S"
    ),
    Group_Colour = factor(Group_Colour, levels = c("Control", "PD_Mid", "PD_Late", "LRRK2 G2019S"))
  )

 p_RIN <- ggplot(p_RIN_data, aes(Group, RIN, colour=Group_Colour))+
  facet_wrap(~Dataset, scales = "free_x")+
  stat_halfeye(
  aes(fill = Group_Colour),
  alpha = 0.5,
  # adjust = 0.5,
  height = 0.8,
  justification = -0.2,
  scale = 0.4,
  point_interval = NULL
)+
   geom_boxplot(
     position = position_nudge(x = -0.2),
     width = 0.2,
     outlier.color = NA,
     alpha = 0.5, 
     outlier.shape = NA,
   )+geom_point(size=1, alpha=0.6)+
   scale_color_manual(values = colours) +
   scale_fill_manual(values = colours)+
   guides(fill='none', colour='none')+
   labs(y="RIN", x="")+
 theme_jb()
 # +
   # coord_flip()

p_RIN


comparisons <- list(c("Control", "PD"), c("Control", "PD"), c("Control", "LRRK2-G2019S")) 

stat.test <-  compare_means(RIN ~ Group, p_RIN_data, group.by = "Dataset")

stat.test 

ggboxplot(data = p_RIN_data,
          x = "Group",y = "RIN",
          color = "Group",
          palette = "RIN",
          add = "jitter",
          facet.by = "Dataset", outlier.shape = NA, scales="free_x")+ 
  stat_pvalue_manual(
  stat.test, x = "Dataset", y.position = 9,
  label = "p.signif",
  position = position_dodge(0.8)
)


# Age ----

p_age_data <- rbind(
  Braak_3_4_metadata %>% distinct(subject_id,.keep_all = T) %>% select(age_at_death, Group) %>% mutate(Dataset = "Mid-stage"),
  Braak_6_metadata %>% distinct(subject_id,.keep_all = T) %>% select(age_at_death, Group) %>% mutate(Dataset = "Late-stage"),
  LRRK2_metadata %>% distinct(Case,.keep_all = T) %>% select(Age, Group) %>% mutate(Dataset = "LRRK2 G2019S") %>% dplyr::rename(age_at_death=Age)
) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("Mid-stage", "Late-stage", "LRRK2 G2019S")),
    
    # Create new variable to differentiate PD across datasets
    Group_Colour = case_when(
      Group == "Control" ~ "Control",
      Group == "PD" & Dataset == "Mid-stage" ~ "PD_Mid",
      Group == "PD" & Dataset == "Late-stage" ~ "PD_Late",
      Group == "LRRK2-G2019S" ~ "LRRK2 G2019S"
    ),
    Group_Colour = factor(Group_Colour, levels = c("Control", "PD_Mid", "PD_Late", "LRRK2 G2019S"))
  )

p_age <- ggplot(p_age_data, aes(Group, age_at_death, colour=Group_Colour))+
  facet_wrap(~Dataset, scales = "free_x")+
  stat_halfeye(
    aes(fill = Group_Colour),
    alpha = 0.5,
    # adjust = 0.5,
    height = 0.8,
    justification = -0.2,
    scale = 0.4,
    point_interval = NULL
  )+
  geom_boxplot(
    position = position_nudge(x = -0.2),
    width = 0.2,
    outlier.color = NA,
    alpha = 0.5, 
    outlier.shape = NA,
  )+geom_point(size=1, alpha=0.6)+
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours)+
  labs(y="Age(Years)", x="")+
guides(fill='none', colour='none')+
theme_jb()+  
  scale_x_discrete(expand = c(1,0))
# +
  # coord_flip()

p_age

stat.test <-  compare_means(age_at_death ~ Group, p_age_data, group.by = "Dataset")

stat.test

# PMI -----

p_PMI_data <- rbind(
  Braak_3_4_metadata %>% distinct(subject_id,.keep_all = T) %>% select(duration_pmi, Group) %>% mutate(Dataset = "Mid-stage"),
  Braak_6_metadata %>% distinct(subject_id,.keep_all = T) %>% select(duration_pmi, Group) %>% mutate(Dataset = "Late-stage")
  # No current LRRK2 PMI data
  # LRRK2_metadata %>% distinct(Case,.keep_all = T) %>% select(Age, Group) %>% mutate(Dataset = "LRRK2 G2019S") %>% rename(age_at_death=Age)
) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("Mid-stage", "Late-stage")),
    
    # Create new variable to differentiate PD across datasets
    Group_Colour = case_when(
      Group == "Control" ~ "Control",
      Group == "PD" & Dataset == "Mid-stage" ~ "PD_Mid",
      Group == "PD" & Dataset == "Late-stage" ~ "PD_Late"
      # ,
      # Group == "LRRK2-G2019S" ~ "LRRK2 G2019S"
    ),
    Group_Colour = factor(Group_Colour, levels = c("Control", "PD_Mid", "PD_Late"
                                                   # , "LRRK2 G2019S"
                                                   ))
  ) %>%  mutate(duration_pmi=as.numeric(duration_pmi)) %>% drop_na(duration_pmi)

p_PMI <- ggplot(p_PMI_data, aes(Group, duration_pmi, colour=Group_Colour))+
  facet_wrap(~Dataset, scales = "free_x")+
  stat_halfeye(
    aes(fill = Group_Colour),
    alpha = 0.5,
    # adjust = 0.5,
    height = 0.8,
    justification = -0.2,
    scale = 0.4,
    point_interval = NULL
  )+
  geom_boxplot(
    position = position_nudge(x = -0.2),
    width = 0.2,
    outlier.color = NA,
    alpha = 0.5, 
    outlier.shape = NA,
  )+geom_point(size=1, alpha=0.6)+
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours)+
  labs(y="PMI(Hours)", x="")+
  guides(fill='none', colour='none')+
  theme_jb()
# +
  # coord_flip()

p_PMI

stat.test <-  compare_means(duration_pmi ~ Group, p_PMI_data, 
                            group.by = "Dataset")

stat.test

## Categorical Variables: ----

# Sex ----

p_sex_data <- rbind(
  Braak_3_4_metadata %>% distinct(subject_id,.keep_all = T) %>% select(sex, Group) %>% mutate(Dataset = "Mid-stage"),
  Braak_6_metadata %>% distinct(subject_id,.keep_all = T) %>% select(sex, Group) %>% mutate(Dataset = "Late-stage"),
  LRRK2_metadata %>% distinct(Case,.keep_all = T) %>% select(Sex, Group) %>% mutate(Dataset = "LRRK2 G2019S") %>% dplyr::rename(sex=Sex) %>%
    mutate(sex=case_when(sex=="F" ~ "Female",
              sex=="M" ~ "Male"))
) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("Mid-stage", "Late-stage", "LRRK2 G2019S")),
    
    # Create new variable to differentiate PD across datasets
    Group_Colour = case_when(
      Group == "Control" ~ "Control",
      Group == "PD" & Dataset == "Mid-stage" ~ "PD_Mid",
      Group == "PD" & Dataset == "Late-stage" ~ "PD_Late",
      Group == "LRRK2-G2019S" ~ "LRRK2 G2019S"
    ),
    Group_Colour = factor(Group_Colour, levels = c("Control", "PD_Mid", "PD_Late", "LRRK2 G2019S"
    ))
  ) 



p_sex <-ggplot(p_sex_data, aes(x = Group, fill = sex)) +
  geom_bar(position = "fill", colour = "black") +  # Use 'fill' to normalize per bar
  facet_wrap(~Dataset, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show % instead of 0–1
  labs(x = "", y = "", fill = "Sex") +
  scale_fill_grey(na.value = "white")+
  theme(legend.position = "bottom",
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, 'cm'))
p_sex

stat.test <-  compare_means(sex ~ Group, p_sex_data, group.by = "Dataset", method = "kruskal.test")

stat.test

# Braak NFT -----

p_nft_data <- rbind(
  Braak_3_4_metadata %>% distinct(subject_id,.keep_all = T) %>% select(path_braak_nft, Group) %>% mutate(Dataset = "Mid-stage"),
  Braak_6_metadata %>% distinct(subject_id,.keep_all = T) %>% select(path_braak_nft, Group) %>% mutate(Dataset = "Late-stage")
  # ,
  
  # No current LRRK2 data for NFTs
  
  # LRRK2_metadata %>% distinct(Case,.keep_all = T) %>% select(Sex, Group) %>% mutate(Dataset = "LRRK2 G2019S") %>% rename(sex=Sex) %>%
  #   mutate(sex=case_when(sex=="F" ~ "Female",
  #                        sex=="M" ~ "Male"))
) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("Mid-stage", "Late-stage"
                                         # , "LRRK2 G2019S"
                                         )),
    
    # Create new variable to differentiate PD across datasets
    Group_Colour = case_when(
      Group == "Control" ~ "Control",
      Group == "PD" & Dataset == "Mid-stage" ~ "PD_Mid",
      Group == "PD" & Dataset == "Late-stage" ~ "PD_Late"
      # ,
      # Group == "LRRK2-G2019S" ~ "LRRK2 G2019S"
    ),
    Group_Colour = factor(Group_Colour, levels = c("Control", "PD_Mid", "PD_Late"
                                                   # , "LRRK2 G2019S"
    ))
  ) %>% drop_na(path_braak_nft) %>%
  mutate(path_braak_nft=factor(path_braak_nft)) %>%
  mutate(path_braak_nft=factor(path_braak_nft, levels = rev(levels(path_braak_nft))))


p_nft <-ggplot(p_nft_data, aes(x = Group, fill = path_braak_nft)
                               ) +
  geom_bar(position = "fill", colour = "black") +  # Use 'fill' to normalize per bar
  facet_wrap(~Dataset, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show % instead of 0–1
  labs(y = "Percentage")+
  labs(x = "", y = "", legend = "")+
  scale_fill_manual(values = MetBrewer::met.brewer(name = "Egypt", n = 7))+
  labs(x = "", y = "", fill = "Braak Stage (Tau)")+
  guides(fill = guide_legend(reverse = TRUE))+
  theme(legend.position = "bottom",
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.3, 'cm'))
# +
  # scale_fill_grey(na.value = "white") 

p_nft

stat.test <-  compare_means(path_braak_nft ~ Group, p_nft_data, group.by = "Dataset", method = "kruskal.test")

stat.test

# Amyloid - thal ----

p_thal_data <- rbind(
  Braak_3_4_metadata %>% distinct(subject_id,.keep_all = T) %>% select(path_thal, Group) %>% mutate(Dataset = "Mid-stage"),
  Braak_6_metadata %>% distinct(subject_id,.keep_all = T) %>% select(path_thal, Group) %>% mutate(Dataset = "Late-stage")
  # ,
  
  # No current LRRK2 data for NFTs
  
  # LRRK2_metadata %>% distinct(Case,.keep_all = T) %>% select(Sex, Group) %>% mutate(Dataset = "LRRK2 G2019S") %>% rename(sex=Sex) %>%
  #   mutate(sex=case_when(sex=="F" ~ "Female",
  #                        sex=="M" ~ "Male"))
) %>%
  mutate(
    Dataset = factor(Dataset, levels = c("Mid-stage", "Late-stage"
                                         # , "LRRK2 G2019S"
    )),
    
    # Create new variable to differentiate PD across datasets
    Group_Colour = case_when(
      Group == "Control" ~ "Control",
      Group == "PD" & Dataset == "Mid-stage" ~ "PD_Mid",
      Group == "PD" & Dataset == "Late-stage" ~ "PD_Late"
      # ,
      # Group == "LRRK2-G2019S" ~ "LRRK2 G2019S"
    ),
    Group_Colour = factor(Group_Colour, levels = c("Control", "PD_Mid", "PD_Late"
                                                   # , "LRRK2 G2019S"
    ))
  ) %>% drop_na(path_thal) %>%
  mutate(path_thal=factor(path_thal)) %>%
  mutate(path_thal=factor(path_thal, levels = rev(levels(path_thal))))


p_thal <-ggplot(p_thal_data, aes(x = Group, fill = path_thal)
) +
  geom_bar(position = "fill", colour = "black") +  # Use 'fill' to normalize per bar
  facet_wrap(~Dataset, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show % instead of 0–1
  labs(y = "Percentage")+
  labs(x = "", y = "", legend = "")+
  scale_fill_manual(values = MetBrewer::met.brewer(name = "Lakota", n = 7))+
  labs(x = "", y = "", fill = "Thal Stage (Amyloid)")+
  guides(fill = guide_legend(reverse = TRUE))+
  theme(legend.position = "bottom")
# +
# scale_fill_grey(na.value = "white") 

p_thal

stat.test <-  compare_means(path_thal ~ Group, p_thal_data, group.by = "Dataset", method = "kruskal.test")

stat.test



# Combining plots:

(p_RIN|p_age|p_PMI)/(p_sex|p_nft)


ggsave(
  plot = (p_RIN | p_age | p_PMI) / (p_sex | p_nft) +
    patchwork::plot_layout(heights = c(3, 2)) +
    patchwork::plot_annotation(tag_levels = "a", tag_prefix = "", tag_suffix = "", 
                               theme = theme(plot.tag = element_text(face = "bold", size = 14))),
  filename = file.path(figure_out_path, "Supplementary_figure_4_demographics.png"),
  width = 8000,
  height = 5250,
  dpi = 600,
  units = 'px'
)


