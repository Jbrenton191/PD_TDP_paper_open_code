# Figure 4 - Cell culture and Western blot result graphs

library(magrittr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(rstatix)

# 0. Setup ----
main_path<-here::here() # main project level i.e TDP43_paper_open_code level

data_path<-file.path(main_path, "Data")

figure_out_path<-file.path(main_path, "Plots/Main_Figures/Figure_4")

# create figure output directory if not present
if (!dir.exists(figure_out_path)) {
  dir.create(figure_out_path, recursive = T)
}

extrafont::font_import(paths = data_path, 
                       pattern = 'Roboto', prompt = F)

theme_jb <- function(text_size=12) {
  theme_bw(base_family = "Roboto", base_size = text_size) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.y = element_text(vjust = 0.6),
      panel.spacing = unit(0.1, "lines")
    )
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1. Cell Culture Averages Graph ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# load FUS data
FUS_ICC_culture_average <- read_csv(file.path(data_path, "LRRK2_MM/FUS ICC culture average.csv"))

FUS_ICC_culture_average$Protein<-"FUS"
# load tdp data
TDP_ICC_Culture_average <- read_csv(file.path(data_path, "LRRK2_MM/TDP ICC Culture average.csv"))

TDP_ICC_Culture_average$Protein<-"TDP-43"

# Merge both into a tibble
Merged_Proteins_ICC_Culture_df<-rbind(FUS_ICC_culture_average,TDP_ICC_Culture_average) %>% pivot_longer(cols = -Protein, names_to = "Group", values_to = "Intensity") %>%
  mutate(Group=case_when(Group=="WT" ~ "WT",
                           Group=="GKI" ~ "LRRK2\np.G2019S KI")) %>% 
  mutate(Group=factor(Group, levels=c("WT", "LRRK2\np.G2019S KI")))

# add comparisons for ggpubr t test function in graph
my_comparisons <- list( c("WT", "LRRK2\np.G20129S KI"))
# perform stat test
stat.test <- Merged_Proteins_ICC_Culture_df %>% 
  group_by(Protein) %>%                             # Group by Protein
  rstatix::t_test(Intensity ~ Group, ref.group = "WT") %>% add_significance()


p1<-ggbarplot(Merged_Proteins_ICC_Culture_df, x = "Group", y = "Intensity", 
                fill="Group", facet.by = "Protein", add = "mean") +
  stat_summary(
    fun.data = mean_se,  # Calculate mean and standard error
    geom = "errorbar",   # Add error bars
    width = 0.25,       # Error bar width
    color = "black",
    alpha=0.7
  )+
  geom_jitter(
    width = 0.075, height=0.5,
    size = 1.8, alpha = 0.5, color="black") +
  facet_wrap(~Protein, ncol = 1)+
  stat_pvalue_manual( # add significance test to plot
    stat.test, 
    label = "p.signif", 
    y.position = c(180,165), size = 3.5
  )+ 
  scale_fill_manual(values = c("#C3B3A1", "#125C5D"))+
  scale_colour_manual(values = c("#C3B3A1", "#125C5D"))+
  labs(x = "", 
       y = "Fluorescence Intensity (a.u)") +
  # theme_jb()+
  theme_bw(base_family = "Roboto")+
  theme(
    text=element_text(size=12),
    strip.text = element_text(size=13),
    strip.background = element_rect(fill = scales::alpha("grey", 0.2)),
    legend.position = "none",
    
    panel.spacing = unit(2.5, "lines"),
    legend.box.margin = margin(-50, 0, 0, 0) 
  )

# save graph
ggsave(
  plot=p1,
  filename = file.path(figure_out_path, "Figure_4d_quantification.png"),
  width = 800,
  height = 2500,
  dpi = 350,
  units = 'px'
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 2. Western Blot Averages Graph ----

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load FUS WB data
Fus_Histone_WB_protein <- read_csv(file.path(data_path, "LRRK2_MM/Fus_Histone WB protein.csv"))
# Label Protein so can merge
Fus_Histone_WB_protein$Protein<-"FUS"
# Load TDP-43 WB data
TDP_Histone_WB_protein <- read_csv(file.path(data_path, "LRRK2_MM/TDP_Histone WB protein.csv"))
# Label Protein so can merge
TDP_Histone_WB_protein$Protein<-"TDP-43"
# rename groups
Merged_Proteins_WB_df<-rbind(Fus_Histone_WB_protein,TDP_Histone_WB_protein) %>% pivot_longer(cols = -Protein, names_to = "Group", values_to = "Intensity") %>%
  mutate(Group=case_when(Group=="WT" ~ "WT",
                         Group=="GKI" ~ "LRRK2\np.G20129S KI")) %>% 
  mutate(Group=factor(Group, levels=c("WT", "LRRK2\np.G20129S KI")))

# perform stat test
stat.test <- Merged_Proteins_WB_df %>% 
  group_by(Protein) %>%                             # Group by Protein
  rstatix::t_test(Intensity ~ Group, ref.group = "WT") %>% add_significance()

# graph
p2<-ggbarplot(Merged_Proteins_WB_df, x = "Group", y = "Intensity", 
            fill="Group", facet.by = "Protein", add = "mean") +
  stat_summary(
    fun.data = mean_se,  # Calculate mean and standard error
    geom = "errorbar",   # Add error bars
    width = 0.25,       # Error bar width
    color = "black",
    alpha=0.7
  )+
  geom_jitter(
    width = 0.02, height=0.02,
    size = 1.8, alpha = 0.5, color="black") +

  facet_wrap(~Protein, ncol = 1)+
  stat_pvalue_manual(
    stat.test, 
    label = "p.signif", 
    y.position = c(2.8,2.35), size = 3.5
  )+ scale_fill_manual(values = c("#C3B3A1", "#125C5D"))+
  scale_colour_manual(values = c("#C3B3A1", "#125C5D"))+
  labs(x = "", 
       y = "Protein of Interest/\nHistone") +
  # theme_jb+
  theme_bw(base_family = "Roboto")+
  theme(
    text=element_text(size=12),
    strip.text = element_text(size=13),
    strip.background = element_rect(fill = scales::alpha("grey", 0.2)),
    legend.position = "none",
    
    panel.spacing = unit(2.5, "lines"),
    legend.box.margin = margin(-50, 0, 0, 0) 
  )

# save graph
ggsave(
  plot=p2,
  filename = file.path(figure_out_path, "Figure_4e_quantification.png"),
  width = 800,
  height = 2500,
  dpi = 350,
  units = 'px'
)