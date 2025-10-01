# Pathways ----  
# Only doing GO pathways that aren't reduced

library(magrittr)
library(patchwork)
library(tidyverse)
library(foreach)
library(doParallel)

# Change to local if needed:
project_dir<-here::here() # if running from project home
figure_out_path<-file.path(project_dir,"Plots/Figures/Extended_data")

# create figure output directory if not present
if (!dir.exists(figure_out_path)) {
  dir.create(figure_out_path, recursive = T)
}

analysis_colours<-c("#0c5394", "#e06666")


stage_colours<-c("#91D1C2", "#4B8D86")

Braak_stage_order<-c("Mid-stage", "Late-stage")

analysis_colours<-c("#0c5394", "#e06666")

Braak_stage_order<-c("Mid-stage", "Late-stage")

extrafont::font_import(paths = file.path(project_dir, "Data"), 
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

# Load Wood/Mid Stage DEG pathways:
# Mid Stage    
Mid_DEG_GO_merge <-
  readRDS(
    file.path(project_dir, "Data/Mid_stage/dge_formula_res.rds")
  )        

Mid_DEG_GO_merge_paths<-Mid_DEG_GO_merge@compareClusterResult %>% filter(
  !tissue=="Collapse-All Brain") %>% distinct(Description) %>% pull()

mid_degs_paths<-tibble(Pathways=Mid_DEG_GO_merge_paths, Analysis="Differential Expression")    

Mid_DS_GO_merge <-
  readRDS(
    file.path(project_dir, "Data/Mid_stage/Wood_Leaf_final_formula_GO_merge.rds")
  )

Mid_DS_GO_merge_paths<-Mid_DS_GO_merge@compareClusterResult %>% distinct(Description) %>% pull()

mid_ds_paths<-tibble(Pathways=Mid_DS_GO_merge_paths, Analysis="Differential Splicing")    

# Combine

mid_overlap_paths<-rbind(mid_ds_paths, mid_degs_paths)
mid_overlap_paths$Braak_stage<-"Mid-stage"

# Late Stage
Late_DEG_GO_merge <-
  readRDS(
    file.path(project_dir, "Data/Late_stage/dge_formula_GO.rds")
  )

Late_DEG_GO_merge_paths<-Late_DEG_GO_merge@compareClusterResult %>% filter(
  !tissue=="Collapse") %>% distinct(Description) %>% pull()

Late_degs_paths<-tibble(Pathways=Late_DEG_GO_merge_paths, Analysis="Differential Expression")   

Late_DS_GO_merge <-
  readRDS(
    file.path(project_dir, "Data/Late_stage/Hardy_Leaf_final_formula_GO_merge.rds")
  )

Late_DS_GO_merge_paths<-Late_DS_GO_merge@compareClusterResult %>% distinct(Description) %>% pull()    

Late_ds_paths<-tibble(Pathways=Late_DS_GO_merge_paths, Analysis="Differential Splicing")   

# Combine

Late_overlap_paths<-rbind(Late_ds_paths, Late_degs_paths)
Late_overlap_paths$Braak_stage<-"Late-stage"


both_overlap<-rbind(mid_overlap_paths, Late_overlap_paths)

both_overlap %<>% mutate(Braak_stage=factor(Braak_stage, levels=Braak_stage_order))

# Number sig pathways plot ----
number_sig_pathways<-both_overlap %>% 
  dplyr::group_by(Braak_stage) %>% 
  dplyr::count(Analysis, .drop = F) %>% 
  ggplot(aes(x = Braak_stage, 
             y=Analysis, fill = n)) +
  geom_tile(color = "black", show.legend = F, linewidth = 0.25, alpha=0.8) +
  geom_text(aes(label = n, color = n > 80), size = 5, show.legend = F, fontface = "bold") +
  scale_x_discrete(expand = expansion()) +
  scale_y_discrete(expand = expansion(), limits = rev) +
  scale_fill_distiller(palette = "Greys", direction = 2 ) +
  scale_color_manual(values = c("black", "white")) +
  # facet_wrap(~Braak_stage, shrink = T)+
  scale_x_discrete(position = "top", expand = expansion()) + # Shift x-axis above
  coord_equal() +
  # labs(x = "Analysis", y = "Dataset",  
  #      title = "Number of Significant Genes") +
  labs(x = "", y = "",  
       title = "") +
  annotate(
    geom = "rect",
    xmin = 1:length(Braak_stage_order) - 0.5,
    xmax = 1:length(Braak_stage_order) + 0.5,
    ymin = length(Braak_stage_order) + 0.5,  # Adjust position above the plot
    ymax = length(Braak_stage_order) + 0.85,  # Height of the boxes
    fill = stage_colours,
    color = "black"
  ) +
  annotate(
    geom = "text",
    x = 1:length(Braak_stage_order),
    y = length(Braak_stage_order)+0.675,  # Centered inside the boxes
    label = Braak_stage_order,
    color = "black",
    size=3.25
  )+
  theme_jb()+
  theme(  axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          text = element_text(size=15),
          # axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5))+
  theme(text =  element_text(family = "Robot"))




both_overlap<-both_overlap %>%
  mutate(Analysis=factor(Analysis),
         Braak_stage=factor(Braak_stage)) 



all_regulations <- levels(both_overlap$Analysis)
all_comparisons <- levels(both_overlap$Braak_stage)
palette = setNames(analysis_colours, c("Differential Expression", "Differential Splicing"))

# Get percentage ovlerlap of pathways
ds_genes<-both_overlap %>% filter(Analysis=="Differential Splicing") %>% distinct(Pathways, .keep_all=T)      
deg_genes<-both_overlap %>% filter(Analysis=="Differential Expression") %>% distinct(Pathways, .keep_all=T) 
all_genes<-both_overlap %>% distinct(Pathways, .keep_all=T) 
# get percentage of overlap from total
length(ds_genes$Pathways[ds_genes$Pathways %in% deg_genes$Pathways])/length(all_genes$Pathways)*100




both_overlap_ij <- both_overlap %>% 
  # dplyr::filter(sig_reg == regulation, test == comparison) %>% 
  dplyr::select(Analysis, Pathways, Braak_stage) %>% 
  group_by(Analysis, Braak_stage) %>%
  distinct(Pathways) %>% 
  ungroup() %>%
  dplyr::mutate(DEG = 1, Braak_stage = droplevels(Braak_stage))

ylab<-"Analysis"

Pathway_list <- both_overlap_ij %>% 
  filter(Braak_stage=="Mid-stage") %>%
  dplyr::arrange(Analysis) %>% 
  dplyr::pull(Pathways) %>% 
  unique()
xlab <-paste0("Pathways (n = ", length(Pathway_list), ")")

path_plot1 <- both_overlap_ij %>% filter(Braak_stage=="Mid-stage") %>%
  dplyr::mutate(Pathways = droplevels(factor(Pathways, levels = Pathway_list))) %>% 
  ggplot(aes(x = Pathways, y = DEG, fill = Analysis)) +
  geom_col(width = 1
  ) +
  scale_x_discrete(expand = expansion(0)) +
  scale_y_continuous(expand = expansion(0))+
  scale_fill_manual(values = palette) +
  facet_grid(cols = vars(Braak_stage), rows = vars(Analysis), 
             scales = "free_x",
             # drop = T, 
             # switch = "y"
  ) +
  labs(x = xlab, y = ylab)+
  theme_jb()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "top",
        strip.text = element_text(size = 12),
        strip.text.y = element_text(size = 14, angle = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(vjust = 0.6),
        axis.title = element_text(size = 12))+
  theme(
    strip.background = element_rect(fill = stage_colours[1]),
    strip.text.y = element_blank() ,
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.45, "lines")
  )+
  theme(text =  element_text(family = "Roboto"))

Pathway_list <- both_overlap_ij %>% 
  filter(Braak_stage=="Late-stage") %>%
  dplyr::arrange(Analysis) %>% 
  dplyr::pull(Pathways) %>% 
  unique()
xlab <-paste0("Pathways (n = ", length(Pathway_list), ")")

path_plot2 <- both_overlap_ij %>% filter(Braak_stage=="Late-stage") %>%
  dplyr::mutate(Pathways = droplevels(factor(Pathways, levels = Pathway_list))) %>% 
  ggplot(aes(x = Pathways, y = DEG, fill = Analysis)) +
  geom_col(width = 1
  ) +
  scale_x_discrete(expand = expansion(0)) +
  scale_y_continuous(expand = expansion(0))+
  scale_fill_manual(values = palette) +
  facet_grid(cols = vars(Braak_stage), rows = vars(Analysis), 
             scales = "free_x",
             # drop = T, 
             # switch = "y"
  ) +
  labs(x = xlab, y = ylab)+
  theme_jb()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        # Increase vertical spacing
        legend.position = "top",
        strip.text = element_text(size = 12),
        strip.text.y = element_text(size = 14, angle = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(vjust = 0.6),
        axis.title = element_text(size = 12))+
  theme(
    strip.background = element_rect(fill = stage_colours[2]),
    strip.text.y = element_blank() ,
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.45, "lines")
    #  panel.spacing.y = unit(2, "lines")
  )+
  theme(text =  element_text(family = "Roboto"))

ggsave(
  plot = ((path_plot1 / path_plot2) | 
            number_sig_pathways) + plot_layout(ncol = 2, widths = c(6, 1)), 
  file.path(figure_out_path, "Ext_data_fig_1.png"),
  width = 10000, height = 3500, dpi = 600, units = 'px'
) # changed scale to 80% in slides/pdf



