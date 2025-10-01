# Generating graphs of numbers of genes and DEGs for Hardy and Wood projects

library(magrittr)
library(patchwork)
library(tidyverse)
library(foreach)
library(doParallel)

# Change to local if needed:
project_dir<-here::here() # if running from project home # main project level i.e TDP43_paper_open_code level
figure_out_path<-file.path(project_dir, "Plots/Main_Figures/Figure_1")

# create figure output directory if not present
if (!dir.exists(figure_out_path)) {
  dir.create(figure_out_path, recursive = T)
}

analysis_colours<-c("#0c5394", "#e06666")

stage_colours<-c("#91D1C2", "#4B8D86")

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


# 1. Loading DGE and Splicing Data: ********************************************************************************************

## 1.1 Braak stage 3 / 4 - Team Wood Data:

### DEG Data:
#### Load differential gene expression results
dge_all <- read_csv(file.path(project_dir, "Data/Mid_stage/dge_all.csv"))

#### Select all DEGs from this - no collapsed analyses so don't double count. Keep significant genes as a vector.
dge_sig<-dge_all %>% filter(
  !tissue=="Collapse-All Brain", !tissue=="Collapse-Cortical", !tissue=="Collapse-Subcortical",
                            adj.P.Val<0.05) %>% distinct(SYMBOL) %>% pull()

#### Add in analysis type and change back to a tibble
wood_degs<-tibble(Genes=dge_sig, Analysis="Differential Expression")

#### Load splicing data:
wood_full_results<-readRDS(file = file.path(project_dir, "Data/Mid_stage/full_results_list_for_default_params.RDS"))

#### Load only significant results and gene names of those:
wood_sig_splice<-wood_full_results$all_results$significant_clusters_0.05_filter %>% 
  distinct(genes) %>% pull()

#### Add in analysis type and change back to a tibble
wood_sig_splice<-tibble(Genes=wood_sig_splice, Analysis="Differential Splicing")

### Combine the two analyses

wood_overlap_deg_ds<-rbind(wood_sig_splice, wood_degs)

#### set dataset name/braak stage
wood_overlap_deg_ds$Braak_stage<-"Mid-stage"

## 1.2 Braak stage 6 - Team Hardy Data:

### DEG Data:
#### Load differential gene expression results

hardy_deg_res<- read_csv(file.path(project_dir, "Data/Late_stage/dge_limma.csv"))

#### Select all DEGs from this - no collapsed analyses so don't double count. Keep significant genes as a vector.
dge_sig<-hardy_deg_res %>% filter(
  !tissue=="Collapse",
  adj.P.Val<0.05) %>% distinct(gene_name) %>% pull()

#### Add in analysis type and change back to a tibble
hardy_degs<-tibble(Genes=dge_sig, Analysis="Differential Expression")

#### Load splicing data:
hardy_full_results<-readRDS(file = file.path(project_dir, "Data/Late_stage/full_results_list_for_default_params.RDS"))

#### Load only significant results and gene names of those:
hardy_sig_splice<-hardy_full_results$all_results$significant_clusters_0.05_filter %>% 
  distinct(genes) %>% pull()

#### Add in analysis type and change back to a tibble
hardy_sig_splice<-tibble(Genes=hardy_sig_splice, Analysis="Differential Splicing")

### Combine the two analyses
hardy_overlap_deg_ds<-rbind(hardy_sig_splice, hardy_degs)

#### set dataset name/braak stage
hardy_overlap_deg_ds$Braak_stage<-"Late-stage"


## 1.3 Combine the two Braak stage datasets 
both_overlap<-rbind(hardy_overlap_deg_ds, wood_overlap_deg_ds)
both_overlap %<>% mutate(Braak_stage=factor(Braak_stage, levels=Braak_stage_order))

# 2. Common DEG plot ---- ********************************************************************************************

# Factor analyses and stagesL
both_overlap<-both_overlap %>%
  mutate(Analysis=factor(Analysis),
         Braak_stage=factor(Braak_stage)) 

# Find number of levels of each:
all_regulations <- levels(both_overlap$Analysis)
all_comparisons <- levels(both_overlap$Braak_stage)

# set colours
palette = setNames(analysis_colours, c("Differential Expression", "Differential Splicing"))

# create binary plot:

# list of genes that don't overlap between DEG and DS analyses in each stage
both_overlap_ij <- both_overlap %>% 
      dplyr::select(Analysis, Genes, Braak_stage) %>% 
      group_by(Analysis, Braak_stage) %>%
      distinct(Genes) %>% 
      ungroup() %>%
      dplyr::mutate(DEG = 1, Braak_stage = droplevels(Braak_stage))

    ylab<-"Analysis"

# create full list of genes that are either DE or DS in mid-stage dataset
    gene_list <- both_overlap_ij %>% 
     filter(Braak_stage=="Mid-stage") %>%
      dplyr::arrange(Analysis) %>% 
      dplyr::pull(Genes) %>% 
      unique()
    xlab <-paste0("Genes (n = ", length(gene_list), ")")

# binary plot for mid-stage:    
    dge_plot1 <- both_overlap_ij %>% filter(Braak_stage=="Mid-stage") %>%
      dplyr::mutate(Genes = droplevels(factor(Genes, levels = gene_list))) %>% 
      ggplot(aes(x = Genes, y = DEG, fill = Analysis)) +
      geom_col(width = 1
               ) +
      scale_x_discrete(expand = expansion(0)) +
      scale_y_continuous(expand = expansion(0))+
      scale_fill_manual(values = palette) +
      facet_grid(cols = vars(Braak_stage), rows = vars(Analysis), 
                 scales = "free_x",
                 drop = T, switch = "y") +
      labs(x = xlab, y = ylab) +
      theme_jb()+
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            legend.position = "top",
            strip.text = element_text(size = 28),
            strip.text.y = element_text(size = 28, angle = 0),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_text(vjust = 0.6),
      )+
      theme(
        strip.background = element_rect(fill = stage_colours[1]),
        strip.text.y = element_blank() ,
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.45, "lines")
      )+
      theme(text =  element_text(family = "Roboto", size=24))
    
# create full list of genes that are either DE or DS in late-stage dataset
    gene_list <- both_overlap_ij %>% 
      filter(Braak_stage=="Late-stage") %>%
      dplyr::arrange(Analysis) %>% 
      dplyr::pull(Genes) %>% 
      unique()
    xlab <-paste0("Genes (n = ", length(gene_list), ")")

    # binary plot for late stage:       
      dge_plot2 <- both_overlap_ij %>% filter(Braak_stage=="Late-stage") %>%
      dplyr::mutate(Genes = droplevels(factor(Genes, levels = gene_list))) %>% 
        ggplot(aes(x = Genes, y = DEG, fill = Analysis)) +
        geom_col(width = 1
        ) +
        scale_x_discrete(expand = expansion(0)) +
        scale_y_continuous(expand = expansion(0))+
        scale_fill_manual(values = palette) +
        facet_grid(cols = vars(Braak_stage), rows = vars(Analysis), 
                   scales = "free_x"
        )+
        labs(x = xlab, y = ylab)+
        theme_jb()+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              legend.position = "top",
              strip.text = element_text(size = 28),
              strip.text.y = element_text(size = 28, angle = 0),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title.y = element_text(vjust = 0.6),
        )+
        theme(
          strip.background = element_rect(fill = stage_colours[2]),
          strip.text.y = element_blank() ,
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.45, "lines")
        )+
        theme(text =  element_text(family = "Roboto", size=24))
      
      
# Number sig genes by stage ----

## mid stage

        number_sig_genes_mid<-
        both_overlap %>% 
        filter(Braak_stage=="Mid-stage") %>%
        mutate(Braak_stage=fct_drop(Braak_stage)) %>%
        dplyr::group_by(Braak_stage) %>% 
        dplyr::count(Analysis, .drop = F) %>% 
        ggplot(aes(x = Braak_stage, 
                   y=Analysis, fill = n)) +
        geom_tile(color = "black", show.legend = F, linewidth = 0.25, alpha=0.8) +
        geom_text(aes(label = n, color = n > 80), size = 10, show.legend = F, fontface = "bold") +
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
          xmin = 0.5,
          xmax =   1.5,
          ymin = 2.5,  # Adjust position above the plot
          ymax = 3,  # Height of the boxes
          fill = stage_colours[1],
          color = "black"
        ) +
        annotate(
          geom = "text",
          x = 1,
          y = 2.75,  # Centered inside the boxes
          label = Braak_stage_order[1],
          color = "black",
          size=8
        )+
        theme_jb()+
        theme(  axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                text = element_text(size=42),
                # axis.text.x = element_blank(),
                axis.ticks = element_blank(),
                plot.title = element_text(hjust = 0.5))+
        theme(text =  element_text(family = "Roboto"))
      
## late-stage
 
        number_sig_genes_late<-
        both_overlap %>% 
        filter(Braak_stage=="Late-stage") %>%
        mutate(Braak_stage=fct_drop(Braak_stage)) %>%
        dplyr::group_by(Braak_stage) %>% 
        dplyr::count(Analysis, .drop = F) %>% 
        ggplot(aes(x = Braak_stage, 
                   y=Analysis, fill = n)) +
        geom_tile(color = "black", show.legend = F, linewidth = 0.25, alpha=0.8) +
        geom_text(aes(label = n, color = n > 80), size = 10, show.legend = F, fontface = "bold") +
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
          xmin = 0.5,
          xmax =   1.5,
          ymin = 2.5,  # Adjust position above the plot
          ymax = 3,  # Height of the boxes
          fill = stage_colours[2],
          color = "black"
        ) +
        annotate(
          geom = "text",
          x = 1,
          y = 2.75,  # Centered inside the boxes
          label = Braak_stage_order[2],
          color = "black",
          size=8
        )+
        theme_jb()+
        theme(  axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                text = element_text(size=42),
                # axis.text.x = element_blank(),
                axis.ticks = element_blank(),
                plot.title = element_text(hjust = 0.5))+
        theme(text =  element_text(family = "Roboto"))
      
    
# Combining and saving the binary plots
      ggsave(
        plot = (
          ((dge_plot1 + theme(plot.margin = unit(c(0,0,0,0), "pt"))) / 
            (dge_plot2+ theme(plot.margin = unit(c(35,0,0,0), "pt")))) + plot_layout(heights = c(1, 1))
         
          |
            (number_sig_genes_mid / number_sig_genes_late)+plot_layout(heights = c(3,3))
        ) + plot_layout(ncol = 2, widths = c(6, 1)),
        filename = file.path(figure_out_path, "Figure_1b.png"),
        width = 10000, height = 4000, dpi = 330, units = 'px'
      ) 
      
      
# The size of the number plots will be increased the in inkscape to match lines box sizes - and not to disrupt the size of the binary plots.