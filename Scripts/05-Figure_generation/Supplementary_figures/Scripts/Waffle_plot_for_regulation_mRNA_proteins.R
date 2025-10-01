library(dplyr)
library(ggplot2)
library(ggfittext)
# 0. Setup ----
figure_out_path<-"/home/jbrenton/TDP43_PD_Paper_Figures/Drafting_figures_code/Supplementary_figures/"

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



# load rosa ma sig genes
ma_paper_sig_spliced<- readxl::read_excel("/home/jbrenton/TDP43_PD_Paper_Figures/Drafting_figures_code/Figure_5/Data/41586_2022_4424_MOESM5_ESM(2).xlsx")

ma_paper_sig_spliced_genes<-ma_paper_sig_spliced %>% select(Genes) %>% distinct()


# Waffle config
n_cols <- 20


# Prepare data and assign x/y per group
df <- ranked_df %>%
  filter(is_tdpkd == TRUE) %>%
  mutate(in_ma_sig=case_when(
    Genes %in% ma_paper_sig_spliced_genes$Genes ~ T,
    .default=F)) %>%
  mutate(tdp_kd_regulation = case_when(
    DE_Log2FC > 0 & prot_Log2FC > 0 ~ "DE_up / Prot_up",
    DE_Log2FC > 0 & prot_Log2FC < 0 ~ "DE_up / Prot_down",
    DE_Log2FC < 0 & prot_Log2FC > 0 ~ "DE_down / Prot_up",
    DE_Log2FC < 0 & prot_Log2FC < 0 ~ "DE_down / Prot_down"
  )) %>%
  arrange(tdp_kd_regulation, Genes) %>%
  group_by(tdp_kd_regulation) %>%
  mutate(
    group_index = row_number() - 1,
    group_size = n(),
    group_rows = ceiling(group_size / n_cols),
    x = group_index %% n_cols,
    y = -1 * (group_index %/% n_cols)  # flip Y to go top-down
  ) %>%
  ungroup()

# Assign row block offsets for each group
group_blocks <- df %>%
  count(tdp_kd_regulation, name = "group_size") %>%
  mutate(group_rows = ceiling(group_size / n_cols)) %>%
  mutate(row_offset = cumsum(lag(group_rows, default = 0)))

# Join row offsets
df <- df %>%
  left_join(group_blocks, by = "tdp_kd_regulation") %>%
  mutate(y = y - row_offset)

# Plot
p1<-ggplot(df, aes(x = x, y = y, fill = tdp_kd_regulation, label = Genes, group=in_ma_sig)) +
  geom_tile(color = "white", size = 0.5, alpha=0.8) +
  geom_tile(data = df %>% filter(in_ma_sig),
            aes(x = x, y = y),
            fill = NA, color = "black", linewidth = 1)+
  geom_fit_text(min.size = 2, 
                # grow = TRUE,
                reflow = TRUE) +
  # theme_void() +

  labs(fill = "TDP-KD Regulation")+
  scale_fill_manual(labels = c(
    "Downregulated mRNA,\nDownregulated Protein",
    "Downregulated mRNA,\nUpregulated Protein",
  
    "Upregulated mRNA,\nDownregulated Protein",
    "Upregulated mRNA,\nUpregulated Protein"
  ),values=
                      c("#0571B0", "#9B5EA3", "#836FA8", "#CA0020")
  )+
  # theme_jb()+
  theme_void()+
  theme(legend.position = "bottom",
        legend.key.size = unit(0.001, 'cm'),
        legend.title = element_blank(),
        legend.margin = margin(t = -30, b = 5, unit = "pt")
        )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 4),
      keyheight = unit(0.2, "cm"),
      keywidth = unit(0.5, "cm")
    )
  )
   


ggsave(
  plot=p1,
  filename = file.path(figure_out_path, "Waffle_name_plot_tdp_prot_mrna_reg_supp.png"),
  width = 4000,
  height = 3700,
  dpi = 350,
  units = 'px'
)

