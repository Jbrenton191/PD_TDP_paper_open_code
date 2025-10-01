## _________________________________________________
##
## Covariate identifcation of Hardy snRNA-seq data
##
## Aim: to develop a complete covariate identification for the Hardy snRNA-seq
## dataset. Focus are on the sources of variation to identify which covariates
## need to be controlled for in downstream analyses.
##
## Author: Jonathan Brenton
##
## Contributors: Aine Fairbrother-Browne, Guillermo Rocamora Pérez
##
## Date Created:
##
## Copyright (c) Author, year
##
## Email:
## _________________________________________________
##
## Notes:
##
## Changelog:
##
## _________________________________________________

# --- 0. Setup -----------------------------------------------------------------

## Import libraries
library(tidyverse)
library(data.table)
library(tidyr)
library(magrittr)
library(vroom)
library(DESeq2)
library(tibble)
library(variancePartition)
library(parallel)
library(factoextra)
library(parallel)
library(ggsci)
library(Hmisc)
library(janitor)
library(caret)
library(scales)
library(patchwork)

options(dplyr.summarise.inform = FALSE)
options(lifecycle_verbosity = "warning")
options(readr.show_progress = F)
options(readr.show_col_types = F)

## G: Set main path
main_path <- "/home/grocamora/RytenLab-Research/18-DGE_limma_dream/"

## G: Source helper files
source(file.path(main_path, "R/run_stat_test_for_all_cols_get_res.R"))

## G: Additional paths
alignment_files_path <- file.path(main_path, "data/Hardy/alignment_run/multiqc_data/")
additional_seq_qc_files_path <- file.path(main_path, "data/Hardy/multiqc_qc/multiqc_data/")
metadata_path <- file.path(main_path, "metadata/Hardy/updated_metadata_w_PMI.csv")
txi.salmon_path <- file.path(main_path, "data/Hardy/txi.salmon.rds")

results_path <- file.path(main_path, "results/Hardy_NewCovs/")

multiqc_fastqc_path <- file.path(alignment_files_path, "multiqc_fastqc.txt")
multiqc_fastp_path <- file.path(alignment_files_path, "mqc_fastp_filtered_reads_plot_1.txt")
multiqc_salmon_path <- file.path(alignment_files_path, "multiqc_salmon.txt")
multiqc_star_path <- file.path(alignment_files_path, "multiqc_star.txt")
multiqc_rseqc_path <- file.path(alignment_files_path, "multiqc_rseqc_read_distribution.txt")

multiqc_picard_rna_seq_path <- file.path(additional_seq_qc_files_path, "multiqc_picard_RnaSeqMetrics.txt")
multiqc_picard_gc_bias_path <- file.path(additional_seq_qc_files_path, "multiqc_picard_gcbias.txt")
multiqc_picard_insert_path <- file.path(file.path(additional_seq_qc_files_path, "multiqc_picard_insertSize.txt"))
multiqc_picard_alignment_path <- file.path(additional_seq_qc_files_path, "multiqc_picard_AlignmentSummaryMetrics.txt")
multiqc_qualimap_alignment_path <- file.path(additional_seq_qc_files_path, "qualimap_rnaseq_genome_results.txt")
multiqc_genomic_origin_path <- file.path(additional_seq_qc_files_path, "mqc_qualimap_genomic_origin_1.txt")

varPart_path <- file.path(results_path, "variance_partition_output.rds")
pca_path <- file.path(results_path, "pca_output.rds")
meta_path <- file.path(results_path, "all_brain_areas_covariates.rds")
colinearDisplot_path <- file.path(results_path, "colinearDistPlot.rds")
meta_pre_path <- file.path(results_path, "metadata_without_removal_of_correlated_terms.rds")
varPartCorPlots_path <- file.path(results_path, "varPartCorPlots.rds")

dir.create(results_path, recursive = T, showWarnings = F)

############################################################################## #
# 1. Load multiqc data & metadata ----

## 1.1 Load original run multiqc ----
### fastqc
fastqc <- readr::read_delim(multiqc_fastqc_path)
fastqc <- fastqc[grep(pattern = "trimmed_1", x = fastqc$Sample), ]
fastqc$bxp_id <- gsub(x = fastqc$Sample, pattern = "(^BX.*_[0-9]+)_S.*", replacement = "\\1")

### fastp
fastp <- readr::read_delim(multiqc_fastp_path)
fastp$bxp_id <- gsub(x = fastp$Sample, pattern = "(^BX.*_[0-9]+)_S.*", replacement = "\\1")
fastp <- dplyr::select(fastp, -Sample)

### Salmon
salmon <- readr::read_delim(multiqc_salmon_path)
salmon$bxp_id <- gsub(x = salmon$Sample, pattern = "(^BX.*_[0-9]+)_S.*", replacement = "\\1")
salmon <- dplyr::select(salmon, -Sample)

### STAR
star <- readr::read_delim(multiqc_star_path)
star$bxp_id <- gsub(x = star$Sample, pattern = "(^BX.*_[0-9]+)_S.*", replacement = "\\1")
star <- dplyr::select(star, -Sample)

### rseqc
rseqc_rd <- readr::read_delim(multiqc_rseqc_path)
rseqc_rd$bxp_id <- gsub(x = rseqc_rd$Sample, pattern = "(^BX.*_[0-9]+)_S.*", replacement = "\\1")
rseqc_rd <- dplyr::select(rseqc_rd, -Sample)

## 1.2 Load additional QC run multiqc data ----
### picard
picard_rna_seq <- readr::read_delim(multiqc_picard_rna_seq_path)
picard_rna_seq$bxp_id <- gsub(x = picard_rna_seq$Sample, pattern = "(^BX.*_[0-9]*)_.*", replacement = "\\1")

picard_gc_bias <- readr::read_delim(multiqc_picard_gc_bias_path)
picard_gc_bias$bxp_id <- gsub(x = picard_gc_bias$Sample, pattern = "(^BX.*_[0-9]*)_.*", replacement = "\\1")

picard_insert <- readr::read_delim(multiqc_picard_insert_path)
picard_insert$bxp_id <- gsub(x = picard_insert$Sample, pattern = "(^BX.*_[0-9]*)_.*", replacement = "\\1")

picard_alignment <- readr::read_delim(multiqc_picard_alignment_path)
picard_alignment$bxp_id <- gsub(x = picard_alignment$Sample, pattern = "(^BX.*_[0-9]*)_.*", replacement = "\\1")

picard <- picard_rna_seq %>%
  dplyr::left_join(picard_gc_bias, by = "bxp_id") %>%
  dplyr::left_join(picard_insert, by = "bxp_id") %>%
  dplyr::left_join(picard_alignment, by = "bxp_id") %>%
  dplyr::select(bxp_id, where(is.numeric))

### qualimap
qualimap_alignment <- readr::read_delim(multiqc_qualimap_alignment_path)
qualimap_alignment$bxp_id <- gsub(x = qualimap_alignment$Sample, pattern = "(^BX.*_[0-9]*)_.*", replacement = "\\1")

qualimap_genomic_origin <- readr::read_delim(multiqc_genomic_origin_path)
qualimap_genomic_origin$bxp_id <- gsub(x = qualimap_genomic_origin$Sample, pattern = "(^BX.*_[0-9]*)_.*", replacement = "\\1")

qualimap <- dplyr::left_join(qualimap_alignment, qualimap_genomic_origin, by = "bxp_id")
qualimap <- dplyr::select(qualimap, bxp_id, where(is.numeric))

## 1.3 Merging ----
run_qc <- fastqc %>%
  dplyr::left_join(fastp, by = "bxp_id") %>%
  dplyr::left_join(star, by = "bxp_id") %>%
  dplyr::left_join(salmon, by = "bxp_id") %>%
  dplyr::left_join(rseqc_rd, by = "bxp_id") %>%
  dplyr::select(bxp_id, where(is.numeric))

alignment_qc <- dplyr::left_join(picard, qualimap, by = "bxp_id")
full_qc <- dplyr::left_join(alignment_qc, run_qc, by = "bxp_id") %>%
  dplyr::select(-any_of(c("PF_ALIGNED_BASES.y", "total_reads.y"))) %>%
  dplyr::rename(any_of(c(PF_ALIGNED_BASES = "PF_ALIGNED_BASES.x",
                         total_reads = "total_reads.x"))) %>%
  dplyr::mutate(bxp_id_full = bxp_id,
                bxp_id = as.numeric(sub(x = bxp_id, pattern = "^BX.*_([0-9].*)$", replacement = "\\1")))

## 1.4 Load metadata ----
sample_metadata <- readr::read_delim(metadata_path)
full_metadata_seq_sample <- sample_metadata %>% dplyr::left_join(full_qc, by = "bxp_id")

meta <- full_metadata_seq_sample %>%
  dplyr::select(-c(order,
                   bxp_id,
                   duration_type,
                   sampled_for_rin,
                   selected_for_rna_extraction,
                   all_regions_sampled,
                   brain_type,
                   sample_shared_with_asap_groups)) %>%
  dplyr::select(-c(dementia,
                   path_autopsy_main_dx,
                   apoe,
                   thal,
                   b_b,
                   cerad,
                   a,
                   b,
                   c,
                   abc,
                   lb_b,
                   mc_keith,
                   artag,
                   path_ad_level,
                   path_dg))
meta <- meta %>% dplyr::mutate(duration = ifelse(is.na(duration), 0, duration))

############################################################################## #
# 2. Aine's script ----

## --- 2.1 Load and wrangle data -----------------------------------------------
### A: if a variable has 1 NA, impute it with median or most common value check
### NA counts in each col
na.count = data.frame(sapply(meta, function(y) sum(length(which(is.na(y)))))) %>%
  tibble::rownames_to_column("col") %>%
  `colnames<-`(c("col", "na.count")) %>%
  tibble::tibble()

### A: get cols with a single NA
cols.with.1.na = na.count %>%
  dplyr::filter(na.count==1) %>%
  dplyr::pull(col)

message(paste("Imputing", length(cols.with.1.na), "meta cols with a single NA value."))

### A: if a column has a single NA, impute it, as this is not likely to have an
### appreciable effect on whether the covariate is deemed important
for(col in cols.with.1.na){

  ### A: if the col is numeric, fill with median
  if(is.numeric(meta[[col]])){
    meta[col][is.na(meta[col])] = median(meta[[col]], na.rm=TRUE)

    ### A: if the col is not numeric, fill with most common val
  } else{
    meta[col][is.na(meta[col])] = sort(table(meta[col]), decreasing=TRUE) %>% names() %>% .[1]
  }
}

### J: one sample with no ath PMI data so will impute
meta$path_PMI_hours[is.na(meta$path_PMI_hours)]<-median(meta$path_PMI_hours, na.rm=TRUE)
saveRDS(meta, meta_pre_path)

## --- 2.2 Normalise counts using DESeq2 ---------------------------------------
### A: running this normalisation procedure as recommended in the
### variancePartition vignette

### J: Using DESeq code as using tximport code is for isoform level analysis
txi.salmon <- readRDS(txi.salmon_path)

### G: Test for identical columns
identical(colnames(txi.salmon$counts), meta$bxp_id_full)
unordered_rows <- which(!colnames(txi.salmon$counts) == meta$bxp_id_full)

#### G: A total of four rows are not ordered between txi.salmon$counts and
#### metadata. Metadata needs to be sorted
meta <- meta[sapply(meta$bxp_id_full, function(x) grep(x, colnames(txi.salmon$counts))), ]

stopifnot(identical(colnames(txi.salmon$counts), meta$bxp_id_full))

### A: create DESeq2 object from pseudobulk count matrix and metadata
dds = DESeq2::DESeqDataSetFromTximport(txi = txi.salmon,
                                       colData = meta,
                                       design = ~ 1)

### A: estimate library size correction scaling factors
dds = DESeq2::estimateSizeFactors(dds)

### A: identify genes that pass expression cutoff - then filter out genes that
### are extremely lowly expressed here, the standard cut-off set is fpm>1 in 50%
### or more of samples
isexpr = rowSums(DESeq2::fpm(dds)>1) >= 0.5 * ncol(dds)

### A: compute log2 Fragments Per Million
quantlog.counts = log2(DESeq2::fpm(dds)[isexpr,] + 1)

## --- 2.3 Deal with co-linearity ---------------------------------------------------------------------------------------------------
### 2.3a filter the meta table to remove unusable variables ----

#### A: get vars that can be included in the variancePartition model deal with
#### NAs, deal with no variance, make sure discrete/categoricals are set to char
#### and continuous set to numeric
meta.varpar.in = meta %>%
  tibble::column_to_rownames("bxp_id_full") %>%
  # A: remove meta variables with NA values - variancePartition can't handle these
  dplyr::select_if(~ !any(is.na(.))) %>%
  # A: remove meta variables with no variance, i.e. they have the same value for all samples
  dplyr::select(where(~n_distinct(.) > 1)) %>%
  # A: converting numeric/logical categorical cols to character so that variancePartition interprets them as categorical
  dplyr::mutate(bxp_plate = as.character(bxp_plate)) %>%
  # A: ensure cols that should be numeric are numeric - value magnitude in these vars is meaningful
  dplyr::mutate(rin_bxp = as.numeric(rin_bxp),
                dod = as.numeric(dod)) %>%
  # G: Remove other columns
  dplyr::select(-case_id)

### 2.3b deal with colinearity involving numeric variables ----
#### A: get names of numeric (num) cols
num.cols = meta.varpar.in %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

#### A: get names of categorical (cat) cols
cat.cols = meta.varpar.in %>%
  dplyr::select(where(is.character) | where(is.factor), -any_of(c("group", "Group", "sample_id"))) %>%
  colnames()

#### A: calculate correlation matrix
cor.mat = meta.varpar.in %>%
  dplyr::select(all_of(num.cols)) %>%
  # center and scale before calculating correlations
  dplyr::mutate(across(where(is.numeric), .fns = ~scale(., center = T, scale = T))) %>%
  tibble::tibble() %>%
  # convert to correlation matrix
  as.matrix() %>%
  cor(., method = "spearman")

#### A: use the R package 'caret' to find and remove co-linear variables
caret.detected.colins = cor.mat %>%
  # A: pass correlation matrix to caret, and calculate correlations, re-calculating each time a co-linear one is removed (selected with exact=T)
  caret::findCorrelation(., cutoff = 0.7, names=TRUE, exact=T, verbose=TRUE)

#### A: remove the colinear variables identified by caret from the main meta dataframe
meta.varpar.in = meta.varpar.in %>%
  dplyr::select(-any_of(c(
    # remove caret-identified
    caret.detected.colins
  )))

### 2.3c now deal with colinearity involving non-numeric variables ----

#### A: numeric vs. categorical and categorical vs. categorical to do this,
#### utilise two different stat tests chi-squared test for categorical vs.
#### categorical k-s chi-squared for numeric vs. categorical get names of
#### numeric (num) cols
#### - this will be reduced now that numeric colinearity has been dealt with
num.cols = meta.varpar.in %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

#### A: chi-squared between categoricals
chi.out = meta.varpar.in %>%
  dplyr::select(all_of(cat.cols)) %>%
  run_stat_test_for_all_cols_get_res(in_df=., stat_test="chisq", colnames1 = colnames(.), colnames(.)) %>%
  dplyr::arrange(p)

#### A: kruskal-wallace chi-square between numerics and categoricals
kruskal.out = meta.varpar.in %>%
  run_stat_test_for_all_cols_get_res(in_df=., stat_test="kruskal", colnames1 = num.cols, colnames2 = cat.cols) %>%
  dplyr::filter(var1 %in% num.cols) %>%
  dplyr::filter(var2 %in% cat.cols) %>%
  dplyr::arrange(p)

all = dplyr::bind_rows(
  chi.out,
  kruskal.out
) %>%
  dplyr::arrange(p) %>%
  dplyr::mutate(stat = abs(stat))

#### A: define Ns for P-value adjustment
n_numeric = length(colnames(cor.mat))
n_categorical = length(unique(c(all$var1, all$var2)))

message("Make manual decision on the cat vs. cat/ cat vs. num colinearity...")

cor_table <- all %>%
  dplyr::mutate(n_tests = case_when(
    test_run == "chisq" ~ n_categorical*n_categorical,
    test_run == "kruskal" ~ n_categorical*n_numeric,
  )) %>%
  dplyr::mutate(is_sig = case_when(
    p<(0.05/n_tests) ~ p,
    TRUE ~ NA
  )) %>%
  dplyr::filter(!is.na(is_sig))

cor_table %>% print()

# #### G: There seems to be a considerable correlation colinearity between
# #### Seq_batch, bxp_plate and many other categorical variables (17 in total with
# #### p < 1e-15). For the moment, I will remove "Seq_batch", "bxp_plate",
# #### "direct_cause_death" and "resequenced"
#
# # filtered_vars <- c("Seq_batch", "bxp_plate", "direct_cause_death", "resequenced")
# filtered_vars <- c("bxp_plate", "resequenced", "sequencing_location", "downsampled",
#                    "Too Many N", "direct_cause_death", "STRAND_BALANCE", "PF_READS_IMPROPER_PAIRS",
#                    "insertion_length", "PCT_R2_TRANSCRIPT_STRAND_READS", "AVG_POS_3PRIME_SOFTCLIP_LENGTH",
#                    "frag_length_sd")

#### G: I developed an algorithm to remove as many covariates as possible:
####
#### i) Select only the covariates with p-values <= 1e-15 (can be modified).
####
#### ii) Remove all covariates that are colinear with the most common covariate
#### in the table.
####
#### iii) Repeat "ii" until no more covariates are left in the colinear table.
collinear_p_limit <- 1e-15
collinear_filter <- "bottom"
tmp_cor_table <- cor_table %>%
  rowwise() %>%
  dplyr::mutate(test = paste0(sort(c(var1, var2)), collapse = "-")) %>%
  dplyr::distinct(test, .keep_all = T) %>%
  dplyr::filter(p <= collinear_p_limit)
filtered_vars <- c()
top_vars <- c(tmp_cor_table$var1, tmp_cor_table$var2) %>% table() %>% sort(decreasing = T)

while(length(top_vars) > 0){
  if(collinear_filter == "bottom"){
    filtered_vars_iter <- tmp_cor_table %>%
      dplyr::filter(var1 == names(top_vars)[1] | var2 == names(top_vars)[1]) %>%
      dplyr::select(var1, var2) %>%
      tidyr::pivot_longer(c(var1, var2)) %>%
      dplyr::filter(value != names(top_vars)[1]) %>%
      dplyr::pull(value) %>%
      unique()
    filtered_vars <- c(filtered_vars, filtered_vars_iter)
  }else{
    filtered_vars <- c(filtered_vars, names(top_vars)[1])
  }

  tmp_cor_table <- tmp_cor_table %>% dplyr::filter(!var1 %in% filtered_vars & !var2 %in% filtered_vars)
  top_vars <- c(tmp_cor_table$var1, tmp_cor_table$var2) %>% table() %>% sort(decreasing = T)
}

meta.varpar.in = meta.varpar.in %>%
  dplyr::select(-all_of(filtered_vars))

message("After colinearity detection, ", ncol(meta.varpar.in), " variables remain")

### 2.3d Visualise collinearity ----

#### A: plot colinearity for numeric vars - heatmap
numnum.heatmap = cor.mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("var1") %>%
  tidyr::pivot_longer(2:ncol(.), names_to="var2", values_to="corr") %>%
  dplyr::mutate(label = case_when(
    corr>=0.7 ~ "*",
    TRUE ~ ""
  )) %>%
  ggplot(data=., aes(x=reorder(var1, corr), y=reorder(var2, corr), fill=abs(corr))) +
  geom_tile() +
  labs(fill="Rho", x="", y="") +
  geom_text(aes(label=label), hjust=0.5, vjust=0.75, size=4) +
  scale_fill_distiller(palette = "Spectral")

#### A: plot colinearity for categorical vs categorical vars - heatmap
catcat.heatmap = all %>%
  dplyr::filter(test_run=="chisq") %>%
  dplyr::mutate(label = case_when(
    (p<(0.05/n_categorical*n_categorical)) ~ "**",
    (p<0.05) & (p>(0.05/n_categorical*n_categorical)) ~ "*",
    TRUE ~ ""
  )) %>%
  ggplot(data=., aes(x=reorder(var1, p), y=reorder(var2, p), fill=stat)) +
  geom_tile() +
  labs(fill="Chi-squared", x="", y="") +
  geom_text(aes(label=label), hjust=0.5, vjust=0.75, size=4) +
  scale_fill_distiller(palette = "Spectral") +
  facet_grid(~test_run)

#### A: plot colinearity for categorical vs numeric vars - heatmap
catnum.heatmap = all %>%
  dplyr::filter(test_run=="kruskal") %>%
  dplyr::mutate(label = case_when(
    (p<(0.05/(n_categorical*n_numeric))) ~ "**",
    (p<0.05) & (p>(0.05/(n_categorical*n_numeric))) ~ "*",
    TRUE ~ ""
  )) %>%
  ggplot(data=., aes(x=reorder(var1, p), y=reorder(var2, p), fill=stat)) +
  geom_tile() +
  labs(fill="K-W chi-squared", x="", y="") +
  geom_text(aes(label=label), hjust=0.5, vjust=0.75, size=4) +
  scale_fill_distiller(palette = "Spectral") +
  facet_grid(~test_run)

#### A: plot colinearity of numeric vars - distribution
colin.distro = cor.mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("var1") %>%
  tidyr::pivot_longer(2:ncol(.), names_to="var2", values_to="corr") %>%
  ggplot(data=., aes(x=abs(corr))) +
  geom_histogram(bins = 35) +
  geom_vline(xintercept=0.7, linetype=2, alpha=0.75) +
  labs(x="Rho", y="Count")

# make joint colinearity fig and export
colin.fig = ((numnum.heatmap + colin.distro) / (catcat.heatmap + catnum.heatmap)) +
  plot_layout(heights=c(1,0.35))
colin.fig %>% saveRDS(colinearDisplot_path)

## --- 2.4 Implement variancePartition ----------------------------------------------------------------------------------------------
### 2.4a running variancePartition on the non-colinear variables ----

#### A: scale and center variables prior to running variancePartition
meta.varpar.in = meta.varpar.in %>%
  dplyr::mutate(across(where(is.numeric), .fns = ~scale(., center = T, scale = T))) %>%
  tibble::as_tibble(rownames = "bxp_id_full")

#### J: error using reformulate because of "side (R or L)" so:
meta.varpar.in %<>% dplyr::rename('brain_hem'="side (R or L)")

#### J: need to fix special characters and spaces so reformulate can work
colnames(meta.varpar.in)<-gsub(pattern = "%", replacement = "percent_", colnames(meta.varpar.in))
colnames(meta.varpar.in)<-gsub(pattern = " ", replacement = "_", colnames(meta.varpar.in))

#### A: save pre-variancePartition metadata (meta.varpar.in)
meta.varpar.in %>% saveRDS(meta_path)

#### A: generate design formula to feed into variancePartition
vp.formula = meta.varpar.in %>%
  # A: get categorical covariates and wrap them in discrete (1|) notation
  dplyr::select(where(is.character)) %>%
  # A: remove patient and group, as we don't want these included in the covariates
  dplyr::select(-any_of(c("sample_id", "Group", "group", "bxp_id_full"))) %>%
  colnames() %>%
  paste0("(1|", ., ")") %>%
  # A: get numerical covariates and keep them in regular continuous notation
  c(., meta.varpar.in %>%
      dplyr::select(where(is.numeric)) %>%
      colnames()) %>%
  reformulate()

message(paste(
  "Running variancePartition with the following formula:",
  paste0(vp.formula, collapse = " ")
))

caret::findLinearCombos(as.matrix(meta.varpar.in %>% dplyr::select(where(is.numeric))))

#### A: run variance partition
#### G: Add multiprocessing
param <- BiocParallel::SnowParam(32, "FORK", progressbar = TRUE)
varPart = variancePartition::fitExtractVarPartModel(quantlog.counts,
                                                    vp.formula,
                                                    meta.varpar.in,
                                                    BPPARAM = param)
varPart %>%
  variancePartition::sortCols() %>%
  variancePartition::plotVarPart(., label.angle=60)

# save object
varPart %>% saveRDS(varPart_path)

############################################################################## #
# 3. Implement PCA ----
res.pca <- prcomp(t(quantlog.counts), scale = T)
res.pca %>% saveRDS(pca_path)

### 3.1 Visualizarion of covariate assessment ----

#### A: Function to read in .rds outputs from 1b and generate useful
#### visualisations for covariate assessment this function generates a 5th .rds
#### object (02-varPartCorPlots_*) which contains visualisations that are then
#### displayed in the accompanying 02-..rmd

varPart <- readRDS(varPart_path)
pca <- readRDS(pca_path)
meta <- readRDS(meta_path)
colinearDistplot <- readRDS(colinearDisplot_path)
meta_pre <- readRDS(meta_pre_path)

#### A: plot varPart
varPart.plot = varPart %>%
  variancePartition::sortCols() %>%
  variancePartition::plotVarPart(., label.angle=60)

#### A: get residual summary stats
varpart = varPart %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  tibble::tibble() %>%
  tidyr::pivot_longer(2:ncol(.)) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(med=median(value),
                   max=max(value)) %>%
  dplyr::filter(name!="Residuals")

#### G: Highest PC to extract
highest_pc <- which(summary(pca)$importance[3, ] > 0.85)[1]

#### A: spearman rho values
cors.r.numericVars = pca %>%
  .[["x"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("bxp_id_full") %>%
  tibble::tibble() %>%
  .[, seq(highest_pc+1)] %>%
  dplyr::left_join(x=.,
                   y=meta %>% dplyr::select(bxp_id_full, where(is.numeric)),
                   by="bxp_id_full") %>%
  dplyr::select(-any_of(c("bxp_id_full", "sample_id"))) %>%
  dplyr::mutate(across(c(where(is.character)), as.factor),
                across(c(where(is.factor)), as.numeric)) %>%
  as.matrix() %>%
  Hmisc::rcorr(., type = "spearman") %>%
  .$r %>%
  as.data.frame() %>%
  dplyr::select(matches("^PC\\d+")) %>%
  # J: need to get rid of some of the meta variables from picard that being with PCT
  dplyr::select(-contains("PCT")) %>%
  tibble::rownames_to_column("meta_var") %>%
  dplyr::filter(!grepl("^PC\\d+", meta_var)) %>%
  tidyr::pivot_longer(2:ncol(.)) %>%
  dplyr::mutate(name = as.factor(as.numeric(gsub("PC", "", name)))) %>%
  dplyr::group_by(meta_var) %>%
  dplyr::mutate(max.PC.cor = max(value),
                mean.PC.cor = mean(value)) %>%
  dplyr::ungroup()

#### A: P-values
cors.p.numericVars = pca %>%
  .[["x"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("bxp_id_full") %>%
  tibble::tibble() %>%
  .[, seq(highest_pc+1)] %>%
  dplyr::left_join(x=.,
                   y=meta %>% dplyr::select(bxp_id_full, where(is.numeric)),
                   by="bxp_id_full") %>%
  dplyr::select(-any_of(c("bxp_id_full", "sample_id"))) %>%
  dplyr::mutate(across(c(where(is.character)), as.factor),
                across(c(where(is.factor)), as.numeric)) %>%
  as.matrix() %>%
  Hmisc::rcorr(., type = "spearman") %>%
  .$P %>%
  as.data.frame() %>%
  dplyr::select(matches("^PC\\d+")) %>%
  # J: need to get rid of some of the meta variables from picard that being with PCT
  dplyr::select(-contains("PCT")) %>%
  tibble::rownames_to_column("meta_var") %>%
  dplyr::filter(!grepl("^PC\\d+", meta_var)) %>%
  tidyr::pivot_longer(2:ncol(.)) %>%
  dplyr::mutate(name = as.factor(as.numeric(gsub("PC", "", name))))

#### A: bind p and stat
cors.num = dplyr::left_join(
  x=cors.p.numericVars %>% dplyr::rename(p=value),
  y=cors.r.numericVars %>% dplyr::rename(stat=value),
  by=c("meta_var", "name")
) %>%
  dplyr::mutate(var_type = "continuous") %>%
  dplyr::group_by(meta_var) %>%
  dplyr::mutate(min.pvalue = min(p)) %>%
  dplyr::relocate(var_type, .after = last_col())

#### A: get PC/categorical correlations using K-W Chi-squared
#### G: added some regular expression to avoid issues with columns that are not
#### PCs.
cors.cat = pca %>%
  .[["x"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("bxp_id_full") %>%
  tibble::tibble() %>%
  .[, seq(highest_pc+1)] %>%
  dplyr::left_join(x=.,
                   y=meta %>% dplyr::select(bxp_id_full, !where(is.numeric)),
                   by="bxp_id_full") %>%
  dplyr::select(-any_of(c("bxp_id_full", "sample_id", "group", "Group"))) %>%
  # A: run K-W test
  run_stat_test_for_all_cols_get_res(in_df=., stat_test="kruskal") %>%
  # A: filter so var1 contains the numerics, var2 contains the categorical
  dplyr::filter(grepl("^PC\\d+", var1),
                !grepl("^PC\\d+", var2)) %>%
  dplyr::select(-test_run) %>%
  dplyr::rename(name=var1, meta_var=var2) %>%
  dplyr::group_by(meta_var) %>%
  dplyr::mutate(max.PC.cor = max(stat),
                mean.PC.cor = mean(stat),
                min.pvalue = min(p)) %>%
  dplyr::mutate(name = as.factor(as.numeric(gsub("PC", "", name)))) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(meta_var, name) %>%
  dplyr::mutate(var_type="categorical")

#### A: generate dataframe containing covariate, PC axis, rho & P-value
#### corresponding to PC-covariate, variance explained by each PC first, binding
#### cont (spearman) and cat (K-S) data to get df containing stat assessment for
#### all putative covariates
cors.all = dplyr::bind_rows(
  cors.num,
  cors.cat
) %>%
  dplyr::mutate(name = as.factor(name)) %>%
  dplyr::mutate(
    padj = p.adjust(p, method="fdr"),
    is.sig = case_when(
      padj<0.05 ~ stat,
      TRUE ~ NA_integer_
    )) %>%
  # A: add variance explained
  dplyr::left_join(
    x=.,
    y=summary(pca) %>%
      .$importance %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("name") %>%
      janitor::clean_names() %>%
      dplyr::mutate(name = as.factor(gsub("PC", "", name))),
    by="name"
  ) %>%
  # A: calculate weighted correlation - weights PC-covariate correlation coefficient by the % variance explained by that PC
  dplyr::mutate(variance_weighted_correlation = proportion_of_variance*(stat**2)) %>%
  # A: add adjusted P-value - adjust for number of putative covariates
  dplyr::group_by(name) %>%
  dplyr::mutate(padj = p.adjust(p, "fdr")) %>%
  dplyr::ungroup()


#### A: generate PC-covariate heatmaps

#### A: plot PC-cov correlations
pca.cor.num = cors.all %>%
  dplyr::mutate(name = as.numeric(as.character(name))) %>%
  dplyr::filter(name %in% c(1:20)) %>%
  dplyr::filter(var_type=="continuous") %>%
  dplyr::arrange(name) %>%
  ggplot(data=., aes(x=name, y=reorder(meta_var, max.PC.cor), fill=stat)) +
  geom_tile() +
  geom_text(aes(label = round(is.sig, 2)), color = "black", size = 1.75) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000",
                       limits = c(-1, 1),
                       name = "Spearman's rho") +
  xlab("PC") +
  ylab("") +
  theme(axis.text.x = element_text(angle=0, vjust = 0.5, hjust=0.5),
        legend.position = "right") +
  scale_x_continuous(breaks=c(1:20)) +
  facet_wrap(~var_type, scales = "free")

pca.cor.cat = cors.all %>%
  dplyr::mutate(name = as.numeric(as.character(name))) %>%
  dplyr::filter(name %in% c(1:20)) %>%
  dplyr::filter(var_type=="categorical") %>%
  dplyr::arrange(name) %>%
  ggplot(data=., aes(x=name, y=reorder(meta_var, max.PC.cor), fill=stat)) +
  geom_tile() +
  geom_text(aes(label = ifelse(p<0.05, round(p, 2), "")), color = "black", size = 1.75) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000",
                       name = "K-W chi-squared") +
  xlab("PC") +
  ylab("") +
  theme(axis.text.x = element_text(angle=0, vjust = 0.5, hjust=0.5),
        legend.position = "right") +
  scale_x_continuous(breaks=c(1:20)) +
  facet_wrap(~var_type, scales = "free")

#### A: use patchwork to get cateorical and numeric/cont heatmaps in one figure
pca.cor.fig = (pca.cor.num / pca.cor.cat)

#### A: plot prop of variance explained
prop.variance.pca = summary(pca) %>%
  .$importance %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PC") %>%
  janitor::clean_names() %>%
  dplyr::mutate(pc = as.numeric(gsub("PC", "", pc))) %>%
  dplyr::filter(pc<=20) %>%
  {
    ggplot(data=.) +
      geom_col(aes(x=pc, y=proportion_of_variance)) +
      geom_line(aes(x=pc, y=cumulative_proportion), linetype=2) +
      xlab("Principle component #") +
      ylab("Proportion of variance") +
      scale_x_continuous(breaks = round(seq(min(.$pc), max(.$pc), by = 1),1)) +
      theme(axis.text.x = element_text(angle=0)) +
      ylim(0,1)
    # geom_text(aes(x=pc, y=proportion_of_variance+0.025, label=paste0(round(proportion_of_variance, 3)*100, "%")), size=2.5)
  }

final.obj = list(
  varPart.plot, # 1. variancePartition distribution plot
  pca.cor.fig, # 2. PCA-covariate correlation heatmap
  prop.variance.pca, # 3. Prop variance explained by each PC
  colinearDistplot, # 4. distribution of covariate colinearity
  cors.all # 5. PC-meta correlations with variance explained and variance_weighted_correlation
)

final.obj %>% saveRDS(varPartCorPlots_path)

# 4. Final covariates ----
dat <- cors.all %>%
  dplyr::select(meta_var, name, p, padj, rho = stat, variance_weighted_correlation, var_type) %>%
  dplyr::arrange(-variance_weighted_correlation) %>%
  dplyr::group_by(meta_var) %>%
  dplyr::filter(variance_weighted_correlation>0) %>%
  dplyr::summarise(variance_weighted_correlation = max(variance_weighted_correlation),
                   min_p = min(p)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    x=.,
    y=cors.all %>%
      dplyr::select(meta_var, name, p, padj, rho = stat, variance_weighted_correlation, var_type),
    by=c("meta_var", "variance_weighted_correlation")
  ) %>%
  dplyr::rename(pc=name) %>%
  dplyr::arrange(-variance_weighted_correlation) %>%
  # then calculate and join variancePartition metrics
  dplyr::left_join(
    x=.,
    y=varPart.plot$data %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(varPart.q3 = as.numeric(quantile(value)[4]),
                       varPart.max = as.numeric(quantile(value)[5])),
    by=c("meta_var"="variable")
  ) %>%
  dplyr::mutate(dataset = "Wood_Bulk") %>%
  dplyr::mutate(is.sig = dplyr::case_when(
    ((p<=0.05) & (padj>0.05)) ~ "P<0.05",
    p>0.05 ~ "NS",
    padj<=0.05 ~ "FDR<0.05",
    TRUE ~ ""))

selected_covs = dat %>%
  dplyr::group_by(meta_var) %>%
  dplyr::summarise(n.PC.sig = sum((padj<=0.05 & var_type=="categorical") |
                                    (padj<=0.05 & variance_weighted_correlation>=0.025 & var_type=="continuous")),
                   n.varPart.q3 = sum(varPart.q3 >= 3),
                   n.varPart.max = sum(varPart.max >= 75)) %>%
  dplyr::ungroup() %>%
  # A: passes (FDR<0.05 & variance-weighted covariate-PC cut-off AND passes
  # variancePartition q3 cut-off) OR (max cut-off at any point)
  dplyr::filter((n.PC.sig >= 1 & n.varPart.q3 >= 1) | (n.varPart.max >= 1)) %>%
  # dplyr::select(meta_var) %>%
  dplyr::distinct()

print(selected_covs)
