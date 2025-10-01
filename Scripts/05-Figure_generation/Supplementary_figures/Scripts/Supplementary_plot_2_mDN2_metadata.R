library(kableExtra)
library(tidyverse)

metadata<-vroom::vroom(file = "/home/jbrenton/TDP43_PD_Paper_Figures/Drafting_figures_code/Supplementary_figures/Data/LRRK2_mDN_G2019S_dataset2_metadata.tsv")

colnames(metadata)<-str_to_title(colnames(metadata))

kable(metadata %>% mutate(Genotype=str_to_title(Genotype)) %>%
        select(c(-1,-2)) %>% select(Line, Induction, everything()), format = "latex")

# Move into overleaf and edit
