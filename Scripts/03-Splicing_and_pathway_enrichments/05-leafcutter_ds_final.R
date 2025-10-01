library(optparse)
library(tidyverse)
library(leafcutter)

arguments <- parse_args(OptionParser(), positional_arguments = 5)

count_file <- arguments$args[1]
min_samples_per_intron  <- arguments$args[2]
min_samples_per_group<-arguments$args[3]
min_coverage<-arguments$args[4]
working_dir<-arguments$args[5]

# count_file<-"~/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/02-Run_with_updated_covariates_and_outliers/Reclustering_of_introns_outliers_removed/16_outliers/leafcutter_perind_numers.counts.gz"

exon_file<-"/home/MRJonathanBrenton/Hardy_ASAP_bulk/Analysis/nextflow_pd/output_2pass_indv/leafcutter/gencode_LC_exon_file.txt.gz"
base_dir<-"/home/MRJonathanBrenton/Hardy_ASAP_bulk/Analysis/nextflow_pd"

# List brain areas super folders containing the group files - as have multiple group files folders to loop over
gf_dir_list <-
  list.dirs(
    '~/Hardy_ASAP_bulk/Analysis/Hardy_Leafcutter_R_analysis/07-Rerun_of_final_PD_groups_collapsed/Creating_groupfiles',
    recursive = F
  )

# Regina original setting
arguments$opt$output_prefix <- ""
arguments$opt$max_cluster_size <- "40"

# removed these as will be looped thought in command line.
# arguments$opt$min_samples_per_intron  <- "5"
# arguments$opt$min_samples_per_group  <- "3"
# arguments$opt$min_coverage <- "20"


arguments$opt$timeout <- "30"

arguments$opt$num_threads <- "64"
opt <- arguments$opt

# print(arguments$args)
# cat("this is gp dir:", groups_file_dir, sep = " ")


for (i in 1:length(gf_dir_list)) {
  print(paste0("this is the working dir:", working_dir))
  setwd(working_dir)
  groups_file_dir <- gf_dir_list[i]
  BA<-basename(groups_file_dir)
  dir.create(BA)
  # As set working dir to output files should work
  setwd(
    paste0(working_dir,"/", BA)
  )
  
  
  group_file_path<-list.files(path = groups_file_dir, pattern = "_group_file.txt", full.names = T)
  group_file_name<-list.files(path = groups_file_dir, pattern = "_group_file.txt", full.names = T) %>%
    str_replace("/.*/", "") %>%
    str_replace("_group_file.txt", "")
  # 
  # 
  group_file_df <- data_frame(group_file_path = group_file_path,
                              group_file_name = group_file_name)
  
  for(i in 1:nrow(group_file_df)){
    
    print(str_c("Performing leafcutter differential splicing using the group file entitled: ", group_file_df$group_file_name[i]))
    
    leafcutter_cmd <- str_c("Rscript ", base_dir, "/leafcutter_ds_threads_Hardy_version24.11.23.R ",
                            count_file, " ", # path to count file
                            group_file_df$group_file_path[i], # path to group file
                            " --output_prefix ", group_file_df$group_file_name[i], # comparison-specific output prefix
                            " --max_cluster_size ", opt$max_cluster_size,
                            " --min_samples_per_intron ", min_samples_per_intron,
                            " --min_samples_per_group ", min_samples_per_group,
                            " --min_coverage ", min_coverage,
                            " --timeout ", opt$timeout,
                            " --exon_file ", exon_file,
                            " --num_threads ", opt$num_threads
    )
    
    system(command = leafcutter_cmd)
    
  }
  
}
  
