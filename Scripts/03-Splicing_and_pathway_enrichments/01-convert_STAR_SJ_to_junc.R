library(optparse)
library(tidyverse)

arguments <- parse_args(OptionParser(), positional_arguments = 3)

sj_dir_path <- arguments$args[1] %>% str_split(",") %>% unlist()
output_path <- arguments$args[2]
blacklist_path <- arguments$args[3]

# sj_dir_path <- "/home/jbrenton/nextflow_pd/output/STAR/align"
# base_dir<- "/home/jbrenton/nextflow_pd"
# output_path <- "/home/jbrenton/nextflow_pd/output/leafcutter"

convert_STAR_SJ_to_junc <- function(sj_dir_path, output_path, filter_out_blacklist_regions=FALSE, path_to_ENCODE_blacklist = NULL){
  
  library(tidyverse)
  
  paths <- list.files(path = sj_dir_path, pattern = "_SJ.out.tab", full.names = TRUE)
  
  if(filter_out_blacklist_regions == TRUE){
    
    library(GenomicRanges)
    library(rtracklayer)
    
    # Load encode blacklist (https://github.com/Boyle-Lab/Blacklist/tree/master/lists)
    ENCODE_blacklist_hg38 <- rtracklayer::import(path_to_ENCODE_blacklist)
    
  }
  
  for(i in 1:length(paths)){
    
    sample_name <- paths[i] %>%
      str_replace("/.*/", "") %>%
      str_replace("_SJ.out.tab", "")
    
    cat("Loading splice junctions from:", paths[i],"\n")
    
    sj_out <- read_delim(file = paths[i],
                         delim = "\t",
                         col_names = c("chr", "intron_start", "intron_end", "strand", "intron_motif", "in_annotation", "unique_reads_junction", "multi_map_reads_junction", "max_splice_alignment_overhang"),
                         col_types = cols(chr = "c", .default = "d"))
    
    if(filter_out_blacklist_regions == TRUE){
      
      sj_out <- sj_out %>%
        dplyr::mutate(strand = ifelse(strand == 0, "*",
                                      ifelse(strand == 1, "+", "-")),
                      unknown = ".") %>% # .junc files include a sixth empty column. '.' denotes this empty column
        makeGRangesFromDataFrame(.,
                                 keep.extra.columns = TRUE,
                                 seqnames.field = "chr",
                                 start.field = "intron_start",
                                 end.field = "intron_end",
                                 ignore.strand = FALSE)
      
      # Remove junctions that overlap with ENCODE blacklist regions
      overlapped_junctions <- GenomicRanges::findOverlaps(query = ENCODE_blacklist_hg38,
                                                          subject = sj_out,
                                                          ignore.strand = F)
      
      indexes <- subjectHits(overlapped_junctions)
      sj_out <- sj_out[-indexes, ] %>%
        as.data.frame() %>%
        dplyr::rename(chr = seqnames,
                      intron_start = start,
                      intron_end = end)
      
      # Convert to leafcutter format
      sj_out_leafcutter_format <-
        sj_out %>%
        # Remove unwanted chromosomes builds
        dplyr::filter(chr %in% c(str_c("chr", 1:22), "chrX", "chrY", "chrM", "chrMT")) %>%
        dplyr::mutate(strand = str_replace(strand, "\\*", "."),
                      unknown = ".") %>% # .junc files include a sixth empty column. '.' denotes this empty column
        dplyr::select(chr, intron_start, intron_end, unknown, unique_reads_junction, strand)
      
      
    } else{
      
      sj_out_leafcutter_format <-
        sj_out %>%
        # Remove unwanted chromosomes builds
        dplyr::filter(chr %in% c(str_c("chr", 1:22), "chrX", "chrY", "chrM", "chrMT")) %>%
        dplyr::mutate(strand = ifelse(strand == 0, ".",
                                      ifelse(strand == 1, "+", "-")),
                      unknown = ".") %>% # .junc files include a sixth empty column. '.' denotes this empty column
        dplyr::select(chr, intron_start, intron_end, unknown, unique_reads_junction, strand)
      
    }
    
    # change the stop and start into characters to avoid saving with the scientific notation
    write_delim(sj_out_leafcutter_format %>%
                  dplyr::mutate(intron_start = as.integer(intron_start),
                                intron_end = as.integer(intron_end),
                                unique_reads_junction = as.integer(unique_reads_junction)),
                path = str_c(output_path, "/", sample_name, "_SJ_leafcutter.junc"),
                delim = "\t",
                col_names = F)
  }
  
  # write a .txt file with each
  #junc_df <- tibble(junc_file_name = list.files(path = output_path, pattern = "_SJ_leafcutter.junc", full.names = TRUE))
  
  #write_delim(junc_df,
  #            path = str_c(output_path, "/list_juncfiles.txt"),
   #           delim = "\t",
    #          col_names = F)
  
}

convert_STAR_SJ_to_junc(sj_dir_path, output_path, filter_out_blacklist_regions=T, path_to_ENCODE_blacklist = blacklist_path)

junc_df <- tibble(junc_file_name = list.files(path = output_path,
                                              pattern = "_SJ_leafcutter.junc", full.names = TRUE))
 
write_delim(junc_df,
             path = str_c(".", "/list_juncfiles.txt"),
             delim = "\t", col_names = F)
