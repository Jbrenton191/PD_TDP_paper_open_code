# Checking novel sequences against those that are different to the annotated transcript longest open reading frame

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. Setup ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(readxl)
  library(Biostrings)
})


conflicted::conflict_prefer_all(winner = "dplyr", quiet = T)


# set paths
main_path<-here::here() # main project level i.e TDP43_paper_open_code level

data_path<-file.path(main_path, "Data/LRRK2_HS_PSC/mDN1")

results_out_path<-file.path(main_path, "Data/LRRK2_HS_PSC/mDN1")

# load peptide sequences from non-uniprot database matching ions
novel_sequences_donors <-
  readxl::read_xlsx(
    path = file.path(
      data_path,
      "DRio_NovelDonor_DIA-NN_Updated.xlsx"
    ),
    sheet = 2
  ) %>% arrange(Protein.Group)

novel_sequences_acceptor <-
  readxl::read_xlsx(
    path = file.path(
      data_path,
      "DRio_NovelAcceptor_DIA-NN_Updated.xlsx"
    ),
    sheet = 2
  ) %>% arrange(Protein.Group)


# load predicted sequences made from junctions
# donors
junction_derived_novel_donors <- read_csv(file = 
                                            file.path(
                                              data_path,
                                              "Upreg_tdp_tsl_nmd_included_Rio_LRRK2_novel_donors_junction_NT_CDS_AA_sequences_filtered_no_annotation_AA_matches.csv"
                                            )
) %>% mutate(Transcript=str_remove_all(pattern = "_ND$", string=Transcript))

# acceptors
junction_derived_novel_acceptors <- read_csv(file = 
                                               file.path(
                                                 data_path,
                                                 "Upreg_tdp_tsl_nmd_included_Rio_LRRK2_novel_acceptor_junction_NT_CDS_AA_sequences_filtered_no_annotation_AA_matches.csv"
                                               )
) %>% mutate(Transcript=str_remove_all(pattern = "_NA$", string=Transcript))



min_match_run<-2

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2. Search for novel peptide MS sequences in junction derived amino acid sequences that don't match relevant annotated transcripts ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 2.1 Novel donor checks ---- 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

merge_results_df<-c()

# Preprocess: separate Protein.Group rows and transcript field
novel_sequences_donors_sep <- novel_sequences_donors %>%
  separate_rows(Protein.Group, sep = ";") %>%
  mutate(transcript = str_remove_all(pattern = "\\_.*", string = Protein.Group)) %>%
  mutate(`Gene names`=str_remove_all(pattern = ";.*", string = `Gene names`))

# Then loop over unique transcripts
for (sel_transcript in unique(novel_sequences_donors_sep$Protein.Group)) {
  
  filt_novel_donor <- novel_sequences_donors_sep %>%
    filter(Protein.Group == sel_transcript)
  
  
  
  filt_junction_derived_novel_donors <-
    junction_derived_novel_donors[junction_derived_novel_donors$Transcript %in% filt_novel_donor$transcript,]
  
  if (nrow(filt_junction_derived_novel_donors)<1) {
    next
  }
  
  for (i in 1:length(filt_junction_derived_novel_donors$seq_aa_from_first_diff)) {
    
    if (is.na(filt_junction_derived_novel_donors$seq_aa_from_first_diff[i])) {
      next
    }
    
    # skip lines where longest transcripts are the same
    if (filt_junction_derived_novel_donors$seq_aa_from_first_diff[i]=="No difference in longest transcript between normal and novel ones") {
      next
    }
    
    # skip lines where there are only < 2 amino acid differences - potential early stop.
    if (str_length(filt_junction_derived_novel_donors$seq_aa_from_first_diff[i])<min_match_run) {
      next
    }
    
    results_df<-c()
    
    for (j in 1:length(filt_novel_donor$transcript)) {
      
      # compare sequences through blasting
      seq1 <- AAString(filt_novel_donor$Stripped.Sequence[j])
      seq2 <- AAString(filt_junction_derived_novel_donors$AA_seq[i])
      #blast
      alignment <- pairwiseAlignment(pattern = seq1, subject = seq2,
                                     substitutionMatrix = "BLOSUM62",
                                     gapOpening = -10, gapExtension = -0.5,
                                     type = "local")
      
      # results and naming for output:
      transcript_status<-"Novel_donor"
      gene<-gsub(x = filt_novel_donor$`Gene names`[j], pattern = "^(.*);.*", replacement = "\\1" )
      transcript<-filt_novel_donor$transcript[j] # transcript name
      query_peptide_ion<-as.character(seq1) # ms ion seq
      AA_seq<-as.character(seq2) # predicted seq
      precursor_charge<-filt_novel_donor$Precursor.Charge[j]
      modified_sequence<-filt_novel_donor$Modified.Sequence[j] 
      alignment_score<-alignment@score # blast score of novel seq vs ms ion
      full_match_within<-alignedPattern(alignment) == alignedSubject(alignment) # is the MS ion fully contained within the novel transcript amion acid sequence
      novel_junction_dpsi<-filt_junction_derived_novel_donors$novel_junction_dpsi # splicing delta psi of the junction
      NMD<-filt_junction_derived_novel_donors$NMD # NMD prediction
      TSL<-filt_junction_derived_novel_donors$TSL # TSL evidence of the annotated transcript connected to the predicted sequence
      
      
      ms_ion_seq_length<-length(seq1)
      number_matches<-nmatch(alignment)     # Number of matches from blast
      percent_matching<-(number_matches/ms_ion_seq_length)*100 # percent of matching from blast
      num_mismatch<-nmismatch(alignment)  # Number of mismatches
      
      indel<-nindel(alignment)     # Number of inserted/deleted residues (gaps)
      
      inserted_residues_number<-indel@insertion[1]
      deleted_residues_number<-indel@deletion[1]
      
      aligned_query<-as.character(alignedPattern(alignment))  # Part of seq1 that was aligned
      
      
      # annotated sequence alignment++++++++++++++++++++++++++++++++++++++++
      
      
      if (is.na(filt_junction_derived_novel_donors$old_transcript_AA_seq[i])) {
        next
      }
      
      
      seq3 <- AAString(filt_junction_derived_novel_donors$old_transcript_AA_seq[i])
      
      # blasting annotated transcript sequence (wthout the novel junction change) to the MS ion sequence
      alignment2 <- pairwiseAlignment(pattern = seq1, subject = seq3,
                                      substitutionMatrix = "BLOSUM62",
                                      gapOpening = -10, gapExtension = -0.5,
                                      type = "local")
      
      annotated_AA_seq<-as.character(seq3)
      annotated_alignment_score<-alignment2@score
      annotated_full_match_within<-alignedPattern(alignment2) == alignedSubject(alignment2)
      
      annotated_number_matches<-nmatch(alignment2)     # Number of matches
      annotated_percent_matching<-(annotated_number_matches/ms_ion_seq_length)*100
      annotated_num_mismatch<-nmismatch(alignment2)  # Number of mismatches
      
      annotated_indel<-nindel(alignment2)     # Number of inserted/deleted residues (gaps)
      
      annotated_inserted_residues_number<-annotated_indel@insertion[1]
      annotated_deleted_residues_number<-annotated_indel@deletion[1]
      
      annotated_aligned_query<-as.character(alignedPattern(alignment2)) 
      
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #+ Checks between unique sequence and annotated
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      # comparing the blasting of the MS ion sequence to that of the novel transcript and annotated transcript to make sure the MS ion sequence would be found in the novel area
      
      extract_clean_unique_nonmatching <- function(unique_seq, annotated_seq, min_match_run = 5, min_start_pos = 1) {
        aln <- pairwiseAlignment(unique_seq, annotated_seq,
                                 substitutionMatrix = "BLOSUM62",
                                 gapOpening = 5, gapExtension = 2)
        
        # Aligned strings with gaps
        aligned_unique <- strsplit(as.character(pattern(aln)), "")[[1]]
        aligned_annotated <- strsplit(as.character(subject(aln)), "")[[1]]
        
        matches <- aligned_unique == aligned_annotated
        rle_matches <- rle(matches)
        run_starts <- cumsum(c(1, head(rle_matches$lengths, -1)))
        
        converged_later <- any(rle_matches$values == TRUE & cumsum(rle_matches$values == FALSE) > 0)
        converged_later_res <- if (converged_later) {
          "The unique amino acid sequence diverged from the annotated one and then re-aligned."
        } else {
          "The unique amino acid sequence did not re-align with the annotated one after sequences diverged"
        }
        
        # Find re-alignment start
        realignment_start <- NA
        for (i in seq_along(rle_matches$values)) {
          if (rle_matches$values[i] &&
              rle_matches$lengths[i] >= min_match_run &&
              run_starts[i] >= min_start_pos) {
            prior_mismatch <- any(!rle_matches$values[1:(i - 1)])
            if (prior_mismatch) {
              realignment_start <- run_starts[i]
              break
            }
          }
        }
        
        # Identify the mismatch region before re-alignment
        mismatch_indices <- which(!matches)
        if (!is.na(realignment_start)) {
          mismatch_indices <- mismatch_indices[mismatch_indices < realignment_start]
        }
        
        # Map to original unique_seq indices, skip gaps
        u_seq <- strsplit(as.character(unique_seq), "")[[1]]
        u_index <- 1
        clean_indices <- c()
        
        for (i in seq_along(aligned_unique)) {
          if (aligned_unique[i] != "-") {
            if (i %in% mismatch_indices) {
              clean_indices <- c(clean_indices, u_index)
            }
            u_index <- u_index + 1
          }
        }
        
        clean_mismatched <- u_seq[clean_indices]
        clean_mismatched <- clean_mismatched[!clean_mismatched %in% c("X", "*")]
        clean_sequence <- paste(clean_mismatched, collapse = "")
        
        return(list(
          convergence = converged_later_res,
          unique_sequence = clean_sequence
        ))
      }
      
      diff_unique_annot_res<-extract_clean_unique_nonmatching(unique_seq = seq2, annotated_seq = seq3)
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #+ Comparing unique non-matching amino acid sequence from junction with MS peptide sequence
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      seq4 <- AAString(diff_unique_annot_res$unique_sequence)
      
      alignment3 <- pairwiseAlignment(pattern = seq1, subject = seq4,
                                      substitutionMatrix = "BLOSUM62",
                                      gapOpening = -10, gapExtension = -0.5,
                                      type = "local")
      
      unique_sequence_nonmatching_annot_AA_seq<-as.character(seq4)
      unique_sequence_nonmatching_annot_alignment_score<-alignment3@score
      unique_sequence_nonmatching_annot_full_match_within<-alignedPattern(alignment3) == alignedSubject(alignment3)
      
      unique_sequence_nonmatching_annot_number_matches<-nmatch(alignment3)     # Number of matches
      unique_sequence_nonmatching_annot_percent_matching<-(unique_sequence_nonmatching_annot_number_matches/ms_ion_seq_length)*100
      unique_sequence_nonmatching_annot_num_mismatch<-nmismatch(alignment3)  # Number of mismatches
      
      unique_sequence_nonmatching_annot_indel<-nindel(alignment3)     # Number of inserted/deleted residues (gaps)
      
      
      unique_sequence_nonmatching_annot_inserted_residues_number<-unique_sequence_nonmatching_annot_indel@insertion[1]
      unique_sequence_nonmatching_annot_deleted_residues_number<-unique_sequence_nonmatching_annot_indel@deletion[1]
      
      unique_sequence_nonmatching_annot_aligned_query<-as.character(alignedPattern(alignment3)) 
      
      # create alignment scores for both the novel alignment and annotated alignment and delta alignment scores:
      # create normalised alignment score
      alignment_length = nchar(aligned_query)
      # create normalised alignment score for 
      annotated_alignment_length = nchar(annotated_aligned_query)
      
      normalized_score_by_alignment_length = alignment_score / alignment_length
      normalized_annotated_score_by_alignment_length = annotated_alignment_score / annotated_alignment_length
      
      normalized_score_by_full_query_length = alignment_score / ms_ion_seq_length
      normalized_annotated_score_by_full_query_length = annotated_alignment_score / ms_ion_seq_length
      
      delta_alignment_score_by_alignment_length = normalized_score_by_alignment_length - normalized_annotated_score_by_alignment_length
      
      delta_alignment_score_by_full_query_length = normalized_score_by_full_query_length - normalized_annotated_score_by_full_query_length
      
      ratio_alignment_score_by_alignment_length = normalized_score_by_alignment_length / normalized_annotated_score_by_alignment_length
      
      ratio_alignment_score_by_full_query_length = normalized_score_by_full_query_length / normalized_annotated_score_by_full_query_length
      
      
      # General Guidance for Thresholding score_diff
      # Rule of Thumb (Normalized Score Differences):
      #   Score Difference (score_diff)	Interpretation
      # < 0.1	Likely ambiguous — aligns similarly to both.
      # 0.1–0.3	Suggestive, but might need manual review or tie-breaking info.
      # > 0.3	Reasonable confidence it aligns better to the higher-scoring one.
      # > 0.5	High confidence match to the better alignment.
      
      # You normalize by the aligned query length (not the original query) because:
      # The alignment score only reflects the portion of the query that actually aligned.
      # Unaligned regions don’t contribute to the score — so including them in the denominator would artificially deflate the normalized value.
      
      # Normalization Method	Use When...
      # score / aligned_query_length	You want to compare alignment quality per aligned residue (e.g. distinguishing weak vs. strong matches).
      # score / full_query_length	You want to penalize partial matches — e.g., you're looking for full-length homology or want to prefer alignments covering the whole query.
      
      
      # Have included both as we want more full matches across all the query not just very local alignments, so want to check the full query length more to see if it makes sense.
      
      
      # Compare percent matching:
      
      delta_percent_matching=percent_matching-annotated_percent_matching
      
      
      # merge results
      
      result<-tibble(
        transcript_status,
        gene,
        transcript,
        query_peptide_ion,
        AA_seq,
        precursor_charge,
        modified_sequence,
        alignment_score,
        full_match_within,
        novel_junction_dpsi,
        NMD,
        TSL,
        ms_ion_seq_length,
        number_matches,
        percent_matching,
        num_mismatch,
        inserted_residues_number,
        deleted_residues_number,
        aligned_query,
        
        annotated_AA_seq,
        annotated_alignment_score,
        annotated_full_match_within,
        
        annotated_number_matches,     # Number of matches
        annotated_percent_matching,
        annotated_num_mismatch,
        annotated_inserted_residues_number,
        annotated_deleted_residues_number,
        
        annotated_aligned_query,
        convergence=diff_unique_annot_res$convergence,
        unique_sequence_nonmatching_annot_AA_seq,
        unique_sequence_nonmatching_annot_alignment_score,
        unique_sequence_nonmatching_annot_full_match_within,
        
        unique_sequence_nonmatching_annot_number_matches,  # Number of matches
        unique_sequence_nonmatching_annot_percent_matching,
        unique_sequence_nonmatching_annot_num_mismatch,
        
        unique_sequence_nonmatching_annot_inserted_residues_number,
        unique_sequence_nonmatching_annot_deleted_residues_number,
        
        unique_sequence_nonmatching_annot_aligned_query,
        
        
        alignment_length,
        annotated_alignment_length,
        normalized_score_by_alignment_length,
        normalized_annotated_score_by_alignment_length,
        normalized_score_by_full_query_length,
        normalized_annotated_score_by_full_query_length,
        delta_alignment_score_by_alignment_length,
        delta_alignment_score_by_full_query_length,
        ratio_alignment_score_by_alignment_length,
        ratio_alignment_score_by_full_query_length,
        delta_percent_matching
        
      )
      
      results_df<-rbind(results_df, result)
      
    } # end of matching loop
    
    merge_results_df<-rbind(merge_results_df, results_df)
    
  }
  
}

novel_donor_matching<-merge_results_df

write_csv(novel_donor_matching, file = file.path(results_out_path, "Upreg_tdp_tsl_nmd_included_Novel_donor_only_All_MS_seq_matching_to_novel_junction_aa_differences_to_annotated_19.5.25.csv"))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 2.2 Novel acceptor checks -------- 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

merge_results_df<-c()

# Preprocess: separate Protein.Group rows and clean transcript field
novel_sequences_acceptors_sep <- novel_sequences_acceptor %>%
  separate_rows(Protein.Group, sep = ";") %>%
  mutate(transcript = str_remove_all(pattern = "\\_.*", string = Protein.Group)) %>%
  mutate(`Gene names`=str_remove_all(pattern = ";.*", string = `Gene names`))

# Then loop over unique transcripts (already clean)
for (sel_transcript in unique(novel_sequences_acceptors_sep$Protein.Group)) {
  
  filt_novel_acceptor <- novel_sequences_acceptors_sep %>%
    filter(Protein.Group == sel_transcript)
  
  
  # for (sel_transcript in unique(novel_sequences_acceptor$Protein.Group)) {
  #   
  #   filt_novel_acceptor <-
  #     novel_sequences_acceptor %>% filter(Protein.Group %in% sel_transcript) %>%
  #     mutate(transcript = str_remove_all(pattern = "\\_.*", string = Protein.Group))
  
  filt_junction_derived_novel_acceptor <-
    junction_derived_novel_acceptors[junction_derived_novel_acceptors$Transcript %in% filt_novel_acceptor$transcript,]
  
  if (nrow(filt_junction_derived_novel_acceptor)<1) {
    next
  }
  
  
  for (i in 1:length(filt_junction_derived_novel_acceptor$seq_aa_from_first_diff)) {
    
    if (is.na(filt_junction_derived_novel_donors$seq_aa_from_first_diff[i])) {
      next
    }
    
    if (filt_junction_derived_novel_acceptor$seq_aa_from_first_diff[i]=="No difference in longest transcript between normal and novel ones") {
      next
    }
    
    # skip lines where there are only minimum match amount (chosen above) amino acid differences - potential early stop.
    if (str_length(filt_junction_derived_novel_acceptor$seq_aa_from_first_diff[i])<min_match_run) {
      next
    }
    
    results_df<-c()
    
    for (j in 1:length(filt_novel_acceptor$transcript)) {
      
      seq1 <- AAString(filt_novel_acceptor$Stripped.Sequence[j])
      seq2 <- AAString(filt_junction_derived_novel_acceptor$AA_seq[i])
      
      alignment <- pairwiseAlignment(pattern = seq1, subject = seq2,
                                     substitutionMatrix = "BLOSUM62",
                                     gapOpening = -10, gapExtension = -0.5,
                                     type = "local")
      # results:
      transcript_status<-"Novel_acceptor"
      # gene<-filt_novel_acceptor$`Gene names`[j]
      gene<-gsub(x = filt_novel_acceptor$`Gene names`[j], pattern = "^(.*);.*", replacement = "\\1" )
      
      transcript<-filt_novel_acceptor$transcript[j]
      query_peptide_ion<-as.character(seq1)
      AA_seq<-as.character(seq2)
      precursor_charge<-filt_novel_acceptor$Precursor.Charge[j]
      modified_sequence<-filt_novel_acceptor$Modified.Sequence[j] 
      alignment_score<-alignment@score
      full_match_within<-alignedPattern(alignment) == alignedSubject(alignment)
      novel_junction_dpsi<-filt_junction_derived_novel_acceptor$novel_junction_dpsi
      NMD<-filt_junction_derived_novel_acceptor$NMD
      TSL<-filt_junction_derived_novel_acceptor$TSL
      
      ms_ion_seq_length<-length(seq1)
      number_matches<-nmatch(alignment)     # Number of matches
      percent_matching<-(number_matches/ms_ion_seq_length)*100
      num_mismatch<-nmismatch(alignment)  # Number of mismatches
      
      indel<-nindel(alignment)     # Number of inserted/deleted residues (gaps)
      
      inserted_residues_number<-indel@insertion[1]
      deleted_residues_number<-indel@deletion[1]
      
      aligned_query<-as.character(alignedPattern(alignment))  # Part of seq1 that was aligned
      
      
      
      # annotated sequence alignment++++++++++++++++++++++++++++++++++++++++
      
      
      if (is.na(filt_junction_derived_novel_acceptor$old_transcript_AA_seq[i])) {
        next
      }
      
      
      seq3 <- AAString(filt_junction_derived_novel_acceptor$old_transcript_AA_seq[i])
      
      alignment2 <- pairwiseAlignment(pattern = seq1, subject = seq3,
                                      substitutionMatrix = "BLOSUM62",
                                      gapOpening = -10, gapExtension = -0.5,
                                      type = "local")
      
      annotated_AA_seq<-as.character(seq3)
      annotated_alignment_score<-alignment2@score
      annotated_full_match_within<-alignedPattern(alignment2) == alignedSubject(alignment2)
      
      annotated_number_matches<-nmatch(alignment2)     # Number of matches
      annotated_percent_matching<-(annotated_number_matches/ms_ion_seq_length)*100
      annotated_num_mismatch<-nmismatch(alignment2)  # Number of mismatches
      
      annotated_indel<-nindel(alignment2)     # Number of inserted/deleted residues (gaps)
      
      annotated_inserted_residues_number<-annotated_indel@insertion[1]
      annotated_deleted_residues_number<-annotated_indel@deletion[1]
      
      annotated_aligned_query<-as.character(alignedPattern(alignment2)) 
      
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #+ Checks between unique sequence and annotated
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      diff_unique_annot_res<-extract_clean_unique_nonmatching(unique_seq = seq2, annotated_seq = seq3)
      
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #+ Comparing unique non-matching amino acid sequence from junction with MS peptide sequence
      #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      seq4 <- AAString(diff_unique_annot_res$unique_sequence)
      
      alignment3 <- pairwiseAlignment(pattern = seq1, subject = seq4,
                                      substitutionMatrix = "BLOSUM62",
                                      gapOpening = -10, gapExtension = -0.5,
                                      type = "local")
      
      unique_sequence_nonmatching_annot_AA_seq<-as.character(seq4)
      unique_sequence_nonmatching_annot_alignment_score<-alignment3@score
      unique_sequence_nonmatching_annot_full_match_within<-alignedPattern(alignment3) == alignedSubject(alignment3)
      
      unique_sequence_nonmatching_annot_number_matches<-nmatch(alignment3)     # Number of matches
      unique_sequence_nonmatching_annot_percent_matching<-(unique_sequence_nonmatching_annot_number_matches/ms_ion_seq_length)*100
      unique_sequence_nonmatching_annot_num_mismatch<-nmismatch(alignment3)  # Number of mismatches
      
      unique_sequence_nonmatching_annot_indel<-nindel(alignment3)     # Number of inserted/deleted residues (gaps)
      
      
      unique_sequence_nonmatching_annot_inserted_residues_number<-unique_sequence_nonmatching_annot_indel@insertion[1]
      unique_sequence_nonmatching_annot_deleted_residues_number<-unique_sequence_nonmatching_annot_indel@deletion[1]
      
      unique_sequence_nonmatching_annot_aligned_query<-as.character(alignedPattern(alignment3)) 
      
      
      
      # create alignment scores for both the novel alignment and annotated alignment and delta alignment scores:
      # create normalised alignment score
      alignment_length = nchar(aligned_query)
      # create normalised alignment score for 
      annotated_alignment_length = nchar(annotated_aligned_query)
      
      normalized_score_by_alignment_length = alignment_score / alignment_length
      normalized_annotated_score_by_alignment_length = annotated_alignment_score / annotated_alignment_length
      
      normalized_score_by_full_query_length = alignment_score / ms_ion_seq_length
      normalized_annotated_score_by_full_query_length = annotated_alignment_score / ms_ion_seq_length
      
      delta_alignment_score_by_alignment_length = normalized_score_by_alignment_length - normalized_annotated_score_by_alignment_length
      
      delta_alignment_score_by_full_query_length = normalized_score_by_full_query_length - normalized_annotated_score_by_full_query_length
      
      ratio_alignment_score_by_alignment_length = normalized_score_by_alignment_length / normalized_annotated_score_by_alignment_length
      
      ratio_alignment_score_by_full_query_length = normalized_score_by_full_query_length / normalized_annotated_score_by_full_query_length
      
      
      # General Guidance for Thresholding score_diff
      # Rule of Thumb (Normalized Score Differences):
      #   Score Difference (score_diff)	Interpretation
      # < 0.1	Likely ambiguous — aligns similarly to both.
      # 0.1–0.3	Suggestive, but might need manual review or tie-breaking info.
      # > 0.3	Reasonable confidence it aligns better to the higher-scoring one.
      # > 0.5	High confidence match to the better alignment.
      
      # You normalize by the aligned query length (not the original query) because:
      # The alignment score only reflects the portion of the query that actually aligned.
      # Unaligned regions don’t contribute to the score — so including them in the denominator would artificially deflate the normalized value.
      
      # Normalization Method	Use When...
      # score / aligned_query_length	You want to compare alignment quality per aligned residue (e.g. distinguishing weak vs. strong matches).
      # score / full_query_length	You want to penalize partial matches — e.g., you're looking for full-length homology or want to prefer alignments covering the whole query.
      
      
      # Have included both as we want more full matches across all the query not just very local alignments, so want to check the full query length more to see if it makes sense.
      
      
      # Compare percent matching:
      
      delta_percent_matching=percent_matching-annotated_percent_matching
      
      
      result<-tibble(
        transcript_status,
        gene,
        transcript,
        query_peptide_ion,
        AA_seq,
        precursor_charge,
        modified_sequence,
        alignment_score,
        full_match_within,
        novel_junction_dpsi,
        NMD,
        TSL,
        ms_ion_seq_length,
        number_matches,
        percent_matching,
        num_mismatch,
        inserted_residues_number,
        deleted_residues_number,
        aligned_query,
        
        annotated_AA_seq,
        annotated_alignment_score,
        annotated_full_match_within,
        
        annotated_number_matches,     # Number of matches
        annotated_percent_matching,
        annotated_num_mismatch,
        annotated_inserted_residues_number,
        annotated_deleted_residues_number,
        
        annotated_aligned_query,
        convergence=diff_unique_annot_res$convergence,
        unique_sequence_nonmatching_annot_AA_seq,
        unique_sequence_nonmatching_annot_alignment_score,
        unique_sequence_nonmatching_annot_full_match_within,
        
        unique_sequence_nonmatching_annot_number_matches,  # Number of matches
        unique_sequence_nonmatching_annot_percent_matching,
        unique_sequence_nonmatching_annot_num_mismatch,
        
        unique_sequence_nonmatching_annot_inserted_residues_number,
        unique_sequence_nonmatching_annot_deleted_residues_number,
        
        unique_sequence_nonmatching_annot_aligned_query,
        
        
        alignment_length,
        annotated_alignment_length,
        normalized_score_by_alignment_length,
        normalized_annotated_score_by_alignment_length,
        normalized_score_by_full_query_length,
        normalized_annotated_score_by_full_query_length,
        delta_alignment_score_by_alignment_length,
        delta_alignment_score_by_full_query_length,
        ratio_alignment_score_by_alignment_length,
        ratio_alignment_score_by_full_query_length,
        
        delta_percent_matching
        
      )
      
      results_df<-rbind(results_df, result)
      
    } # end of matching loop
    
    
    merge_results_df<-rbind(merge_results_df, results_df)
    
  }
  
}

novel_acceptor_matching<-merge_results_df

write_csv(novel_acceptor_matching, file = file.path(results_out_path, "Upreg_tdp_tsl_nmd_included_Novel_acceptor_only_All_MS_seq_matching_to_novel_junction_aa_differences_to_annotated_19.5.25.csv"))


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3. Save results ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

novel_donor_and_acceptor_matching_results<-rbind(novel_donor_matching,
                                                 novel_acceptor_matching) %>%
  select(gene, transcript, transcript_status, everything()) %>% arrange(desc(alignment_score), desc(percent_matching))

write_csv(novel_donor_and_acceptor_matching_results, file = file.path(results_out_path, "Upreg_tdp_tsl_nmd_included_ALL_MS_seq_matching_to_novel_junction_aa_differences_to_annotated_7.5.25.csv"))



