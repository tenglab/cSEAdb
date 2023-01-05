#------------------------------------------------------------------------
#functions of track plot
#------------------------------------------------------------------------

# function of create track dataframe for plot based on se_region
create_se_track_data<- function (se_region,se_file,ce_matrix_file,bin_size=100) {
  # create se bin grange
  se_tmp <- se_file[which(se_file$V4==se_region),]
  se_gr <- GRanges(seqnames=se_tmp$V1,
                   ranges=IRanges(start = as.integer(se_tmp$V2), end = as.integer(se_tmp$V3)),
                   ce_name=se_tmp$V4)
  se_gr_bins <- unlist(tile(se_gr, width=bin_size))

  # create ce grange
  ce_df <- data.frame(ce_name=unique(ce_matrix_file$merge_e_name[which(ce_matrix_file$se_name==se_region)]))
  ce_df <- separate(ce_df,ce_name,c("chr","start","end"),remove = F)

  ce_gr <- GRanges(seqnames=ce_df$chr,
                   ranges=IRanges(start = as.integer(ce_df$start), end = as.integer(ce_df$end)),
                   ce_name=ce_df$ce_name)

  # create se bin count for 60 cell line
  out_df <- data.frame()
  for (c in 1:length(cell_list)){
    #print(c)
    in_path <- paste0("results_v7/final_tables_for_function/bw_subset/",cell_list[c],"_norm_sub.bw")
    bw_gr <- import(in_path,which=se_gr_bins,
                    as="GRanges")


    # sum up score to each bins
    overlap_bin <- findOverlaps(se_gr_bins,bw_gr)

    agg <- aggregate(bw_gr, overlap_bin, score=sum(score),drop=F)
    se_bin_df <- data.frame(se_gr_bins)

    se_bin_df$score <- agg$score

    # add se names
    se_bin_df$se_name <- se_tmp$V4

    # add ce names (for color group)
    se_bin_df$ce_name <- NA
    overlap_ce <- findOverlaps(se_gr_bins,ce_gr)
    se_bin_df$ce_name[queryHits(overlap_ce)] <- mcols(ce_gr)$ce_name[subjectHits(overlap_ce)]

    # replace NA ce_name with dummy ce_name
    se_bin_df$ce_name[which(is.na(se_bin_df$ce_name))] <- "dummy"

    # add x-axis name
    se_bin_df$x_axis <- paste0("x_",seq(1:nrow(se_bin_df)))
    # add cell line name
    se_bin_df$Cell_line <- cell_list[c]

    # make final plot_df
    se_bin_df <- se_bin_df[,-c(1:5)]

    out_df <- rbind(out_df,se_bin_df)

  }
  return(out_df)
}

# Transfer SE specific string to dataframe
string2df <- function(spec_string,se_cell_spec_file) {
  ce_list_tmp <- unlist(strsplit(spec_string,"\\|"))
  out_df <- data.frame()
  for (i in 1:length(ce_list_tmp)){
    ce_tmp <- ce_list_tmp[i]
    ce_name_tmp <- unlist(strsplit(unlist(strsplit(ce_tmp,":")),","))[1]
    spec_name_tmp <- unlist(strsplit(unlist(strsplit(ce_tmp,":")),","))[-1]

    tmp_df <- data.frame(spec_type = spec_name_tmp,
                         ce_name=rep(ce_name_tmp,length(spec_name_tmp)))

    out_df <- rbind(out_df,tmp_df)
  }
  return(out_df)
}
