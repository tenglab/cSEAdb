setwd("~/Projects/super_enhancer/se_data_portal/cSEAdb/")
library(GenomicRanges)
library(rtracklayer)
library(RColorBrewer)
library(tidyr)
library(tidyverse)
library(gtable)
library(ggh4x)
library(ggplot2)
library(grid)

se_db <- readRDS("data/cSEAdb.rds")

#------------------------------------------------------------------------
# create bigwig track plot df for ggplot
#------------------------------------------------------------------------

se_region <- "chr5_151039481_151102681"

# examples: chr5_151039481_151102681, chr1_109927087_110035367, chr1_154931426_155027337, chr13_113177058_113210511,
# double check: se_region <- "chr13_113177058_113210511"
se_region <- "chr13_113177058_113210511"
se_test <- create_se_track_data(se_region,se_ce,ce_matrix,100)
track_plot <- se_plot(se_test,se_cell_spec)

pdf(file="results_v7/final_tables_for_function/test_output/track_plot_example_4_old_mix.pdf",width = 20,height=30)
grid.draw(track_plot)
#plot_out
dev.off()

#------------------------------------------------------------------------
# function of create track dataframe for plot based on se_region
#------------------------------------------------------------------------

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


#------------------------------------------------------------------------
# se plot function
#------------------------------------------------------------------------
se_plot <- function(se_plot_df,se_cell_spec_file) {
  se_plot_df$Cell_line[which(se_plot_df$Cell_line=="786-0")] <- "X786-0"
  # extract se
  se_target <- unique(se_plot_df$se_name)
  se_tmp <- ce_matrix[which(ce_matrix$se_name==se_target),]
  se_cell_spec_tmp <- se_cell_spec_file[which(se_cell_spec_file$se_name==se_target),]

  # add cancer type
  plot_df <- merge(se_plot_df,meta_data[,c(1,2)],by="Cell_line")
  x_axis_sort <- unique(plot_df$x_axis)[order(nchar(unique(plot_df$x_axis)), unique(plot_df$x_axis))]
  plot_df$x_axis <- factor(plot_df$x_axis,levels = x_axis_sort)
  plot_df$Cancer <- factor(plot_df$Cancer,levels = unique(meta_data$Cancer))
  plot_df$Cell_line <- factor(plot_df$Cell_line,levels = unique(meta_data$Cell_line))

  # identify cell and cancer specific
  cell_gain_tmp <- string2df(se_cell_spec_tmp$gain_ce,se_cell_spec_tmp)
  cell_loss_tmp <- string2df(se_cell_spec_tmp$loss_ce,se_cell_spec_tmp)

  # add cell spec group
  plot_df$cell_group <- NA
  if(nrow(cell_gain_tmp) > 0 & nrow(cell_loss_tmp) > 0){
    plot_df$cell_group[which(paste0(plot_df$Cell_line,"_",plot_df$ce_name) %in%
                               paste0(cell_gain_tmp$spec_type,"_",cell_gain_tmp$ce_name))] <- "cell_gain"
    plot_df$cell_group[which(paste0(plot_df$Cell_line,"_",plot_df$ce_name) %in%
                               paste0(cell_loss_tmp$spec_type,"_",cell_loss_tmp$ce_name))] <- "cell_loss"
    plot_df$cell_group[is.na(plot_df$cell_group)] <- "none"
  } else if (nrow(cell_gain_tmp) > 0 & nrow(cell_loss_tmp) == 0){
    plot_df$cell_group[which(paste0(plot_df$Cell_line,"_",plot_df$ce_name) %in%
                               paste0(cell_gain_tmp$spec_type,"_",cell_gain_tmp$ce_name))] <- "cell_gain"
    plot_df$cell_group[is.na(plot_df$cell_group)] <- "none"
  } else if (nrow(cell_gain_tmp) == 0 & nrow(cell_loss_tmp) > 0) {
    plot_df$cell_group[which(paste0(plot_df$Cell_line,"_",plot_df$ce_name) %in%
                               paste0(cell_loss_tmp$spec_type,"_",cell_loss_tmp$ce_name))] <- "cell_loss"
    plot_df$cell_group[is.na(plot_df$cell_group)] <- "none"
  } else if (nrow(cell_gain_tmp) == 0 & nrow(cell_loss_tmp) == 0) {
    plot_df$cell_group[is.na(plot_df$cell_group)] <- "none"
  }

  #------------------------------------
  # cell specific plot
  #------------------------------------
  # add rectange around spec loss
  if (nrow(cell_loss_tmp) > 0) {
    spec_cell_loss <- unique(plot_df[which(plot_df$cell_group=="cell_loss"),])

    tmp_df <- spec_cell_loss[,c("Cell_line","ce_name","Cancer","x_axis")]
    tmp_df$idx <- as.integer(gsub("x_","",tmp_df$x_axis))
    ce_loss <- unique(spec_cell_loss$ce_name)
    rect_df <- data.frame()
    for (i in 1:length(ce_loss)) {
      rect_df_tmp <- tmp_df[which(tmp_df$ce_name==ce_loss[i]),]
      rect_df_tmp$idx_min <- min(rect_df_tmp$idx)
      rect_df_tmp$idx_max <- max(rect_df_tmp$idx)
      rect_df_tmp <- unique(rect_df_tmp[,-c(4,5)])
      rect_df <- rbind(rect_df,rect_df_tmp)
    }

    plot_out <- ggplot(plot_df,aes(x=x_axis, y=score,fill=cell_group))+
      geom_bar(stat="identity")+
      geom_rect(data=rect_df,
                aes(xmin = idx_min, xmax = idx_max, ymin = 0, ymax = Inf),
                alpha = 0, color="blue",inherit.aes = F,size=0.4)+
      facet_nested(Cancer+Cell_line~.,switch="y")+
      theme_classic()+
      ggtitle(se_target)+
      theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            strip.background = element_rect(colour=alpha("black",0.1),
                                            fill=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text.y.left = element_text(angle = 0),
            legend.position = "none",
            panel.spacing = unit(0, "lines"),
            plot.title = element_text(size = 20, face = "bold"),
            strip.text = element_text(size = 10))+
      scale_fill_manual(values=c(cell_gain="red",
                                 none="grey50",
                                 cell_loss="blue"))
  } else {
    plot_out <- ggplot(plot_df,aes(x=x_axis, y=score,fill=cell_group))+
      geom_bar(stat="identity")+
      facet_nested(Cancer+Cell_line~.,switch="y")+
      theme_classic()+
      ggtitle(se_target)+
      theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            strip.background = element_rect(colour=alpha("black",0.1),
                                            fill=NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text.y.left = element_text(angle = 0),
            legend.position = "none",
            panel.spacing = unit(0, "lines"),
            plot.title = element_text(size = 20, face = "bold"),
            strip.text = element_text(size = 10))+
      scale_fill_manual(values=c(cell_gain="red",
                                 none="grey50",
                                 cell_loss="blue"))
  }


  # pdf(file="results_v5/final_tables_for_function/test_output/test_plot_new_2.pdf",width = 15,height=30)
  # grid.draw(g)
  # #plot_out
  # dev.off()

  #ylab("mean of CE number")+
  #xlab("cutoff")
  #-----------------------------------------
  # add strip fort color for cancer and cell
  #-----------------------------------------
  g <- ggplot_gtable(ggplot_build(plot_out))

  strips <- which(grepl('strip-l', g$layout$name))
  # change facet width
  for (i in seq_along(strips)) {
    g$grobs[[strips[i]]]$widths <-unit(c(10.5,4), "cm")

  }
  # change margins
  g$widths[1] <- unit(0,"cm")
  g$widths[5] <- unit(14.8,"cm")

  # map strips number with cell and cancer name
  spec_cell_name_gain <- as.character(unique(plot_df$Cell_line[which(plot_df$cell_group=="cell_gain")]))
  spec_cell_name_loss <- as.character(unique(plot_df$Cell_line[which(plot_df$cell_group=="cell_loss")]))

  if (length(spec_cell_name_gain) > 0 & length(spec_cell_name_loss) > 0) {
    spec_common <- intersect(spec_cell_name_gain,spec_cell_name_loss)

    strips_df <- data.frame(strips_n=strips,
                            cell_cancer_name=c(unique(meta_data$Cancer),
                                               unique(meta_data$Cell_line)))
    strips_df$strp_col <- NA
    strips_df$strp_col[which(strips_df$cell_cancer_name %in% spec_cell_name_gain &
                               !strips_df$cell_cancer_name %in% spec_common)] <-"red"
    strips_df$strp_col[which(strips_df$cell_cancer_name %in% spec_cell_name_loss &
                               !strips_df$cell_cancer_name %in% spec_common)] <-"blue"
    strips_df$strp_col[which(strips_df$cell_cancer_name %in% spec_common)] <-"purple"
    strips_df$strp_col[is.na(strips_df$strp_col)] <- "black"
  } else if (length(spec_cell_name_gain) > 0 & length(spec_cell_name_loss) == 0) {
    strips_df <- data.frame(strips_n=strips,
                            cell_cancer_name=c(unique(meta_data$Cancer),
                                               unique(meta_data$Cell_line)))
    strips_df$strp_col <- NA
    strips_df$strp_col[which(strips_df$cell_cancer_name %in% spec_cell_name_gain)] <-"red"
    strips_df$strp_col[is.na(strips_df$strp_col)] <- "black"
  } else if (length(spec_cell_name_gain) == 0 & length(spec_cell_name_loss) > 0) {
    strips_df <- data.frame(strips_n=strips,
                            cell_cancer_name=c(unique(meta_data$Cancer),
                                               unique(meta_data$Cell_line)))
    strips_df$strp_col <- NA
    strips_df$strp_col[which(strips_df$cell_cancer_name %in% spec_cell_name_loss)] <-"blue"
    strips_df$strp_col[is.na(strips_df$strp_col)] <- "black"
  } else if (length(spec_cell_name_gain) == 0 & length(spec_cell_name_loss) == 0){
    strips_df <- data.frame(strips_n=strips,
                            cell_cancer_name=c(unique(meta_data$Cancer),
                                               unique(meta_data$Cell_line)))
    strips_df$strp_col <- "black"
  }




  for (i in seq_along(strips)) {
    l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- strips_df$strp_col[i]
  }

  return(g)


}


#------------------------------------
# functions
#------------------------------------

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

