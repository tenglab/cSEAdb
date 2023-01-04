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
