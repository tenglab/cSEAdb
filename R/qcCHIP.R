#' Guide the selection of VAF, DP, or population metric through permutation analysis
#'
#' A function of generating plots of permutation consistency based on different setting of VAF, DP, or population metric.
#'
#' @details
#' This function will split the whole input file to N groups.
#' A permutation will be used to compare the results from subseting sample and whole sample.
#' Finally, a summary file and a permutation consistency plot will be generated to guide the selection of VAF, DP, or population metric.
#'
#' @param input_df a data frame of text annotated file with proper column names. (detailed in instruction)
#' @param out_path path of output directory
#' @param permut_metrics which metric want to permute.
#'                       VAF: minimum VAF setting;
#'                       DP: minimum read depth;
#'                       population: maximum percentage of sample size allowed for one variant.
#' @param permut_n number of permutations.
#' @param sub_sample_n integers of subset sample groups.
#'                     For example: c(2,5) means the function will split the whole sample into 2 and 5 groups, each group has same size.
#' @param metric_min minimum value to start with.
#' @param metric_max maximum value to end with.
#' @param metric_step the step of metric changes.
#' @param remove_tmp_dir if remove temporary directory.
#' @param core number of cores used for parallel computing.
#' @param blacklist_f a bed file contain blacklist regions to exclude.
#' @param show_info if print processing massage.
#'
#' @return
#' A text file with comparison results of subsample and whole sample will be saved and written to out_path.
#' A pdf file including 3 permutation consistency plots will be written to out_path.
#' A temp directory contains all permutation and whole sample output files of metric setting for all sub-sample approaches (if remove_tmp_dir==F).
#'
#' @import data.table
#' @import GenomicRanges
#' @import ggplot2
#' @import ggpubr
#' @import parallel
#'
#' @export
#' @examples
#' # input file
#' in_df <- fread("data/demo_input.txt")
#' bl_f <- fread("blacklist.bed")
#'
#' # permutation on VAF
#' qcCHIP(input_df,
#'       out_path = "out_path",
#'       metric_min = 0,
#'       metric_step = 0.005,
#'       metric_max=0.1,
#'       permut_metrics = "VAF")
#'
#' # permutation on DP with blacklist region
#' qcCHIP(input_df,
#'       out_path = "out_path",
#'       metric_min = 0,
#'       metric_step = 5,
#'       metric_max=50,
#'       permut_metrics = "DP",
#'       blacklist_f = bl_f)
#'
#' # permutation on population with blacklist region and do not print info
#' qcCHIP(input_df,
#'       out_path = "out_path",
#'       metric_min = 0.05,
#'       metric_step = 0.005,
#'       metric_max=0.2,
#'       permut_metrics = "population",
#'       blacklist_f = bl_f,
#'       show_info=F)

qcCHIP <- function(input_df,
                   blacklist_f=F,
                   out_path,
                   permut_metrics="VAF",
                   permut_n=20,
                   sub_sample_n=c(2,5,10),
                   metric_min=0,
                   metric_max=0.1,
                   metric_step=0.05,
                   remove_tmp_dir=T,
                   core=2,
                   show_info=T) {

  ###########
  #create output directors
  #######
  if(!dir.exists(out_path)) {
    dir.create(out_path,recursive=T)
  }

  # tmp director
  tmp_dir <- paste0(out_path,"/temp")
  if(!dir.exists(tmp_dir)) {
    dir.create(tmp_dir,recursive=T)
  }

  # whole sample file director
  whole_sample_dir <- paste0(tmp_dir,"/whole_sample")
  if(!dir.exists(whole_sample_dir)) {
    dir.create(whole_sample_dir,recursive=T)
  }

  #####################
  # main: permutation and generate summary df for plot
  #####################

  seting_list <- seq(metric_min,metric_max,metric_step)
  out_seting_name <- paste0(permut_metrics,gsub("0\\.","_",seting_list))
  summary_tmp <- data.frame()
  for (m in 1:length(seting_list)) {
    if (show_info==T) {
      print(paste0(permut_metrics,": ",seting_list[m]))
    }

    #--------
    # whole sample
    #--------
    # create director
    out_name <- paste0(whole_sample_dir,"/",out_seting_name[m],".txt")

    if (permut_metrics=="VAF") {
      out_tmp <- CHIPfilter(input=input_df,
                                      VAF_min = seting_list[m],
                                      blacklist_f = blacklist_f,
                            info=F)

    } else if (permut_metrics=="population") {
      out_tmp <- CHIPfilter(input=input_df,
                                      max_percent = seting_list[m],
                                      blacklist_f = blacklist_f,
                            info=F)
    } else if (permut_metrics=="DP") {
      out_tmp <- CHIPfilter(input=input_df,
                                      DP_min = seting_list[m],
                                      blacklist_f = blacklist_f,
                            info=F)
    }

    write.table(out_tmp,out_name,sep="\t",quote=F,row.names=F)

    #--------
    # subsample permutation
    #--------
    n_total <- length(unique(input_df$SampleID))
    # sampleID_df <- data.frame(sampleID=unique(input_df$SampleID),
    #                           seq=seq(1,n_total))

    sub_size <- c(ceiling(n_total/sub_sample_n))

    for (s in 1:length(sub_size)) {
      n_sub <- sub_size[s]
      if (show_info==T) {
        print(sub_size[s])
      }

      subsample_dir <- paste0(tmp_dir,"/",sub_size[s],"/",out_seting_name[m])
      if(!dir.exists(subsample_dir)) {
        dir.create(subsample_dir,recursive=T)
      }

      # loop permutation
      for (p in 1:permut_n) {
        if (show_info==T) {
          if (p%%10==0) {
            print(paste0("permutation:",p))
          }
        }

        sub_out_name <- paste0(subsample_dir,"/","permut_",p,".txt")

        # random assgin sample ids
        n_split <- sub_sample_n[s]
        group_labels <- rep(1:n_split, each = n_sub)
        random_groups <- sample(group_labels,size=n_total)
        tmp1 <- data.frame(Element = unique(input_df$SampleID), Group = random_groups)

        g_name <- unique(group_labels)
        out_all <- data.frame()

        split_list <- list()
        for (g in 1:length(g_name)) {
          split_list[[g]] <- input_df[which(input_df$SampleID %in% tmp1$Element[which(tmp1$Group==g_name[g])]),]
        }

        # run CHIPfilter paralle
        if (permut_metrics=="VAF") {
          out_sub <- do.call(rbind,mclapply(split_list,
                                            CHIPfilter,
                                            VAF_min = seting_list[m],
                                            blacklist_f = blacklist_f,
                                            info=F,
                                            mc.cores = core))

        } else if (permut_metrics=="population") {
          out_sub <- do.call(rbind,mclapply(split_list,
                                            CHIPfilter,
                                            max_percent = seting_list[m],
                                            blacklist_f = blacklist_f,
                                            info=F,
                                            mc.cores = core))
        } else if (permut_metrics=="DP") {
          out_sub <- do.call(rbind,mclapply(split_list,
                                            CHIPfilter,
                                            DP_min = seting_list[m],
                                            blacklist_f = blacklist_f,
                                            info=F,
                                            mc.cores = core))
        }

        write.table(out_sub,sub_out_name,sep="\t",quote=F,row.names=F)

        #-------------------
        # compare each permutation with whole sample
        # save the common counts in a dataframe
        #-------------------
        # union
        total_whole_sub <- unique(c(out_tmp$mut_sample,out_sub$mut_sample))

        # common
        common_whole_sub<- unique(intersect(out_tmp$mut_sample,out_sub$mut_sample))

        # only in all
        whole_only <- out_tmp$mut_sample[!(out_tmp$mut_sample %in% common_whole_sub)]

        # only in shuffle
        sub_only <- out_sub$mut_sample[!(out_sub$mut_sample %in% common_whole_sub)]

        # summary data frame
        summary_tmp_tmp <- data.frame(metric_name=permut_metrics,
                                      metric_setting=seting_list[m],
                              group_size=sub_sample_n[s],
                              permut_index=p,
                              var_n_whole=length(unique(out_tmp$mut_sample)),
                              var_n_sub=length(unique(out_sub$mut_sample)),
                              union_n=length(total_whole_sub),
                              common_n=c(length(common_whole_sub)),
                              whole_only=c(length(whole_only)),
                              sub_only=c(length(sub_only)),
                              common_whole=c(round(length(common_whole_sub)/length(unique(out_tmp$mut_sample)),
                                                 digits = 3)),
                              common_sub=c(round(length(common_whole_sub)/length(unique(out_sub$mut_sample)),
                                                 digits = 3)))

        summary_tmp <- rbind(summary_tmp,summary_tmp_tmp)

      }
    }
  }

  write.table(summary_tmp,paste0(out_path,"/",permut_metrics,"_compare_summary.txt"),sep="\t",quote=F,row.names=F)
  if (remove_tmp_dir==T) {
    unlink(tmp_dir, recursive = TRUE)
  }

  ##################
  # plot generation
  ##################

  plot_tmp <- summary_tmp
  plot_tmp$group_size <- as.factor(plot_tmp$group_size)
  plot_tmp$average <- (plot_tmp$common_sub+plot_tmp$common_whole)/2

  plot_mean <- aggregate(.~metric_setting+group_size,plot_tmp[,-c(1,4)],FUN=mean)

  # average
  p0 <- ggplot(plot_tmp,aes(x=metric_setting,y=average,color=group_size))+
    geom_point(data=plot_mean,aes(x=metric_setting,y=average,color=group_size),
               size=1)+
    geom_smooth(se=T)+
    #geom_line()+
    theme_bw()+
    theme(strip.text = element_text(size=15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab(unique(plot_tmp$metric_name))+
    ylab("Permutation consistency")+
    ggtitle("Comparision with whole and sub sample size (average)")

  # common in whole sample
  p1 <- ggplot(plot_tmp,aes(x=metric_setting,y=common_whole,color=group_size))+
    geom_point(data=plot_mean,aes(x=metric_setting,y=common_whole,color=group_size),
               size=1)+
    geom_smooth()+
    #geom_line()+
    theme_bw()+
    theme(strip.text = element_text(size=15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab(unique(plot_tmp$metric_name))+
    ylab("Permutation consistency")+
    ggtitle("Comparision with whole and sub sample size (whole sample)")

  # common in sub-sample
  p2 <- ggplot(plot_tmp,aes(x=metric_setting,y=common_sub,color=group_size))+
    geom_point(data=plot_mean,aes(x=metric_setting,y=common_sub,color=group_size),
               size=1)+
    geom_smooth()+
    #geom_line()+
    theme_bw()+
    theme(strip.text = element_text(size=15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab(unique(plot_tmp$metric_name))+
    ylab("Permutation consistency")+
    ggtitle("Comparision with whole and sub sample size (sub-sample)")


  p_out <- ggarrange(p0,p1,p2,common.legend = T)

  pdf(paste0(out_path,"/",permut_metrics,"_compare_summary.pdf"),width = 12,height=8)
  print(p_out)
  dev.off()

  out_list <- list(summary_df=summary_tmp,
                   figs=p_out)
  return(out_list)


}











