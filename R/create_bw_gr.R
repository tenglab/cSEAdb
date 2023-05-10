#' Create coverage GRanges from bigwig file with specific genome ranges
#'
#' A function to create a coverage Granges of bigwig file
#'
#' @details
#' This function will return a list of coverage granges of multiple bigwig file at specific regions
#'
#' @param plot_region a string of region of interest
#' @param bw_path a vector of bigwig file path
#' @param se_spec_table se specifity table generated from "search_db" function
#'
#' @return
#' A list of bigwig Granges for each bigwig file and plot region gr
#'
#' @import GenomicRanges
#' @import rtracklayer
#' @import tidyr
#'
#' @examples
#' # bigwig file path
#' bw_file_path <- c("merged_normalized_bw/786-0_normalized.bw","merged_normalized_bw/A549_normalized.bw")
#' plot_region <- "chr1_1109053_1176954"
#'
#' bw_gr <- create_bw_gr(plot_region,bw_file_path)


create_bw_gr<- function (plot_region,bw_path,se_spec_table) {
  #---------------
  # make plot gr
  #---------------
  plot_region_tmp <- unlist(strsplit(plot_region,"_|-|:"))

  # merge with CE and SE to get the longest region
  se_tmp <-separate(unique(se_spec_table[,c("se_name","query")]),
                    se_name,c("chr","start","end"),remove=F)

  ce_tmp <- separate(unique(se_spec_table[,c("spec_ce","query")]),
                     spec_ce,c("chr","start","end"),remove=F)

  plot_gr_tmp <- GRanges(c(plot_region_tmp[1],se_tmp$chr,ce_tmp$chr),
                         IRanges(c(as.integer(plot_region_tmp[2]),as.integer(se_tmp$start),as.integer(ce_tmp$start)),
                                 c(as.integer(plot_region_tmp[3]),as.integer(se_tmp$end),as.integer(ce_tmp$end))))

  plot_gr <- GenomicRanges::reduce(plot_gr_tmp)

  #----------------
  # make plot list
  #----------------
  bw_gr <- list()
  for (a in 1:length(bw_path)) {
    #print(a)
    in_path <- bw_path[a]
    cell <- gsub("_normalized.bw","",basename(in_path))
    bw_gr_tmp <- import(in_path, format="BigWig",
                        which=plot_gr,
                        as="GRanges")

    bw_gr[[a]] <-  bw_gr_tmp
    names(bw_gr)[a] <- cell
  }

  output <- list(bw_gr=bw_gr,
                 plot_gr=plot_gr)
  return(output)
}
