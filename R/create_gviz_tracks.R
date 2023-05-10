#' Create Gviz plot tracks for specific regions
#'
#' A function to create Gviz data tracks and annotation tracks with specific CEs
#'
#' @details
#' This function will create different Gviz data tracks based on the selected cell and
#' annotation tracks of specific CEs and SEs
#'
#' @param bw_gr_list list of bigwig Granges generated from "create_be_gr" function
#' @param se_spec_table se specifity table generated from "search_db" function
#' @param cell a comman separated string of which cell to plot.
#'              Examples:"cell1,cell2,cell3","specific","specific,cell1,cell2".
#'             (Default: specific. Cell have specific CEs within the plot region.)
#' @param plot_region a string indicate region of interest. (Example:"chr1_1109053_1176954"or "chr1:112388326-112392055")
#' @param txdb a txdb file for gene annotation track. (Default: "none" No gene annotation track)
#'
#' @return
#' A output list with different type of gviz tracks:
#'     htrack: major data tracks with highlighted specific CEs.
#'     atrack_ce: annotation track indicating CE locations.
#'     atrack_se: annotation track indicating SE locations overlap with plot region.
#'     atrack_query: annotation track indicating plot region location.
#'     gtrack: gene annotation track.
#'     plot_region: vector of plot region information (chr,start,end).
#'     axistrack: genome axis track.
#'
#' @import GenomicRanges
#' @import rtracklayer
#' @import Gviz
#' @import data.table
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#'
#' @examples
#' # region of interests
#' plot_region <- "chr1_1109053_1176954"
#'
#' # search db to find overlap specific SEs
#' search_db_table <- search_db(plot_region,query_type="se_region",cSEAdb)
#'
#' # create bw_gr list
#' bw_gr <- create_bw_gr(plot_region,bigWigs$V1)
#'
#' # default: only plot specific cells with no gene track
#' tracks_default <- create_gviz_tracks(bw_gr_list=bw_gr,
#'                              se_spec_table=search_db_table,
#'                              plot_region=plot_region)
#'
#' # load txdb file
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' # plot selected cells within bw_gr and gene track from txdb object
#' tracks <- create_gviz_tracks(bw_gr_list=bw_gr,
#'                              se_spec_table=search_db_table,
#'                              cell="MCF7,A549,K-562,specific",
#'                              plot_region=plot_region,
#'                              txdb=TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' # plot tracks using Gviz
#' plotTracks(c(plot_out$htrack,plot_out$atrack_ce,
#'              plot_out$atrack_se,plot_out$gtrack,plot_out$atrack_query,
#'              plot_out$axistrack),
#'            chr=plot_out$plot_region[1],
#'            from=as.integer(plot_out$plot_region[2]),
#'            to=as.integer(plot_out$plot_region[3]))


create_gviz_tracks <- function(bw_gr_list,se_spec_table,cell="specific",plot_region,txdb="none") {
  #-------------------------
  # plot region vector
  #-------------------------
  plot_region_v <- unlist(strsplit(plot_region,"_|:|-"))
  #-------------------------
  # higlight specific gain or loss ce region
  #-------------------------
  spec_cell_table <- se_spec_table[which(se_spec_table$object_type=="cell" & se_spec_table$number_spec_ce!=0),]
  ce_hl_2 <- separate(unique(spec_cell_table[,c("spec_ce","spec_object","specifity")]),
                      spec_ce,c("chr","start","end"),remove=F)

  se_hl <- separate(unique(spec_cell_table[,c("se_name","query")]),
                    se_name,c("chr","start","end"),remove=F)
  # cell have gain ce
  ce_hl_gain <- ce_hl_2[which(ce_hl_2$specifity=="gain"),]
  spec_gain <- paste(ce_hl_gain$spec_object,collapse=",")
  spec_gain <- unique(unlist(strsplit(spec_gain,",")))

  # cell have loss ce
  ce_hl_loss <- ce_hl_2[which(ce_hl_2$specifity=="loss"),]
  spec_loss <- paste(ce_hl_loss$spec_object,collapse=",")
  spec_loss <- unique(unlist(strsplit(spec_loss,",")))


  #--------------------------
  # annotation tracks
  #--------------------------
  # axistrack
  axistrack <- GenomeAxisTrack(chr=plot_region_v[1],
                               from=as.integer(plot_region_v[2]),
                               to=as.integer(plot_region_v[3]),
                               scale=0.5,
                               labelPos="below")


  # gene annotation track
  if (is.character(txdb)) {
    gtrack <- NULL
  } else {
    gtrack <- GeneRegionTrack(txdb,
                              chromosome=plot_region_v[1],
                              start=as.integer(plot_region_v[2]),
                              end=as.integer(plot_region_v[3]),
                              collapseTranscripts="meta",
                              transcriptAnnotation = "symbol",
                              name="hg38")
  }


  # ce annotation track
  atrack_ce <- AnnotationTrack(start = ce_hl_2$start,
                               end=ce_hl_2$end,
                               chromosome = ce_hl_2$chr,
                               name = "spec_ce",
                               shape = "box",
                               fill="grey70")

  # se annotation track
  atrack_se <- AnnotationTrack(start = se_hl$start,
                               end=se_hl$end,
                               chromosome = se_hl$chr,
                               name = "hit_se",
                               shape = "box",
                               fill="grey50")

  # plot region annotation track
  atrack_query <- AnnotationTrack(start = plot_region_v[2],
                               end=plot_region_v[3],
                               chromosome = plot_region_v[1],
                               name = "query",
                               fill="orange")

  # highlighted datatrack of cell to plot
  # default is all cells have specific

  spec_cell_tmp <- unlist(strsplit(cell,","))
  if ("specific" %in% spec_cell_tmp) {
    spec_cell <- unique(c(spec_gain,spec_loss,spec_cell_tmp[spec_cell_tmp!="specific"]))
  } else {
    spec_cell <- spec_cell_tmp
  }

  htrack <- list()
  for(i in 1:length(spec_cell)){
    #print(i)
    if (spec_cell[i] =="X786-0") {
      spec_cell[i] <- "786-0"
    }
    track_cell_tmp <- bw_gr_list[[which(names(bw_gr_list)==spec_cell[i])]]
    trackList_tmp <- DataTrack(range=track_cell_tmp,type="histogram",
                               name=spec_cell[i],
                               col.histogram="grey60",
                               fill.histogram="grey60",
                               ylim=c(0,200))


    # highlight gain as red, loss as blue
    ce_hl_gain_tmp <- ce_hl_gain[which(ce_hl_gain$spec_object %like% spec_cell[i]),]
    ce_hl_loss_tmp <- ce_hl_loss[which(ce_hl_loss$spec_object %like% spec_cell[i]),]

    if (nrow(ce_hl_gain_tmp)==0 & nrow(ce_hl_loss_tmp)==0) {
      ht_tmp <- trackList_tmp
    } else {
      ht_tmp <- HighlightTrack(trackList = trackList_tmp,
                               start = c(ce_hl_gain_tmp$start,ce_hl_loss_tmp$start),
                               end = c(ce_hl_gain_tmp$end,ce_hl_loss_tmp$end),
                               chromosome = c(ce_hl_gain_tmp$chr,ce_hl_loss_tmp$chr),
                               fill=c(rep("red",nrow(ce_hl_gain_tmp)),
                                      rep("blue",nrow(ce_hl_loss_tmp))),
                               alpha=0.4)
    }


    htrack[[i]] <- ht_tmp
  }


  output <- list(htrack=htrack,
                 atrack_ce=atrack_ce,
                 atrack_se=atrack_se,
                 atrack_query=atrack_query,
                 gtrack=gtrack,
                 plot_region=plot_region_v,
                 axistrack=axistrack)


  return(output)
}
