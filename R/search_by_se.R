#' Search database with query and return SE specificity
#'
#' A function to return the overlapped SE regions with a SE query string or SE bed file
#'
#' @details
#' This function will return the SE regions in the database which are overlapped
#' with a query SE region.
#'
#' @param query a query to search.
#' Can be SE regions (with format of "chr_start_end"), SE bed file, cell name,
#' cancer name, or tissue name.(The full list of cell name, cancer name and tissue name
#' available in this database can be found in the description file)
#'
#' @param query_type a character parameter indicates which type of query it is:
#' "se_region", "se_bed","cancer", "cell", or "tissue". (Default: se_region)
#'
#' @param se_db the cSEAdb object loaded from the package
#'
#' @return
#' if query_type is "se_region" or "se_bed":
#' A list of SE regions and there cell and cancer specificity will return.
#'
#' if query_type is "cancer" or "cell":
#' A list of SE regions which only specific in "cancer" or "cell" will return.
#'
#' if query_type is "tissue":
#' A list of SE regions which specific in "cancer" and "cell" within the tissue will return.
#'
#' @import tidyr
#'
#' @export
#' @examples
#'
#' # search with SE string
#' hit_se <- search_db(c("chr1_121000_123250","chr5_111111_222222"),query_type="se_region")
#'
#' # search with SE bed file
#' hit_se <- search_db("se_query_list.bed",query_type="se_bed")
#'
#' # search with cell name
#' hit_se <- search_db(c("HL-60","MCF7"),query_type="cell")
#'
#' # search with cancer name
#' hit_se <- search_db(c("Amelanotic melanoma","Prostate carcinoma"),query_type="cancer")
#'
#' # search with SE bed file
#' hit_se <- search_db(c("Colon","Breast"),query_type="tissue")

search_db <- function(query, query_type="se_region",se_db) {
  #----------------------------------
  # query is se_region or se_bed
  #----------------------------------
  if (query_type %in% c("se_region","se_bed")){
    # make query SE gr
    if (query_type=="se_region") {
      se_region_df <- as.data.frame(query)
      se_region_df <- se_region_df %>%
        separate(query,c("V1","V2","V3"))

      q_gr <- GRanges(
        seqnames=se_region_df[,1],
        ranges=IRanges(as.integer(se_region_df[,2]),
                       as.integer(se_region_df[,3]))
      )
    } else if (query_type=="se_bed"){
      se_bed <- read.table(query,sep="\t",header=F)
      q_gr <- GRanges(
        seqnames=se_bed[,1],
        ranges=IRanges(as.integer(se_bed[,2]),as.integer(se_bed[,3]))
      )
    }

    # create db gr
    db_gr <- GRanges(
      seqnames=se_db$se_region$V1,
      ranges=IRanges(as.integer(se_db$se_region$V2),
                     as.integer(se_db$se_region$V3),
                     names = se_db$se_region$V4)
    )

    #------------------------------
    # overlap query and db
    overlap <- findOverlaps(q_gr,db_gr)
    db_hit <- names(db_gr[subjectHits(overlap),])

    # retrieve cell and cancer specific information
    cell_spec <- se_db$se_cell_specificity[which(se_db$se_cell_specificity$se_name %in% db_hit),]
    cancer_spec <- se_db$se_cancer_specificity[which(se_db$se_cancer_specificity$se_name %in% db_hit),]

    #------------------------------
    # create summary dataframe
    # cell gain
    cell_gain_summary <- cell_spec[,-c(5:7)] %>%
      separate_rows(gain_ce,sep="\\|")
    cell_gain_summary <- cell_gain_summary %>%
      separate(gain_ce,c("ce_name","spec_cell"),":")

    # cell loss
    cell_loss_summary <- cell_spec[,-c(2:4)] %>%
      separate_rows(loss_ce,sep="\\|")
    cell_loss_summary <- cell_loss_summary %>%
      separate(loss_ce,c("ce_name","spec_cell"),":")

    # make final cell summary dataframe
    colnames(cell_gain_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cell")
    cell_gain_summary$spec_type <- ifelse(cell_gain_summary$cell_specific=="yes","gain","none")

    colnames(cell_loss_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cell")
    cell_loss_summary$spec_type <- ifelse(cell_loss_summary$cell_specific=="yes","loss","none")

    cell_spec_summary <- rbind(cell_gain_summary,cell_loss_summary)

    #-----------------------------
    # cancer gain
    cancer_gain_summary <- cancer_spec[,-c(5:7)] %>%
      separate_rows(gain_ce,sep="\\|")
    cancer_gain_summary <- cancer_gain_summary %>%
      separate(gain_ce,c("ce_name","spec_cancer"),":")

    # cancer loss
    cancer_loss_summary <- cancer_spec[,-c(2:4)] %>%
      separate_rows(loss_ce,sep="\\|")
    cancer_loss_summary <- cancer_loss_summary %>%
      separate(loss_ce,c("ce_name","spec_cell"),":")

    # make final cancer summary dataframe
    colnames(cancer_gain_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cancer")
    cancer_gain_summary$spec_type <- ifelse(cancer_gain_summary$cell_specific=="yes","gain","none")

    colnames(cancer_loss_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cancer")
    cancer_loss_summary$spec_type <- ifelse(cancer_loss_summary$cell_specific=="yes","loss","none")

    cancer_spec_summary <- rbind(cancer_gain_summary,cancer_loss_summary)

    # final output
    out_list <- list(cell_spec_summary=cell_spec_summary,
                     cancer_spec_summary=cancer_spec_summary)
  }
  #----------------------------------
  # query is cell
  #----------------------------------
  else if (query_type == "cell") {
    n_cell <- length(query)
    cell_spec <- data.frame()
    for (c in 1:n_cell) {
      # retrieve cell and cancer specific information
      cell_tmp <- se_db$se_cell_specificity[which(se_db$se_cell_specificity$gain_ce %like% query[c] |
                                                    se_db$se_cell_specificity$loss_ce %like% query[c]),]
      cell_spec <- rbind(cell_spec,cell_tmp)
    }

    #------------------------------
    # create summary dataframe
    # cell gain
    cell_gain_summary <- cell_spec[,-c(5:7)] %>%
      separate_rows(gain_ce,sep="\\|")
    cell_gain_summary <- cell_gain_summary %>%
      separate(gain_ce,c("ce_name","spec_cell"),":")

    # cell loss
    cell_loss_summary <- cell_spec[,-c(2:4)] %>%
      separate_rows(loss_ce,sep="\\|")
    cell_loss_summary <- cell_loss_summary %>%
      separate(loss_ce,c("ce_name","spec_cell"),":")

    # make final cell summary dataframe
    colnames(cell_gain_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cell")
    cell_gain_summary$spec_type <- ifelse(cell_gain_summary$cell_specific=="yes","gain","none")

    colnames(cell_loss_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cell")
    cell_loss_summary$spec_type <- ifelse(cell_loss_summary$cell_specific=="yes","loss","none")

    cell_spec_tmp <- rbind(cell_gain_summary,cell_loss_summary)

    # only keep query cancer
    n_cell <- length(query)
    cell_spec_summary <- data.frame()
    for (c in 1:n_cancer) {
      # retrieve cell and cancer specific information
      cell_tmp <- cell_spec_tmp[which(cell_spec_tmp$cell %like% query[c]),]
      cell_spec_summary <- rbind(cell_spec_summary,cell_tmp)
    }

    # final output
    out_list <- list(cell_spec_summary=cell_spec_summary)
  }
  #----------------------------------
  # query is cancer
  #----------------------------------
  else if (query_type == "cancer") {
    n_cancer <- length(query)
    cancer_spec <- data.frame()
    for (c in 1:n_cancer) {
      # retrieve cell and cancer specific information
      cancer_tmp <- se_db$se_cancer_specificity[which(se_db$se_cancer_specificity$gain_ce %like% query[c] |
                                                    se_db$se_cancer_specificity$loss_ce %like% query[c]),]
      cancer_spec <- rbind(cancer_spec,cancer_tmp)
    }

    #-----------------------------
    # cancer gain
    cancer_gain_summary <- cancer_spec[,-c(5:7)] %>%
      separate_rows(gain_ce,sep="\\|")
    cancer_gain_summary <- cancer_gain_summary %>%
      separate(gain_ce,c("ce_name","spec_cancer"),":")

    # cancer loss
    cancer_loss_summary <- cancer_spec[,-c(2:4)] %>%
      separate_rows(loss_ce,sep="\\|")
    cancer_loss_summary <- cancer_loss_summary %>%
      separate(loss_ce,c("ce_name","spec_cell"),":")

    # make final cancer summary dataframe
    colnames(cancer_gain_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cancer")
    cancer_gain_summary$spec_type <- ifelse(cancer_gain_summary$cell_specific=="yes","gain","none")

    colnames(cancer_loss_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cancer")
    cancer_loss_summary$spec_type <- ifelse(cancer_loss_summary$cell_specific=="yes","loss","none")

    cancer_spec_tmp <- rbind(cancer_gain_summary,cancer_loss_summary)

    # only keep query cancer
    n_cancer <- length(query)
    cancer_spec_summary <- data.frame()
    for (c in 1:n_cancer) {
      # retrieve cell and cancer specific information
      cancer_tmp <- cancer_spec_tmp[which(cancer_spec_tmp$cancer %like% query[c]),]
      cancer_spec_summary <- rbind(cancer_spec_summary,cancer_tmp)
    }

    # final output
    out_list <- list(cancer_spec_summary=cancer_spec_summary)
  }
  #----------------------------------
  # query is tissue
  #----------------------------------
  else if (query_type == "tissue") {
    # extract all cell and cancer in the query tissue
    query_meta <- se_db$cell_cancer_meta[which(se_db$cell_cancer_meta$Tissue %in% query),]
    cell_q <- unique(query_meta$Cell_line)
    cancer_q <- unique(query_meta$Cancer)

    n_tissue <- length(query)
    #------------------------------------
    # cell
    #------------------------------------
    n_cell <- length(cell_q)
    cell_spec <- data.frame()
    for (c in 1:n_cell) {
      # retrieve cell and cancer specific information
      cell_tmp <- se_db$se_cell_specificity[which(se_db$se_cell_specificity$gain_ce %like% cell_q[c] |
                                                    se_db$se_cell_specificity$loss_ce %like% cell_q[c]),]
      cell_spec <- rbind(cell_spec,cell_tmp)
    }

    #------------------------------
    # create summary dataframe
    # cell gain
    cell_gain_summary <- cell_spec[,-c(5:7)] %>%
      separate_rows(gain_ce,sep="\\|")
    cell_gain_summary <- cell_gain_summary %>%
      separate(gain_ce,c("ce_name","spec_cell"),":")

    # cell loss
    cell_loss_summary <- cell_spec[,-c(2:4)] %>%
      separate_rows(loss_ce,sep="\\|")
    cell_loss_summary <- cell_loss_summary %>%
      separate(loss_ce,c("ce_name","spec_cell"),":")

    # make final cell summary dataframe
    colnames(cell_gain_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cell")
    cell_gain_summary$spec_type <- ifelse(cell_gain_summary$cell_specific=="yes","gain","none")

    colnames(cell_loss_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cell")
    cell_loss_summary$spec_type <- ifelse(cell_loss_summary$cell_specific=="yes","loss","none")

    cell_spec_tmp <- rbind(cell_gain_summary,cell_loss_summary)

    # only keep query cancer
    cell_spec_summary <- data.frame()
    for (c in 1:n_cancer) {
      # retrieve cell and cancer specific information
      cell_tmp <- cell_spec_tmp[which(cell_spec_tmp$cell %like% cell_q[c]),]
      cell_spec_summary <- rbind(cell_spec_summary,cell_tmp)
    }

    #------------------------------------
    # cancer
    #------------------------------------
    n_cancer <- length(cancer_q)
    cancer_spec <- data.frame()
    for (c in 1:n_cancer) {
      # retrieve cell and cancer specific information
      cancer_tmp <- se_db$se_cancer_specificity[which(se_db$se_cancer_specificity$gain_ce %like% cancer_q[c] |
                                                        se_db$se_cancer_specificity$loss_ce %like% cancer_q[c]),]
      cancer_spec <- rbind(cancer_spec,cancer_tmp)
    }

    #-----------------------------
    # cancer gain
    cancer_gain_summary <- cancer_spec[,-c(5:7)] %>%
      separate_rows(gain_ce,sep="\\|")
    cancer_gain_summary <- cancer_gain_summary %>%
      separate(gain_ce,c("ce_name","spec_cancer"),":")

    # cancer loss
    cancer_loss_summary <- cancer_spec[,-c(2:4)] %>%
      separate_rows(loss_ce,sep="\\|")
    cancer_loss_summary <- cancer_loss_summary %>%
      separate(loss_ce,c("ce_name","spec_cell"),":")

    # make final cancer summary dataframe
    colnames(cancer_gain_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cancer")
    cancer_gain_summary$spec_type <- ifelse(cancer_gain_summary$cell_specific=="yes","gain","none")

    colnames(cancer_loss_summary) <- c("se_name","cell_specific","n_ce_spec","ce_name","cancer")
    cancer_loss_summary$spec_type <- ifelse(cancer_loss_summary$cell_specific=="yes","loss","none")

    cancer_spec_tmp <- rbind(cancer_gain_summary,cancer_loss_summary)

    # only keep query cancer
    cancer_spec_summary <- data.frame()
    for (c in 1:n_cancer) {
      # retrieve cell and cancer specific information
      cancer_tmp <- cancer_spec_tmp[which(cancer_spec_tmp$cancer %like% cancer_q[c]),]
      cancer_spec_summary <- rbind(cancer_spec_summary,cancer_tmp)
    }

    # final output
    out_list <- list(query_meta=query_meta,
                     cell_spec_summary=cell_spec_summary,
                     cancer_spec_summary=cancer_spec_summary)
  }

  return(out_list)
}











