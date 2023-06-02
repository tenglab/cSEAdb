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
#' "se_region", "se_bed", "gene", "cancer", or "cell". (Default: se_region)
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
#' @import GenomicRanges
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
#' # search with gene symbol
#' hit_se <- search_db(c("MYC","FOXA2"),query_type="gene")
#'
#' # search with cell name
#' hit_se <- search_db(c("HL-60","MCF7"),query_type="cell")
#'
#' # search with cancer name
#' hit_se <- search_db(c("Amelanotic melanoma","Prostate carcinoma"),query_type="cancer")
#'

search_db <- function(query, query_type="se_region",se_db,gene_region=F) {
  # se specificity table
  se_specificity <- se_db$se_specificity
  gene_promotor <- se_db$gene_promotor
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
                       as.integer(se_region_df[,3]),
                       names=paste(se_region_df[,1],se_region_df[,2],se_region_df[,3],sep="_"))
      )
    } else if (query_type=="se_bed"){
      se_bed <- read.table(query,sep="\t",header=F)
      q_gr <- GRanges(
        seqnames=se_bed[,1],
        ranges=IRanges(as.integer(se_bed[,2]),
                       as.integer(se_bed[,3]),
                       names=paste(se_bed[,1],se_bed[,2],se_bed[,3],sep="_"))
      )
    }

    # create db gr
    db_gr <- GRanges(
      seqnames=se_db$se_bed$chr,
      ranges=IRanges(as.integer(se_db$se_bed$start),
                     as.integer(se_db$se_bed$end),
                     names = se_db$se_bed$se_name)
    )

    #------------------------------
    # overlap query and db
    overlap <- findOverlaps(q_gr,db_gr)
    query_hit_df <- data.frame(query=names(q_gr[queryHits(overlap),]),
                               se_name=names(db_gr[subjectHits(overlap),]))

    # retrieve cell and cancer specific information
    out_table <- merge(query_hit_df,se_specificity,by="se_name")

  }
  #----------------------------------
  # query is gene
  #----------------------------------
  else if (query_type=="gene") {
    gene_df <- gene_promotor[which(gene_promotor$V4 %in% query),]

    # gene +- 500k region
    q_gr <- GRanges(
      seqnames=gene_df[,1],
      ranges=IRanges(as.integer(gene_df[,2]-500000),
                     as.integer(gene_df[,3]+500000),
                     names=gene_df[,4])
    )

    # create db gr
    db_gr <- GRanges(
      seqnames=se_db$se_bed$chr,
      ranges=IRanges(as.integer(se_db$se_bed$start),
                     as.integer(se_db$se_bed$end),
                     names = se_db$se_bed$se_name)
    )

    #------------------------------
    # overlap query and db
    overlap <- findOverlaps(q_gr,db_gr)
    query_hit_df <- data.frame(query=names(q_gr[queryHits(overlap),]),
                               se_name=names(db_gr[subjectHits(overlap),]))

    # retrieve cell and cancer specific information
    out_table <- merge(query_hit_df,se_specificity,by="se_name")
  }
  #----------------------------------
  # query is cell
  #----------------------------------
  else if (query_type == "cell") {
    n_cell <- length(query)
    cell_spec <- data.frame()
    for (c in 1:n_cell) {
      # retrieve cell and cancer specific information
      cell_tmp <- se_specificity[which(se_specificity$spec_object %like% query[c]),]
      cell_spec <- rbind(cell_spec,cell_tmp)
    }
    out_table <- unique(cell_spec)
    }
  #----------------------------------
  # query is cancer
  #----------------------------------
  else if (query_type == "cancer") {
    n_cancer <- length(query)
    cancer_spec <- data.frame()
    for (c in 1:n_cancer) {
      # retrieve cell and cancer specific information
      cancer_name <- gsub(" |,|\\-","\\.",query[c])
      cancer_tmp <- se_specificity[which(se_specificity$spec_object %like% cancer_name),]
      cancer_spec <- rbind(cancer_spec,cancer_tmp)
    }
    out_table <- unique(cancer_spec)
    }

  return(out_table)
}

