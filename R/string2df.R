#' Transfer SE specific files to dataframe
#'
#' A function to filter super-enhancers (SEs)
#'
#' @details
#' This function firstly remove SEs whose width are within lower and upper 1 percent.
#' Then filter out SEs in the range of super-enhancer blacklist.
#' Finally, remained SEs from different samples are merged to the longest representative SE.
#'
#' @param se_in merged SE region bed file of all samples.
#' @param bl_file super-enhancer blacklist bed file download from ENCODE (ENCFF356LFX).
#' @param custom_range a vector of extra customized blacklist to ignore.
#' Format: c(chr:start-end, chr:start-end, ...).
#'
#' @return
#' se_filtered_no_bl: filtered SE datasets with merged SE names,
#' se_filtered_in_bl: SE datasets identified in the blacklist
#' se_merged_meta: merged SE metadata
#'
#' @import data.table
#'
#' @export
#' @examples
#' # run without blacklist file
#' se_list <- SEfilter(se_in)
#'
#' # run with blacklist file
#' string2df <- SEfilter(se_in,bl_file,has_bl_file=TRUE)
#'
#'
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
