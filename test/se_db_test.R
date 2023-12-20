setwd("~/Projects/super_enhancer/se_data_portal/cSEAdb/")
# library(GenomicRanges)
# library(rtracklayer)
# library(Gviz)
# library(tidyr)
# library(data.table)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)

devtools::document()
devtools::load_all()

library(devtools)
devtools::install_github("https://github.com/tenglab/cSEAdb.git")

# load database
cSEAdb <- readRDS("data/cSEAdb.rds")

library(cSEAdb)
#-----------------------------------
# test search_db
#-----------------------------------

# test SE_region
test_se_r <- search_db(c("chr1_112388326_112392055"),
                     query_type="se_region",cSEAdb)


# test SE bed
test_se_b <- search_db("inst/extdata/se_query_example.bed",
                     query_type="se_bed",cSEAdb)

# test gene
test_g<- search_db(query=c("MYC","FOXA2"),
                        query_type="gene",cSEAdb)

# test cell
test_cell <- search_db(c("HL-60"),
                       query_type="cell",cSEAdb)

# test cancer
test_cancer<- search_db(query=c("Lung non small cell carcinoma"),
                       query_type="cancer",cSEAdb)


#-----------------------------------
# track plot
#-----------------------------------


# create a bigwig list
bigWigs  <- read.table("../bw_path.txt",sep="\t",header=F)


#test region "chr1:1156819-1178981","chr1_100574859_101152255"

################
#plot_region <- "chr5:1156819-1178981"
plot_region <- "chr1:1156819-1178981"
search_db_table <- search_db(plot_region,
                             query_type="se_region",cSEAdb)

usr_bigwig <- "inst/extdata/dummy_coverage.bw"
bw_gr <- create_bw_gr(plot_region,c(bigWigs$V1),se_spec_table=search_db_table)

# for (i in 1:60) {
#   export.bw(bw_gr$bw_gr[[i]],paste0("inst/extdata/example_bw_2/",names(bw_gr$bw_gr)[[i]],"_dummy.bw"))
# }



plot_out_2 <- create_gviz_tracks(bw_gr_list=bw_gr$bw_gr,
                               cell="specific",
                               se_spec_table=search_db_table,
                               plot_region=plot_region,
                               txdb=TxDb.Hsapiens.UCSC.hg38.knownGene)



pdf("test/test.pdf",width = 6,
    height = 10)
plotTracks(c(plot_out_2$htrack,
             plot_out_2$atrack_ce,
             plot_out_2$atrack_se,
             plot_out_2$atrack_query,
             plot_out_2$axistrack,
             plot_out_2$gtrack),
           size=5,
           cex.title = 0.5)
dev.off()
################




plot_out <- create_gviz_tracks(bw_gr_list=bw_gr$bw_gr,
                               cell="specific",
                               se_spec_table=search_db_table,
                               plot_region=plot_region,
                               txdb="none")


pdf("test_1.pdf",width = 6,
    height = 9)
plotTracks(c(plot_out$htrack,
             plot_out$atrack_ce,
             plot_out$atrack_se,
             plot_out$atrack_query,
             plot_out$axistrack))
dev.off()

pdf("test_1.pdf",width = 6,
    height = 10)
plotTracks(c(plot_out_2$htrack,
             plot_out_2$atrack_ce,
             plot_out_2$atrack_se,
             plot_out_2$atrack_query,
             plot_out_2$axistrack,
             plot_out_2$gtrack),
           cex.title = 0.4)
dev.off()




