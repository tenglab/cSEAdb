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

# load database
cSEAdb <- readRDS("data/cSEAdb.rds")

#-----------------------------------
# test search_db
#-----------------------------------

# test SE_region
test_se_r <- search_db(c("chr1_100584859_100673000","chr1_112388326_112392055"),
                     query_type="se_region",cSEAdb)

# test SE bed
test_se_b <- search_db("inst/extdata/se_query_example.bed",
                     query_type="se_bed",cSEAdb)

# test cell
test_cell <- search_db(c("HL-60","MCF7"),
                       query_type="cell",cSEAdb)

# test cancer
test_cancer<- search_db( c("Lung_non_small_cell_carcinoma","Melanoma"),
                       query_type="cancer",cSEAdb)

# test tissue
test_tissue<- search_db(c("Ovarian","Prostate","Non-Small Cell Lung"),
                        query_type="tissue",cSEAdb)

#-----------------------------------
# test plot function
#-----------------------------------
#



