setwd("~/Projects/super_enhancer/se_data_portal/cSEAdb/")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(data.table)
library(RColorBrewer)
library(chromoMap)

cSEAdb <- readRDS("data/cSEAdb.rds")

# chrom background
chrom_size <- read.table("data/hg38_chrom_loc_w_centromere.txt",sep="\t",header=F)

#------------------------
# all specific SE
#------------------------
cell_sepc <- cSEAdb$se_cancer_specificity[,c(1,2,5)]

cell_sepc$group <- NA
cell_sepc$group[which(cell_sepc$cell_specific_gain=="yes" &
                        cell_sepc$cell_specific_loss=="no")] <- "gain"
cell_sepc$group[which(cell_sepc$cell_specific_gain=="no" &
                        cell_sepc$cell_specific_loss=="yes")] <- "loss"
cell_sepc$group[which(cell_sepc$cell_specific_gain=="yes" &
                        cell_sepc$cell_specific_loss=="yes")] <- "both"
cell_sepc$group[which(cell_sepc$cell_specific_gain=="no" &
                        cell_sepc$cell_specific_loss=="no")] <- "none"

# make se_df for plot

se_spec_df <- separate(data = cell_sepc[,-c(2,3)],se_name, into=c("V1","V2","V3"), sep='_',remove = F)

chromoMap(list(chrom_size),list(se_spec_df),
          #n_win.factor = 2,
          # chr_width = 20,
          chr_length = 6,
          chr_color = "white",
          data_based_color_map = T,
          data_type = "categorical",
          legend = T,
          lg_x = 200,
          lg_y = 500,
          text_font_size = c(16),
          data_colors = list(c("firebrick3","green3","gold2","lightskyblue")))







#----------------------
# make SE annotation file
#----------------------
# cell specific example:MCF7

cell_sepc <- cSEAdb$se_cell_specificity

# extract all MCF7 spec SE
se_spec_gain <- cell_sepc$se_name[which(cell_sepc$gain_ce %like% "MCF7")]
se_spec_loss <- cell_sepc$se_name[which(cell_sepc$loss_ce %like% "MCF7")]

# se have gain and loss
se_gain_loss <- intersect(se_spec_gain,se_spec_loss)
se_gain_only <- se_spec_gain[!(se_spec_gain %in% se_gain_loss)]
se_loss_only <- se_spec_loss[!(se_spec_loss %in% se_gain_loss)]

# make se_df for plot
se_spec_df <- data.frame(se_name=c(se_gain_only,se_loss_only,se_gain_loss),
                         group=c(rep("gain",length(se_gain_only)),
                                 rep("loss",length(se_loss_only)),
                                 rep("both",length(se_gain_loss))))

se_spec_df_2 <- separate(data = se_spec_df,se_name, into=c("V1","V2","V3"), sep='_',remove = F)

chromoMap(list(chrom_size),list(se_spec_df_2),
          n_win.factor = 2,
          # chr_width = 20,
          # chr_length = 5,
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("firebrick3","lightskyblue","green3")))




# reformat
se_df <- se_df[,c(4,1,2,3)]

se_test <- head(se_df)

chromoMap(list(chrom_size),list(se_df),
          n_win.factor = 3)


chromoMap(list(chrom_size),list(se_test))



