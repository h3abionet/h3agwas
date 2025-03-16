#!/usr/bin/env Rscript
# 12.12.2018
# Daniel Schreyer
# Imputation quality score (r2-values) and SNP count plotted against MAF bin

#----------------------------------------------------------------
# Generates a scatterplot: x = MAF bin,
#                          y - axis left: SNPs count, y - axis right: mean r-squared
# Input: info files of one/multiple dataset/s
# e.g.: Rscript r2_freq_hist.R -i panelname==<path.to.info.file>,panel2==<path.to.info.2>,...
#       -o r2_freq_hist.png
#----------------------------------------------------------------


# required packages
library(ggplot2)
library(dplyr)
library(optparse)
library(data.table)
library(ggsci)

# Takes paths to Input .info file and Output .png file as options
option_list <- list(
  make_option(c("-i", "--infos"), action="store", default = "${infos}", type = 'character',
              help = "Imputation .info file of each reference panel"),
  make_option(c("-i", "--ref"), action="store", default = "${ref_panels}", type = 'character',
              help = "Imputation .info file of each reference panel"),
  make_option(c("-o", "--output"), action="store", default = "${plot_out}", type = 'character',
              help = "Output .png file")
)
args <- parse_args(OptionParser(option_list = option_list))

# read in info files of both reference panels
input <- as.character(args[1])
inputs <- unlist(strsplit(input,","))
inp <- as.character(args[2])
inputed <- unlist(strsplit(inp,","))

panels <- list()
for(n in inputed){
  idx_name <- which(inputed==n)  # Get the index of name in names
  file <- inputs[idx_name]
  panels[paste0(n)] <- file
}

# read in .info files of each reference panel and merge them together in one table
i <- 1
for(file in panels){
    name <- names(panels)[i]
    panel <- fread(as.character(file), sep = "\\t", header = T, select = c("SNP","MAF","Rsq","Genotyped"))
    panel <- panel %>% mutate(R_Panel = paste0(name))
    if(i > 1){
        full <- rbind(full, panel)
    }else{full <- panel}
    i <- i+1
}

# filter out non-imputed SNPs and missing Rsq values
Imputed <- filter(full, Genotyped == "Imputed")
Imputed <- filter(Imputed, Rsq != "-" | !is.na(Rsq) | MAF != "-" | is.na(MAF))

# calculate mean Rsq and frequency of both reference panel
# MAF are rounded to 2 decimal places <- bin
Imputed <- Imputed %>% mutate( MAF = round(as.numeric(MAF), 2)) %>% group_by(R_Panel, MAF) %>%
    summarise(Rsq_mean = mean(as.numeric(Rsq)), N = n())

#calculate maximum SNP count to define second y axis
max <- max(Imputed\$N)

# Plot Minor Allele Frequency(MAF) bin (x), r-squared (y-left) / SNPs count (y-right)
p <- ggplot(Imputed, aes(x = MAF, color = R_Panel)) +
    geom_point(aes(y = Rsq_mean*max, shape = "Rsq"), size = 1) +
    geom_line(aes(y = Rsq_mean*max)) +
    geom_point(aes(y = N, shape = "Frq"), size = 2) +
    geom_line(aes(y = N))+
    scale_y_continuous(sec.axis = sec_axis(~./max, name = "mean IQS per MAF bin",
    breaks = c(seq(0, 1, 0.2))), name = "SNPs count",
    labels = scales::comma) + labs(colour = "") +
    geom_point(aes(x = MAF, y = Rsq_mean*max, shape = "Rsq", color = R_Panel)
    , inherit.aes = F, size = 1) +
    geom_point(aes(y = N, shape = "Frq"), size = 2) +
    scale_color_npg(name = "Reference Panel") +
    scale_shape_manual(name = "", breaks = c("Rsq", "Frq"),
    labels = c("mean IQS per MAF bin","SNPs count" ), values = c(18,4)) +
    theme_bw()

# save plot as .png file
ggsave(file = as.character(args[3]) ,plot = p, width = 8, height = 5, units = "in")