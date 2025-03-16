#!/usr/bin/env Rscript
# 25.10.2018
# Daniel Schreyer
# SNP counts vs Mean Imputation Quality Score ####

# required packages
library(ggplot2)
library(dplyr)
library(optparse)
library(data.table)
library(ggsci)

# takes input files as arguments
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
panels

# read in .info files of each reference panel and merge them together in one table
i <- 1
for(file in panels){
  name <- inputed[i]
  panel <- fread(as.character(file), sep = "\t", header = T, select = c("SNP", "MAF","Rsq","Genotyped"),
                 stringsAsFactors=F)
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
Imputed <- Imputed %>% mutate( Rsq = round(as.numeric(Rsq), 2)) %>% group_by(R_Panel, Rsq) %>% summarise(Rsq_mean = mean(Rsq), N = n())

#### plot frequency vs r2 ####
r2_frequency_plot <- ggplot(Imputed, aes(x = Rsq_mean, y = N, color = R_Panel)) + 
  geom_point() +
  geom_line() +
  theme_classic() + labs(x = "Mean Imputation Quality Score", y = "SNP Count") + 
  scale_x_continuous(breaks = seq(0, 1, 0.2)) + scale_y_continuous(labels = scales::comma) +
  geom_vline(aes(xintercept = 0.3),colour = "red", show.legend = F) + theme(axis.line = element_line(size = 0.8)) +
  scale_color_npg(name = "Reference Panel") +
  theme(legend.position = "bottom")

# save plot as .png file 
ggsave(filename = as.character(args[3]) ,plot = r2_frequency_plot, width = 8, height = 5, units = "in")
