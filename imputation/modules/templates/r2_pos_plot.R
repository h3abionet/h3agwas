#!/usr/bin/env Rscript
# 25.10.2018
# Daniel Schreyer
# R-squared - Position - Plot ####

# required packages
library(dplyr)
library(ggplot2)
library(optparse)
library(tidyr)
library(data.table)

# takes input files as arguments
option_list <- list(
make_option(c("-i", "--info"), action = "store", default = "${info}", type = 'character',
help = "Imputation .info file"),
make_option(c("-t", "--target"), action = "store", default = "${target}", type = 'character',
help = "Target .frq file"),
make_option(c("-o", "--output"), action = "store", default = "${output}", type = 'character',
help = "Output .png file"),
make_option(c("-m", "--maf"), action="store", default = "${maf_thresh}", type = 'character',
help = "Minor allele frequency threshold in %")
)

args <- parse_args(OptionParser(option_list = option_list))


maf_thresh <- 0
# Read in the allele frequency files
info <- fread(file = as.character(args[1]), sep = "\\t", header = T, select = c("SNP", "Rsq", "Genotyped", "MAF", "ALT_Frq"))
frq <- fread(file = as.character(args[2]), sep = "\\t", header = T, select= c("SNP","CHR","POS"))

# modify SNP ID of frq file
frq <- frq %>% separate("SNP", c("CHR", "Position", "REF", "ALT"), "_" )
frq <- frq %>% mutate(SNP = paste(CHR, POS, REF, ALT, sep = ":"))

# merge tables
full <- merge(frq, info, by = "SNP")

# change chromosome names
full\$CHR <- paste("Chromosome", full\$CHR, sep = " ")

# extract imputed SNPs and filter out SNPs with missing R-squared values
Imputed <- filter(full, Genotyped == "Imputed")
Imputed <- filter(Imputed, Rsq != "-" | !is.na(Rsq))

Genotyped <- filter(full, Genotyped == "Genotyped")

# display only  ~ 50,000 Imputed SNPs in the r2 - position plot
N <- ifelse(nrow(Imputed)> 50000, as.integer(nrow(Imputed)/50000), 1)
Imputed <- Imputed[seq(1, nrow(Imputed),N),]

# display only ~ 1,000 Genotyped SNPs as ticks
N2 <- ifelse(nrow(Genotyped)>1000, as.integer(nrow(Genotyped)/1000), 1)
Genotyped <- Genotyped[seq(1, nrow(Genotyped),N2),]

# filter out MAF below a given value
AF_thresh <- ifelse(!is.na(maf_thresh) & maf_thresh > 0, maf_thresh/100, 0)
# AF_thresh
Imputed <- filter(Imputed, MAF > AF_thresh & MAF != 0)
# Imputed

# categorize MAF levels: extreme rare, moderate rare, rare, moderate, common, extreme common
Imputed <- mutate(Imputed, MAF2 = MAF)

Imputed\$MAF2[Imputed\$MAF2 > 0 & Imputed\$MAF2 <= 0.001]<- "extreme rare (0,0.001]"
Imputed\$MAF2[Imputed\$MAF2 > 0.001 & Imputed\$MAF2 <= 0.01]<- "moderate rare (0.001,0.01]"
Imputed\$MAF2[Imputed\$MAF2 > 0.01 & Imputed\$MAF2 <= 0.02]<- "rare (0.01,0.02]"
Imputed\$MAF2[Imputed\$MAF2 > 0.02 & Imputed\$MAF2 <= 0.05]<- "moderate (0.02,0.05]"
Imputed\$MAF2[Imputed\$MAF2 > 0.05 & Imputed\$MAF2 <= 0.2]<- "common (0.05,0.2]"
Imputed\$MAF2[Imputed\$MAF2 > 0.2 & Imputed\$MAF2 <= 0.5]<- "extreme common (0.2,0.5]"

head(Imputed)
# plot rsquared vs. SNP_position
r2_position_plot <- ggplot(data = Imputed, aes(x = POS, y = Rsq)) +
    geom_point(aes(colour = MAF2)) +
    labs(x = "Position [bp]", y = "Imputation accuracy (r-squared)") +
    scale_y_discrete(breaks = seq(0, 1, 0.2)) + #scale_x_continuous(labels = scales::comma) +
    facet_grid(~CHR) + theme_classic() +
    geom_hline(aes(yintercept = 0.3, colour = "red"), show.legend = F) +
    scale_color_discrete(breaks = c("extreme rare (0,0.001]","moderate rare (0.001,0.01]","rare (0.01,0.02]",
    "moderate (0.02,0.05]","common (0.05,0.2]","extreme common (0.2,0.5]")) +
    geom_rug(data=Genotyped, aes(x = POS), inherit.aes = F)


# save plot as a .png file
ggsave(file = as.character(args[3]), r2_position_plot, width = 10, height = 5.5, units = "in")
