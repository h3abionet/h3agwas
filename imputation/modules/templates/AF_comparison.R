#!/usr/bin/env Rscript
# 12.12.2018
# Daniel Schreyer

#----------------------------------------------------------------
# Correlation between ref allele freq and target allele freq Plot
# Takes Imputation info file and target file as an Input and outputs two differently colored scatter plots
# one plot is colored based on the SNP r-squared values 
# and one based on the difference between the allele frequencies
#----------------------------------------------------------------

# required packages
library(ggplot2)
library(dplyr)
library(optparse)
library(tidyr)

# takes input/output filepaths, r-squared threshold and the number of SNPs in the graph as an option
option_list <- list(
make_option(c("-i", "--info"), action="store", default= "${info}", type='character',
              help = "Imputation .info file"),
make_option(c("-t", "--target"), action="store", default = "${target}", type = 'character',
              help = "Target .frg file"),
make_option(c("-f", "--frq"), action = "store", default = "${frq}", type = 'character',
              help = "Reference Panel .frq file"),
  make_option(c("-r", "--rsq"), action="store", default = 0, type = 'character',
              help = "Filter files R-squared threshold "),
make_option(c("-o", "--outputcolor"), action = "store", default = "${outputcolor}", type = 'character',
              help = "Output Plot 2: SNP color based on r-squared values"),
make_option(c("-s", "--subset"), action = "store", default = 20000, type = 'integer',
help = "Display [] number of SNPs [Default = 20000]")
)
# make_option(c("-o", "--output"), action="store", default = "output", type = 'character',
#               help = "Output Plot: SNP color based on the ref AF and target AF difference"),
args <- parse_args(OptionParser(option_list = option_list))

# read in info and frequency file of imputed sample
info <- read.table(args\$i, header = T, sep ="\\t")
frq <- read.table(args\$t, header = T, sep = "\\t")

# modify SNP ID 
frq <- frq %>% separate("SNP", c("CHR", "Position", "REF", "ALT"), "_")
frq <- frq %>% mutate(SNP = paste(frq\$CHR, frq\$POS, frq\$REF, frq\$ALT, sep = ":")) %>%
  select(c("CHR","POS","SNP"))

# merge tables together and read in frequency file of reference panel
full <- merge(frq, dplyr::select(info, c("SNP","ALT_Frq", "Rsq", "Genotyped")),
              by = "SNP")
full <- merge(full, read.table(args\$f , sep = "\\t", header = T),
by = c("CHR","POS"))

# filter rsquared SNPs above a given threshold and calculate the frequency difference
rsq.thresh <- ifelse(!is.na(args\$r), args\$r, 0)
imputed <- full%>% filter(Genotyped == "Imputed") %>% 
  mutate(diff = abs(ALT_Frq-AF))
imputed <- filter(imputed, Rsq != "-" | !is.na(Rsq)) #%>% 
  #filter(Rsq > rsq.thresh)

# filter every Nth SNP to extract [] SNPs [Default = 20000]
n <- ifelse(nrow(imputed)> args\$s, as.integer(nrow(imputed)/args\$s), 1)
imputed <- imputed[seq(1, nrow(imputed),n),]


# plot the MAF against each other and color it based on the r-squared values
plot.rsq.colored <- ggplot(imputed , aes(x = ALT_Frq, y = AF, color = Rsq)) +
  geom_point(size = 0.9) + theme_classic() + 
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) + 
  labs(x = "RefAllele Frequency (Uploaded Samples)", 
       y = "Ref Allele Frequency (Reference Panel)") +
  geom_text(color = "black", x = 0.2, y = 1.02,
  label = paste("R-squared threshold:", rsq.thresh, sep = " ")) +
  geom_text(color = "black", x = 0.2, y = 0.95, 
            label = paste(nrow(imputed),"SNPs", sep = " ")) #+
  # scale_color_gradient(low = "lightblue", high = "darkblue")


# generate the plot and color the SNPs based on their ref AF and target AF difference
# SNP color for AF difference greater than 0.15 is black 
# plot.diff.colored <- ggplot(filter(imputed, diff <= 0.15),
#                        aes(x = ALT_Frq, y = AF, color = diff)) +
#   geom_point(size = 0.9) + theme_classic() +
#   scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
#   scale_color_gradient(low = "darkblue", high = "lightblue") +
#   scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
#   labs(x = "Ref Allele Frequency (Uploaded Samples)",
#        y = "Ref Allele Frequency (Reference Panel)") +
#   theme(legend.position="none") +
#   geom_text(color = "black", x = 0.2, y = 1.02,
#   label = paste("R-squared threshold:", rsq.thresh, sep = " ")) +
#   geom_text(color = "black", x = 0.2, y = 0.95,
#             label = paste(nrow(imputed), "SNPs", sep = " ")) +
#   geom_point(data = filter(imputed, diff > 0.15), aes(x = ALT_Frq, y = AF),
#              shape = 1, color = "black", size = 0.6)

# save both plots
ggsave(filename = as.character(args[5]), plot = plot.rsq.colored, width = 7, height = 7, units = "in")
# ggsave(filename = as.character(args[5]), plot = plot.diff.colored, width = 7, height = 7, units = "in")
