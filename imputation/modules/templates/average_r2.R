#!/usr/bin/env Rscript
library(dplyr)
library(data.table)

info_file <- fread("${file_infos}", sep="\\t", header=TRUE, data.table=FALSE)
input <- info_file %>%
  filter(!Rsq=="-") %>% select(Rsq) %>% 
  mutate(Rsq=as.numeric(Rsq)) %>% pull(Rsq)
average <- mean(input, na.rm = TRUE)
# Rsquared <- info_file\$Rsq
# input <- as.numeric(Rsquared)
# average <- mean(input, na.rm = TRUE)

write.table(average, file="${meanr2_out}", quote = FALSE, row.names = FALSE)
