#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2017/10/04

library(plyr)
library(ggplot2)
library(reshape2)

superpop.plot.colours <- c("#8A2BE2", "#CC0066", "#6495ED", "#66CDAA", "#A52A2A", "#CDAA7D", "#66CD00", "#7AC5CD", "#CD5B45",
"#CDC8B1", "#00CDCD", "#CD950C", "#8B7500", "#800000 ", "#808000", "#008000", "#800080", "#008080", "#000080")
linetypes <- c("dashed", "solid", "dotted", "dotdash", "longdash", "twodash", "twodash")

data <- read.table("${report}", h = T, sep = '\\t')

colnames(data) <- c("${group}", "(0,0.001]", "(0.001,0.01]", "(0.01,0.02]", "(0.02,0.05]", "(0.05,0.2]", "(0.2,0.5]")
mdata <- melt(data, id = c("${group}"))
mdata\$"${group}" = as.factor(mdata\$"${group}")

p <- ggplot(mdata, aes(x = variable, y = value, fill = $group))+
geom_bar(stat="identity", position=position_dodge()) +
scale_fill_manual(values=superpop.plot.colours) +
ylab("${ylab}") + xlab("${xlab}") +
theme_bw()
ggsave("${plot_by_maf}", height=6, width=10, units='in', dpi=150)
