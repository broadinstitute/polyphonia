#!/usr/bin/env Rscript

# Visualizes plate map with wells colored by total estimated contamination volume
# and with arrows to contaminated wells from proposed potential sources of contamination.

# plate map visualization based on
# https://rstudio-pubs-static.s3.amazonaws.com/427185_abcc25a951c9436680dc6a8fcc471ca9.html

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
input_file_path  <- args[1]
output_file_path <- args[2]
number_rows      <- as.numeric(args[3])
number_columns   <- as.numeric(args[4])

library(ggplot2)
library(plyr)

WIDTH <- 7.5
HEIGHT <- 5

scaling_factor <- min(12/number_columns, 8/number_rows)
WELL_CIRCLE_SIZE <- 12 * scaling_factor
EXPAND_X <- 0.06 * scaling_factor
EXPAND_Y <- 0.08 * scaling_factor

# reads in input table
input_table <- read.table(input_file_path, sep="\t", header=TRUE)

# expands input table to include all wells, including wells not included in input table
letters <- LETTERS[1:number_rows]
numbers <- c(1:number_columns)
well <- paste(rep(letters, each = length(numbers)), numbers, sep = "")
all_wells <- data.frame(well)
input_table_all_wells <- merge(x=all_wells, y=input_table, by="well", all=TRUE)

# calculates sum of estimated contamination volume for each well
# (in case a well is potentially contaminated by multiple samples)
volume_summed_by_well <- aggregate(x=input_table_all_wells$estimated_contamination_volume,
  by=list(input_table_all_wells$well), FUN=sum, na.rm=TRUE)
colnames(volume_summed_by_well) <- c("well","estimated_contamination_volume_sum")

# splits out columns and rows of wells
plate_map <- mutate(volume_summed_by_well,
  Row=as.numeric(match(toupper(substr(well, 1, 1)), LETTERS)),
  Column=as.numeric(substr(well, 2, 5)))
input_table <- mutate(input_table,
  Row=as.numeric(match(toupper(substr(well, 1, 1)), LETTERS)),
  Column=as.numeric(substr(well, 2, 5)))
input_table <- mutate(input_table,
  Row0=as.numeric(match(toupper(substr(contamination_source_well, 1, 1)), LETTERS)),
  Column0=as.numeric(substr(contamination_source_well, 2, 5)))

# label for maximum contamination volume in legend
maximum_contamination_volume <- max(plate_map$estimated_contamination_volume_sum)
maximum_contamination_volume_text <- paste(100*maximum_contamination_volume, "%", sep="")

# generates figures
plate_figure <- ggplot() +
  geom_point(data=plate_map, aes(x=Column, y=Row, fill=estimated_contamination_volume_sum),
    shape=21, size=WELL_CIRCLE_SIZE, colour="black") +
  geom_segment(data=input_table, mapping=aes(x=Column0, y=Row0, xend=Column, yend=Row),
    arrow=arrow(type="open", angle=30)) +
  coord_fixed(ratio=1, expand=TRUE, clip="off") +
  scale_y_reverse(breaks=seq(1, number_rows), labels=LETTERS[1:number_rows], expand=c(EXPAND_Y,EXPAND_Y)) +
  scale_x_continuous(breaks=seq(1, number_columns), position = "top", expand=c(EXPAND_X,EXPAND_X)) +
  xlab("") + ylab("") +
  theme(
    legend.background=element_blank(),
    legend.key = element_blank(),
    axis.ticks = element_blank(),
    legend.position="bottom",
    panel.background=element_blank(),
    panel.border = element_rect(colour="black", fill=NA, size=0.3),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_gradient("Total Estimated Contamination Volume", low = "white", high = "#CC857E",
    breaks=c(0, maximum_contamination_volume), labels=c(0, maximum_contamination_volume_text))

ggsave(paste(output_file_path, ".pdf", sep=""), plate_figure, width=WIDTH, height=HEIGHT)
ggsave(paste(output_file_path, ".jpg", sep=""), plate_figure, width=WIDTH, height=HEIGHT)

# May 20, 2021
# June 9, 2021
# June 14, 2021
