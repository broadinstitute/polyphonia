# Visualizes plate map with wells colored by total estimated contamination volume
# and with arrows to contaminated wells from proposed potential sources of contamination.

# plate map visualization based on
# https://rstudio-pubs-static.s3.amazonaws.com/427185_abcc25a951c9436680dc6a8fcc471ca9.html

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
input_file_path = args[1]
output_file_path = args[2]

library(ggplot2)
library(plyr)

WIDTH <- 7.5
HEIGHT <- 5 #4.25

WELL_CIRCLE_SIZE <- 12

# reads in input table
input_table <- read.table(input_file_path, sep="\t", header=TRUE)

# expands input table to include all wells, including wells not included in input table
well = c(
  "A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","A11","A12",
  "B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12",
  "C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12",
  "D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12",
  "E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12",
  "F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12",
  "G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12",
  "H1","H2","H3","H4","H5","H6","H7","H8","H9","H10","H11","H12")
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
plate_figure <-
ggplot() +
  geom_point(data=plate_map, aes(x=Column, y=Row, fill=estimated_contamination_volume_sum),
    shape=21, size=WELL_CIRCLE_SIZE, colour="black") +
  geom_segment(data=input_table, mapping=aes(x=Column0, y=Row0, xend=Column, yend=Row), arrow=arrow()) +
  coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.5, 12.5), ylim=c(0.5, 8.5)) +
  scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
  scale_x_continuous(breaks=seq(1, 12)) +
  xlab("") + ylab("") +
  theme(
    legend.background=element_blank(),
    legend.key = element_blank(),
    axis.ticks = element_blank(),
    legend.position="bottom",
    panel.background=element_blank(),
    panel.border = element_rect(colour="black", fill=NA, size=0.3)
  ) +
  scale_fill_gradient("Total Estimated Contamination Volume", low = "white", high = "#CC857E",
    breaks=c(0, maximum_contamination_volume), labels=c(0, maximum_contamination_volume_text))

ggsave(paste(output_file_path, ".pdf", sep=""), plate_figure, width=WIDTH, height=HEIGHT)
ggsave(paste(output_file_path, ".jpg", sep=""), plate_figure, width=WIDTH, height=HEIGHT)

# May 20, 2021
# June 9, 2021
