#!/usr/bin/env Rscript

# Visualizes plate map with wells colored by total estimated contamination volume
# and with arrows to contaminated wells from proposed potential sources of contamination.

# plate map visualization based on
# https://rstudio-pubs-static.s3.amazonaws.com/427185_abcc25a951c9436680dc6a8fcc471ca9.html

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
input_file_path  <- args[1]
output_file_path <- args[2]
number_columns   <- as.numeric(args[3])
number_rows      <- as.numeric(args[4])

input_file_path = "/Users/lakras/2021_05_19_contamination_detection_from_minor_alleles/2021_06_01_polishing_software/intermediate0/plate_map.txt_potential_cross_contamination_2x3.txt"
output_file_path = input_file_path
number_columns = 12
number_rows = 8

library(ggplot2)
library(plyr)

WIDTH <- 7.5
HEIGHT <- 5

# plate layouts:
#   6   2×3
#   12  3×4
#   24  4×6
#   48  6×8
#   96  8×12
#   384 16×24
#   1536 32×48
#   3456 48×72

alphaToRowNum <- function(row)
{
  # convert from letter row index to 
  # test: 1 = A, 26 = Z, 27 = AA, 703 = AAA
  # from:
  #   https://stackoverflow.com/a/34537691/2292993
  
  # codes from 
  s = row
  # Uppercase
  s_upper <- toupper(s)
  # Convert string to a vector of single letters
  s_split <- unlist(strsplit(s_upper, split=""))
  # Convert each letter to the corresponding number
  s_number <- sapply(s_split, function(x) {which(LETTERS == x)})
  # Derive the numeric value associated with each letter
  numbers <- 26^((length(s_number)-1):0)
  # Calculate the row number
  row_number <- sum(s_number * numbers)
  return(row_number)
}

rowNumToAlpha <- function(row_num)
{
  # convers a number of rows into the corresponding letter format
  #  32 rows = 'AF'
  #  678     = 'ZB'
  #  729     = 'ABA')
  row_alpha<-""
  
  row_num<-row_num
  while(row_num > 0)
  {
    letters_into_alphabet <- (row_num-1)%%26
    row_alpha <- paste(intToUtf8(letters_into_alphabet + utf8ToInt('A')), row_alpha , sep="")
    row_num <- ((row_num-1)%/%26)
  }
  return(row_alpha)
}

plate_list <- function(max_row, num_columns)
{
  # accepts either letters of last row or number of rows
  #   Ex. for a 1536-well plate (32x48):
  #     plate_list("AF",48)
  #     plate_list(32,48)
  #
  # Can also be used to return the ID of a given
  # 1-indexed well number
  #   Ex. for a 96-well plate (8x12)
  #     wells <- plate_list("H",12)
  #     wells[28] # (slice gives "C4")
  
  if(is.character(max_row))
  {
    num_rows<-alphaToRowNum(max_row)
  }
  else
  {
    num_rows<-max_row
  }
  col_idx <- seq(1, num_columns, by=1)
  seq(from = 1, to = num_rows)
  row_idx<-unlist(lapply(seq(from = 1, to = num_rows),rowNumToAlpha))
  wells<-expand.grid(row=row_idx,col=col_idx)
  index_pairs<-wells[with(wells, order(row, col)), ]
  paste(index_pairs$row,index_pairs$col,sep="")
}

# scales circle sizes and padding to panel edges
scaling_factor <- min(12/number_columns, 8/number_rows)
well_circle_size <- 12 * scaling_factor

expand_x <- 0.06 * scaling_factor
expand_y <- 0.08 * scaling_factor
if(number_columns < number_rows)
{
  # swap
  expand_y <- 0.06 * scaling_factor
  expand_x <- 0.08 * scaling_factor
}

# reads in input table
input_table <- read.table(input_file_path, sep="\t", header=TRUE)

# expands input table to include all wells, including wells not included in input table
well <- plate_list(number_rows, number_columns)
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
    shape=21, size=well_circle_size, colour="black") +
  geom_segment(data=input_table, mapping=aes(x=Column0, y=Row0, xend=Column, yend=Row),
    arrow=arrow(type="open", angle=30)) +
  coord_fixed(ratio=1, expand=TRUE, clip="off") +
  scale_x_continuous(breaks=seq(1, number_columns), position = "top", expand=c(expand_x,expand_x)) +
  scale_y_reverse(breaks=seq(1, number_rows), labels=LETTERS[1:number_rows], expand=c(expand_y,expand_y)) +
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
