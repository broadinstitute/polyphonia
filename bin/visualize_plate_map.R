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
input_file_type  <- tolower(args[5]) # "isnvs" or "contamination" or "contamination_minor" or "contamination_consensus"


library(ggplot2)
library(plyr)

WIDTH <- 7.5
HEIGHT <- 5


# HELPER FUNCTIONS:

# helper function to convert from letter row index to number
# test: 1 = A, 26 = Z, 27 = AA, 703 = AAA
# from https://stackoverflow.com/a/34537691/2292993
row_letters_to_row_number <- function(row)
{
  # convert string to a vector of single uppercase letters
  row <- toupper(row)
  row_letters <- unlist(strsplit(row, split=""))
  
  # convert each letter to corresponding number
  row_numbers <- sapply(row_letters, function(x) {which(LETTERS == x)})
  
  # calculate row number
  numbers <- 26^((length(row_numbers)-1):0)
  row_number <- sum(row_numbers * numbers)
  return(row_number)
}

# helper function to convert row number to corresponding letter format
#  32 rows = 'AF'
#  678     = 'ZB'
#  729     = 'ABA'
row_number_to_row_letters <- function(row_number)
{
  row_letters <- ""
  while(row_number > 0)
  {
    letters_into_alphabet <- (row_number-1)%%26
    row_letters <- paste(intToUtf8(letters_into_alphabet + utf8ToInt('A')), row_letters , sep="")
    row_number <- ((row_number-1)%/%26)
  }
  return(row_letters)
}

# helper function to generate list of wells in plate
#
# accepts either letters of last row or number of rows
# ex. for a 1536-well plate (32x48):
#  plate_list("AF",48)
#  plate_list(32,48)
#
# can also be used to return the ID of a given 1-indexed well number
# ex. for a 96-well plate (8x12):
#  wells <- plate_list("H",12)
#  wells[28] # (slice gives "C4")
plate_list <- function(max_row, num_columns)
{
  if(is.character(max_row))
  {
    num_rows <- row_letters_to_row_number(max_row)
  }
  else
  {
    num_rows <- max_row
  }
  col_idx <- seq(1, num_columns, by=1)
  seq(from=1, to=num_rows)
  row_idx <- unlist(lapply(seq(from=1, to=num_rows), row_number_to_row_letters))
  wells <- expand.grid(row=row_idx, col=col_idx)
  index_pairs <- wells[with(wells, order(row, col)), ]
  paste(index_pairs$row, index_pairs$col, sep="")
}


# depending on plate size, scales circle size, axis label size, arrow head size
# and padding to panel edges
scaling_factor <- min(12/number_columns, 8/number_rows)
well_circle_size <- 12 * scaling_factor
well_circle_line_thickness <- min(0.3 * scaling_factor, 0.5)
arrow_head_length <- 0.63 * scaling_factor
arrow_thickness <- 2 * well_circle_line_thickness

text_in_wells_size <- 4 * scaling_factor
if(scaling_factor < 0.26)
{
  text_in_wells_size <- text_in_wells_size * 1.4 # make the within-well text larger
}

axis_text_size <- max(6, 9 * scaling_factor)
axis_text_n <- 1 # display every value in axis text
if(scaling_factor < 0.2)
{
  axis_text_n <- 2 # display every 2nd value in axis text
}

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

# reads in input table estimated contamination volume
input_table$estimated_contamination_volume <- as.numeric(sub("%", "", input_table$estimated_contamination_volume)) / 100

# orders options for appearance of potential contamination
input_table$appearance_of_potential_contamination <- factor(input_table$appearance_of_potential_contamination,
                                                            levels=c("minor alleles", "minor and consensus-level", "consensus-level"))

# subsets input file to display consensus-level potential contamination only
if(input_file_type == "contamination_consensus")
{
  input_table <- subset(input_table, appearance_of_potential_contamination == "consensus-level")
  input_file_type <- "contamination"
}

# subsets input file to display minor allele or minor-consensus mix potential contamination only
if(input_file_type == "contamination_minor")
{
  input_table <- subset(input_table, appearance_of_potential_contamination == "minor alleles" | appearance_of_potential_contamination == "minor and consensus-level")
  input_file_type <- "contamination"
}

# expands input table to include all wells, including wells not included in input table
well <- plate_list(number_rows, number_columns)
all_wells <- data.frame(well)
input_table$well <- input_table$potential_contaminated_sample_plate_position
plate_map <- merge(x=all_wells, y=input_table, by="well", all=TRUE)

# calculates sum of estimated contamination volume for each well
# (in case a well is potentially contaminated by multiple samples)
if(input_file_type == "contamination")
{
  volume_summed_by_well <- aggregate(x=plate_map$estimated_contamination_volume, by=list(plate_map$well), FUN=sum, na.rm=TRUE)
  colnames(volume_summed_by_well) <- c("well","estimated_contamination_volume_sum")
  volume_summed_by_well$estimated_contamination_volume_sum <- sapply(volume_summed_by_well$estimated_contamination_volume_sum, function(x) min(x,1)) # caps volume at 100%
  plate_map <- volume_summed_by_well
}

# splits out columns and rows of wells
# based on https://stackoverflow.com/questions/9756360/split-character-data-into-numbers-and-letters
plate_map$Row <- as.numeric(lapply(gsub("[[:digit:]]","",plate_map$well), FUN=row_letters_to_row_number)) # retrieves letters, converts to digits
plate_map$Column <- as.numeric(gsub("[^[:digit:]]", "", plate_map$well)) # retrieves digits

input_table$Row <- as.numeric(lapply(gsub("[[:digit:]]","",input_table$well), FUN=row_letters_to_row_number)) # retrieves letters, converts to digits
input_table$Column <- as.numeric(gsub("[^[:digit:]]", "", input_table$well)) # retrieves digits

if(input_file_type == "contamination")
{
  input_table$contamination_source_well <- input_table$potential_contaminating_sample_plate_position
  input_table$Row0 <- as.numeric(lapply(gsub("[[:digit:]]","",input_table$contamination_source_well), FUN=row_letters_to_row_number)) # retrieves letters, converts to digits
  input_table$Column0 <- as.numeric(gsub("[^[:digit:]]", "", input_table$contamination_source_well)) # retrieves digits
  
  # determines amount to jitter--more jitter if more lines of potential cross-contamination
  MIN_JITTER <- 0
  MAX_JITTER <- 0.25
  jitter_distance <- max(MIN_JITTER, MAX_JITTER * min(length(input_table$Row0) / (number_columns * number_rows), 1))
  
  # adds vertical and horizontal jitter values
  input_table$jitter_vertical <- sign(runif(length(input_table$Row0), min=-1, max=1)) * runif(length(input_table$Row0),min=MIN_JITTER,max=jitter_distance)
  input_table$jitter_horizontal <- sign(runif(length(input_table$Row0), min=-1, max=1)) * runif(length(input_table$Row0),min=MIN_JITTER,max=jitter_distance)
  
  # removes horizontal jitter for vertical lines; removes vertical jitter for horizontal lines
  input_table$jitter_horizontal[input_table$Row == input_table$Row0] <- 0
  input_table$jitter_vertical[input_table$Column == input_table$Column0] <- 0
}


# generates visualization of plate map (to overlay contamination or iSNVs)
plate_figure_base <- ggplot() +
  coord_fixed(ratio=1, expand=TRUE, clip="off") +
  scale_x_continuous(
    breaks=seq(1, number_columns, axis_text_n),
    position = "top",
    expand=c(expand_x, expand_x)) +
  scale_y_reverse(
    breaks=seq(1, number_rows, axis_text_n),
    labels=lapply(seq(1, number_rows, axis_text_n), FUN=row_number_to_row_letters),
    expand=c(expand_y, expand_y)) +
  xlab("") + ylab("") +
  theme(
    legend.background=element_blank(),
    legend.key=element_blank(),
    axis.ticks=element_blank(),
    axis.text=element_text(size=axis_text_size),
    legend.position="bottom",
    panel.background=element_blank(),
    panel.border=element_rect(colour="black", fill=NA, size=0.3),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  )


# generates contamination visualization
if(input_file_type == "contamination")
{
  # label for maximum contamination volume in legend
  # maximum_contamination_volume <- max(plate_map$estimated_contamination_volume_sum)
  # if(maximum_contamination_volume == 0)
  # {
  #   maximum_contamination_volume = 0.01
  # }
  # maximum_contamination_volume_text <- paste(signif(100*maximum_contamination_volume, digits=2), "%", sep="")
  minimum_contamination_volume <- 0
  maximum_contamination_volume <- 1
  minimum_contamination_volume_text <- "0%"
  maximum_contamination_volume_text <- "100%"
  
  # whether or not a well is involved in any detected potential contamination (giving or receiving)
  plate_map$well_involved <- FALSE
  for(i in 1:nrow(plate_map))
  {
    plate_map$well_involved[i] <- any(input_table$contamination_source_well==plate_map$well[i]) || any(input_table$well==plate_map$well[i])
  }
  
  # colors plate map by total potential contamination with arrows from
  # potential sources of contamination to potential contaminated wells
  plate_figure_contamination <- plate_figure_base +
    geom_point(
      data=plate_map,
      aes(x=Column, y=Row, fill=estimated_contamination_volume_sum),
      colour=ifelse(plate_map$well_involved, "black", "darkgrey"),
      shape=21, size=well_circle_size, stroke=well_circle_line_thickness) +
    geom_segment(
      data=subset(input_table, appearance_of_potential_contamination == "consensus-level"),
      mapping=aes(x=Column0+jitter_horizontal, y=Row0+jitter_vertical, xend=Column+jitter_horizontal, yend=Row+jitter_vertical,
        color=appearance_of_potential_contamination),
      arrow=arrow(type="open", angle=30, length=unit(arrow_head_length,"cm"))) +
    geom_segment(
      data=subset(input_table, appearance_of_potential_contamination == "minor and consensus-level"),
      mapping=aes(x=Column0+jitter_horizontal, y=Row0+jitter_vertical, xend=Column+jitter_horizontal, yend=Row+jitter_vertical,
                  color=appearance_of_potential_contamination),
      arrow=arrow(type="open", angle=30, length=unit(arrow_head_length,"cm"))) +
    geom_segment(
      data=subset(input_table, appearance_of_potential_contamination == "minor alleles"),
      mapping=aes(x=Column0+jitter_horizontal, y=Row0+jitter_vertical, xend=Column+jitter_horizontal, yend=Row+jitter_vertical,
                  color=appearance_of_potential_contamination),
      arrow=arrow(type="open", angle=30, length=unit(arrow_head_length,"cm"))) +
    guides(size=FALSE) +
    scale_size_manual(values = c("minor alleles"=arrow_thickness*0.6, "minor and consensus-level"=arrow_thickness*1, "consensus-level"=arrow_thickness*1.4)) +
    scale_fill_gradient("Total Estimated Contamination Volume", low="white", high="#CC857E",
      limits=c(minimum_contamination_volume, maximum_contamination_volume),
      breaks=c(minimum_contamination_volume, maximum_contamination_volume),
      labels=c(minimum_contamination_volume_text, maximum_contamination_volume_text)) +
    scale_color_manual("Appearance of Contamination",
                       values=c("minor alleles"="black", "minor and consensus-level"="darkgrey", "consensus-level"="lightgrey")) +
    theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())
  
  ggsave(paste(output_file_path, ".pdf", sep=""), plate_figure_contamination, width=WIDTH, height=HEIGHT)
  ggsave(paste(output_file_path, ".jpg", sep=""), plate_figure_contamination, width=WIDTH, height=HEIGHT)
}


# generates iSNVs visualization
if(input_file_type == "isnvs")
{
  # label for maximum contamination number iSNVs in legend
  maximum_iSNVs <- max(plate_map$iSNVs, na.rm=T)
  if(maximum_iSNVs == 0)
  {
    maximum_iSNVs = 1
  }
  
  # colors plate map by total number iSNVs
  plate_figure_iSNVs <- plate_figure_base +
    geom_point(
      data=plate_map,
      aes(x=Column, y=Row, fill=iSNVs),
      colour=ifelse(is.na(plate_map$iSNVs), "grey", "black"),
      shape=21, size=well_circle_size, stroke=well_circle_line_thickness) +
    geom_text(
      data=plate_map,
      aes(x=Column, y=Row, label=ifelse(is.na(iSNVs), "NA", iSNVs)),
      colour=ifelse(is.na(plate_map$iSNVs), "grey", "black"),
      size=text_in_wells_size) +
    scale_fill_gradient("Number iSNVs", low="white", high="#CC857E", na.value="#ededed",
      limits=c(0, maximum_iSNVs),
      breaks=c(0, maximum_iSNVs),
      labels=c(0, maximum_iSNVs))
  
  ggsave(paste(output_file_path, ".pdf", sep=""), plate_figure_iSNVs, width=WIDTH, height=HEIGHT)
  ggsave(paste(output_file_path, ".jpg", sep=""), plate_figure_iSNVs, width=WIDTH, height=HEIGHT)
}

# May 20, 2021
# June 9, 2021
# June 14, 2021
