#!/usr/bin/env Rscript

# UNNAMED is a classifier of outlier sequences in multiple alignments
# Copyright (C) 2011  Aleksandr Levchuk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

options(width=10000)

SH <- function(cmd){lines <- readLines(p <- pipe(cmd)); close(p); lines}

# Name of this script
script_path <- sub("--file=","", commandArgs()[grep("--file=", commandArgs())])
script_path_list <- strsplit(script_path, '/')[[1]]
script_name <- script_path_list[[length(script_path_list)]]

# Directory where this script is located
script_dir_list  <- script_path_list[1:(length(script_path_list)-1)]
script_dir  <- paste(script_dir_list, sep='/', collapse="/")


# Primary locations
EXPERIMENT_NAME <- SH("ls -1d ./trial-*/ | tail -1 | cut -d/ -f 2")
UPSTREAM <- paste(EXPERIMENT_NAME, "/029-collect-scores-data-out", sep="")
#OUT_DIR <- paste(EXPERIMENT_NAME, "/", script_name, "-data-out", sep="")

#dir.create(OUT_DIR, showWarnings=FALSE)

# How to use this script
usage_and_quit <- function() {
  cat(paste(script_name, " collects accounting information", "\n\n"))
  cat(paste("Usage:", script_name, "<sample-size>", "\n\n"))
  cat(paste(
    "Required Arguments",
    "  <sample-size> ",
    "",
    sep="\n"))
  quit()
}

# Check number of arguments
script_args <- commandArgs(trailingOnly=T)
if(length(script_args) == 0) usage_and_quit()
if(length(script_args) != 1) {
  cat("ERROR: Invalid number of arguments", length(script_args), "\n")
  usage_and_quit()
}

# Get arguments
SAMPLE_SIZE  <- script_args[[1]]


# Porject conventions
SAMPLE_FILES <- SH(paste(
  "ls ", UPSTREAM, "/sample-size", SAMPLE_SIZE, "*-seq-scores", sep=""))


reshape_data <- function(scores_df, rand) {
  # Mark +/- based on the sequence name
  n <- colnames(scores_df)
  scores_df <- cbind(scores_df,
    factor("Positive", levels=c("Positive", "Negative")))
  colnames(scores_df) <- c(n, "Actual")
  if (rand > 0) {
    locations_of_negatives <- grep("_RAND[0-9]", scores_df[,'seqid'])
    scores_df[locations_of_negatives, "Actual"] <- "Negative"
  }

  scores_df
}



print_accounting <- function(scores_df, info) {
    cat(nrow(scores_df), " ",
      info$sample_id, " ",
      sprintf("%.2f", info$injection_type), " ",
      info$method,
      "\n",
    sep="")
}


cat("\t\t\t0%\t5%\t10%\t15%\t20%\t25%\n")

for(sample_file in SAMPLE_FILES) {
  sample_id <- SH(sprintf(
   "ruby -e '\"%s\" =~ /sample-size[0-9]+-id([0-9]+)/; puts $1'",
   sample_file))

  cat(sample_id, "\n", sep="")

  scores_df <- read.table(sample_file)
  colnames(scores_df) <- c('scr', 'rand', 'unused', 'msa', 'method', 'seqid')
  injection_types <- sort(unique(scores_df[,'rand']))
  scoring_methods <- sort(as.character(unique(scores_df[,'method'])))


  for(scoring_method in scoring_methods) {
    local({
      scores_df <- scores_df[scores_df[,'method'] == scoring_method,]

      cat(scoring_method, "\t", sep="")
      for(injection_type in injection_types) {
        local({
          scores_df <- scores_df[scores_df[,'rand'] == injection_type,]

          scores_df <- reshape_data(scores_df, injection_type)

          cat(nrow(scores_df), "\t", sep="")

          #cat("\n")
          #print(scores_df[1:10,])

        })
      }
      cat("\n")

    })
  }
}



