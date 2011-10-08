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

# To run Bash commands
SH <- function(cmd){lines <- readLines(p <- pipe(cmd)); close(p); lines}

# Name of this script
script_path <- sub("--file=","", commandArgs()[grep("--file=", commandArgs())])
script_path_list <- strsplit(script_path, '/')[[1]]
script_name <- script_path_list[[length(script_path_list)]]

# Directory where this script is located
script_dir_list  <- script_path_list[1:(length(script_path_list)-1)]
script_dir  <- paste(script_dir_list, sep='/', collapse="/")

# Arguments
script_args <- commandArgs(trailingOnly=T)

if(length(script_args) != 1) {
       cat(paste(
         paste("Usage:", script_name, "<blast-results-dir>\n"),
         "Required Options",
         "  <blast-results-dir>  Directory with the Blast m8 files",
         "",
         sep="\n"))
       quit()
}

UPSTREAM_DIR <- script_args[[1]]


# Get all blast m8 files
file_list <- SH(paste("ls ", UPSTREAM_DIR, "/*-m8 |",
  " grep -v 'tab-new-'", # filter out new - those are incomplete/failed
  sep=""))

# Read tables
for (f in file_list) {

  print(f)
  blast_m8 <- read.table(f, sep="\t")
  names(blast_m8) <-
    c('query',
      'subject',
      'pid',
      'alignment_length',
      'mismatches',
      'gap_openings',
      'q_start',
      'q_end',
      's_start',
      's_end',
      'e_value',
      'bit_score')

  blast_m8$query <- as.character(blast_m8$query)
  blast_m8$subject <- as.character(blast_m8$subject)

  # Same as all_seqs_in_sample because of self hits 
  all_seqs_in_blast_results <- unique(c(blast_m8$query, blast_m8$subject))
  n <- length(all_seqs_in_blast_results)
  a2a_matrix <- matrix(rep(rep(1.0, n), n), nrow=n, ncol=n)

  rownames(a2a_matrix) <- all_seqs_in_blast_results
  colnames(a2a_matrix) <- all_seqs_in_blast_results

  apply(blast_m8, 1, function(x) {
    a2a_matrix[x["query"], x["subject"]] <<- as.numeric(x["e_value"])
  })

   
  #options(width=10000)
  #print(blast_m8[, c('query', 'subject', 'e_value')])
  #print(a2a_matrix)
  #cat("\n")
 
  # Get counts of e_value >= 0.01
  scores1 <- rowSums(a2a_matrix <= 0.01) + rowSums(t(a2a_matrix) <= 0.01)

  #scores2 <- 
  # apply(a2a_matrix, 1, function(x) length(which(x <= 0.01))) +
  # apply(t(a2a_matrix), 1, function(x) length(which(x <= 0.01)))

  #if(!all(scores1 == scores2)) {
  #  print(scores1)
  #  print(scores2)
  #  cat("Which:\n") 
  #  print(which(!(scores1 == scores2)))
  #  cat("Matrix:\n") 
  #  #rownames(a2a_matrix) <- NULL
  #  #colnames(a2a_matrix) <- NULL
  #  options(width=10000)
  #  print(a2a_matrix)
  #  cat("\n")

  #  ##a2a_matrix_t <- t(a2a_matrix)
  #  #  apply(1:nrows(a2a_matrix), 1,
  #  #    function(r) {
  #  #      smaller_in_col <- a2a_matrix[,r] > a2a_matrix[r,]
  #  #      a2a_matrix[smaller_in_col,] <<- a2a_matrix[,smaller_in_col]
  #  #    })

  #  #print(a2a_matrix)

  #  q()
  #}
  
  
  # Write seq_id + score tab files
  write.table(file=paste(f, "-scores", sep=""), scores1, sep="\t", col.names=F)


}
