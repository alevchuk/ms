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

source("opt/bio3d/R/read.fasta.R")
source("opt/bio3d/R/pairwise.R")
source("opt/bio3d/R/is.gap.R")
source("opt/bio3d/R/identity.R")
matrixFromList <- function(listX) t(sapply(listX, function(x, n) {
  c(x, rep(NA, n))[1:n]}, n = max(sapply(listX, length))))

#CAND_NUMSEQ_MIN <- 30 # Minimum # of sequences for an MSA to be considered
#                      # This adjusts for the fact that this fishing method
#                      # will results in alignments with reduced numseq
NUMSEQ_MIN <- 19 # Minimum # of sequences per output MSA
PID_HIT_MIN <- .50 # Minumum sequence % identity to be considered a hit
MIN_HITS <- 19 # Miniumun # of hits for a seqeunce to be approximated as
               # a member of a tightly knit protein group

MAX_AVR_SEQLEN <- 1000


script.path <- sub("--file=","", commandArgs()[grep("--file=", commandArgs())])
script.args <- commandArgs(trailingOnly=T)

if(length(script.args) != 2) {
       cat(paste(
         "Usage:", script.path, " <source-info> <source-dir>\n",
         "Required Options\n",
         "  <source-info>  List with columns: msa, numseq, numres, avrseqlen\n",
         "  <source-dir>   Directory where all MSAs are located\n",
         "\n",
         sep=""))
       quit()
}

SOURCE_INFO <- script.args[[1]]
SOURCE_DIR <- script.args[[2]]

all <- read.table(SOURCE_INFO)
colnames(all) <- c("msa", "numseq", "numres", "avrseqlen")

#all <- all[all[,"numseq"] >= CAND_NUMSEQ_MIN,]

for (x in 1:nrow(all)) {

  candidate <- as.character(all[x, "msa"])
  
  file_name_tab <- paste(SOURCE_DIR, "/", candidate, sep="")

  msa_table <- read.table(file_name_tab)
  msa <- as.character(msa_table$V2)
  msa <- matrixFromList(strsplit(msa, ""))
  rownames(msa) <- msa_table$V1


  # All-agains-all % identity matrix
  aaa_id <- identity(msa)

  # Step 1/2 (SELECT): Pull out the desired seqences
  hits_per_seq <- 
    apply(aaa_id, 1, function(x) length(which(x >= PID_HIT_MIN)) - 1)
    # -1 to avoid counting self
 
  # Which have at least 19 ID50 mates?
  hits <- which(hits_per_seq >= MIN_HITS)


  aaa_id_selected <- aaa_id[hits, hits]

  if (length(hits) >= NUMSEQ_MIN) {


    #print(aaa_id)
    #print(msa)
    #print(msa[hits,])
    #quit()

    msa <- msa[hits,] 


    # Step 2/2 (FILTER): Remove any non-id50 members
    n <- length(hits)
    triangle <- sapply(1:n, function(x) sapply(1:n, function(y) {x / y > 1}))
    misses_tri <- (aaa_id_selected < PID_HIT_MIN) & triangle
    bad_columns <- sapply(1:n, function(x) T %in% misses_tri[,x])

    hits <- which(!bad_columns)
    aaa_id_selected_and_filtered <- aaa_id_selected[hits, hits]

    #print(hits)
    #print(aaa_id)
    #print(misses_tri)
    #print(bad_columns)
    #print(aaa_id_selected_and_filtered)
    #print(msa)
    #quit()


    if (length(hits) >= NUMSEQ_MIN) {

      #print(aaa_id)
      #print(misses_tri)
      #print(bad_columns)
      #print(aaa_id_selected_and_filtered)
      #print(msa)
      #quit()



      # Determine the average sequence lenght in hits
      num_seq <- nrow(msa[hits,])
      num_residues <- sum(apply(msa[hits,], 1, function(row) {
        length(which(row != "-"))
      }))
      rounded_avr_seq_length <- round(num_residues / num_seq)

      # Don't grab hits that exceed MAX_AVR_SEQLEN because then
      # we will not be able to inject foregin sequences from Uniport
      if (rounded_avr_seq_length <= MAX_AVR_SEQLEN) {

        hit_seq <- apply(msa[hits,], 1, paste, sep="", collapse="")

        sapply(names(hit_seq), function(x) {
          # MSA_NAME	SEQ_NAME	SEQUENCE
          cat("ID50_", candidate, "\t", x, "\t", hit_seq[x], "\n",sep="") 
        })
      }


    }
  }

}
