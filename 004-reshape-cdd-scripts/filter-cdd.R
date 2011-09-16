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

script.path <- sub("--file=","", commandArgs()[grep("--file=", commandArgs())])
script.args <- commandArgs(trailingOnly=T)

if(length(script.args) != 7) {
       cat(paste(
         paste("Usage:", script.path,
         "<all_alignments_dir>", "<msa_sizes_numseq>", "<min_msa_residues>", 
         "<max_msa_residues>", "<min_msa_avr_seqlen>", "<max_msa_avr_seqlen>",
         "<min_msa_numseq>", "<max_msa_numseq>", sep=" "),
         "",
         "Required Options",
         "  <all_alignments_dir> directory will all MSAs in tab format",
         "  <msa_numseq_file>    file that lists all MSAs + their NumSeq size",
         "",
         "  <max_msa_residues>   maximum number of residues (non-gap) in MSA",
         "",
         "  <min_msa_avr_seqlen> minimum average sequence length in MSA",
         "  <max_msa_avr_seqlen> maximum -//-",
         "",
         "  <min_msa_numseq>     minimum number of sequences (NumSeq) in MSA",
         "  <max_msa_numseq>     maximum -//-",
         "", sep="\n"))
       quit()
}

ALL_ALIGNMENTS_DIR <- script.args[[1]]
MSA_SIZES_NUMSEQ_FILE <- script.args[[2]]

MIN_MSA_RESIDUES <- NA
MAX_MSA_RESIDUES <- as.numeric(script.args[[3]])

MIN_MSA_AVR_SEQLEN <- as.numeric(script.args[[4]])
MAX_MSA_AVR_SEQLEN <- as.numeric(script.args[[5]])

MIN_MSA_NUMSEQ <- as.numeric(script.args[[6]])
MAX_MSA_NUMSEQ <-  as.numeric(script.args[[7]])





cat("Counting residues in MSAs and selecting only MSAs that have:", "\n",
  "  * No LESS than", MIN_MSA_NUMSEQ, "sequences\n",
  "  * No LESS than", MIN_MSA_AVR_SEQLEN, "average sequence length\n",
  "  * No MORE than", MAX_MSA_NUMSEQ, "sequences\n",
  "  * No MORE than", MAX_MSA_AVR_SEQLEN, "average sequence length\n",
  "  * No MORE than", sprintf("%dK", round(MAX_MSA_RESIDUES / 1000)),
    "residues total\n",
  file=stderr())


msa_list <- read.delim(MSA_SIZES_NUMSEQ_FILE)
colnames(msa_list) <- c("NumSeq", "AlignmentFile")
msa_list <- msa_list[order(msa_list[,"NumSeq"]),]

# Filter #1
msa_list <-
  msa_list[msa_list[,"NumSeq"] >= MIN_MSA_NUMSEQ, ]

invisible(
  apply(msa_list, 1, function(x) {
    num_seq <- as.numeric(x[["NumSeq"]]) # apply casts msa_list to char :(
    alignmet_file <- x[["AlignmentFile"]]

    msa <- read.delim(paste(ALL_ALIGNMENTS_DIR, alignmet_file, sep="/"))
    colnames(msa) <- c("SeqID", "Sequence")

    num_residues <- sum(sapply(msa[,"Sequence"], function (x) {
      sum(strsplit(as.character(x), "")[[1]] != "-")
    }))

    rounded_avr_seq_length <- round(num_residues / num_seq)

    filter_it_out <- F

    # Filter #2
    if (rounded_avr_seq_length < MIN_MSA_AVR_SEQLEN)
      filter_it_out <- T

    # Filter #3
    if (rounded_avr_seq_length > MAX_MSA_AVR_SEQLEN)
      filter_it_out <- T

    # Filter #4
    if (num_residues > MAX_MSA_RESIDUES)
      filter_it_out <- T

    if (!filter_it_out)
      cat(alignmet_file, "\t", num_seq, "\t", num_residues, "\t",
        rounded_avr_seq_length, "\n")
  })
)

#length <- MSA_NUMRES / MSA_NUMSEQ
#length_in_filename <- sprintf("%04d", round(length))
#seq_file <-
#  paste(UNIPROT_PREFIX, "/", "length", length_in_filename, "-tab", sep="")
#population <- readLines(seq_file)
