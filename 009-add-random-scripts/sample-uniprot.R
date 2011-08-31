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

if(length(script.args) != 5) {
       cat(paste(
         paste("Usage:", script.path,
         "<src-prefix> <seed> <num-seq> <num-res> <rand-ratio>\n"),
         "Required Options",
         "  <src-prefix>  Place of the Uniprot sequences, binned by length",
         "  <seed>        For the Random number generator",
         "  <num-seq>     How many sequences are in the original MSA?",
         "  <num-res>     How many residues are in the original MSA?",
         "  <rand-ratio>  How much randomness would you like (0.05, 0.1, ..)?",
         "",
         sep="\n"))
       quit()
}

UNIPROT_PREFIX <- script.args[[1]]
SEED <- as.numeric(script.args[[2]])
NUM_SEQ <- as.numeric(script.args[[3]])
NUM_RES <- as.numeric(script.args[[4]])
RAND_RATIO <- as.numeric(script.args[[5]])

formula_a <- function(n, r) {
 # n is num-seq in the original MSA
 # r is the desired random ratio

 # This formula finds the number of random sequences needed to achive the
 # actual random ratio A that is bound to the desired r in the following way:
 #   r - 0.05 < A <= r

 # PRCONDITION: random sequences always have the length of the average in MSA

 (r * n) / (1 - r) # Derived from: r - 0.05 < x / (n + x) <= r
}


rounded_avr_seq_length <- round(NUM_RES / NUM_SEQ)
sample_size <- formula_a(NUM_SEQ, RAND_RATIO)

length_in_filename <- sprintf("%04d", rounded_avr_seq_length)
seq_file <-
  paste(UNIPROT_PREFIX, "/", "length", length_in_filename, "-tab", sep="")
population <- readLines(seq_file)

set.seed(SEED)
writeLines(sample(population, sample_size))
