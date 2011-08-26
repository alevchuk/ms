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

if(length(script.args) != 3) {
       cat(paste(
         paste("Usage:", script.path, "<size> <source>\n"),
         "Required Options",
         "  <seed>     Seed. Any number (e.g. 1). For sampling reproducibly",
         "  <size>     Size of the disired sample",
         "  <source>   File where lines constitute the population",
         "",
         sep="\n"))
       quit()
}

seed <- script.args[[1]]
size <- script.args[[2]]
source_file <- script.args[[3]]

set.seed(seed)
writeLines(sample(readLines(source_file), size))
