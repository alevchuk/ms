#!/usr/bin/env Rscript

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
