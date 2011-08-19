#!/usr/bin/env Rscript

script.path <- sub("--file=","", commandArgs()[grep("--file=", commandArgs())])
script.args <- commandArgs(trailingOnly=T)

if(length(script.args) != 2) {
       cat(paste(
         paste("Usage:", script.path, "<size> <source>\n"),
         "Required Options",
         "  <size>     Size of the disired sample",
         "  <source>   File where lines constitute the population",
         "",
         sep="\n"))
       quit()
}

size <- script.args[[1]]
source_file <- script.args[[2]]

writeLines(sample(readLines(source_file), size))
