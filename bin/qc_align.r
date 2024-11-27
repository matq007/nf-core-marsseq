#!/usr/bin/env Rscript
# Adapted source code from
# https://tanaylab.github.io/old_resources/pages/672.html

args = commandArgs(trailingOnly = TRUE)

if (length(args) == 1 && args[1] == "--version") {
    message("v1.0")
    quit()
} else if (length(args) == 2) {
    read_sam = args[1]
    read_qc = args[2]
} else {
    stop("usage: Rscript qc_align.r [read.sam] [read_qc.txt]")
}

flag_counts = table(read.table(
    pipe(paste("grep -v ^@ ", read_sam, " | cut -f2", sep = ""))
)[, 1])

stats = read.table(read_qc, header = T)
stats$mapped = flag_counts["0"] + flag_counts["16"]

write.table(file = paste0("_", read_qc),
    stats, sep = "\t", col.names = T, row.names = F, quote = F)
