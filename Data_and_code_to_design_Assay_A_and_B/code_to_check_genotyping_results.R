## Verify that the norovirus GI sequences from GenBank actually are norovirus GI

## Import and check typing results
x <- read.table("genotyping_results.txt", sep = "\t", header = TRUE)
all(x$BLAST == "Caliciviridae Norovirus GI")
nrow(x)
