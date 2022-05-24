## Code to exclude sequences with long N stretches
library(Biostrings)

## Import target sequences
seq <- readDNAStringSet(
    "Norovirus_GI_sequences_from_GenBank.fasta", format = "fasta"
)

## Pattern to search for
pattern <- DNAString(paste(rep("N", 30), collapse = ""))

## Identify sequences to exclude
exclude <- vmatchPattern(pattern, seq)

exclude <- as.data.frame(exclude)

exclude <- unique(exclude$group)

exclude

## Subset sequences and save to file
seq <- seq[-exclude]

length(seq)

