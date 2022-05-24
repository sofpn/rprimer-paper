## Code to reproduce design of Assay B

library(rprimer)

## Import target alignment and place relevant masks
origAln <- Biostrings::readDNAMultipleAlignment("alignment.fasta", format = "fasta")
aln <- origAln
Biostrings::colmask(aln, invert = TRUE) <- c(3800:5300, 5500:7000)

## Make consensus profile
consensus <- consensusProfile(aln, ambiguityThreshold = 0.01)
plotData(consensus)

## Design  primers
primerCandidates <- designOligos(consensus,
                           lengthPrimer = 25:40,
                           tmPrimer = c(60, 75),
                           gcClampPrimer = FALSE,
                           designStrategyPrimer = "mixed",
                           probe = FALSE
)

## Design assays
assayCandidates <- designAssays(
    primerCandidates, length = c(300, 1000), tmDifferencePrimers = 2
)

# I chose the following one from the Shiny app
r <- "CATGTTGCCAACCCAACCRTTRTACA"
f <- "CTTCACAGGTGAACAGCATAAAYCAYTGG"

finalPrimerF <- primerCandidates[
    primerCandidates$iupacSequence == f, ]

finalPrimerR <- primerCandidates[
    primerCandidates$iupacSequenceRc == r, ]

finalAssay <- assayCandidates[
    assayCandidates$iupacSequenceFwd == f,
]

finalAssay <- finalAssay[
    finalAssay$iupacSequenceRev == r,
]

# Indicate amplicon region
plotData(consensus, highlight = c(finalAssay$start, finalAssay$end))

# View oligo binding regions

fwd <- consensus[
    consensus$position >= finalAssay$startFwd &
        consensus$position <= finalAssay$endFwd,
]

rev <- consensus[
    consensus$position >= finalAssay$startRev &
        consensus$position <= finalAssay$endRev,
]

plotData(fwd, type = "nucleotide")
plotData(rev, type = "nucleotide", rc = TRUE)

# Check match
fMatch <- checkMatch(finalPrimerF, aln)

plotData(fMatch)

rMatch <- checkMatch(finalPrimerR, aln)

plotData(rMatch)
