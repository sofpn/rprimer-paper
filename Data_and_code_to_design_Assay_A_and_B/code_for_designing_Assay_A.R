## Code to reproduce design of Assay A

library(rprimer)

## Import target alignment
aln <- Biostrings::readDNAMultipleAlignment("alignment.fasta", format = "fasta")

## Make consensus profile
consensus <- consensusProfile(aln, ambiguityThreshold = 0.1)

## Design  primers and probes
oligoCandidates <- designOligos(consensus,
                          lengthPrimer = 18:26,
                          maxDegeneracyPrimer = 2,
                          gcClampPrimer = FALSE,
                          lengthProbe = 16:26,
                          maxDegeneracyProbe = 2
)

## Design assays
assayCandidates <- designAssays(oligoCandidates, tmDifferencePrimers = 2)
nrow(assayCandidates)

# View score, and based on that, filter on score
boxplot(assayCandidates$score)
assayCandidates <- assayCandidates[assayCandidates$score < 10, ]
nrow(assayCandidates)

# View all candidates
#View(as.data.frame(assayCandidates))

# Make plots
plotData(
    consensus,
    highlight = c(min(assayCandidates$start), max(assayCandidates$end))
)
plotData(oligoCandidates)
plotData(assayCandidates)

## All 24 assays are basically the same, with a few nt difference
## Upon some further investigation (especially off target binding
## for the rev primer!), I chose the following one: (from Shiny app)

r <- "CGTCCTTAGACGCCATCATCATTTAC"
f <- "GCCATGTTCCGCTGGATG"
p <- "CGRTCTCCTGTCCACA"

assayCandidates <- assayCandidates[assayCandidates$iupacSequenceRev == r, ]
assayCandidates <- assayCandidates[assayCandidates$iupacSequenceFwd == f, ]
assayCandidates <- assayCandidates[assayCandidates$iupacSequenceRcPr == p, ]

finalAssay <- assayCandidates[1, ]

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

pr <- consensus[
    consensus$position >= finalAssay$startPr &
        consensus$position <= finalAssay$endPr,
]

plotData(fwd, type = "nucleotide")
plotData(rev, type = "nucleotide", rc = TRUE)
plotData(pr, type = "nucleotide", rc = TRUE)

## Check match
finalPrimerF <- oligoCandidates[
    oligoCandidates$iupacSequence == f, ]
fMatch <- checkMatch(finalPrimerF, aln)

finalPrimerR <- oligoCandidates[
    oligoCandidates$iupacSequenceRc == r, ]
rMatch <- checkMatch(finalPrimerR, aln)

finalProbe <- oligoCandidates[
    oligoCandidates$iupacSequenceRc == p, ]
pMatch <- checkMatch(finalProbe, aln)

plotData(fMatch[1, ])
plotData(pMatch[1, ])
plotData(rMatch[1, ])

