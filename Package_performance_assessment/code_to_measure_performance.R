## Code for measuring package performance
library(rprimer)
library(microbenchmark)
library(ggplot2)

## Custom functions ============================================================

makeMicrobenchmarkAssayA <- function(aln, times = 20) {

    ## Generate data
    consensus <- consensusProfile(aln, ambiguityThreshold = 0.1)
    oligoCandidates <- designOligos(consensus,
                                    lengthPrimer = 18:26,
                                    maxDegeneracyPrimer = 2,
                                    gcClampPrimer = FALSE,
                                    lengthProbe = 16:26,
                                    maxDegeneracyProbe = 2
    )

    ## Measure time
    timeConsensus <- microbenchmark(consensusProfile(
        aln, ambiguityThreshold = 0.1), unit = "s", times = times
    )

    timeOligos <-  microbenchmark(designOligos(consensus,
                                               lengthPrimer = 18:26,
                                               maxDegeneracyPrimer = 2,
                                               gcClampPrimer = FALSE,
                                               lengthProbe = 16:26,
                                               maxDegeneracyProbe = 2
    ), unit = "s", times = times)

    timeAssays <-  microbenchmark(designAssays(
        oligoCandidates, tmDifferencePrimers = 2
    ), unit = "s", times = times)

    timeCheckMatch <- microbenchmark(checkMatch(
        oligoCandidates, aln
    ), unit = "s", times = times)

    list(
        "timeConsensus" = timeConsensus,
        "timeOligos" = timeOligos,
        "timeAssays" = timeAssays,
        "timeCheckMatch" = timeCheckMatch
    )
}

makeMicrobenchmarkAssayB <- function(aln, times = 20) {

    ## Generate data
    consensus <- consensusProfile(aln, ambiguityThreshold = 0.01)
    primerCandidates <- designOligos(consensus,
                                     lengthPrimer = 25:40,
                                     tmPrimer = c(60, 75),
                                     gcClampPrimer = FALSE,
                                     designStrategyPrimer = "mixed",
                                     probe = FALSE
    )

    ## Measure time

    timeConsensus <- microbenchmark(consensusProfile(
        aln, ambiguityThreshold = 0.01), unit = "s", times = times
    )

    timeOligos <-  microbenchmark(designOligos(consensus,
                                               lengthPrimer = 25:40,
                                               tmPrimer = c(60, 75),
                                               gcClampPrimer = FALSE,
                                               designStrategyPrimer = "mixed",
                                               probe = FALSE
    ), unit = "s", times = times)

    timeAssays <-  microbenchmark(designAssays(
        primerCandidates, length = c(300, 1000), tmDifferencePrimers = 2
    ), unit = "s", times = times)

    timeCheckMatch <- microbenchmark(checkMatch(
        primerCandidates, aln
    ), unit = "s", times = times)

    list(
        "timeConsensus" = timeConsensus,
        "timeOligos" = timeOligos,
        "timeAssays" = timeAssays,
        "timeCheckMatch" = timeCheckMatch
    )
}

prepareDataset <- function(x, label) {
    if (is.list(x)) {
        x <- do.call("rbind", x)
    }
    levels(x$expr) <- c(
        "consensusProfile", "designOligos", "designAssays", "checkMatch"
    )
    x <- as.data.frame(x)
    rownames(x) <- NULL
    x$time <- x$time/1e9 ## time unit is in ns in the data frame format
    Assay <- rep(label, nrow(x))
    Type <- c(
        rep("Workflow", nrow(x[x$expr != "checkMatch", ])),
        rep("Addidional checks", nrow(x[x$expr == "checkMatch", ]))
    )
    x <- cbind(Type, Assay, x)
    levels(x$Type) <- c("Workflow", "Additional checks")
    x
}

plotMicrobenchmark <- function(x) {

    x$expr <- factor(x$expr, levels = c(
        "consensusProfile", "designOligos", "designAssays", "checkMatch")
    )

    ggplot(x, aes(x = expr, y = time)) +
        facet_wrap(~ Type, scales = "free") +
        geom_violin(aes(fill = Assay)) +
        scale_y_log10() +
        ylab("Time (s)") +
        xlab("") +
        theme_bw() +
        theme(
            axis.text.x = element_text(face = "italic"),
            strip.text.x = element_blank(),
            text = element_text(size = 18)
        )
}

## Import target alignment =====================================================

aln <- Biostrings::readDNAMultipleAlignment(
    "alignment.fasta", format = "fasta"
)
maskedAln <- aln
Biostrings::colmask(maskedAln, invert = TRUE) <- c(3800:5300, 5500:7000)

## Measure and plot ============================================================

# Set number of times
t <- 20

assayA <- makeMicrobenchmarkAssayA(aln, times = t)
assayAdf <- prepareDataset(assayA, "A")

assayB <- makeMicrobenchmarkAssayB(maskedAln, times = t)
assayBdf <- prepareDataset(assayB, "B")

both <- rbind(assayAdf, assayBdf)

mbPlot <- plotMicrobenchmark(both)
mbPlot

## Save result output ==========================================================

#write.table(both, "microbenchmark_results.txt", sep = "\t", header = TRUE,
#row.names = FALSE, quote = FALSE)

# ggsave("performance_plot.png", mbPlot, dpi = 300) #, width = 85, height = 40, units = "mm")

## Measure time to generate oligos with high degeneracy ========================

consensusA <- consensusProfile(aln, ambiguityThreshold = 0.1)
timeOligosA <-  microbenchmark(designOligos(consensusA,
                                            lengthPrimer = 18:26,
                                            maxDegeneracyPrimer = 64,
                                            gcClampPrimer = FALSE,
                                            lengthProbe = 16:26,
                                            maxDegeneracyProbe = 64
), unit = "s", times = t)

consensusB <- consensusProfile(maskedAln, ambiguityThreshold = 0.01)
timeOligosB <-  microbenchmark(designOligos(consensusB,
                                            lengthPrimer = 25:40,
                                            tmPrimer = c(60, 75),
                                            gcClampPrimer = FALSE,
                                            designStrategyPrimer = "mixed",
                                            maxDegeneracyPrimer = 64,
                                            probe = FALSE
), unit = "s", times = t)

timeOligosA
timeOligosB
