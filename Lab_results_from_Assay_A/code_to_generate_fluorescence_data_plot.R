## Code to plot amplification curves

library(ggplot2)
library(reshape2)

## Import fluorescence data
# Use clipboard to import qpcr_fluorescence_data, for some reason it doesn't
# work otherwise
x <- read.table("clipboard", header = TRUE, sep = "\t")
x <- x[, -1]

## Import and add sample names to fluorescence data
sampleNames <- read.table("qpcr_sample_names.txt", header = TRUE, sep = "\t")

sampleNames$Sample.Name <- trimws(sampleNames$Sample.Name, "both", " ")
sampleNames$Gene.Name <- trimws(sampleNames$Gene.Name, "both", " ")
id <- paste0(sampleNames$Sample.Name, "_", sampleNames$Gene.Name)
names(x) <- id

## Keep only clinical samples and NTCs
x <- x[!grepl("Sample|EC", names(x))]

## Convert to long format
y <- melt(x)

## Add important variables
Cycle <- rep(seq_len(45), length(unique(y$variable)))
Assay <- rep(c("Assay A", "Reference"), each = nrow(y) / 2)
Replicate <- rep(rep(seq_len(2), each = 45), length(unique(y$variable)))
Sample <- gsub("_Ny.1|ISO.1|S|_ISO|_Ny|.1$", "", y$variable)
df <- cbind.data.frame(
  Assay, Cycle, Replicate, Sample, "Fluorescence" = y$value
)
df$Assay <- factor(df$Assay, levels = c("Assay A", "Reference"))
df$Sample <- factor(
  df$Sample,
  levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "NTC")
)

## Make plot
pcrPlot <- ggplot(
  data = df, aes(
    x = Cycle, y = Fluorescence, group = interaction(Sample, Replicate)
  )
) +
  facet_wrap(~ Assay) +
  geom_line(aes(color = Sample)) +
  theme_bw() + theme(text = element_text(size = 18))

pcrPlot

## Save plot

# ggsave("qPCR.png", pcrPlot)
