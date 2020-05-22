# Process raw sequences into phyloseq object for analyses

# Load dada2 and prep ####
library(dada2); packageVersion("dada2")
library(vegan)
library(ggplot2)
library(dplyr)

# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "D:/GitHub/Seagrass Diversity/seagrass-diversity/1912KMI-0009/Final_Processed" # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "Quality_Filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- list.files(path)
fastqs <- fns[grepl("R1.fastq.gz$", fns)] # CHANGE if different file extensions or to target only certain sequences
rm(fns)
list.files()
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

# visualize a couple of fwd read quality profiles to help you decide reasonable filtration parameters
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Filter and trim ####
filtFs <- file.path(path, "Quality_Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "Quality_Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# learn error rates ####
if(!file.exists("./Environment/Run 2/errF.Rds")){
  errF <- learnErrors(filtFs, multithread=FALSE, MAX_CONSIST = 10, nbases = 1e8)
  saveRDS(errF, file = "./Environment/Run 2/errF.Rds")
} else {
  errF <- readRDS(file = "./Environment/Run 2/errF.Rds")
}
if(!file.exists("./Environment/Run 2/errR.Rds")){
  errR <- learnErrors(filtRs, multithread=FALSE, MAX_CONSIST = 10, nbases = 1e8)
  saveRDS(errR, file = "./Environment/Run 2/errR.Rds")
} else {
  errR <- readRDS(file = "./Environment/Run 2/errR.Rds")
}
# sanity check
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Loop taking care of dereplication, sample inference and merging paired reads
dadaFs <- vector("list", length(sample.names))
dadaRs <- vector("list", length(sample.names))
mergers <- vector("list", length(sample.names))
names(dadaFs) <- sample.names
names(dadaRs) <- sample.names
names(mergers) <- sample.names
for (sam in sample.names){
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]], verbose = TRUE)
  derepR <- derepFastq(filtRs[[sam]], verbose = TRUE)
  dadaFs[[sam]] <- dada(derepF, err=errF, multithread = FALSE)
  dadaRs[[sam]] <- dada(derepR, err=errR, multithread = FALSE)
  mergers[[sam]] <- mergePairs(dadaFs[[sam]], derepF, dadaRs[[sam]], derepR, verbose = TRUE)
}
head(mergers[[1]])

# Make a sequence table ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]

# Remove Chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seq2tab)


# Track Reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$filter.loss = (track[,1]-track[,2])/track[,1]
write.csv(track, file = "../Output/Run2_read_counts_at_each_step.csv", row.names = TRUE)

###

rm(dadaFs)
rm(dadaRs)
rm(derepF)
rm(derepR)
rm(errF)
rm(errR)
rm(mergers)
rm(out)
rm(seqtab)


# remove contaminants using controls ####
library(decontam)

# import metadata ####
meta = read.csv("./Run 1 metadata.csv")[,1:5]
row.names(meta) <- NULL
row.names(meta) <- meta$SampleID
row.names(seqtab.nochim)

# reorder metadata
meta = meta[order(row.names(meta)),]

# Find controlsamples (extraction negatives) ####
meta$controls <- meta$Location == "Blank"

# find contaminants
contams = isContaminant(seqtab.nochim, neg = meta$controls, normalize = TRUE)
table(contams$contaminant)
write.csv(contams, file = "../Output/Run 2 likely_contaminants.csv", row.names = TRUE)

# remove them
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$controls == FALSE,]
meta = meta[meta$controls == FALSE,]

# Assign Taxonomy ####
taxa <- assignTaxonomy(seqtab.nochim, "../tax/silva_nr_v138_train_set.fa.gz", multithread=FALSE)
# taxa <- addSpecies(taxa, "../tax/silva_species_assignment_v138.fa.gz")

write.csv(as.data.frame(seqtab.nochim), file = "../Output/Run 2 SeqTable_no-chimera_no-contams.csv", row.names = TRUE, quote = FALSE)
saveRDS(seqtab.nochim, file = "../Output/Run 2 clean_dada2_seqtable.RDS")
saveRDS(taxa, file = "../Output/Run 2 Silva_Taxonomy_from_dada2.RDS")

seqtab.nochim <- readRDS(file = "../Output/Run 1 clean_dada2_seqtable.RDS")
taxa <- readRDS(file = "../Output/Run 1 Silva_Taxonomy_from_dada2.RDS")
######################

# Hand off to Phyloseq ####

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))


names(meta)
plot_richness(ps, x="Location", measures = "Shannon")

saveRDS(ps, file = "../Output/Run 1 clean_phyloseq_object.RDS")
ps = readRDS(file = "../Output/Run 1 clean_phyloseq_object.RDS")


