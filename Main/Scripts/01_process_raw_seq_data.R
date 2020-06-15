# Process raw sequences into phyloseq object for analyses

# Load dada2 and prep ####
library(dada2); packageVersion("dada2")
library(vegan)
library(ggplot2)
library(dplyr)


# File parsing - For this, we will use only the forward illumina reads
path <- "D:/GitHub/Seagrass Diversity/seagrass-diversity/1912KMI-0008/Final_Processed" # CHANGE to the directory containing your demultiplexed fastq files
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
if(!file.exists("./Environment/Run 1/errF.Rds")){
  errF <- learnErrors(filtFs, multithread=FALSE, MAX_CONSIST = 10, nbases = 1e8)
  saveRDS(errF, file = "./Environment/Run 1/errF.Rds")
} else {
  errF <- readRDS(file = "./Environment/Run 1/errF.Rds")
}
if(!file.exists("./Environment/Run 1/errR.Rds")){
  errR <- learnErrors(filtRs, multithread=FALSE, MAX_CONSIST = 10, nbases = 1e8)
  saveRDS(errR, file = "./Environment/Run 1/errR.Rds")
} else {
  errR <- readRDS(file = "./Environment/Run 1/errR.Rds")
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

# Filtering according to length as we are only looking at V4 region of 16S which is about 250-256 bps long and most of our sequences lie within this region.
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]


# Remove Chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)


# Track Reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "lengthFiltered", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$filter.loss = (track[,1]-track[,2])/track[,1]
write.csv(track, file = "../Output/Run 1/Run 1 read_counts_at_each_step.csv", row.names = TRUE)


# Removing objects not required ahead to conserve RAM because some of the next steps are memory intensive
rm(dadaFs)
rm(dadaRs)
rm(derepF)
rm(derepR)
rm(errF)
rm(errR)
rm(mergers)
rm(out)
rm(seqtab)
rm(seqtab2)


# remove contaminants using controls ####
library(decontam)

# import metadata ####
meta = read.csv("./Run 1 metadata.csv", stringsAsFactors = FALSE)[,1:5]
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
write.csv(contams, file = "../Output/Run 1/Run 1 likely_contaminants.csv", row.names = TRUE)

# remove them
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$controls == FALSE,]
meta = meta[meta$controls == FALSE,]

saveRDS(meta, file = "../Output/Run 1/Run 1 final_meta.RDS")


# Assign Taxonomy ####
taxa <- assignTaxonomy(seqtab.nochim, "../tax/silva_nr_v138_train_set.fa.gz", multithread=FALSE)
# taxa <- addSpecies(taxa, "../tax/silva_species_assignment_v138.fa.gz") # Species level not required for our purposes and also memory intensive

write.csv(as.data.frame(seqtab.nochim), file = "../Output/Run 1/Run 1 SeqTable_no-chimera_no-contams.csv", row.names = TRUE, quote = FALSE)
saveRDS(seqtab.nochim, file = "../Output/Run 1/Run 1 clean_dada2_seqtable.RDS")
saveRDS(taxa, file = "../Output/Run 1/Run 1 Silva_Taxonomy_from_dada2.RDS")

seqtab.nochim <- readRDS(file = "../Output/Run 1/Run 1 clean_dada2_seqtable.RDS")
taxa <- readRDS(file = "../Output/Run 1/Run 1 Silva_Taxonomy_from_dada2.RDS")

# Removing Mitochondria and Chloroplast as we are only analysing bacteria

# Remove "mitochondria" taxa
is.mitochon <- taxa[,"Family"] %in% "Mitochondria"
taxa <- taxa[!is.mitochon,]
seqtab.nomitochon <- seqtab.nochim[,!is.mitochon]

# Remove "Chloroplast" taxa
is.chloro <- taxa[,"Order"] %in% "Chloroplast"
taxa <- taxa[!is.chloro,]
seqtab.nochloro <- seqtab.nomitochon[,!is.chloro]

saveRDS(seqtab.nochloro, file = "../Output/Run 1/Run 1 final_clean_dada2_seqtable.RDS")
saveRDS(taxa, file = "../Output/Run 1/Run 1 final_Silva_Taxonomy_from_dada2.RDS")

# Adding Read Counts for new steps
nochim <- rowSums(seqtab.nochim)
nomitochon <- rowSums(seqtab.nomitochon)
nochloro <- rowSums(seqtab.nochloro)
netloss <- (nochim-nochloro)/nochim
track_decontam <- data.frame(nochim, nomitochon, nochloro, netloss)
write.csv(track_decontam, file = "../Output/Run 1/Run 1 read_counts_at_each_step_after_decontam.csv", row.names = TRUE)


######################
### Merging runs ###
rm(list = ls())

# Import Run 1 data
seqtab.nochloro.Run1 <- readRDS(file = "../Output/Run 1/Run 1 final_clean_dada2_seqtable.RDS")
taxa.Run1 <- readRDS(file = "../Output/Run 1/Run 1 final_Silva_Taxonomy_from_dada2.RDS")
meta.Run1 <- readRDS(file = "../Output/Run 1/Run 1 final_meta.RDS")

# Import Run 2 data
seqtab.nochloro.Run2 <- readRDS(file = "../Output/Run 2/Run 2 final_clean_dada2_seqtable.RDS")
taxa.Run2 <- readRDS(file = "../Output/Run 2/Run 2 final_Silva_Taxonomy_from_dada2.RDS")
meta.Run2 <- readRDS(file = "../Output/Run 2/Run 2 final_meta.RDS")

# Merging
seqtab.final <- mergeSequenceTables(table1 = seqtab.nochloro.Run1, table2 = seqtab.nochloro.Run2, repeats = "error", orderBy = NULL, tryRC = TRUE)

taxa.final <- assignTaxonomy(seqtab.final, "../tax/silva_nr_v138_train_set.fa.gz", multithread=FALSE)
taxa.final <- addSpecies(taxa.final, "../tax/silva_species_assignment_v138.fa.gz")

meta.final <- merge(meta.Run1, meta.Run2, by=c("SampleID", "Location", "Species", "Structure.DNA.Extracted.from", "GPS.Coordinates", "controls"), all=TRUE)
meta.final <- meta.final[,-6]
row.names(meta.final) <- NULL
row.names(meta.final) <- meta.final$SampleID

#write.csv(as.data.frame(seqtab.final), file = "../Output/Merged/SeqTable_no-chimera_no-contams_no-mitochon_no-chloro.csv", row.names = TRUE, quote = FALSE)
#saveRDS(seqtab.final, file = "../Output/Merged/clean_dada2_seqtable.RDS")
#saveRDS(taxa.final, file = "../Output/Merged/Silva_Taxonomy_from_dada2.RDS")
#saveRDS(meta.final, file = "../Output/Merged/meta.RDS")

write.csv(t(seqtab.final), file = "../Output/Merged/clean_dada2_seqtable.csv", row.names = TRUE)
write.csv(taxa.final, file = "../Output/Merged/Silva_Taxonomy_from_dada2.csv", row.names = TRUE)
write.csv(meta.final, file = "../Output/Merged/meta.csv", row.names = TRUE)

seqtab.final <- readRDS(file = "../Output/Merged/clean_dada2_seqtable.RDS")
taxa.final <- readRDS(file = "../Output/Merged/Silva_Taxonomy_from_dada2.RDS")
meta.final <- readRDS(file = "../Output/Merged/meta.RDS")
######################
# Hand off to Phyloseq ####

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

ps <- phyloseq(otu_table(seqtab.final, taxa_are_rows=FALSE), 
               sample_data(meta.final), 
               tax_table(taxa.final))

# Changing actual sequences to custom ASV IDs while actual sequences stored in refseq()
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#saveRDS(ps, file = "../Output/Merged/clean_phyloseq_object.RDS")
ps = readRDS(file = "../Output/Merged/clean_phyloseq_object.RDS")
