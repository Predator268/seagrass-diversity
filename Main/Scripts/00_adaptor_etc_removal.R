# Adaptor, Primer, Pad, and Linker Removal

# Raw data is used as input to the 01_Cutadapt_Adaptor_etc_Removal_Script.R script to remove adaptors, primers, paddings, and linkers from sequences.
# This outputs N_filtered folder, Adaptor_Removed folder and Adaptor_Primer_Removed folder, all inside the raw data folder.
# The folders are fed to the script in the above order and the final folder is the one used for analysis.

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")


# Declaring path to raw data
path <- "../1912KMI-0008"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))


# Declaring adaptor and primmers and pad+linker
FWD_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"  
REV_adaptor <- "AATGATACGGCGACCACCGAGATCTACAC"

FWD_primer <- "GTGCCAGCMGCCGCGGTAA"  
REV_primer <- "GGACTACHVGGGTWTCTAAT"

FWD_pl <- "TATGGTAATTGT"  
REV_pl <- "AGTCAGTCAGCC"

# Getting all orientations of each sequence
allOrients <- function(sequence) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(sequence)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD_adaptor.orients <- allOrients(FWD_adaptor)
REV_adaptor.orients <- allOrients(REV_adaptor)
FWD_adaptor.orients

FWD_primer.orients <- allOrients(FWD_primer)
REV_primer.orients <- allOrients(REV_primer)
FWD_primer.orients

FWD_pl.orients <- allOrients(FWD_pl)
REV_pl.orients <- allOrients(REV_pl)
FWD_pl.orients


# Pre-filtering to remove ambiguous bases (Ns)
fnFs.filtN <- file.path(path, "N_filtered", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "N_filtered", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)


### Removal Process ###

#Using cutadapt for trimming
cutadapt <- "C:/Users/User/AppData/Roaming/Python/Python38/Scripts/cutadapt"
system2(cutadapt, args = "--version")

## Adaptor Removal ##
# Counting presence of adaptors before removal
sequenceHits <- function(sequence, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(sequence, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.filtN[[5]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.filtN[[5]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.filtN[[5]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.filtN[[5]]))

# Removing adaptors
path.cut_adaptor <- file.path(path, "Adaptor_Removed")
if(!dir.exists(path.cut_adaptor)) dir.create(path.cut_adaptor)
fnFs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnFs))
fnRs.cut_adaptor <- file.path(path.cut_adaptor, basename(fnRs))

FWD_adaptor.RC <- dada2:::rc(FWD_adaptor)
REV_adaptor.RC <- dada2:::rc(REV_adaptor)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_adaptor, "-a", REV_adaptor.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_adaptor, "-A", FWD_adaptor.RC) 

# Run Cutadapt
outputStatsAdaptor <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "-o", fnFs.cut_adaptor[i], "-p", fnRs.cut_adaptor[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }
)
cat(outputStatsAdaptor, file="./Output/cutadapt_adaptor_trimming_stats.txt", sep="\n", append = FALSE)

# Counting adaptors after removal
rbind(FWD_adaptor.ForwardReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[5]]), 
      FWD_adaptor.ReverseReads = sapply(FWD_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[5]]), 
      REV_adaptor.ForwardReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnFs.cut_adaptor[[5]]), 
      REV_adaptor.ReverseReads = sapply(REV_adaptor.orients, sequenceHits, fn = fnRs.cut_adaptor[[5]]))

## Primer Removal
#Counting Primers before removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_adaptor[[5]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_adaptor[[5]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_adaptor[[5]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_adaptor[[5]]))

# Removing primers
path.cut_primer <- file.path(path, "Primer_Adaptor_Removed")
if(!dir.exists(path.cut_primer)) dir.create(path.cut_primer)
fnFs.cut_primer <- file.path(path.cut_primer, basename(fnFs))
fnRs.cut_primer <- file.path(path.cut_primer, basename(fnRs))

FWD_primer.RC <- dada2:::rc(FWD_primer)
REV_primer.RC <- dada2:::rc(REV_primer)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_primer, "-a", REV_primer.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_primer, "-A", FWD_primer.RC) 

# Run Cutadapt
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "-o", fnFs.cut_primer[i], "-p", fnRs.cut_primer[i], # output files
                               fnFs.cut_adaptor[i], fnRs.cut_adaptor[i])) # input files
  }
)
cat(outputStatsPrimer, file="./Output/cutadapt_primer_trimming_stats.txt", sep="\n", append = FALSE)

# Counting primers after removal
rbind(FWD_primer.ForwardReads = sapply(FWD_primer.orients, sequenceHits, fn = fnFs.cut_primer[[5]]), 
      FWD_primer.ReverseReads = sapply(FWD_primer.orients, sequenceHits, fn = fnRs.cut_primer[[5]]), 
      REV_primer.ForwardReads = sapply(REV_primer.orients, sequenceHits, fn = fnFs.cut_primer[[5]]), 
      REV_primer.ReverseReads = sapply(REV_primer.orients, sequenceHits, fn = fnRs.cut_primer[[5]]))

## Padding and Linker removal
# Counting padding + linker before removal
rbind(FWD_pl.ForwardReads = sapply(FWD_pl.orients, sequenceHits, fn = fnFs.cut_primer[[5]]), 
      FWD_pl.ReverseReads = sapply(FWD_pl.orients, sequenceHits, fn = fnRs.cut_primer[[5]]), 
      REV_pl.ForwardReads = sapply(REV_pl.orients, sequenceHits, fn = fnFs.cut_primer[[5]]), 
      REV_pl.ReverseReads = sapply(REV_pl.orients, sequenceHits, fn = fnRs.cut_primer[[5]]))

# Removing pad and linker
path.cut_pl <- file.path(path, "Final_Processed")
if(!dir.exists(path.cut_pl)) dir.create(path.cut_pl)
fnFs.cut_pl <- file.path(path.cut_pl, basename(fnFs))
fnRs.cut_pl <- file.path(path.cut_pl, basename(fnRs))

FWD_pl.RC <- dada2:::rc(FWD_pl)
REV_pl.RC <- dada2:::rc(REV_pl)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD_pl, "-a", REV_pl.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV_pl, "-A", FWD_pl.RC) 
# Run Cutadapt
#system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
#                           "-m", 1,
#                           "-o", fnFs.cut_pl[5], "-p", fnRs.cut_pl[5], # output files
#                           fnFs.cut_primer[5], fnRs.cut_primer[5])) # input files
outputStatsPrimer <- capture.output(
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-m", 1,
                               "-o", fnFs.cut_pl[i], "-p", fnRs.cut_pl[i], # output files
                               fnFs.cut_primer[i], fnRs.cut_primer[i])) # input files
  }
)
cat(outputStatsPL, file="./Output/cutadapt_pl_trimming_stats.txt", sep="\n", append = FALSE)

# Counting pad and linker after removal
rbind(FWD_pl.ForwardReads = sapply(FWD_pl.orients, sequenceHits, fn = fnFs.cut_pl[[5]]), 
      FWD_pl.ReverseReads = sapply(FWD_pl.orients, sequenceHits, fn = fnRs.cut_pl[[5]]), 
      REV_pl.ForwardReads = sapply(REV_pl.orients, sequenceHits, fn = fnFs.cut_pl[[5]]), 
      REV_pl.ReverseReads = sapply(REV_pl.orients, sequenceHits, fn = fnRs.cut_pl[[5]]))


# Forward and reverse fastq filenames have the format:
cutFs_raw <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs_raw <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))
cut_raw <- sort(list.files(path, pattern = "fastq.gz", full.names = TRUE))

path.filtN <- file.path(path, "N_filtered")
cutFs_filtN <- sort(list.files(path.filtN, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs_filtN <- sort(list.files(path.filtN, pattern = "_R2.fastq.gz", full.names = TRUE))
cut_filtN <- sort(list.files(path.filtN, pattern = "fastq.gz", full.names = TRUE))

cutFs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs_adaptor <- sort(list.files(path.cut_adaptor, pattern = "_R2.fastq.gz", full.names = TRUE))
cut_adaptor <- sort(list.files(path.cut_adaptor, pattern = "fastq.gz", full.names = TRUE))

cutFs_primer <- sort(list.files(path.cut_primer, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs_primer <- sort(list.files(path.cut_primer, pattern = "_R2.fastq.gz", full.names = TRUE))
cut_primer <- sort(list.files(path.cut_primer, pattern = "fastq.gz", full.names = TRUE))

cutFs_pl <- sort(list.files(path.cut_pl, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs_pl <- sort(list.files(path.cut_pl, pattern = "_R2.fastq.gz", full.names = TRUE))
cut_pl <- sort(list.files(path.cut_pl, pattern = "fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs_pl, get.sample.name))
head(sample.names)


## Quality Profiles ##
# Raw
plotQualityProfile(cutFs_raw[6:11])
plotQualityProfile(cutRs_raw[6:11])

# Pre-filtered
plotQualityProfile(cutFs_filtN[6:11])
plotQualityProfile(cutRs_filtN[6:11])

# Adaptor Removed
plotQualityProfile(cutFs_adaptor[6:11])
plotQualityProfile(cutRs_adaptor[6:11])

# Primer Removed
plotQualityProfile(cutFs_primer[6:11])
plotQualityProfile(cutRs_primer[6:11])

# Pad+Linker Removed (Final)
plotQualityProfile(cutFs_pl[6:11])
plotQualityProfile(cutRs_pl[6:11])


## Comparing files and sample properties for each step
library(seqTools)
packageVersion("seqTools")
library(microseq)
packageVersion("microseq")

samples <- sort(list.files(path, pattern = "fastq.gz", full.names = FALSE))

# Commented code could be used to count and compare the file size of each sample file during each step
#raw_size <- numeric(length(cut_raw))
#prefiltered_size <- numeric(length(cut_filtN))
#adaptorremoved_size <- numeric(length(cut_adaptor))
#primerremoved_size <- numeric(length(cut_primer))
#padlinkerremoved_size <- numeric(length(cut_pl))
#fractionOfRawDataInFinal_size <- numeric(length(samples))

#for (j in 1:length(samples)){
#  raw_size[j] <- file.size(cut_raw[j])
#  prefiltered_size[j] <- file.size(cut_filtN[j])
#  adaptorremoved_size[j] <- file.size(cut_adaptor[j])
#  primerremoved_size[j] <- file.size(cut_primer[j])
#  padlinkerremoved_size[j] <- file.size(cut_pl[j])
#  fractionOfRawDataInFinal_size[j] <- padlinkerremoved_size[j] / raw_size[j]
#}

#file_info <- data.frame(samples, raw_size, prefiltered_size, adaptorremoved_size, primerremoved_size, padlinkerremoved_size, fractionOfRawDataInFinal_size, stringsAsFactors = FALSE)

# Sample Summaries of Trimming
#raw_reads <- nReads(fastqq(cut_raw[1:3]))
#prefiltered_reads <- nReads(fastqq(cut_filtN[1:3]))
#adaptorremoved_reads <- nReads(fastqq(cut_adaptor[1:3]))
#primerremoved_reads <- nReads(fastqq(cut_primer[1:3]))
#padlinkerremoved_reads <- nReads(fastqq(cut_pl[1:3]))
#fractionOfRawDataInFinal_reads <- padlinkerremoved_reads / raw_reads

#sample_summaries <- data.frame(samples[1:3], raw_reads, prefiltered_reads, adaptorremoved_reads, primerremoved_reads, padlinkerremoved_reads, fractionOfRawDataInFinal_reads,raw_seqs, prefiltered_seqs, adaptorremoved_seqs, primerremoved_seqs, padlinkerremoved_seqs, fractionOfRawDataInFinal_seqs, stringsAsFactors = FALSE)

# Comparing number of reads for each sample during each step

raw_seqs <- numeric(length(cut_raw))
prefiltered_seqs <- numeric(length(cut_raw))
adaptorremoved_seqs <- numeric(length(cut_raw))
primerremoved_seqs <- numeric(length(cut_raw))
padlinkerremoved_seqs <- numeric(length(cut_raw))
fractionOfRawDataInFinal_seqs <- padlinkerremoved_seqs / raw_seqs

for (j in 1:length(samples)){
  raw_seqs[j] <- nrow(readFastq(cut_raw[j]))
  prefiltered_seqs[j] <- nrow(readFastq(cut_filtN[j]))
  adaptorremoved_seqs[j] <- nrow(readFastq(cut_adaptor[j]))
  primerremoved_seqs[j] <- nrow(readFastq(cut_primer[j]))
  padlinkerremoved_seqs[j] <- nrow(readFastq(cut_pl[j]))
}
fractionOfRawDataInFinal_seqs <- padlinkerremoved_seqs / raw_seqs


sample_summaries <- data.frame(samples, raw_seqs, prefiltered_seqs, adaptorremoved_seqs, primerremoved_seqs, padlinkerremoved_seqs, fractionOfRawDataInFinal_seqs, stringsAsFactors = FALSE)
write.csv(sample_summaries, file = "Number of Sequence Summary - Run1.csv")
summary(sample_summaries)
hist(sample_summaries$raw_seqs)
hist(sample_summaries$padlinkerremoved_seqs)
