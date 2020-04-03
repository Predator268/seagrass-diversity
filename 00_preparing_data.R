library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
path <- "D:/GitHub/Seagrass Diversity/seagrass-diversity/data_raw/Run 1/1912KMI-0008"  ## CHANGE ME to the directory containing the fastq files.

list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

fnFsT <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRsT <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

# Identifying Primers
FWD <- "GTGCCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME...

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Pre-filtering to remove Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
pre_filter <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
pre_filter
head(pre_filter)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[2]]))

cutadapt <- "C:/users/user/appdata/roaming/python/python37/site-packages/cutadapt/" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

fnFsT.filtN <- file.path(path, "filtN/trimmed", basename(fnFsT))
fnRsT.filtN <- file.path(path, "filtN/trimmed", basename(fnRsT))
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFsT.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRsT.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFsT.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRsT.filtN[[1]]))
