#Dada2 + decontamination

library(dada2)

#path (chnage paths based on your file directory)

trimmed_path <- "trimmed"
filtered_path <- "filtered"
output_path <- "results/dada2"

dir.create(filtered_path, showWarnings = FALSE)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

#List trimmed fastq files
fnFs <- sort(list.files(trimmed_path, pattern = "_R1_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(trimmed_path, pattern = "_R2_trimmed.fastq.gz", full.names = TRUE))

#extract sample name

get.sample.name <- function(x) {
  sub("_R[12]_trimmed.fastq.gz", "", basename(x))
}
sample.names <- sapply(fnFs, get.sample.name)

#filtere
filtFs <- file.path(filtered_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_path, paste0(sample.names, "_R_filt.fastq.gz"))

cat("Running filterAndTrim...\n")

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen=c(95,95),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

write.csv(out, file.path(output_path, "filtering_summary.csv"))


#learn errors
cat("Learning error rates...\n")

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

#dada2 inference
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)



#Merge paired reads

cat("Merging paired reads...\n")

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)


seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE
)

write.csv(seqtab.nochim,
          file = file.path(output_path, "seqtab_nochim.csv"))


#SIMPLE DECONTAMINATION


metadata <- read.csv("metadata.csv",
                     stringsAsFactors = FALSE,
                     check.names = TRUE)

#only seqtab nochim
metadata <- metadata[metadata$Sample.ID %in% rownames(seqtab.nochim), ]

cat("Running decontamination...\n")

seqtab.clean <- seqtab.nochim

# for my sample EB represent extraction blanks(rep1, rep2, rep3)
eb.reps <- rownames(seqtab.clean)[grepl("^(EB|neg)", rownames(seqtab.clean))]


if (length(eb.reps) == 0) {
  stop("No EB samples found")
}


#combine all blank replication
EB.combined <- apply(
  seqtab.clean[eb.reps, , drop = FALSE],
  2,
  max
)

#for remove conatamination from all samples
for (asv in colnames(seqtab.clean)) {

  background <- EB.combined[asv]

  seqtab.clean[, asv] <- seqtab.clean[, asv] - background

  seqtab.clean[, asv][seqtab.clean[, asv] < 0] <- 0
}

write.csv(
  seqtab.clean,
  file = file.path(output_path, "seqtab_decontaminated_EBcombined.csv")
)

cat("DONE")
