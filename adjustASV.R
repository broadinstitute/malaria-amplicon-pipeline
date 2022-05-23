#!/bin/r env

library(seqinr)
library(data.table)
library(argparse)
library(Biostrings)
library(parallel)
library(doMC)

parser <- ArgumentParser()
parser$add_argument("-s", "--seqtab", 
                    help="Path to input")
parser$add_argument("-ref", "--reference",
                    help="Path to reference fasta sequences")
parser$add_argument("-o", "--output",
                    help="Path to output for corrected ASV list")

args <- parser$parse_args()
path_to_refseq <- args$reference

if (file.exists(path_to_refseq)) {
  ref <- toupper(sapply(read.fasta(path_to_refseq),c2s))
} else {
  stop("Reference file not found!")
}

if (!is.null(args$seqtab)) {
  seqfile <- args$seqtab
  if (file.exists(seqfile)) {
    seqtab <- as.matrix(fread(seqfile), rownames=1)
  } else {
    stop(paste("ASV sequence table file",seqtab,"not found!"))
  }
} else {
  stop("Sequence table file (--seqtab) is required")
}

sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
seqs <- as.character(colnames(seqtab))

registerDoMC(detectCores())
df <- foreach(i=1:length(seqs), .combine = "rbind") %dopar% {
  map <- pairwiseAlignment(ref, seqs[i],substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = TRUE)
  tar = ref[which.max(map)]
  seq <- strsplit(seqs[i],"NNNNNNNNNN")[[1]]
  aln <- pairwiseAlignment(seq[1:2], tar, substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE, type = 'overlap')
  con <- compareStrings(consensusString(aln[1]),consensusString(aln[2]))
  overlap <- unlist(gregexpr("[[:alpha:]]",con))
  if (overlap == -1) {
    N = (nchar(seq[1])+nchar(seq[2])) - nchar(tar)
    stkN <- paste0(rep('N',abs(N)),collapse = '')
    correctedASV <- paste0(seq[1],stkN,seq[2])
  } else {
    N = length(overlap)
    correctedASV <- paste0(seq[1],substr(seq[2],(N+1),nchar(seq[2])))
  }
  if (nchar(correctedASV) != nchar(tar)) {
    N = NA
    correctedASV = NA
  }
  data.frame(target = names(tar),
             ASV = seqs[i],
             correctedASV = correctedASV,
             overlap = N)
}
write.table(df, file = args$output, sep = "\t", quote = FALSE, row.names = FALSE)
seqfile_corrected <- paste0(dirname(seqfile),"/seqtab_corrected.tsv")
colnames(seqtab) <- as.character(df$correctedASV)
write.table(seqtab, file = seqfile_corrected, sep = "\t", quote = FALSE, row.names = FALSE)