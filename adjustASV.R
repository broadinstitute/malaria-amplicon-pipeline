#!/bin/r env

# Script Created by Ruchit Panchal 
# Updated by Jason T. Mohabir & Angela Early 

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

#path_to_refseq <- "pf3d7_ref_updated_v3.fasta"

if (file.exists(path_to_refseq)) {
  ref <- toupper(sapply(read.fasta(path_to_refseq),c2s))
} else {
  stop("Reference file not found!")
}

if (!is.null(args$seqtab)) {
  seqfile <- args$seqtab
  #seqfile <- "seqtab_nop.tsv"
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
df <- foreach(i=1:length(seqs),.combine = "rbind") %dopar% {

  # Figure out the amplicon of origin for the ASV 
  map <- pairwiseAlignment(ref,seqs[i],substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = TRUE)
  tar = ref[which.max(map)]
  
  # Split the ASVs at Ns 
  seq <- strsplit(seqs[i],"NNNNNNNNNN")[[1]]
  
  # Generate alignment of seq sequences to find the overlap (overlap < 12bp)
  aln <- pairwiseAlignment(seq[1:2], tar, substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE, type = 'overlap')
  
  # Generate the overlap sequences 
  con <- compareStrings(consensusString(aln[1]),consensusString(aln[2]))
  overlap <- unlist(gregexpr("[[:alpha:]]",con)) # Check the question 
  
  # Check if the ASVs overlap (merge)
  if (overlap == -1) {
    
    # Concatenation 
    N = abs((nchar(seq[1])+nchar(seq[2])) - nchar(tar))
    stkN <- paste0(rep('N',N),collapse = '')
    correctedASV <- paste0(seq[1],stkN,seq[2])
    status_flag <- 'concatenated'
    
  } 
  else {
    
    # Merging  
    N = length(overlap)
    correctedASV <- paste0(seq[1],substr(seq[2],(N+1),nchar(seq[2])))
    status_flag <- 'merged'
    
  }
  
  # Check correctedASV is the correct length 
  #if (nchar(correctedASV) != nchar(tar)) {
  #  # Concatenating 
  #  N = NA
  #  correctedASV = NA
  #}
  
  # Check if alignment contains unexpected characters such as '?' 
  if (grepl("\\?",con)) {
    correctedASV <- NA
    status_flag <- 'overlapMismatch'
  }
  
  data.frame(ix = i, 
             target = names(tar),
             uncorrectedASV = seqs[i],
             correctedASV = correctedASV,
             overlap = length(overlap),
             correctASVLength = nchar(correctedASV),
             targetLength = nchar(tar),
             correctASVLength_match_targetLength = nchar(correctedASV) == nchar(tar),
             statusFlag = status_flag,
             consensusCompare = con)
}

write.table(df, file = args$output, sep = "\t", quote = FALSE, row.names = FALSE)

seqfile_corrected <- paste0(dirname(seqfile),"/seqtab_corrected.tsv")

# NOTE: Duplicated ASVs do not have read support merged 
colnames(seqtab) <- as.character(df$correctedASV)
filt_seqtab <- seqtab[,!is.na(colnames(seqtab))]
filt_seqtab <- filt_seqtab[,!duplicated(colnames(filt_seqtab))]

write.table(filt_seqtab, file = seqfile_corrected, sep = "\t", quote = FALSE, row.names = TRUE)

save(list = ls(all.names = TRUE), file = "currentEnvironment_adjustASV.RData")


