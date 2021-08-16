#!/bin/r env

library(dada2)
library(limma)
library(argparse)
library(Biostrings)
library(data.table)
library(stringr)
library(seqinr)
library(parallel)
library(doMC)

seq_align <- function(seqs_df, path_to_ref, overlap = TRUE, parallel = TRUE)
{
  align_df <- data.frame()
  if (file.exists(path_to_ref)) {
    ref <- read.fasta(path_to_ref)
    ref_str <- toupper(sapply(ref, c2s))
  } else {
    stop(paste("File",path_to_ref,"not found!"))
  }
  
  sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
  if (!overlap) {
    split <- paste0(rep("N",10), collapse = "")
    seq_split <- strsplit2(x = seqs_df[,1], split = split)
    seq_all <- data.frame(sequence=paste0(seq_split[,1],seq_split[,2]), hapid = seqs_df[,2])
  } else {
    seq_all <- data.frame(sequence=seqs_df[,1], hapid = seqs_df[,2])
  }
  if (parallel) {
    registerDoMC(detectCores())
    align_df <- foreach(seq_1=1:length(seq_all$sequence), .combine = "rbind") %dopar% {
      aln <- pairwiseAlignment(ref_str, seq_all$sequence[seq_1], substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
      num <- which.max(score(aln))
      patt <- c(alignedPattern(aln[num]), alignedSubject(aln[num]))
      dist <- adist(as.character(patt)[1],as.character(patt)[2])
      ind <- sum(str_count(as.character(patt),"-"))
      data.frame(hapid = seq_all$hapid[seq_1], 
                 hapseq = as.character(patt)[2], 
                 refseq = as.character(patt)[1], 
                 refid = names(patt)[1], 
                 aln_score = score(aln[num]), 
                 snv_dist = (dist - ind),
                 indel_dist = ind)
    }
  } else {
    for (seq_1 in 1:length(seq_all$sequence)) {
      aln <- pairwiseAlignment(ref_str, seq_all$sequence[seq_1], substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
      num <- which.max(score(aln))
      patt <- c(alignedPattern(aln[num]), alignedSubject(aln[num]))
      dist <- adist(as.character(patt)[1],as.character(patt)[2])
      ind <- sum(str_count(as.character(patt),"-"))
      df <- data.frame(hapid = seq_all$hapid[seq_1], 
                       hapseq = as.character(patt)[2], 
                       refseq = as.character(patt)[1], 
                       refid = names(patt)[1], 
                       aln_score = score(aln[num]), 
                       snv_dist = (dist - ind),
                       indel_dist = ind)
      align_df <- rbind(align_df,df)
    }
  }
  return(align_df)
}

parser <- ArgumentParser()
parser$add_argument("-s", "--seqtab", 
                    help="Path to input")
parser$add_argument("-ref", "--reference",
                    help="Path to reference fasta sequences")
parser$add_argument("-ref2", "--reference2")
parser$add_argument("-b","--bimera", help="ASV File with identifed bimeras")
parser$add_argument("-o", "--output", 
                    help="Path to output file")
parser$add_argument("--fasta", action='store_true', help="Write ASV sequences separately into fasta file")
parser$add_argument("-snv", "--snv_filter", 
                    help="Path to file for filtering ASVs based on edit distance")
parser$add_argument("--indel_filter", 
                    help="Specify proportion of ASV length (between 0 and 1) to target length for filtering based on indels")
parser$add_argument("--strain", default="3D7", help="Name of Specific strain to map to. Defaults to 3D7")
parser$add_argument("--strain2", help="Name of second strain if mapping to 2 different strains")
parser$add_argument("--parallel", action='store_true', help="Enable parallel processing")



args <- parser$parse_args()
if (!is.null(args$reference)) {
  path_to_refseq <- args$reference
  if (!is.null(args$strain)) {
    strains <- args$strain
  } else {
    stop("Name of target strain (--strain) is required")
  }
  if (!is.null(args$reference2)) {
    path_to_refseq <- c(args$reference,args$reference2)
    if (!is.null(args$strain2)) {
      strains <- c(args$strain,args$strain2)
    } else {
      stop("Name of second target strain (--strain2) required if --reference2 is given")
    }
  }
} else {
  stop("Reference fasta file with target sequences (--reference) is required")
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

if (!is.null(args$output)) {
  output <- args$output
  } else {
    stop("Output filename not provided!")
  }

if (exists("seqtab"))
{
  seqs <- colnames(seqtab)
  nsample=nrow(seqtab)
  hapid <- paste0("ASV",1:length(seqs))
  # DataFrame for aligning to truth set
  seqs_df <- data.frame(sequence = seqs, hapid = hapid)
  # Change colnames of ASV from sequences to ASV ids
  seqtab_haps <- seqtab
  colnames(seqtab_haps) <- hapid
  ## ASV summary table
  total_reads <- apply(seqtab_haps,2,sum)
  total_samples <- apply(seqtab_haps,2,function(x) sum(x != 0))
  asvdf <- data.frame(hapid = seqs_df$hapid,
                      total_reads = total_reads,
                      total_samples = total_samples,
                      strain = "N")
  asvdf$hapid <- as.character(asvdf$hapid)
  asvdf$strain <- as.character(asvdf$strain)

} else {
  stop("cannot find DADA2 sequence table!")
}

for (p in 1:length(path_to_refseq)) {
  refasta <- read.fasta(path_to_refseq[p])
  refseq <- toupper(sapply(refasta,c2s))
  amplicons <- as.character(names(refseq))

  # Alignment with RefSet
  align_df <- seq_align(seqs_df = seqs_df, path_to_ref = path_to_refseq[p], overlap = TRUE, parallel = args$parallel)
  align_df$refid <- as.character(align_df$refid)
  align_df$hapid <- as.character(align_df$hapid)

  ## Map truthset onto ASV summary table based on exact and inexact matches to truth set
  df <- align_df[,c(1,4,6,7)]
  colnames(df) <- c("hapid", paste0("refid_",strains[p]), paste0("snv_dist_from_",strains[p]),paste0("indel_dist_from_",strains[p]))
  asvdf <- merge(asvdf, df, by = "hapid", sort = FALSE)
  asvdf$strain[(align_df[,6] == 0 & align_df[,7] == 0)] <- as.character(strains[p])
}

if (!is.null(args$snv_filter)) {
  filter_file <- args$snv_filter
  if (file.exists(filter_file)) {
    VariantCounts <- fread(filter_file)
    asvdf$snv_filter <- NA
    # For Sliding edit distance and length based filter
    for (i in 1:nrow(asvdf)) {
      refid <- asvdf[i,paste0("refid_",strains[1])]
      hapdist <- asvdf[i,paste0("snv_dist_from_",strains[1])]
      tardist <- VariantCounts[c(VariantCounts[,1] == refid),2]
      if (hapdist <= tardist) {
        sfil <- "PASS"
      } else {
        sfil <- "FAIL"
      }
      asvdf$snv_filter[i] <- sfil
    }
  } else {
    warning(paste("File",filter_file,"not found!. Skipping SNV based filtering.."))
  }
} else {
  print("SNV filter file not given. Skipping SNV based filtering..")
}

if (!is.null(args$indel_filter)) {
  if (file.exists(args$indel_filter)) {
    Indelcutoff <- fread(args$indel_filter)
  } else {
    Indelcutoff <- as.numeric(args$indel_filter)
  }
  if (is.na(Indelcutoff)) {
    warning("INDEL based filter threshold argument not valid. Using default proportion of +-90%..")
    Indelcutoff <- 0.90
  }
  asvdf$indel_filter <- NA
  refseq <- toupper(sapply(read.fasta(path_to_refseq[1]),c2s))
  for (i in 1:nrow(asvdf)) {
    haplen <- nchar(as.character(seqs_df$sequence[i]))
    refid <- asvdf[i,paste0("refid_",strains[1])]
    if (length(Indelcutoff) != 1) {
      ind <- Indelcutoff[c(Indelcutoff[,1] == refid),2]
    } else {
      ind <- Indelcutoff
    }
    tarlen <- nchar(refseq[refid])
    lprop <- (haplen/tarlen)
    if (lprop < ind | lprop > (1/ind)) {
      ifil <- "FAIL"
    } else {
      ifil <- "PASS"
    }
    asvdf$indel_filter[i] <- ifil
  }
} else {
  print("INDEL based filter threshold not given. Skipping length filtering..")
}

if (!is.null(args$bimera)) {
  if (file.exists(args$bimera)) {
    bimera <- fread(args$bimera)
    seqs_df <- merge(seqs_df,bimera, by = "sequence", all = TRUE, sort = FALSE)
    asvdf <- merge(asvdf,seqs_df[,2:3], by = "hapid", all = TRUE, sort = FALSE)
  } else {
    warning(paste("File",args$bimera,"not found. Skipping bimera flag.."))
  }
} else {
  print("Bimeric ASV file not given. Skipping bimera flag..")
}

write.table(asvdf, file = output, sep = "\t", quote = FALSE, row.names = FALSE)

if (args$fasta) {
  write.fasta(sapply(seqs, s2c), names = hapid, file.out = paste0(dirname(output),"/ASVSeqs.fasta"), nbchar = 600)
}
