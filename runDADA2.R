#!/bin/R env

<<<<<<< HEAD
library(dada2)
library(limma)
library(argparse)
library(data.table)
=======
library(argparse)
>>>>>>> ab1b0ba044df3bf30a708387c47c9c0788a39602

# Custom filtering, denoising parameters (if not default) can be provided as a separate config file?

parser <- ArgumentParser()
parser$add_argument("-p", "--path_to_meta", help="Path to input meta file listing fastqs (required)")
parser$add_argument("-c", "--class", help="Class specifying 'parasite' or 'vector' (required if '--default' is specified)")
parser$add_argument("-d", "--dir", help="Working directory path for writing all dada2 output files")
parser$add_argument("-o", "--output_filename", help="output tab-separated filename (required)")
parser$add_argument("-s", "--save_run", help="save Run as R workspace image")
parser$add_argument("-ee", "--maxEE",
                    help="Maximum expected errors for filtering forward and reverse read")
parser$add_argument("-tR", "--trimRight",
                    help="Length for trimming from right for both forward and reverse reads")
parser$add_argument("-mL", "--minLen", type="integer",
                    help="minimum length required for reads on both end. Shorter reads are discarded")
parser$add_argument("-tQ", "--truncQ",
                    help="truncate reads to first occurence of truncQ. All filtered reads have quality >= truncQ")
parser$add_argument("-id", "--matchIDs", type="integer",
                    help="match ids on fastqs to make sure reads on forward and reverse end are in same order")
parser$add_argument("-mC", "--max_consist", type="integer",
                    help="Maximum cycles for error model until consistency. If no convergence, error values at max_consist cycle are used")
parser$add_argument("-wA", "--omega_a", type="double",
                    help="P-value threshold in sample inference for forming a new partition")
parser$add_argument("-jC", "--justConcatenate", type="integer",
                    help="Specify whether ASVs need to be concatinated with Ns instead of merging")
parser$add_argument("--bimera", action='store_true', help="Optionally output list of sequences identified as bimeras")
args <- parser$parse_args()

<<<<<<< HEAD
=======
library(dada2)
library(limma)
library(data.table)

>>>>>>> ab1b0ba044df3bf30a708387c47c9c0788a39602
# Universal parameters
work_dir <- args$dir
path_to_meta <- args$path_to_meta
if (file.exists(path_to_meta)) {
  metafile <- fread(path_to_meta, sep = "\t", header=FALSE)
  sample.names <- metafile$V1
  fnFs <- metafile$V2
  fnRs <- metafile$V3
} else {
  stop(paste("metafile",path_to_meta,"not found!"))
}

# obtain/initialize Parameters
# (Universal) Parameters
randomize=TRUE
selfConsist=TRUE
filter = TRUE
matchIDs <- args$matchIDs
# Parameters for merging
justConcatenate <- args$justConcatenate

# DADA2 and Filtering parameters
maxEE <- args$maxEE
trimRight <- args$trimRight
minLen <- args$minLen
truncQ <- args$truncQ
max_consist <- args$max_consist
omega_a <- args$omega_a

if (is.null(matchIDs)||matchIDs == '') {
  matchIDs = TRUE
} else {
  matchIDs = as.logical(as.numeric(args$matchIDs))
}
if (is.null(trimRight)||trimRight == '') {
  trimRight = c(0,0)
} else {
  trimRight <- as.numeric(strsplit(trimRight,',')[[1]])
}
if (is.null(truncQ)||truncQ == '') {
  truncQ = c(5,5)
} else {
  truncQ <- as.numeric(strsplit(truncQ,',')[[1]])
}


if (args$class == "parasite") {
  # Parameters for filtering
  if (is.null(maxEE)||maxEE == '') {
    maxEE = c(5,5)
  } else {
    maxEE <- as.numeric(strsplit(maxEE,',')[[1]])
  }
  if (is.null(minLen)||minLen == '') {
    minLen=30
  } else {
    minLen = as.numeric(minLen)
  }

  # Parameters for Denoising
  if (is.null(max_consist)||max_consist == '') {
    max_consist=10
  } else {
    max_consist = as.numeric(max_consist)
  }
  if (is.null(omega_a)||omega_a == '') {
    omega_a=1e-120
  } else {
    omega_a = as.numeric(omega_a)
  }

  #Parameters for merging
  if (is.null(justConcatenate)||justConcatenate == '') {
    justConcatenate = FALSE
  } else {
    justConcatenate = as.logical(as.numeric(justConcatenate))
  }
  #justConcatenate=FALSE

} else if (args$class == "vector") {
  # Parameters for filtering
  if (is.null(maxEE)||maxEE == '') {
    maxEE = c(2,2)
  } else {
    maxEE <- as.numeric(strsplit(maxEE,',')[[1]])
  }
  if (is.null(minLen)||minLen == '') {
    minLen=75
  } else {
    minLen = as.numeric(minLen)
  }

  # Parameters for Denoising
  if (is.null(max_consist)||max_consist == '') {
    max_consist=20
  } else {
    max_consist = as.numeric(max_consist)
  }
  if (is.null(omega_a)||omega_a == '') {
    omega_a=1e-40
  } else {
    omega_a = as.numeric(omega_a)
  }

  #Parameters for merging
  if (is.null(justConcatenate)||justConcatenate == '') {
    justConcatenate = TRUE
  } else {
    justConcatenate = as.logical(as.numeric(justConcatenate))
  }
  #justConcatenate=TRUE
} else {
  stop("Please provide valid option for the '--class' argument")
}

#Output parameters
if (dirname(args$output_filename) != ".") {
  output_filename <- args$output_filename
  } else {
    output_filename <- paste0(work_dir,"/",args$output_filename)
  }


#Datatable to summarize parmeters
parameter_df <- data.frame(maxEE=maxEE,
		trimRight=trimRight,
  		minLen=minLen,
		truncQ=truncQ,
  		matchIDs=matchIDs,
		max_consist=max_consist,
		randomize=randomize,
		selfConsist=selfConsist,
		OMEGA_A=omega_a,
  		justConcatenate=justConcatenate)

print(parameter_df)

# List files and sample names

if (length(fnFs) == 0 || length(fnFs) != length(fnRs)) {
	stop("fastq files incomplete or not found")
}

# Plot Quality profiles before filering
png(paste0(work_dir,"/qualityF.png"), height = 800, width = 700)
try(print(plotQualityProfile(fnFs[1:2])), silent = TRUE)
dev.off()
png(paste0(work_dir,"/qualityR.png"),height = 800, width = 700)
try(print(plotQualityProfile(fnRs[1:2])), silent = TRUE)
dev.off()

# Create paths for filtered fastq
filtFs <- file.path(work_dir, "filtered", paste0(sample.names, "_filt_R1.fastq.gz"))
filtRs <- file.path(work_dir, "filtered", paste0(sample.names, "_filt_R2.fastq.gz"))
<<<<<<< HEAD

print("DEBUG in runDADA2.R =======")

save(list = ls(all.names = TRUE), file = "currentEnvironment_runDADA2.RData")

#print("samples.names",sample.names)
#print("names(filtFs)",names(filtFs))

=======
>>>>>>> ab1b0ba044df3bf30a708387c47c9c0788a39602
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter read
if (filter == TRUE) {
	print("filtering samples...")
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
            maxN=0, maxEE=maxEE, trimRight=trimRight, truncQ=truncQ, minLen=minLen,
<<<<<<< HEAD
            rm.phix=TRUE, compress=TRUE, multithread=10, verbose=TRUE,
=======
            rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose=TRUE,
>>>>>>> ab1b0ba044df3bf30a708387c47c9c0788a39602
            matchIDs=matchIDs)
	print("filtering done!")
} else {
	print("skipping filter except mandatory removal of N's... ")
	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=c(0,0), maxN=0, rm.phix=TRUE,
<<<<<<< HEAD
            compress=TRUE, multithread=10, verbose=TRUE, matchIDs=matchIDs)
=======
            compress=TRUE, multithread=TRUE, verbose=TRUE, matchIDs=matchIDs)
>>>>>>> ab1b0ba044df3bf30a708387c47c9c0788a39602
}

# Report and Correct for samples with zero reads after filter
zeros <- row.names(out)[out[,2] == 0]
write.table(zeros, paste0(work_dir,"/zeroReadSamples.txt"), sep = "\t", quote = FALSE)
filtFs <- filtFs[out[,2] != 0]
filtRs <- filtRs[out[,2] != 0]
sample.names <- sample.names[out[,2] != 0]

# Update Out table
out <- out[(out[,2] != 0),]

#Compute the error model
print("starting error model learning for forward reads...")
<<<<<<< HEAD
errF <- learnErrors(filtFs, multithread=10, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)
print("starting error model learning for reverse reads...")
errR <- learnErrors(filtRs, multithread=10, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)
=======
errF <- learnErrors(filtFs, multithread=TRUE, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)
print("starting error model learning for reverse reads...")
errR <- learnErrors(filtRs, multithread=TRUE, verbose=2, randomize=randomize, MAX_CONSIST=max_consist)
>>>>>>> ab1b0ba044df3bf30a708387c47c9c0788a39602

#Plot the Errors
png(paste0(work_dir,"/errF.png"), height = 800, width = 700)
try(print(plotErrors(errF, nominalQ=TRUE)), silent = TRUE)
dev.off()
png(paste0(work_dir,"/errR.png"), height = 800, width = 700)
try(print(plotErrors(errR, nominalQ=TRUE)), silent = TRUE)
dev.off()

#DeReplicate the reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Run core DADA2 algorithm
print("starting dada2 for forward reads...")
<<<<<<< HEAD
dadaFs <- dada(derepFs, err=errF, selfConsist=selfConsist, multithread=10, verbose=TRUE, OMEGA_A=omega_a)
print("starting dada2 for reverse reads...")
dadaRs <- dada(derepRs, err=errR, selfConsist=selfConsist, multithread=10, verbose=TRUE, OMEGA_A=omega_a)
=======
dadaFs <- dada(derepFs, err=errF, selfConsist=selfConsist, multithread=TRUE, verbose=TRUE, OMEGA_A=omega_a)
print("starting dada2 for reverse reads...")
dadaRs <- dada(derepRs, err=errR, selfConsist=selfConsist, multithread=TRUE, verbose=TRUE, OMEGA_A=omega_a)
>>>>>>> ab1b0ba044df3bf30a708387c47c9c0788a39602

# Merge reads
print("merging paird ends...")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate=justConcatenate, trimOverhang = TRUE)

#Generate sequence table
print("generating sequence table...")
seqtab <- makeSequenceTable(mergers)
print("Number of sequences in table")
print(dim(seqtab))
# Inspect distribution of sequence lengths
print(table(nchar(getSequences(seqtab))))

#Remove Chimeras
if(args$bimera) {
  print("identifying bimeric sequences...")
<<<<<<< HEAD
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=10, verbose=TRUE)
=======
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
>>>>>>> ab1b0ba044df3bf30a708387c47c9c0788a39602
  print("Number of non-bimeric sequences:")
  print(dim(seqtab.nochim)[2])
  print("Percentage of reads which are non-bimeric:")
  print(sum(seqtab.nochim)/sum(seqtab))
  bimeras <- !(colnames(seqtab) %in% colnames(seqtab.nochim))
  write.table(data.frame(sequence = colnames(seqtab), bimera = bimeras), file=paste0(work_dir,"/ASVBimeras.txt"),
    quote=FALSE, sep="\t", row.names=FALSE)
} else {
  print("skipping Bimera identification..")
}

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# sink summary from stdout to a file
sink(paste0(work_dir,"/reads_summary.txt"))
print(track)
#close sink
sink()

#Show the barplot of length distribution
png(paste0(work_dir,"/sequences_barplot.png"), height = 800, width = 700)
print(barplot(table(nchar(getSequences(seqtab)))))
dev.off()

#Generate output: sequence table to a tsv
write.table(seqtab, file=output_filename, quote = FALSE, sep = "\t")

# Save Run as R workspace image (Optional)
if (is.null(args$save_run)||args$save_run == '') {
    print("--save_run not found or empty. skip saving Rdata image to a file")
  } else {
    rm(fnFs,fnRs,filtFs,filtRs,derepFs,derepRs,errF,errR)
    save.image(paste0(work_dir,"/",args$save_run))
  }
