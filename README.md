# malaria-amplicon-pipeline
Contact in case of issues running the tool:  
Ruchit Panchal  
rpanchal@broadinstitute.org

## Pre-Requisites
### Install Anaconda3
The online documentation on how to install Anaconda 3 is given here:  
https://docs.anaconda.com/anaconda/install/linux/  
Follow your Operating System specific instructions on how to install Anaconda3

### Create conda environment for running the tool
Use the ```environment.yml``` file to create a conda virtual environment
```
conda env create --file environment.yml -p /path/to/env/<name-of-environment>/
```
To activate the conda environment
```
source activate <name-of-environment>
```
A detail description on creating a conda environment is given here: 
https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file

## Flowchart
The flowchart of the pipeline as well as some (active) modifications:
https://lucid.app/lucidchart/invitations/accept/inv_a81a5b95-5516-4d2e-a305-8743ceb503d6

## Inputs :
The tool accepts inputs either independently or combined in the form of a JSON file. The list of inputs are given in the file ```inputs.json```

#### Required Inputs :
1. Path_to_meta : The workflow depends on a list of fastq file names tab separated and first column representing sample ID and successive columns (2 and 3) for full path to the forward and reverse end fastq files. This JSON input is required and must point to the text file having the list of raw fastqs.
2. Class : This input is required and currently can accept only two entries, specific to lowercase ; parasite and vector 

#### Optional JSON Inputs :
1. run_dir : This must point to the working directory. If left blank (with empty quotes “”), it will take the parent directory where the metafile with a list to fastq files are and all successive files will be written into this directory.
2. tool : Path to where the workflow scripts are. This path is used in order to call other dependent modules. This is generally different from the working directory to keep things clean and more organized. If blank (“”) it will default to the current central location where the pipeline is on the server.
3. preprocess : This is a numeric flag. If numeric 1, the script will go ahead with the pre-processing step. If numeric 0, it will skip preprocessing (adapter removal and QC). Default is 1, if kept blank (“”)
4. remove_primer : This is a numeric flag. If numeric 1, the script will remove primers from the data based on the two forward and reverse fasta files as primer sequences. If numeric 0, it will skip the primer removal step.

### Creating Metafile (list of fastqs) :

This is a short python script for writing all fastqs into a tab-separated text file that the pipeline uses for processing. It may be run manually first to initiate the main workflow and later on the tool would automate the process of making next step metafiles.

```
usage: python create_meta.py [-h] --path_to_fq PATH_TO_FQ --output_file OUTPUT_FILE
                      --pattern_fw PATTERN_FW --pattern_rv PATTERN_RV

optional arguments:
  -h, --help            show this help message and exit
  --path_to_fq PATH_TO_FQ
                        Path to fastq files
  --output_file OUTPUT_FILE
                        Path to output meta file
  --pattern_fw PATTERN_FW
                        pattern followed in forward end fastq for naming (*_L001_R1_001.fastq.gz)
  --pattern_rv PATTERN_RV
                        pattern followed in reverse end fastq for naming (*_L001_R2_001.fastq.gz)
```

**_NOTE_** : _It is highly recommended to create the fastq metafile in the working directory of the new run. By default the pipeline considers the parent directory of the metafile as the standard working directory and keep the analysis organized. Successive files and directories are written into the working directory._

## Run the tool
If using the JSON input file, the main script can be initiated as given
```
python AmpliconPipeline.py --json inputs.json
```
If however one wishes to supply the inputs independently without requiring a JSON input file, one could do so in which case each input would have its own flag. The entire list of input flags with their description are given below.
```
usage: AmpliconPipeline.py [-h] [--json JSON] [--path_to_meta PATH_TO_META]
                           [--skip_preprocess] [--keep_primers] [--pr1 PR1]
                           [--pr2 PR2] [--skip_dada2] [--Class CLASS]
                           [--maxEE MAXEE] [--trimRight TRIMRIGHT]
                           [--minLen MINLEN] [--truncQ TRUNCQ]
                           [--max_consist MAX_CONSIST] [--omegaA OMEGAA]
                           [--justConcatenate JUSTCONCATENATE]
                           [--saveRdata SAVERDATA] [--iseq]
                           [--reference REFERENCE] [--demux_by_amp]
                           [--overlap_pr1 OVERLAP_PR1]
                           [--overlap_pr2 OVERLAP_PR2]

optional arguments:
  -h, --help            show this help message and exit
  --json JSON           Path to json inputs
  --path_to_meta PATH_TO_META
                        Path to input fastq files
  --skip_preprocess     Mention if preprocessing is not needed
  --keep_primers        Skip primer removal step
  --pr1 PR1             Path to forward primers FASTA file
  --pr2 PR2             Path to reverse primers FASTA file
  --skip_dada2          Mention if DADA2 processing is not needed
  --Class CLASS         Specify Analysis class. Accepts one of two:
                        parasite/vector
  --maxEE MAXEE         Maximum Expected errors (dada2 filtering argument)
  --trimRight TRIMRIGHT
                        Hard trim number of bases at 5` end (dada2 filtering
                        argument)
  --minLen MINLEN       Minimum length filter (dada2 filtering argument)
  --truncQ TRUNCQ       Soft trim bases based on quality (dada2 filtering
                        argument)
  --max_consist MAX_CONSIST
                        Number of cycles for consistency in error model (dada2
                        argument)
  --omegaA OMEGAA       p-value for the partitioning algorithm (dada2
                        argument)
  --justConcatenate JUSTCONCATENATE
                        whether reads should be concatenated with N's during
                        merge (dada2 argument)
  --saveRdata SAVERDATA
                        Optionally save dada2 part of this run as Rdata object
  --iseq                Specify whether the data is from a paired iseq 2 x 150
                        bp sequencing run
  --reference REFERENCE
                        Path to reference target sequences (If --iseq flag is
                        set)
  --demux_by_amp        Demultiplex reads in each sample by Amplicon (requires
                        pr1 and pr2)
  --overlap_pr1 OVERLAP_PR1
                        Path to forward primers for shorter overlapping
                        targets FASTA file (For iseq run only)
  --overlap_pr2 OVERLAP_PR2
                        Path to reverse primers for shorter overlapping
                        targets FASTA file (For iseq run only)
```
## De-multiplex by amplicon
If the amplicon list is relatively small (eg. < 10 ), it is possible to demultiplex reads by amplicon via ```--demux_by_amp``` option. This will produce separate directories and paired fastq files for each amplicon based on reads matching with primers (primers will be removed from the reads) and individual dada2 runs on each amplicon. Note that for larger amplicon panels the dada2 runs are linear and therefore slower and computationally expensive hence it is recommended for panels fewer than 10 amplicons. The inputs that go along are as follows.
```
--demux_by_amp : Logical Flag. Enables workflow to demultiplex reads by amplicon.
--justConcatenate : This option will now allow a string of 0's and 1's to indicate which amplicons (following the same order of primers provided) are to be merged (> 12 bp) versus which are to be contatenated with N's following the DADA2 step. NOTE that if left with a single 0 or 1 would imply that all amplicons would follow that option. Unequal number of 0's and 1's than the primers would result in error.
```

## Handling iSeq (2x150bp) short reads
There is an additional step added in for reads shorter than an illumina Miseq run (2x250 bp). Since there is a possibility of shorter reads (Eg. Illumina iseq 2x150 bp) not covering the entire amplicon length (0 < x < 12 bp overlap between forward and reverse reads) for a subset of amplicons, one could use the ```--iseq``` option which would essentially split the files into reads corresponding to amplicons which have sufficent overlap versus reads corresponding to amplicons with little to no overlap. This would require the user to provide additional inputs as follows
```
--iseq : Logical flag. Enables workflow to bifurcate files to handle shorter read length which may include uncovered region of an amplicon.
--reference : Path to reference target sequences in FASTA format. Whole genome FASTA is currently not supported.
--overlap_pr1 : A subset of forward primers of targets which conntain sufficent (> 12 bp) overlap between forward and reverse ends.
--overlap_pr2 : A subset of reverse primers of targets which contain sufficient (> 12 bp) overlap between forward and reverse ends. 
```
**_NOTE_** : _Use of ```--demux_by_amp``` option along with ```--iseq``` option is currently not supported._

## Post-DADA2 Filters (optional processing parasite only) :  
There is an additional semi-workflow for Post-processing the obtained Amplicon Sequence Variant (ASV) output from the main workflow. This step is intended to be a follow-up step if the target amplicons are from the Parasite genome and for a larger set (> 10 amplicons). This step assumes ```--demux_by_amp``` is not used. Briefly, this step will map the given ASV sequences to the target amplicons while keeping track of non-matching sequences with the number of nucleotide differences and insertions/deletions. It will then output a table of ASV sequences with the necessary information. Optionally, a FASTA file can be created in addition to the table output, listing the sequences in standard FASTA format. A filter tag can be provided to tag the sequences above certain nucleotide (SNV) and length differences due to INDELs and a bimera column to tag sequences which are bimeric (a hybrid of two sequences). A complete list of inputs given below.  
```
usage: Rscript postProc_dada2.R [-h] [-s SEQTAB] [-ref REFERENCE]
                                         [-ref2 REFERENCE2] [-b BIMERA]
                                         [-o OUTPUT] [--fasta]
                                         [-snv SNV_FILTER]
                                         [--indel_filter INDEL_FILTER]
                                         [--strain STRAIN] [--strain2 STRAIN2]
                                         [--parallel]

optional arguments:
  -h, --help            show this help message and exit
  -s SEQTAB, --seqtab SEQTAB
                        Path to input
  -ref REFERENCE, --reference REFERENCE
                        Path to reference fasta sequences
  -ref2 REFERENCE2, --reference2 REFERENCE2
  -b BIMERA, --bimera BIMERA
                        ASV File with identifed bimeras
  -o OUTPUT, --output OUTPUT
                        Path to output file
  --fasta               Write ASV sequences separately into fasta file
  -snv SNV_FILTER, --snv_filter SNV_FILTER
                        Path to file for filtering ASVs based on edit distance
  --indel_filter INDEL_FILTER
                        Specify proportion of ASV length (between 0 and 1) to
                        target length for filtering based on indels
  --strain STRAIN       Name of Specific strain to map to. Defaults to 3D7
  --strain2 STRAIN2     Name of second strain if mapping to 2 different
                        Strains (Required if --ref2 mentioned)
  --parallel            Enable parallel processing
```
## ASV to CIGAR / Variant Calling (optional processing) :
This Step is also an additional step concerning Parasite target amplicons. This will change the representation of the ASV sequences in a kind of pseudo-CIGAR string format, while also masking homopolymer runs, low complexity runs and filtering out sequences tagged in the previous step. This step assumes Post-DADA2 step is used as it requires mapped ASV table as one of the inputs.

#### General idea:
1. Parse DADA2 pipeline outputs to get ASVs → amplicon target.  
   a. Fasta file of ASV sequences.  
   b. ASV → amplicon table from the previous step.  
   c. seqtab tsv file from DADA2.  
2. Build multi-fasta file for each amplicon target containing one or more ASVs
3. Run Muscle on each fasta file to generate alignment file (*.msa)
4. Parse alignments per amplicon target, masking on polyN homopolymer runs (and optionally DUST-masker low complexity sequences)
5. Output seqtab read count table, but the columns are amplicon,CIGAR instead of ASV sequence
6. Optional output ASV to amplicon + CIGAR string table (--asv_to_cigar)

An ASV matching perfectly to the 3D7 reference is currently indicated by “.”  
A complete list of inputs for running this step is given below
```
usage: ASV_to_CIGAR.py [options] fasta table alignments out

Convert ASVs from DADA2 pipeline to pseudo-CIGAR strings.

positional arguments:
  fasta             	Fasta file of ASV sequences from DADA2 pipeline
  table             	ASV table from DADA2 pipeline
  seqtab            	DADA2 seqtab tsv file
  out               	Output seqtab tsv file with amplicon/variant counts

optional arguments:
  -h, --help        	show this help message and exit
  --asv_to_cigar ASV_TO_CIGAR
                    	Output file for ASV -> CIGAR string table
  -a ALIGNMENTS, --alignments ALIGNMENTS
                    	Directory to store ASV alignment files (default: alignments)
  -p POLYN, --polyN POLYN
                    	Mask homopolymer runs length >= polyN (default: 5; disabled < 2)
  -r MIN_READS, --min_reads MIN_READS
                    	Minimum total reads to include ASV (default: 0, disabled)
  -n MIN_SAMPLES, --min_samples MIN_SAMPLES
                    	Minimum samples to include ASV (default: 0, disabled)
  -f, --include_failed  INCLUDE ASVs that failed post-DADA2 filters (default: False)
  -b, --exclude_bimeras
                    	EXCLUDE ASVs that DADA2 flagged as bimeras (default: False)
  -s MAX_SNV_DIST, --max_snv_dist MAX_SNV_DIST
                    	Maximum SNV distance to include ASV (default: -1, disabled)
  -i MAX_INDEL_DIST, --max_indel_dist MAX_INDEL_DIST
                    	Maximum indel distance to include ASV (default: -1, disabled)
  -d AMP_DB, --amp_db AMP_DB
                    	Amplicon sequence fasta file (default: /gsap/garage-protistvector/ampseq_data/AmpSeQC/amplicons_noprimers.fasta)
  -m AMP_MASK, --amp_mask AMP_MASK
                    	Amplicon low complexity mask info (default: None, disabled)
  -v, --verbose     	Increase verbosity

(C)2021 Broad Institute
```

