# malaria-amplicon-pipeline

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

## Inputs :
The tool accepts inputs either independently or combined in the form of a JSON file. The list of inputs are given in the file ```inputs.json```

#### Required Inputs :
1. Path_to_meta : The workflow depends on a list of fastq file names tab separated and first column representing sample ID and successive columns (2 and 3) for full path to the forward and reverse end fastq files. This JSON input is required and must point to the text file having the list of raw fastqs.
2. Class : This input is required and currently can accept only two entries, specific to lowercase ; parasite and vector 

#### Optional Inputs :
1. run_dir : This must point to the working directory. If left blank (with empty quotes “”), it will take the parent directory where the metafile with a list to fastq files are and all successive files will be written into this directory.
2. tool : Path to where the workflow scripts are. This path is used in order to call other dependent modules. This is generally different from the working directory to keep things clean and more organized. If blank (“”) it will default to the current central location where the pipeline is on the server.
3. preprocess : This is a numeric flag. If numeric 1, the script will go ahead with the pre-processing step. If numeric 0, it will skip preprocessing (adapter removal and QC). Default is 1, if kept blank (“”)
remove_primer : This is a numeric flag. If numeric 1, the script will remove primers from the data based on the two forward and reverse fasta files as primer sequences. If numeric 0, it will skip the primer removal step.

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

