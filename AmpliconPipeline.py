#!/usr/bin/env python

import argparse
import subprocess
import threading
import pandas
import multiprocessing
from Bio import Align
import sys
import os
import glob
import json

def create_meta(path_to_fq,output_file,pattern_fw,pattern_rv):
	proc = subprocess.Popen(['python', os.path.join(path,'create_meta.py'),  '--path_to_fq', path_to_fq,
		'--output_file', output_file, '--pattern_fw', pattern_fw, '--pattern_rv',
		pattern_rv], stdout=sys.stdout, stderr=sys.stderr)
	proc.wait()
	return()

def preprocess(sampleid,fileF,fileR, qvalue = 5, length = 20):
	if os.path.isfile(fileF) and os.path.isfile(fileR):
		proc = subprocess.Popen(['trim_galore', '--paired', '--gzip', '--quality', str(qvalue), '--length', str(length),
			'--output_dir', os.path.join(run_dir,"preprocess_fq"), '--basename', sampleid, fileF, fileR],
			stdout=sys.stdout, stderr=sys.stderr)
		proc.wait()
	else:
		sys.exit('Pre-process halted : one or both of the fastq files not found! Exiting..')
	return()

def Runfastqc(Filelist):
	for file in glob.glob(Filelist):
		if os.path.isfile(file):
			proc = subprocess.Popen(['fastqc', file],stdout=sys.stdout, stderr=sys.stderr)
			proc.wait()
		else:
			print('Warning : Fastq file %s not found! Moving to next..' % file)
	return(True)

def MakeFasta(dataFrame,columnNumber,outputFile):
	ids = dataFrame.iloc[:,0]
	seqs = dataFrame.iloc[:,columnNumber]
	file = open(outputFile,'w+')
	for i in range(len(seqs)):
		header = '>' + ids[i] + ';' + 'min_overlap=' + str(len(seqs[i]))
		file.write(header + "\n" + "^" + seqs[i] + "\n")
	file.close()
	return(outputFile)

## add functionality to have diff Fq1 and Fq2 that match primer info ##
def trim_primer(sampleid,fileF,fileR,pr1,pr2,prefix,keep_untrim=False):
	if os.path.isfile(fileF) and os.path.isfile(fileR):
		if keep_untrim:
			cmd = ['cutadapt', '-g', ('file:'+pr1), '-G', ('file:'+pr2),
			'-o', os.path.join(run_dir,"prim_fq",sampleid+"_"+prefix+"_1.fq.gz"), 
			'-p', os.path.join(run_dir,"prim_fq",sampleid+"_"+prefix+"_2.fq.gz"),
			'--untrimmed-output', os.path.join(run_dir,"prim_fq",sampleid+"_temp_1.fq.gz"),
			'--untrimmed-paired-output', os.path.join(run_dir,"prim_fq",sampleid+"_temp_2.fq.gz"),
			'--pair-adapters','--action=trim',
			fileF, fileR]
		else:
			cmd = ['cutadapt', '-g', ('file:'+pr1), '-G', ('file:'+pr2),
			'-o', os.path.join(run_dir,"prim_fq",sampleid+"_"+prefix+"_1.fq.gz"), 
			'-p', os.path.join(run_dir,"prim_fq",sampleid+"_"+prefix+"_2.fq.gz"),
			'--pair-adapters','--discard-untrimmed','--action=trim',
			fileF, fileR]
		proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
		proc.wait()
	else:
		sys.exit('Pre-process halted : one or both of the fastq files not found! Exiting..')
	return()

def merge_seqtab(path_op,path_nop):

	print('path_op',path_op)
	print('path_nop',path_nop)

	if os.path.isfile(path_op) and os.path.isfile(path_nop):
		seqtab_op = pandas.read_csv(path_op, sep = "\t")
		seqtab_nop = pandas.read_csv(path_nop, sep = "\t")
		seqtab = pandas.concat([seqtab_op,seqtab_nop],axis=1)
	else:
		sys.exit('iseq overlapping and/or non-overlapping dada2 tables not found! Exiting..')

	return(seqtab)

def main():
	global run_dir
	global path
	path = os.path.dirname(__file__)
	print('path',path) 

	parser = argparse.ArgumentParser()
	parser.add_argument('--json', help="Path to json inputs")
	parser.add_argument('--path_to_meta', help="Path to a text file that lists paths to input fastq files (one sample per line")
	parser.add_argument('--skip_preprocess', action="store_true", help="Mention if preprocessing is not needed")
	parser.add_argument('--keep_primers', action="store_true", help="Skip primer removal step")
	parser.add_argument('--pr1', help="Path to forward primers FASTA file")
	parser.add_argument('--pr2', help="Path to reverse primers FASTA file")
	parser.add_argument('--skip_QC', action="store_true", help="Skip FastQC for intermediate FastQs")
	parser.add_argument('--skip_dada2', action="store_true", help="Mention if DADA2 processing is not needed")
	parser.add_argument('--Class', help="Specify Analysis class. Accepts one of two: parasite/vector")
	parser.add_argument('--maxEE', help="Maximum Expected errors (dada2 filtering argument)")
	parser.add_argument('--trimRight', help="Hard trim number of bases at 5` end (dada2 filtering argument)")
	parser.add_argument('--minLen', help="Minimum length filter (dada2 filtering argument)")
	parser.add_argument('--truncQ', help="Soft trim bases based on quality (dada2 filtering argument)")
	parser.add_argument('--max_consist', help="Number of cycles for consistency in error model (dada2 argument)")
	parser.add_argument('--omegaA', help="p-value for the partitioning algorithm (dada2 argument)")
	parser.add_argument('--justConcatenate', help="whether reads should be concatenated with N's during merge (dada2 argument)")
	parser.add_argument('--saveRdata', help="Optionally save dada2 part of this run as Rdata object")
	parser.add_argument('--iseq', action="store_true", help="Specify whether the data is from a paired iseq 2 x 150 bp sequencing run")
	parser.add_argument('--reference', help="Path to reference target sequences (If --iseq flag is set)")
	parser.add_argument('--demux_by_amp', action="store_true", help="Demultiplex reads in each sample by Amplicon (requires pr1 and pr2)")
	parser.add_argument('--overlap_pr1', help="Path to forward primers for shorter overlapping targets FASTA file (For iseq run only)")
	parser.add_argument('--overlap_pr2', help="Path to reverse primers for shorter overlapping targets FASTA file (For iseq run only)")
	parser.add_argument('--inputs', help="Path to a TXT file contnaining amplicon List, primer seqs, fragment length, concat/merge option")

	# Converting args parser for compatibility with JSON
	args = parser.parse_args()
	argparse_dict = vars(args)
	
	print("Current args:", argparse_dict)

	if argparse_dict['json'] is None:
		print('NOTE : JSON input file not provided. Checking for individual arguments')
	elif os.path.isfile(argparse_dict['json']):
		jsonfile = open(argparse_dict['json'],'r')
		jsoninputs = json.load(jsonfile)
		jsonfile.close()

		argparse_dict.update(jsoninputs)

		print("Updated JSON args:", argparse_dict)

	else:
		print('NOTE : JSON input file provided but not found. Using individual arguments')

	if argparse_dict['path_to_meta'] is not None:
		path_to_meta = argparse_dict['path_to_meta']
		if os.path.isfile(path_to_meta):
			run_dir = os.path.abspath(os.path.join(path_to_meta, os.pardir))
			sys.stdout = open((run_dir + "/stdout.txt"),"a")
			sys.stderr = open((run_dir + "/stderr.txt"),"a")
		else:
			sys.exit('Execution halted : Metafile with sample list not found! Exiting..')
	else:
		sys.exit('Execution halted : JSON inputs not provided and --path_to_meta not found! Exiting..')

#### Reading Inputs File ####
	
	if argparse_dict['inputs'] is not None:
		inputs = argparse_dict['inputs']
		if os.path.isfile(inputs):
			inputsDF = pandas.read_csv(inputs, sep="\t")
		else:
			sys.exit('NOTE : inputs file not found!')
	else:
		argparse_dict['keep_primers'] = True
		print('NOTE : inputs file not provided. This will skip primer removal and demux_by_amp options')

#### Adapter removal steps ####

	if argparse_dict['skip_preprocess']:
		print("skipping Preprocess step..")
		pass
	else:
		print("Now running Preprocess..")
		if not os.path.exists(os.path.join(run_dir,"preprocess_fq")):
			os.mkdir(os.path.join(run_dir,"preprocess_fq"))
		else:
			print("Directory %s already exists.." % (os.path.join(run_dir,"preprocess_fq")))
			## If preprocess Fastq directory exsists, consider err/warning message ##
		meta = open(path_to_meta,'r')
		samples = meta.readlines()
		p = multiprocessing.Pool()
		for sample in samples:
			slist = sample.split()
			p.apply_async(preprocess, args=(slist[0],slist[1],slist[2]))
		p.close()
		p.join()

		create_meta(os.path.join(run_dir,"preprocess_fq"),os.path.join(run_dir,"preprocess_meta.txt"),
			pattern_fw="*_val_1.fq.gz", pattern_rv="*_val_2.fq.gz")
		path_to_meta = os.path.join(run_dir,"preprocess_meta.txt")

#### QC on Preprocessed FQs ####

	if argparse_dict['skip_QC']:
		print("skipping FastQC on Preprocessed files..")
		pass
	else:
		print("Now running FastQC on Preprocessed files..")
		QC = Runfastqc(os.path.join(run_dir,"preprocess_fq","*.fq.gz"))
		if QC:
			cmd = ['multiqc', '-o', os.path.join(run_dir,"preprocess_fq"), os.path.join(run_dir,"preprocess_fq")]
			proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
			proc.wait()


#### Primer removal steps ####

	if argparse_dict['keep_primers']:
		print("skipping Primer removal step..")
		pass
	else:
		## Make Primer files from input ##
		pr1 = MakeFasta(inputsDF,1,os.path.join(run_dir,"PrimersF.fasta"))
		pr2 = MakeFasta(inputsDF,2,os.path.join(run_dir,"PrimersR.fasta"))
		
		## Make Primer trimming directory ##
		if os.path.isfile(pr1) and os.path.isfile(pr2):
			print("Now running Primer removal..")	
			if not os.path.exists(os.path.join(run_dir,"prim_fq")):
				os.mkdir(os.path.join(run_dir,"prim_fq"))
			else:
				print("Directory %s already exists.." % (os.path.join(run_dir,"prim_fq")))
			
			## Run Primer Removal ##
			meta = open(path_to_meta,'r')
			samples = meta.readlines()
			if argparse_dict['iseq']:
				## Subset primers that do not need concatenation ##
				concat = inputsDF.iloc[:,-1]
				inputsDF_overlap = inputsDF[~concat]
				overlap_pr1 = MakeFasta(inputsDF_overlap,1,os.path.join(run_dir,"OverlappingPrimersF.fasta"))
				overlap_pr2 = MakeFasta(inputsDF_overlap,2,os.path.join(run_dir,"OverlappingPrimersR.fasta"))

				# Trim primers off Overlapping short targets and write them to different file
				p = multiprocessing.Pool()
				for sample in samples:
					slist = sample.split()
					p.apply_async(trim_primer, args=(slist[0],slist[1],slist[2],overlap_pr1,overlap_pr2,"iseq_op",True))
				p.close()
				p.join()
				
				# Metafile for trimmed overlapping target reads
				create_meta(os.path.join(run_dir,"prim_fq"),os.path.join(run_dir,"iseq_op_prim_meta.txt"),
					pattern_fw="*_iseq_op_1.fq.gz", pattern_rv="*_iseq_op_2.fq.gz")
				# Metafile for un-trimmed non-op target reads
				create_meta(os.path.join(run_dir,"prim_fq"),os.path.join(run_dir,"iseq_temp_meta.txt"),
					pattern_fw="*_temp_1.fq.gz", pattern_rv="*_temp_2.fq.gz")
				temp_meta = open(os.path.join(run_dir,"iseq_temp_meta.txt"),'r')
				samples = temp_meta.readlines()

				# Trim primers off second subset of non-op long targets 
				p = multiprocessing.Pool()
				for sample in samples:
					slist = sample.split()
					p.apply_async(trim_primer, args=(slist[0],slist[1],slist[2],pr1,pr2,"iseq_nop"))
				p.close()
				p.join()
				# Metafile for trimmed non-op target reads
				create_meta(os.path.join(run_dir,"prim_fq"),os.path.join(run_dir,"iseq_nop_prim_meta.txt"),
					pattern_fw="*_iseq_nop_1.fq.gz", pattern_rv="*_iseq_nop_2.fq.gz")
				temp_meta.close()
			## Make demux_by_amp prior to iseq ##
			elif argparse_dict['demux_by_amp']:
				p = multiprocessing.Pool()
				for sample in samples:
					slist = sample.split()
					p.apply_async(trim_primer, args=(slist[0],slist[1],slist[2],pr1,pr2,"{name}"))
				p.close()
				p.join()
			else:	
				p = multiprocessing.Pool()
				for sample in samples:
					slist = sample.split()
					p.apply_async(trim_primer, args=(slist[0],slist[1],slist[2],pr1,pr2,"prim"))
				p.close()
				p.join()

				create_meta(os.path.join(run_dir,"prim_fq"),os.path.join(run_dir,"prim_meta.txt"),
					pattern_fw="*_prim_1.fq.gz", pattern_rv="*_prim_2.fq.gz")
				path_to_meta = os.path.join(run_dir,"prim_meta.txt")
			meta.close()
		else:
			sys.exit("Either set of primer files are missing! Exiting..")

	print("Done with primer removal")

#### DADA2 Steps ####

	if argparse_dict['skip_dada2']:
		print("Skipping dada2 processing..")
		pass
	else:
		if argparse_dict['Class'] is not None:
			#print('NOTE : --Class argument found. This overrides any json input provided')
			Class = argparse_dict['Class']
		else:
			sys.exit('Execution halted : JSON inputs not provided and --Class not found! Exiting..')

		if argparse_dict['maxEE'] is not None:
			#print('NOTE : --maxEE argument found. This overrides any json input provided')
			maxEE = argparse_dict['maxEE']
		else:
			maxEE = ''
			print('NOTE : JSON inputs not provided and --maxEE argument not found. Default will be used depending on the Class')

		if argparse_dict['trimRight'] is not None:
			#print('NOTE : --trimRight argument found. This overrides any json input provided')
			trimRight = argparse_dict['trimRight']
		else:
			trimRight = ''
			print('NOTE : JSON inputs not provided and --trimRight argument not found. Default will be used depending on the Class')

		if argparse_dict['minLen'] is not None:
			minLen = argparse_dict['minLen']
		else:
			minLen = ''
			print('NOTE : JSON inputs not provided and --minLen argument not found. Default will be used depending on the Class')

		if argparse_dict['truncQ'] is not None:
			truncQ = argparse_dict['truncQ']
		else:
			truncQ = ''
			print('NOTE : JSON inputs not provided and --truncQ argument not found. Default will be used depending on the Class')

		if argparse_dict['max_consist'] is not None:
			max_consist = argparse_dict['max_consist']
		else:
			max_consist = ''
			print('NOTE : JSON inputs not provided and --max_consist argument not found. Default will be used depending on the Class')

		if argparse_dict['omegaA'] is not None:
			omegaA = argparse_dict['omegaA']
		else:
			omegaA = ''
			print('NOTE : JSON inputs not provided and --omegaA argument not found. Default will be used depending on Class')

		if argparse_dict['saveRdata'] is not None:
			saveRdata = argparse_dict['saveRdata']
		else:
			saveRdata = ''
			print('NOTE : JSON inputs not provided and --saveRdata argument not found. By Default, DADA2 run is not saved as an Rdata object')
		
		## Check on relevance ##
		if argparse_dict['iseq']:
			print('NOTE : with --iseq enabled, --justConcatenate is irrelevant and ignored')
		elif argparse_dict['justConcatenate'] is not None:
			justConcatenate = argparse_dict['justConcatenate']
		else:
			justConcatenate = ''
			print('NOTE : JSON inputs not provided and --justConcatenate argument not found. Default will be used depending on Class')

		print("Now running DADA2..")

		if not os.path.exists(os.path.join(run_dir,"run_dada2")):
			os.mkdir(os.path.join(run_dir,"run_dada2"))
		else:
			print("Directory %s already exists.." % (os.path.join(run_dir,"run_dada2")))

		## Check for short iseq reads option

		if argparse_dict['iseq']:
			
			if not os.path.exists(os.path.join(run_dir,"run_dada2","dada2_op")):
				os.mkdir(os.path.join(run_dir,"run_dada2","dada2_op"))
			else:
				print("Directory %s already exists.." % (os.path.join(run_dir,"run_dada2","dada2_op")))

			print("Running DADA2 on overlapping targets")
			# Run DADA2 on op targets
			cmdOp = ['Rscript', os.path.join(path,'runDADA2.R'), '-p', os.path.join(run_dir,"iseq_op_prim_meta.txt"),
			'-d', os.path.join(run_dir,'run_dada2','dada2_op'),
			'-o', 'seqtab_op.tsv', '-c', Class, '-ee', str(maxEE), '-tR', str(trimRight), '-mL', str(minLen), '-tQ', str(truncQ),
			'-mC', str(max_consist), '-wA', str(omegaA), '-jC', str(0), '-s', saveRdata, '--bimera']
			procOp = subprocess.Popen(cmdOp, stdout=sys.stdout, stderr=sys.stderr)
			procOp.wait()

			seqtab_op = os.path.join(run_dir,'run_dada2','dada2_op','seqtab_op.tsv')

			if not os.path.exists(os.path.join(run_dir,"run_dada2","dada2_nop")):
				os.mkdir(os.path.join(run_dir,"run_dada2","dada2_nop"))
			else:
				print("Directory %s already exists.." % (os.path.join(run_dir,"run_dada2","dada2_nop")))

			print("Running DADA2 on non-op targets")
			# Run DADA2 on non-op targets
			cmdNOp = ['Rscript', os.path.join(path,'runDADA2.R'), '-p', os.path.join(run_dir,"iseq_nop_prim_meta.txt"),
			'-d', os.path.join(run_dir,'run_dada2','dada2_nop'),
			'-o', 'seqtab_nop.tsv', '-c', Class, '-ee', str(maxEE), '-tR', str(trimRight), '-mL', str(minLen), '-tQ', str(truncQ),
			'-mC', str(max_consist), '-wA', str(omegaA), '-jC', str(1), '-s', saveRdata, '--bimera']
			procNOp = subprocess.Popen(cmdNOp, stdout=sys.stdout, stderr=sys.stderr)
			procNOp.wait()

			# ASV modification block for non-op targets
			seqtab_nop = os.path.join(run_dir,'run_dada2','dada2_nop','seqtab_nop.tsv')
			if argparse_dict['reference'] is not None:

				print('--reference given %s' % (str(argparse_dict['reference'])))
				if os.path.isfile(argparse_dict['reference']):
					adjASV = ['Rscript', os.path.join(path, 'adjustASV.R'), '-s', seqtab_nop, '-ref', str(argparse_dict['reference']),
					'-o', os.path.join(run_dir,'run_dada2','dada2_nop','correctedASV.txt')]

					procASV = subprocess.Popen(adjASV, stdout=sys.stdout, stderr=sys.stderr)
					procASV.wait()
					seqtab_corrected = os.path.join(run_dir,'run_dada2','dada2_nop','seqtab_corrected.tsv')
					seqtab = merge_seqtab(seqtab_op,seqtab_corrected)
				else:
					print('--reference file not found. skipping ASV correction..')
					seqtab = merge_seqtab(seqtab_op,seqtab_nop)
			else:
				print('--reference not given. Skipping ASV correction..')
				seqtab = merge_seqtab(seqtab_op,seqtab_nop)

			# Merge two ASV tables
			seqtab.to_csv(os.path.join(run_dir,'run_dada2','seqtab_iseq.tsv'), sep = "\t")

		## Check for demux by amplicon option
		elif argparse_dict['demux_by_amp']:
			# Loop DADA2 over amplicons
			n = 0
			amplist = open(pr1,'r').readlines()
			amplicons = list(filter(lambda x:'>' in x, amplist))
			if len(str(justConcatenate)) == len(amplicons):
				jc = str(justConcatenate)
			else:
				jc = '0'*len(amplist)
			for amplicon in amplicons:
				a2 = amplicon.rstrip().split('>')[1]
				if not os.path.exists(os.path.join(run_dir,"run_dada2",a2.split(';')[0])):
					os.mkdir(os.path.join(run_dir,"run_dada2",a2.split(';')[0]))
				else:
					print("Directory %s already exists.." % (os.path.join(run_dir,"run_dada2",a2.split(';')[0])))
				path_to_meta = os.path.join(run_dir,a2.split(';')[0]+"_meta.txt")
				create_meta(os.path.join(run_dir,"prim_fq"),path_to_meta,
					pattern_fw="*_"+a2+"_1.fq.gz", pattern_rv="*_"+a2+"_2.fq.gz")
				cmd = ['Rscript', os.path.join(path,'runDADA2.R'), '-p', path_to_meta, '-d', os.path.join(run_dir,"run_dada2",a2.split(';')[0]),
				'-o', a2.split(';')[0]+'_seqtab.tsv', '-c', Class, '-ee', str(maxEE), '-tR', str(trimRight), '-mL', str(minLen), '-tQ', str(truncQ),
				'-mC', str(max_consist), '-wA', str(omegaA), '-jC', jc[n], '-s', saveRdata, '--bimera']
				proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
				proc.wait()
				print(n)
				print(jc[n])
				n += 1

		else:
			cmd = ['Rscript', os.path.join(path,'runDADA2.R'), '-p', path_to_meta, '-d', os.path.join(run_dir,'run_dada2'),
			'-o', 'seqtab.tsv', '-c', Class, '-ee', str(maxEE), '-tR', str(trimRight), '-mL', str(minLen), '-tQ', str(truncQ),
			'-mC', str(max_consist), '-wA', str(omegaA), '-jC', str(justConcatenate), '-s', saveRdata, '--bimera']
			proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
			proc.wait()
		
		print('DADA2 step complete!')
	

	return()

if __name__ == "__main__":
	main()
