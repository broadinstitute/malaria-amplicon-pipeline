#!/usr/bin/env python

import argparse
import subprocess
import threading
import pandas
import multiprocessing
from Bio import Align
import sys
import os
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

	parser = argparse.ArgumentParser()
	parser.add_argument('--json', help="Path to json inputs")
	parser.add_argument('--path_to_meta', help="Path to input fastq files")
	parser.add_argument('--skip_preprocess', action="store_true", help="Mention if preprocessing is not needed")
	parser.add_argument('--keep_primers', action="store_true", help="Skip primer removal step")
	parser.add_argument('--pr1', help="Path to forward primers FASTA file")
	parser.add_argument('--pr2', help="Path to reverse primers FASTA file")
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

	args = parser.parse_args()
	if args.json is None:
		print('NOTE : JSON input file not provided. Checking for individual arguments')
	elif os.path.isfile(args.json):
		jsonfile = open(args.json,'r')
		jsoninputs = json.load(jsonfile)
		jsonfile.close()
	else:
		print('NOTE : JSON input file provided but not found. Using individual arguments')

	if args.path_to_meta is not None:
		print('NOTE : --path_to_meta argument found. This overrides any json input provided')
		path_to_meta = args.path_to_meta
	elif args.json is not None:
		path_to_meta = jsoninputs['path_to_meta']
	else:
		sys.exit('Execution halted : JSON inputs not provided and --path_to_meta not found! Exiting..')

	if os.path.isfile(path_to_meta):
		run_dir = os.path.abspath(os.path.join(path_to_meta, os.pardir))
		sys.stdout = open((run_dir + "/stdout.txt"),"a")
		sys.stderr = open((run_dir + "/stderr.txt"),"a")
	else:
		sys.exit('Execution halted : Metafile with sample list not found! Exiting..')

#### Adapter removal steps ####

	if args.skip_preprocess:
		print("skipping Preprocess step..")
		pass
	else:
		print("Now running Preprocess..")
		if not os.path.exists(os.path.join(run_dir,"preprocess_fq")):
			os.mkdir(os.path.join(run_dir,"preprocess_fq"))
		else:
			print("Directory %s already exists.." % (os.path.join(run_dir,"preprocess_fq")))

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

#### Primer removal steps ####

	if args.keep_primers:
		print("skipping Primer removal step..")
		pass
	else:
		if args.pr1 is not None:
			print('NOTE : --pr1 argument found. This overrides any json input provided')
			pr1 = args.pr1
			mark = True
		elif args.json is not None:
			pr1 = jsoninputs['pr1']
			mark = True
		else:
			print('NOTE : JSON inputs not provided and --pr1 argument not found. Skipping primer removal..')
			mark = False

		if args.pr2 is not None:
			print('NOTE : --pr2 argument found. This overrides any json input provided')
			pr2 = args.pr2
			mark = True
		elif args.json is not None:
			pr2 = jsoninputs['pr2']
			mark = True
		else:
			print('NOTE : JSON inputs not provided and --pr2 argument not found. Skipping primer removal..')
			mark = False

		if mark:
			print("Now running Primer removal..")
			if not os.path.exists(os.path.join(run_dir,"prim_fq")):
				os.mkdir(os.path.join(run_dir,"prim_fq"))
			else:
				print("Directory %s already exists.." % (os.path.join(run_dir,"prim_fq")))
			
			meta = open(path_to_meta,'r')
			samples = meta.readlines()
			if args.iseq:
				if args.overlap_pr1 is not None:
					print('NOTE : --overlap_pr1 argument found. This overrides any json input provided')
					overlap_pr1 = args.overlap_pr1
				elif args.json is not None:
					overlap_pr1 = jsoninputs['overlap_pr1']
				else:
					sys.exit('Execution halted: iseq flag set and overlapping Fwd primer not found')

				if args.overlap_pr2 is not None:
					print('NOTE : --overlap_pr2 argument found. This overrides any json input provided')
					overlap_pr2 = args.overlap_pr2
				elif args.json is not None:
					overlap_pr2 = jsoninputs['overlap_pr2']
				else:
					sys.exit('Execution halted: iseq flag set and overlapping Rev primer not found')

				# Trim primers off Overlapping short targets and demux them to different file
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
			elif args.demux_by_amp:
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

#### DADA2 Steps ####

	if args.skip_dada2:
		print("Skipping dada2 processing..")
		pass
	else:
		if args.Class is not None:
			print('NOTE : --Class argument found. This overrides any json input provided')
			Class = args.Class
		elif args.json is not None:
			Class = jsoninputs['Class']
		else:
			sys.exit('Execution halted : JSON inputs not provided and --Class not found! Exiting..')

		if args.maxEE is not None:
			print('NOTE : --maxEE argument found. This overrides any json input provided')
			maxEE = args.maxEE
		elif args.json is not None:
			maxEE = jsoninputs['maxEE']
		else:
			maxEE = ''
			print('NOTE : JSON inputs not provided and --maxEE argument not found. Default will be used depending on the Class')

		if args.trimRight is not None:
			print('NOTE : --trimRight argument found. This overrides any json input provided')
			trimRight = args.trimRight
		elif args.json is not None:
			trimRight = jsoninputs['trimRight']
		else:
			trimRight = ''
			print('NOTE : JSON inputs not provided and --trimRight argument not found. Default will be used depending on the Class')

		if args.minLen is not None:
			print('NOTE : --minLen argument found. This overrides any json input provided')
			minLen = args.minLen
		elif args.json is not None:
			minLen = jsoninputs['minLen']
		else:
			minLen = ''
			print('NOTE : JSON inputs not provided and --minLen argument not found. Default will be used depending on the Class')

		if args.truncQ is not None:
			print('NOTE : --truncQ argument found. This overrides any json input provided')
			truncQ = args.truncQ
		elif args.json is not None:
			truncQ = jsoninputs['truncQ']
		else:
			truncQ = ''
			print('NOTE : JSON inputs not provided and --truncQ argument not found. Default will be used depending on the Class')

		if args.max_consist is not None:
			print('NOTE : --max_consist argument found. This overrides any json input provided')
			max_consist = args.max_consist
		elif args.json is not None:
			max_consist = jsoninputs['max_consist']
		else:
			max_consist = ''
			print('NOTE : JSON inputs not provided and --max_consist argument not found. Default will be used depending on the Class')

		if args.omegaA is not None:
			print('NOTE : --omegaA argument found. This overrides any json input provided')
			omegaA = args.omegaA
		elif args.json is not None:
			omegaA = jsoninputs['omegaA']
		else:
			omegaA = ''
			print('NOTE : JSON inputs not provided and --omegaA argument not found. Default will be used depending on Class')

		if args.saveRdata is not None:
			print('NOTE : --saveRdata argument found. This overrides any json input provided')
			saveRdata = args.saveRdata
		elif args.json is not None:
			saveRdata = jsoninputs['saveRdata']
		else:
			saveRdata = ''
			print('NOTE : JSON inputs not provided and --saveRdata argument not found. By Default, DADA2 run is not saved as an Rdata object')

		if args.iseq:
			print('NOTE : with --iseq enabled, --justConcatenate is irrelavant and ignored')
		elif args.justConcatenate is not None:
			print('NOTE : --justConcatenate argument found. This overrides any json input provided')
			justConcatenate = args.justConcatenate
		elif args.json is not None:
			justConcatenate = jsoninputs['justConcatenate']
		else:
			justConcatenate = ''
			print('NOTE : JSON inputs not provided and --justConcatenate argument not found. Default will be used depending on Class')

		print("Now running DADA2..")
		if not os.path.exists(os.path.join(run_dir,"run_dada2")):
			os.mkdir(os.path.join(run_dir,"run_dada2"))
		else:
			print("Directory %s already exists.." % (os.path.join(run_dir,"run_dada2")))

		## Check for short iseq reads option
		if args.iseq:
			
			if not os.path.exists(os.path.join(run_dir,"run_dada2","dada2_op")):
				os.mkdir(os.path.join(run_dir,"run_dada2","dada2_op"))
			else:
				print("Directory %s already exists.." % (os.path.join(run_dir,"run_dada2","dada2_op")))

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

			# Run DADA2 on non-op targets
			cmdNOp = ['Rscript', os.path.join(path,'runDADA2.R'), '-p', os.path.join(run_dir,"iseq_nop_prim_meta.txt"),
			'-d', os.path.join(run_dir,'run_dada2','dada2_nop'),
			'-o', 'seqtab_nop.tsv', '-c', Class, '-ee', str(maxEE), '-tR', str(trimRight), '-mL', str(minLen), '-tQ', str(truncQ),
			'-mC', str(max_consist), '-wA', str(omegaA), '-jC', str(1), '-s', saveRdata, '--bimera']
			procNOp = subprocess.Popen(cmdNOp, stdout=sys.stdout, stderr=sys.stderr)
			procNOp.wait()

			# ASV modification block for non-op targets
			seqtab_nop = os.path.join(run_dir,'run_dada2','dada2_nop','seqtab_nop.tsv')
			if args.reference is not None:
				print('--reference given %s' % (str(args.reference)))
				if os.path.isfile(args.reference):
					adjASV = ['Rscript', os.path.join(path, 'adjustASV.R'), '-s', seqtab_nop, '-ref', str(args.reference),
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
		elif args.demux_by_amp:
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