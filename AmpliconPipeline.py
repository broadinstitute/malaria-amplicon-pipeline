#!/usr/bin/env python

import argparse
import subprocess
import threading
import multiprocessing
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

def trim_primer(sampleid,fileF,fileR,pr1,pr2):
	if os.path.isfile(fileF) and os.path.isfile(fileR):
		proc = subprocess.Popen(['cutadapt', '-g', ('file:'+pr1), '-G', ('file:'+pr2),
			'-o', os.path.join(run_dir,"prim_fq",sampleid+"_prim_1.fq.gz"), 
			'-p', os.path.join(run_dir,"prim_fq",sampleid+"_prim_2.fq.gz"),
			'--pair-adapters','--discard-untrimmed','--action=trim',
			fileF, fileR], stdout=sys.stdout, stderr=sys.stderr)
		proc.wait()
	else:
		sys.exit('Pre-process halted : one or both of the fastq files not found! Exiting..')
	return()

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
	parser.add_argument('--Class', help="Specify Analysis class. Accepts one of two: parasite/vector")
	parser.add_argument('--maxEE', help="Maximum Expected errors (dada2 filtering argument)")
	parser.add_argument('--trimRight', help="Hard trim number of bases at 5` end (dada2 filtering argument)")
	parser.add_argument('--minLen', help="Minimum length filter (dada2 filtering argument)")
	parser.add_argument('--truncQ', help="Soft trim bases based on quality (dada2 filtering argument)")
	parser.add_argument('--max_consist', help="Number of cycles for consistency in error model (dada2 argument)")
	parser.add_argument('--omegaA', help="p-value for the partitioning algorithm (dada2 argument)")
	parser.add_argument('--justConcatenate', help="whether reads should be concatenated with N's during merge (dada2 argument)")
	parser.add_argument('--saveRdata', help="Optionally save dada2 part of this run as Rdata object")

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

	if args.justConcatenate is not None:
		print('NOTE : --justConcatenate argument found. This overrides any json input provided')
		saveRdata = args.justConcatenate
	elif args.json is not None:
		justConcatenate = jsoninputs['justConcatenate']
	else:
		justConcatenate = ''
		print('NOTE : JSON inputs not provided and --justConcatenate argument not found. Default will be used depending on Class')

	if args.skip_preprocess:
		print("skipping Preprocess step..")
		pass
	else:
		# preprocessing part here
		if not os.path.exists(os.path.join(run_dir,"preprocess_fq")):
			os.mkdir(os.path.join(run_dir,"preprocess_fq"))
		else:
			print("Directory %s already exists.." % (os.path.join(run_dir,"preprocess_fq")))

		print("Now running Preprocess..")
		meta = open(path_to_meta,'r')
		samples = meta.readlines()
		p = multiprocessing.Pool()
		for sample in samples:
			slist = sample.split()
			p.apply_async(preprocess, args=(slist[0],slist[1],slist[2]))
		p.close()
		p.join()
		meta.close()

		create_meta(os.path.join(run_dir,"preprocess_fq"),os.path.join(run_dir,"preprocess_meta.txt"),
			pattern_fw="*_val_1.fq.gz", pattern_rv="*_val_2.fq.gz")
		path_to_meta = os.path.join(run_dir,"preprocess_meta.txt")

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
			# Primer removal steps here
			if not os.path.exists(os.path.join(run_dir,"prim_fq")):
				os.mkdir(os.path.join(run_dir,"prim_fq"))
			else:
				print("Directory %s already exists.." % (os.path.join(run_dir,"prim_fq")))

			print("Now running Primer removal..")
			meta = open(path_to_meta,'r')
			samples = meta.readlines()
			p = multiprocessing.Pool()
			for sample in samples:
				slist = sample.split()
				p.apply_async(trim_primer, args=(slist[0],slist[1],slist[2],pr1,pr2))
			p.close()
			p.join()
			meta.close()

			create_meta(os.path.join(run_dir,"prim_fq"),os.path.join(run_dir,"prim_meta.txt"),
				pattern_fw="*_prim_1.fq.gz", pattern_rv="*_prim_2.fq.gz")
			path_to_meta = os.path.join(run_dir,"prim_meta.txt")


	# Steps after Primer removal
	if not os.path.exists(os.path.join(run_dir,"run_dada2")):
		os.mkdir(os.path.join(run_dir,"run_dada2"))
	else:
		print("Directory %s already exists.." % (os.path.join(run_dir,"run_dada2")))

	print("Now running DADA2..")
	cmd = ['Rscript', os.path.join(path,'runDADA2.R'), '-p', path_to_meta, '-d', os.path.join(run_dir,'run_dada2'),
	'-o', 'seqtab.tsv', '-c', Class, '-ee', str(maxEE), '-tR', str(trimRight), '-mL', str(minLen), '-tQ', str(truncQ),
	'-mC', str(max_consist), '-wA', str(omegaA), '-jC', str(justConcatenate), '-s', saveRdata, '--bimera']
	proc = subprocess.Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
	proc.wait()

	print('DADA2 step complete!')

	return()

if __name__ == "__main__":
	main()