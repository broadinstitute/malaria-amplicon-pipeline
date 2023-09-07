#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import os, fnmatch

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--path_to_fq', required=True, help="Path to fastq files")
	parser.add_argument('--output_file', required=True, help="Path to output meta file")
	parser.add_argument('--pattern_fw', required=True, help="pattern followed in processed forward fastq for naming")
	parser.add_argument('--pattern_rv', required=True, help="pattern followed in processed reverse fastq for naming")

	args = parser.parse_args()

	pathdir = args.path_to_fq
	odir = args.output_file
	meta_df = pd.DataFrame()
	filelist = os.listdir(pathdir)
	pattern_fw = args.pattern_fw
	pattern_rv = args.pattern_rv

	for entry in filelist:
		if fnmatch.fnmatch(entry, pattern_fw):
			sampleid = entry.split(pattern_fw[1:],1)[0]
			ipath_fw = os.path.join(pathdir, sampleid+pattern_fw[1:])
			ipath_rv = os.path.join(pathdir, sampleid+pattern_rv[1:])
			df = pd.DataFrame({'id' : [sampleid], 'ip1' : [ipath_fw], 'ip2' : [ipath_rv]})
			df = df[['id','ip1','ip2']]
			meta_df = meta_df.append(df, ignore_index=True)
	meta_df.to_csv(odir, sep = "\t", header=False, index=False)
	print( "meta file generated at " + odir )
	return()

if __name__ == "__main__":
	main()
