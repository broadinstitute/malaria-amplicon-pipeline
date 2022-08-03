#python3
#Make DADA2 output into FASTA
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#Imput Parser

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('dada2_table', help='Table style output of DADA2, ', type=str)
parser.add_argument('-out','--OutName', help='Output name to append to all output FASTA, including / at the end',default='out' ,type=str)
parser.add_argument('-dir','--DirectoryPath', help='Path of the directory to store the output fastas. Default, current',default='' ,type=str)
args = parser.parse_args()

infile = args.dada2_table
append_name=args.OutName
dir_path=args.DirectoryPath

#Load the output into a table
dada_output=pd.read_csv(infile,sep='\t',index_col=0)

#Output a large FASTA file with all the samples
count_id=0

#Add tuples with (record, read_count) structure for sorting
persample_records={sample:[] for sample in dada_output.index}

records_allASVs=[]
fasta_out=dir_path+'AllASV_'+append_name+'.fasta'

asv_ids_dict={}

for asv in dada_output.columns:
    count_id=count_id+1

    #Make the new record for all ASVs
    idrec='ASV_'+append_name+'_'+str(count_id)
    namerec='ASV_'+append_name+'_'+str(count_id)
    newASV = SeqRecord(Seq(asv),id=idrec,name=namerec, description=namerec)
    records_allASVs.append(newASV)

    #Save the sequence to id equivalence
    asv_ids_dict[asv]=idrec

    #Parse in per-sample basis, where only non-zero readcount asvs are saved
    for sample, read_counts in dada_output[asv].iteritems():
        if read_counts > 0:
            namerec='ASV_'+append_name+'_'+str(count_id)+'_Reads-'+str(read_counts)
            newASV = SeqRecord(Seq(asv),id=idrec,name=namerec, description=namerec)
            persample_records[sample].append((newASV, read_counts))
        else:
            continue

#Output all the asv's fasta
with open(fasta_out, 'w') as output_handle:
    SeqIO.write(records_allASVs, output_handle, "fasta")

#Sort the per sample records by read count, and output their per sample FASTA
for sample in persample_records:
    persample_records[sample].sort(key=lambda rec: rec[1], reverse=True)

    sample_name=sample.replace(' ', '-')
    fasta_out=dir_path+sample_name+'_'+append_name+'.fasta'

    records=[rec[0] for rec in persample_records[sample]]
    with open(fasta_out, 'w') as output_handle:
        SeqIO.write(records, output_handle, "fasta")


#Rename the reads table where you replace the ASVs with their ID.
dada_output_renamed=dada_output.rename(columns=asv_ids_dict)

#Add a heat map to the output
#Transpose the matrix
results_trans=dada_output.T
normalized_sum=results_trans/results_trans.sum()

fig_size_y=int(0.25*len(results_trans))
fig_size_x=len(dada_output)
fig=plt.figure(figsize=(fig_size_x,fig_size_y))
sns.set(font_scale = 1)
try:
    ax = sns.heatmap(normalized_sum, annot=results_trans, yticklabels=False, fmt="d", linewidths=.1)
    fig.savefig(append_name+'_heatmapASVs.png')

except ValueError:
    print("Too many samples to produce a plot, inspect manually")

#Output the read counts table
dada_output_renamed.to_csv(dir_path+'ReadCountsperASV_'+append_name+'.tsv', sep='\t')
