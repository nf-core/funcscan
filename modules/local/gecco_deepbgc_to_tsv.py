#!/usr/bin/env python3

# Import libraries
from Bio import SeqIO
import pandas as pd
import argparse
import os
import re

# Initialize parser
parser = argparse.ArgumentParser(prog = 'combgc', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description=('''\
    .............................................................................
                                    *comBGC*
    .............................................................................
                This tool parses the results of BGC prediction tools
            For detailed usage documentation please refer to <github_repo>
    .............................................................................'''),
                                epilog='''Your comBGC-summary file is in the output folder. Enjoy your day!''',
                                add_help=True)
# Input options
parser.add_argument("--input_gecco", dest="gecco", nargs='?', help="Enter the path to the folder that contains the different tool's output files in sub-folders named by sample name.",
                    type=str, default='./gecco/')
parser.add_argument("--input_deepbgc", dest="deepbgc", nargs='?', help="Enter the path to the folder that contains the different tool's output files in sub-folders named by sample name.",
                    type=str, default='./deepbgc/')

# get command line arguments
args = parser.parse_args()

# assign input arguments to variables
input_gecco = args.gecco
input_deepbgc = args.deepbgc

########################
# GECCO FUNCTIONS
########################
# retrieve InterPro ID and contig length from gbk file
def getInterProID(gbk):
    with open(gbk, 'r') as fp:
        # get InterPro IDs
        ip_id = []
        id_pattern = 'InterPro\:(.*)\"'
        lines = fp.readlines()
        for row in lines:
            if row.find("InterPro:") != -1:
                new_id = re.search(id_pattern, row).group(1)
                ip_id.append(new_id)
        ipid_str = ';'.join(map(str, ip_id))
    return(ipid_str)

def walk_gecco_path(path):
    gbk_lst = []
    tsv_lst = []
    for dirpath, subdirs, files in os.walk(path):
        gbk_sub = []
        for file in files:
            if (file.endswith(('.gbk'))):
                gbk_sub.append(os.path.join(dirpath, file))
            if (file.endswith(('clusters.tsv'))):
                tsv_lst.append(os.path.join(dirpath, file))
        gbk_lst.append(gbk_sub)
    # remove empty lists (may be caused by .DS_Store)
    gbk_lst = [x for x in gbk_lst if x]
    sample_names = [x.split('/')[-1] for x in tsv_lst]
    sample_names = [x.split('.')[0] for x in sample_names]
    # create a dictionary of samplename, result.tsv and contig.gbf
    bgc_dict = {}
    for i in range(0, len(gbk_lst)):
        bgc_dict[sample_names[i]] = [tsv_lst[i],gbk_lst[i]]
    return bgc_dict

# Retrieve contig length from input fasta
# transform fasta to dataframe with three columns: contig_id, sequence, sequence length
# def fasta2table(fasta):
#     #read the fasta with SeqIO
#     fasta_seq = SeqIO.parse(open(fasta), 'fasta')
#     #initiate the dataframe containing contig ids sequence and sequence length in three columns
#     fasta_df = pd.DataFrame(columns=['contig_id', 'sequence', 'length'])
#     #append contig information to df
#     for contig in fasta_seq:
#         contig_id, sequence = contig.id, str(contig.seq)
#         fasta_df = fasta_df.append({'contig_id':contig_id, 'sequence':sequence, 'length':len(sequence)}, ignore_index=True)
#     return fasta_df

# Function: Determine on which side of the contig the BGC stops. Region of 50 base pairs on each end counts as truncated
# def edge_position_gecco(df, startcol, endcol, contig, fasta):
#     fasta_df = fasta2table(fasta)
#     contig_length = fasta_df['length'][fasta_df['contig_id']==contig]
#     contig_idx = df.index[df['sequence_id']==contig]
#     bgc_start = int(df[startcol].loc[contig_idx]) # BGC start position in contig from output tsv
#     bgc_end = int(df[endcol].loc[contig_idx]) # BGC stop position in contig from output tsv

#     try:
#         if bgc_start <= 50 and bgc_end >= contig_length - 50: # Contig edge on both sides
#             side = "left, right"
#         elif bgc_start > 50 and bgc_end < contig_length - 50: # No contig edge
#             side = "none"
#         elif bgc_start <= 50 and bgc_end < contig_length - 50: # Contig edge on the left
#             side = "left"
#         elif bgc_start > 50 and bgc_end >= contig_length - 50: # Contig edge on the right
#             side = "right"
#     except KeyError: # Return contig IDs for which no entry in the prokka file was found
#         side = "NA"
#     return side

def gecco_workflow(gecco_path):
    # GECCO output columns that can be mapped (comBGC:GECCO)
    unused_cols = ['max_p', 'alkaloid_probability', 'polyketide_probability', 'ripp_probability', 'saccharide_probability', 'terpene_probability', 'nrp_probability']
    map_dict = {'sequence_id':'Contig_ID', 'bgc_id':'bgc_id', 'type':'Product_class', 'average_p':'BGC_probability', 'start':'BGC_start', 'end':'BGC_end', 'domains':'PFAM_ID', 'proteins':'CDS_ID'}
    summary_cols = ['Sample_ID', 'Prediction_tool', 'Contig_ID', 'Product_class', 'BGC_probability', 'BGC_complete', 'BGC_start', 'BGC_end', 'BGC_length', 'CDS_ID', 'CDS_count', 'PFAM_ID', 'MIBiG_ID', 'InterPro_ID']

    # Dict of Sample names (key) and paths to [tsv, [all contigs.gbk]]
    gecco_dict = walk_gecco_path(gecco_path)
    # Initiate dataframe
    gecco_out = pd.DataFrame(columns=summary_cols)
    # add information by sample in gecco_dict
    for i in gecco_dict.keys():
        gecco_df = pd.read_csv(gecco_dict[i][0], sep='\t').drop(unused_cols, axis=1).rename(columns=map_dict)
        # Add column 'Sample_ID'
        gecco_df['Sample_ID'] = i
        # Add column 'BGC_length'
        gecco_df['BGC_length'] = gecco_df['BGC_end']-gecco_df['BGC_start']
        # Add column 'CDS_count' (number of contigs in 'Annotation_ID')
        gecco_df['CDS_count'] = [len(gecco_df['CDS_ID'].iloc[i].split(';')) for i in range(0, gecco_df.shape[0])]
        # Add column 'Prediction_tool'
        gecco_df['Prediction_tool'] = 'GECCO'
        # Add column 'InterPro_ID'
        for gtf_path in gecco_dict[i][1]:
            bgc_id = gtf_path.split('/')[-1][0:-4]
            gecco_df.loc[gecco_df['bgc_id']==bgc_id, 'InterPro_ID'] = getInterProID(gtf_path)
        # Add empty columns which have no output in GECCO
        gecco_df['BGC_complete'] = ''
        gecco_df['MIBiG_ID'] = ''
        gecco_out = pd.concat([gecco_out, gecco_df])
    # sort columns in the specified order (summary_cols)
    gecco_out = gecco_out[summary_cols]
    gecco_out.to_csv('comBGC_GECCO_summary.tsv', sep='\t', index=False)
    return gecco_out

########################
# DEEPBGC FUNCTIONS
########################
def deepbgc_initiate(deepbgc_path):
    # Go over every subdir in main dir
    # Append all files with bgc.tsv ending in a list of lists 
    deepbgc_list = []
    for root, dirs, files  in os.walk(deepbgc_path):
        path_tsv = os.path.abspath(deepbgc_path)
        for file in files:
            if file.endswith(".bgc.tsv"):
                deepbgc_list.append(os.path.abspath(os.path.join(root, file)))
    # Grab sample names from the list of lists
    sample_names = [os.path.basename(x).rsplit( ".bgc", 1 )[ 0 ] for x in deepbgc_list]
    # Create a dictionary for sample name and sample tsv
    deepbgc_dict = {}
    for i in range(0, len(deepbgc_list)):
        deepbgc_dict[sample_names[i]] = [deepbgc_list[i]]
    return deepbgc_dict

def deepbgc_workflow(deepbgc_path):
    # Prepare input and output columns
    deepbgc_map_dict = {'sequence_id':'Contig_ID', 
                    'detector':'Prediction_tool', 
                	'nucl_start':'BGC_start', 
                	'nucl_end':'BGC_end', 
                	'nucl_length':'BGC_length',
                 	'num_proteins':'CDS_count',
                  	'deepbgc_score':'BGC_probability',
					'product_class':'Product_class',
                   	'protein_ids':'CDS_ID',
                    'pfam_ids':'PFAM_ID'}
    deepbgc_unused_cols = ['detector_version', 'detector_label', 'bgc_candidate_id', 'num_domains', 'num_bio_domains', 'product_activity', 'antibacterial', 'cytotoxic', 'inhibitor', 'antifungal', 'Alkaloid', 'NRP', 'Other', 'Polyketide', 'RiPP', 'Saccharide', 'Terpene', 'bio_pfam_ids']
    deepbgc_sum_cols = ['Sample_ID', 'Prediction_tool', 'Contig_ID', 'Product_class', 'BGC_probability', 'BGC_complete', 'BGC_start', 'BGC_end', 'BGC_length', 'CDS_ID', 'CDS_count', 'PFAM_ID', 'MIBiG_ID', 'InterPro_ID']
    # Grab the deepbgc dict
    dict = deepbgc_initiate(deepbgc_path)
    # initiate dataframe
    deepbgc_out = pd.DataFrame(columns=deepbgc_sum_cols)
    # Add information by sample in dict
    for i in dict.keys():
        # Open the file-rename headers according to dict-remove unused columns
        deepbgc_df = pd.read_csv(dict[i][0], sep='\t').drop(deepbgc_unused_cols, axis=1).rename(columns=deepbgc_map_dict)
        # Add column 'Sample_ID'
        deepbgc_df['Sample_ID'] = i
        # Concatenate df to out w/o common index column (e.g. sample_id) due to duplicate row names
        deepbgc_out = pd.concat([deepbgc_out, deepbgc_df], ignore_index=True, sort=False)
    # Reordre the columns according to sum_cols
    deepbgc_out = deepbgc_out[deepbgc_sum_cols]
    # Save it to a csv
    deepbgc_out.to_csv('comBGC_deepbgc_summary.tsv', sep='\t', index=False)
    return deepbgc_out                 

########################
# MAIN
########################
if __name__ == "__main__":
    gecco_workflow(input_gecco)
    deepbgc_workflow(input_deepbgc)