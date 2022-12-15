#!/usr/bin/env python3

# Import libraries
from Bio import SeqIO
import pandas as pd
import argparse
import os
import re

# Initialize parser
parser = argparse.ArgumentParser(prog = 'comBGC', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description=('''\
    .............................................................................
                                    *comBGC*
    .............................................................................
                This tool parses the results of BGC prediction tools
            For detailed usage documentation please refer to <github_repo>
    .............................................................................'''),
                                epilog='''Your comBGC summary file is in the specified output folder. Enjoy your day!''',
                                add_help=True)
# Input options
parser.add_argument("--input_antismash", dest="antismash", nargs='?', help="Enter the path to the folder that contains the antiSMASH output in subfolders named by sample name.")
parser.add_argument("--input_deepbgc", dest="deepbgc", nargs='?', help="Enter the path to the folder that contains the DeepBGC output files in subfolders named by sample name.",
                    type=str, default='./deepbgc/')
parser.add_argument("--input_gecco", dest="gecco", nargs='?', help="Enter the path to the folder that contains the GECCO output files in subfolders named by sample name.",
                    type=str, default='./gecco/')

# get command line arguments
args = parser.parse_args()

# assign input arguments to variables
input_antismash = args.antismash
input_gecco = args.gecco
input_deepbgc = args.deepbgc

########################
# ANTISMASH FUNCTIONS
########################

# Extract MIBiG IDs from the knownclusterblast txt file
def parse_knownclusterblast(kcb_file):
    with open(kcb_file) as infile:
        hits = 0
        MIBiG_IDs = []

        for line in infile:
            if line == "Significant hits: \n" and not hits:
                hits = 1
            elif line == "\n" and hits:
                break
            elif line != "Significant hits: \n" and hits:
                MIBiG_ID = re.search("(BGC\d+)", line).group(1)
                MIBiG_IDs.append(MIBiG_ID)
    return MIBiG_IDs

# Create summary file:
# - Iterate over list of sample folders with antiSMASH output.
# - Open summary GBK and parses the information.
# - Extract the knownclusterblast output from the antiSMASH folder.
# - Store everything into a data frame and write to TSV file.
def antismash_workflow(antismash_path):
    summary_cols = ['Sample_ID', 'Prediction_tool', 'Contig_ID', 'Product_class', 'BGC_probability', 'BGC_complete', 'BGC_start', 'BGC_end', 'BGC_length', 'CDS_ID', 'CDS_count', 'PFAM_domains', 'MIBiG_ID', 'InterPro_ID']
    antismash_out = pd.DataFrame(columns=summary_cols)

    for sample in os.listdir(antismash_path):
        print("  - " + sample)

        Sample_ID = sample # Assuming folder name equals sample name

        sample_path = "/".join([antismash_path.rstrip("/"), sample]) + "/"
        CDS_ID = []
        CDS_count = 0

        kcb_path = sample_path + "knownclusterblast/"
        kcb_files = []

        if "knownclusterblast" in os.listdir(sample_path):
            kcb_files = [file for file in os.listdir(kcb_path) if file.startswith("c") and file.endswith(".txt")]

        gbk_path = sample_path + sample + ".gbk"
        with open(gbk_path) as gbk:
            for record in SeqIO.parse(gbk, "genbank"):
                cluster_num = 1
                antismash_out_line = {}
                Contig_ID = record.id
                Product_class = ""
                BGC_complete = ""
                BGC_start = ""
                BGC_end = ""
                BGC_length = ""
                PFAM_domains = []
                MIBiG_ID = "NA"

                for feature in record.features:

                    # Use the first protocluster feature from the contig record file to extract all the infos
                    if feature.type == "protocluster":

                        if antismash_out_line: # If we are having more than 1 BGC per contig, reset the output line for new BGC. Assuming that BGCs do not overlap.
                            antismash_out_line = {
                                'Sample_ID': Sample_ID,
                                'Prediction_tool': "antiSMASH",
                                'Contig_ID': Contig_ID,
                                'Product_class': Product_class,
                                'BGC_probability': "NA",
                                'BGC_complete': BGC_complete,
                                'BGC_start': BGC_start,
                                'BGC_end': BGC_end,
                                'BGC_length': BGC_length,
                                'CDS_ID': ";".join(CDS_ID),
                                'CDS_count': CDS_count,
                                'PFAM_domains': ";".join(PFAM_domains),
                                'MIBiG_ID': MIBiG_ID,
                                'InterPro_ID': "NA"
                            }
                            antismash_out_line = pd.DataFrame([antismash_out_line])
                            antismash_out = pd.concat([antismash_out, antismash_out_line], ignore_index=True)
                            antismash_out_line = {}
                            CDS_ID = []
                            CDS_count = 0
                            PFAM_domains = []

                        Product_class = feature.qualifiers["product"]
                        for i in range(len(Product_class)):
                            Product_class[i] = Product_class[i][0].upper() + Product_class[i][1:] # Make first letters uppercase, e.g. lassopeptide -> Lassopeptide

                        if feature.qualifiers["contig_edge"] == ['True']:
                            BGC_complete = "No"
                        elif feature.qualifiers["contig_edge"] == ['False']:
                            BGC_complete = "Yes"

                        BGC_start = feature.location.start + 1 # +1 because zero-based start position
                        BGC_end = feature.location.end
                        BGC_length = feature.location.end - feature.location.start + 1

                        if kcb_files: # If there are elements in the knownclusterblast file list
                            kcb_file = '{}_c{}.txt'.format(record.id, str(cluster_num)) # Check if this filename is among the knownclusterblast files
                            if kcb_file in kcb_files:
                                MIBiG_ID = ";".join(parse_knownclusterblast(kcb_path + kcb_file))
                                if MIBiG_ID == "":
                                    MIBiG_ID = "NA"
                                cluster_num += 1

                    # Count functional CDSs (no pseudogenes)
                    elif feature.type == "CDS" and "translation" in feature.qualifiers.keys() and BGC_start != "": # Make sure not to count pseudogenes (which would have no "translation tag") and count no CDSs before first BGC
                        if feature.location.end <= BGC_end: # Make sure CDS is within the current BGC region
                            CDS_ID.append(feature.qualifiers["locus_tag"][0])
                            CDS_count += 1
                            if "sec_met_domain" in feature.qualifiers.keys():
                                for PFAM_domain in feature.qualifiers["sec_met_domain"]:
                                    PFAM_domain_name = re.search("(.+) \(E-value", PFAM_domain).group(1)
                                    PFAM_domains.append(PFAM_domain_name)

                antismash_out_line = {
                    'Sample_ID': Sample_ID,
                    'Prediction_tool': "antiSMASH",
                    'Contig_ID': Contig_ID,
                    'Product_class': Product_class,
                    'BGC_probability': "NA",
                    'BGC_complete': BGC_complete,
                    'BGC_start': BGC_start,
                    'BGC_end': BGC_end,
                    'BGC_length': BGC_length,
                    'CDS_ID': ";".join(CDS_ID),
                    'CDS_count': CDS_count,
                    'PFAM_domains': ";".join(PFAM_domains),
                    'MIBiG_ID': MIBiG_ID,
                    'InterPro_ID': "NA"
                }

                if BGC_start != "": # Only keep records with BGCs
                    antismash_out_line = pd.DataFrame([antismash_out_line])
                    antismash_out = pd.concat([antismash_out, antismash_out_line], ignore_index=True)
                    CDS_ID = []
                    CDS_count = 0
                    PFAM_domains = []
    antismash_out.to_csv('comBGC_antiSMASH_summary.tsv', sep="\t", index=False)
    return antismash_out

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
                    'pfam_ids':'PFAM_domains'}
    deepbgc_unused_cols = ['detector_version', 'detector_label', 'bgc_candidate_id', 'num_domains', 'num_bio_domains', 'product_activity', 'antibacterial', 'cytotoxic', 'inhibitor', 'antifungal', 'Alkaloid', 'NRP', 'Other', 'Polyketide', 'RiPP', 'Saccharide', 'Terpene', 'bio_pfam_ids']
    deepbgc_sum_cols = ['Sample_ID', 'Prediction_tool', 'Contig_ID', 'Product_class', 'BGC_probability', 'BGC_complete', 'BGC_start', 'BGC_end', 'BGC_length', 'CDS_ID', 'CDS_count', 'PFAM_domains', 'MIBiG_ID', 'InterPro_ID']
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

def gecco_workflow(gecco_path):
    # GECCO output columns that can be mapped (comBGC:GECCO)
    unused_cols = ['max_p', 'alkaloid_probability', 'polyketide_probability', 'ripp_probability', 'saccharide_probability', 'terpene_probability', 'nrp_probability']
    map_dict = {'sequence_id':'Contig_ID', 'bgc_id':'bgc_id', 'type':'Product_class', 'average_p':'BGC_probability', 'start':'BGC_start', 'end':'BGC_end', 'domains':'PFAM_ID', 'proteins':'CDS_ID'}
    summary_cols = ['Sample_ID', 'Prediction_tool', 'Contig_ID', 'Product_class', 'BGC_probability', 'BGC_complete', 'BGC_start', 'BGC_end', 'BGC_length', 'CDS_ID', 'CDS_count', 'PFAM_domains', 'MIBiG_ID', 'InterPro_ID']

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
# MAIN
########################
if __name__ == "__main__":
    gecco_workflow(input_gecco)
    deepbgc_workflow(input_deepbgc)
    antismash_workflow(input_antismash)
