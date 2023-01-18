#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import argparse
import os
import re

welcome = '''\
                ........................
                    * comBGC v.0.5 *
                ........................
    This tool aggregates the results of BGC prediction tools:
                antiSMASH, deepBGC, and GECCO
           For detailed usage documentation please
             refer to https://nf-co.re/funcscan
    .........................................................'''

# Initialize parser
parser = argparse.ArgumentParser(prog = 'comBGC', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description=(welcome),
                                add_help=True)

# Input options
parser.add_argument("-a", "--antismash", metavar='PATH', dest="antismash", nargs='?', help="path to the folder that contains antiSMASH output",  type=str, default="")
parser.add_argument('-d', '--deepbgc', metavar='PATH', dest="deepbgc", nargs='?', help="path to the folder that contains DeepBGC output", type=str, default="")
parser.add_argument('-g', '--gecco', metavar='PATH', dest="gecco", nargs='?', help="path to the folder that contains GECCO output", type=str, default="")
parser.add_argument('-o', '--outdir', metavar='PATH', dest="outdir", nargs='?', help="directory for comBGC output. Default: current directory", type=str, default=".")
parser.add_argument("-vv", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-v", "--version", help="print version number and exit", action="store_true")

# Get command line arguments
args = parser.parse_args()

# Assign input arguments to variables
input_antismash = args.antismash
input_deepbgc = args.deepbgc
input_gecco = args.gecco
outdir = args.outdir
verbose = args.verbose
version = args.version

if version:
    exit("comBGC 0.5")

# Make sure that at least one input argument is given
if not (input_antismash or input_gecco or input_deepbgc):
    exit("Please specify at least one input directory (--antismash, --deepbgc, --gecco) or see --help")
if not outdir:
    exit("Please specify an output directory (--outdir) or see --help")


########################
# ANTISMASH FUNCTIONS
########################


def parse_knownclusterblast(kcb_file_path):
    '''
    Extract MIBiG IDs from knownclusterblast TXT file.
    '''

    with open(kcb_file_path) as kcb_file:
        hits = 0
        MIBiG_IDs = []

        for line in kcb_file:
            if line == "Significant hits: \n" and not hits:
                hits = 1 # Indicating that the following lines contain relevant information
            elif line == "\n" and hits:
                break
            elif line != "Significant hits: \n" and hits:
                MIBiG_ID = re.search("(BGC\d+)", line).group(1)
                MIBiG_IDs.append(MIBiG_ID)
    return MIBiG_IDs


def antismash_workflow(antismash_path):
    """
    Create data frame with aggregated antiSMASH output:
    - Iterate over list of sample folders with antiSMASH output.
    - Open summary GBK and grab relevant information.
    - Extract the knownclusterblast output from the antiSMASH folder (MIBiG annotations) if present.
    - Return data frame with aggregated info.
    """

    if verbose:
        print("\nParsing antiSMASH files in " + antismash_path + "\n... ", end="")

    antismash_sum_cols = ['Sample_ID', 'Prediction_tool', 'Contig_ID', 'Product_class', 'BGC_probability', 'BGC_complete', 'BGC_start', 'BGC_end', 'BGC_length', 'CDS_ID', 'CDS_count', 'PFAM_domains', 'MIBiG_ID', 'InterPro_ID']
    antismash_out = pd.DataFrame(columns=antismash_sum_cols)

    # Aggregate information
    Sample_ID = antismash_path.rstrip("/").split("/")[-1] # Assuming folder name equals sample name
    CDS_ID = []
    CDS_count = 0
    kcb_path = os.path.join(antismash_path, "knownclusterblast/")
    kcb_files = []

    if "knownclusterblast" in os.listdir(antismash_path):
        kcb_files = [file for file in os.listdir(kcb_path) if file.startswith("c") and file.endswith(".txt")]

    gbk_file = [gbk for gbk in os.listdir(antismash_path) if gbk.endswith(".gbk") and not gbk.endswith("region001.gbk")][0]
    gbk_path = os.path.join(antismash_path, gbk_file)
    with open(gbk_path) as gbk:
        for record in SeqIO.parse(gbk, "genbank"): # GBK records are contigs in this case

            # Initiate variables per contig
            cluster_num        = 1
            antismash_out_line = {}
            Contig_ID          = record.id
            Product_class      = ""
            BGC_complete       = ""
            BGC_start          = ""
            BGC_end            = ""
            BGC_length         = ""
            PFAM_domains       = []
            MIBiG_ID           = "NA"

            for feature in record.features:

                # Extract relevant infos from the first protocluster feature from the contig record
                if feature.type == "protocluster":

                    if antismash_out_line: # If there is more than 1 BGC per contig, reset the output line for new BGC. Assuming that BGCs do not overlap.
                        antismash_out_line = { # Create dictionary of BGC info
                            'Sample_ID'      : Sample_ID,
                            'Prediction_tool': "antiSMASH",
                            'Contig_ID'      : Contig_ID,
                            'Product_class'  : ";".join(Product_class),
                            'BGC_probability': "NA",
                            'BGC_complete'   : BGC_complete,
                            'BGC_start'      : BGC_start,
                            'BGC_end'        : BGC_end,
                            'BGC_length'     : BGC_length,
                            'CDS_ID'         : ";".join(CDS_ID),
                            'CDS_count'      : CDS_count,
                            'PFAM_domains'   : ";".join(PFAM_domains),
                            'MIBiG_ID'       : MIBiG_ID,
                            'InterPro_ID'    : "NA"
                        }
                        antismash_out_line = pd.DataFrame([antismash_out_line])
                        antismash_out = pd.concat([antismash_out, antismash_out_line], ignore_index=True)
                        antismash_out_line = {}

                        # Reset variables per BGC
                        CDS_ID = []
                        CDS_count = 0
                        PFAM_domains = []

                    # Extract all the BGC info
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

                    # If there are knownclusterblast files for the BGC, get MIBiG IDs of their homologs
                    if kcb_files:
                        kcb_file = '{}_c{}.txt'.format(record.id, str(cluster_num)) # Check if this filename is among the knownclusterblast files
                        if kcb_file in kcb_files:
                            MIBiG_IDs = ";".join(parse_knownclusterblast(os.path.join(kcb_path, kcb_file)))
                            if MIBiG_IDs != "":
                                MIBiG_ID = MIBiG_IDs
                            cluster_num += 1

                # Count functional CDSs (no pseudogenes) and get the PFAM annotation
                elif feature.type == "CDS" and "translation" in feature.qualifiers.keys() and BGC_start != "": # Make sure not to count pseudogenes (which would have no "translation tag") and count no CDSs before first BGC
                    if feature.location.end <= BGC_end: # Make sure CDS is within the current BGC region
                        CDS_ID.append(feature.qualifiers["locus_tag"][0])
                        CDS_count += 1
                        if "sec_met_domain" in feature.qualifiers.keys():
                            for PFAM_domain in feature.qualifiers["sec_met_domain"]:
                                PFAM_domain_name = re.search("(.+) \(E-value", PFAM_domain).group(1)
                                PFAM_domains.append(PFAM_domain_name)

            # Create dictionary of BGC info
            antismash_out_line = {
                'Sample_ID'      : Sample_ID,
                'Prediction_tool': "antiSMASH",
                'Contig_ID'      : Contig_ID,
                'Product_class'  : ";".join(Product_class),
                'BGC_probability': "NA",
                'BGC_complete'   : BGC_complete,
                'BGC_start'      : BGC_start,
                'BGC_end'        : BGC_end,
                'BGC_length'     : BGC_length,
                'CDS_ID'         : ";".join(CDS_ID),
                'CDS_count'      : CDS_count,
                'PFAM_domains'   : ";".join(PFAM_domains),
                'MIBiG_ID'       : MIBiG_ID,
                'InterPro_ID'    : "NA"
            }

            if BGC_start != "": # Only keep records with BGCs
                antismash_out_line = pd.DataFrame([antismash_out_line])
                antismash_out = pd.concat([antismash_out, antismash_out_line], ignore_index=True)

                # Reset variables per BGC
                CDS_ID = []
                CDS_count = 0
                PFAM_domains = []

    if verbose:
        print("Done.")
    return antismash_out


########################
# DEEPBGC FUNCTIONS
########################


def deepbgc_initiate(deepbgc_path):
    '''
    Create dictionary with sample name and corresponding path to deepBGC output TSV.
    '''

    # Go over every sample directory in deepbgc_path
    # Append all files ending with bgc.tsv in a list of lists
    deepbgc_list = []
    for file in os.listdir(deepbgc_path):
        if file.endswith(".bgc.tsv"):
            deepbgc_list.append(os.path.join(deepbgc_path, file))

    # Grab sample names from the list of lists
    sample_names = [os.path.basename(tsv).rsplit(".bgc", 1)[0] for tsv in deepbgc_list]

    # Return dictionary with sample name and corresponding TSV
    deepbgc_dict = {}
    for i in range(0, len(deepbgc_list)):
        deepbgc_dict[sample_names[i]] = [deepbgc_list[i]]
    return deepbgc_dict


def deepbgc_workflow(deepbgc_path):
    '''
    Create data frame with aggregated deepBGC output.
    '''

    if verbose:
        print("\nParsing deepBGC files in " + input_deepbgc + "\n... ", end="")

    # Prepare input and output columns
    deepbgc_map_dict = {'sequence_id'  :'Contig_ID',
                        'nucl_start'   :'BGC_start',
                        'nucl_end'     :'BGC_end',
                        'nucl_length'  :'BGC_length',
                        'num_proteins' :'CDS_count',
                        'deepbgc_score':'BGC_probability',
                        'product_class':'Product_class',
                        'protein_ids'  :'CDS_ID',
                        'pfam_ids'     :'PFAM_domains'}
    deepbgc_sum_cols = ['Sample_ID', 'Prediction_tool', 'Contig_ID', 'Product_class', 'BGC_probability', 'BGC_complete', 'BGC_start', 'BGC_end', 'BGC_length', 'CDS_ID', 'CDS_count', 'PFAM_domains', 'MIBiG_ID', 'InterPro_ID']
    deepbgc_unused_cols = ['detector_version', 'detector_label', 'bgc_candidate_id', 'num_domains', 'num_bio_domains', 'product_activity', 'antibacterial', 'cytotoxic', 'inhibitor', 'antifungal', 'Alkaloid', 'NRP', 'Other', 'Polyketide', 'RiPP', 'Saccharide', 'Terpene', 'bio_pfam_ids']

    # Grab deepBGC sample files and paths
    dict = deepbgc_initiate(deepbgc_path)

    # Initiate dataframe
    deepbgc_out = pd.DataFrame(columns=deepbgc_sum_cols)

    # Add information by sample
    for sample in dict.keys():

        # Add relevant deepBGC output columns per BGC
        deepbgc_df = pd.read_csv(dict[sample][0], sep='\t').drop(deepbgc_unused_cols, axis=1).rename(columns=deepbgc_map_dict)
        deepbgc_df['Sample_ID']       = sample
        deepbgc_df['Prediction_tool'] = "deepBGC"
        deepbgc_df['BGC_complete']    = "NA"
        deepbgc_df['MIBiG_ID']        = "NA"
        deepbgc_df['InterPro_ID']     = "NA"

        # Concatenate data frame to out w/o common index column (e.g. sample_id) due to duplicate row names
        deepbgc_out = pd.concat([deepbgc_out, deepbgc_df], ignore_index=True, sort=False)

    # Return data frame with ordered columns
    deepbgc_out = deepbgc_out[deepbgc_sum_cols]
    if verbose:
        print("Done.")
    return deepbgc_out


########################
# GECCO FUNCTIONS
########################


def getInterProID(gbk_path):
    '''
    Retrieve InterPro IDs from GECCO GBK file.
    '''

    with open(gbk_path) as gbk:
        ip_ids = []
        id_pattern = 'InterPro\:(.*)\"'

        for line in gbk:
            if line.find("InterPro:") != -1:
                new_id = re.search(id_pattern, line).group(1)
                ip_ids.append(new_id)
        ipids_str = ';'.join(map(str, ip_ids))
    return(ipids_str)


def walk_gecco_path(gecco_path):
    '''
    Get list of GECCO GBK files (1 file per BGC).
    '''

    gbks = []
    for file in os.listdir(gecco_path):
        if (file.endswith(('.gbk'))):
            gbks.append(os.path.join(gecco_path, file))
        if (file.endswith(('clusters.tsv'))):
            tsv = os.path.join(gecco_path, file)

    sample_name = tsv.split('/')[-1].split('.')[0]

    # Return dictionary of sample name (key) and [result.tsv and [contig.gbk]] (value)
    bgc_dict = { sample_name: [tsv, gbks] }
    return bgc_dict


def gecco_workflow(gecco_path):
    '''
    Create data frame with aggregated GECCO output.
    '''

    if verbose:
        print("\nParsing GECCO files in " + input_gecco + "\n... ", end="")

    # GECCO output columns that can be mapped (comBGC:GECCO)
    map_dict = {'sequence_id':'Contig_ID',
                'bgc_id'     :'bgc_id',
                'type'       :'Product_class',
                'average_p'  :'BGC_probability',
                'start'      :'BGC_start',
                'end'        :'BGC_end',
                'domains'    :'PFAM_domains',
                'proteins'   :'CDS_ID'}
    summary_cols = ['Sample_ID', 'Prediction_tool', 'Contig_ID', 'Product_class', 'BGC_probability', 'BGC_complete', 'BGC_start', 'BGC_end', 'BGC_length', 'CDS_ID', 'CDS_count', 'PFAM_domains', 'MIBiG_ID', 'InterPro_ID']
    unused_cols = ['max_p', 'alkaloid_probability', 'polyketide_probability', 'ripp_probability', 'saccharide_probability', 'terpene_probability', 'nrp_probability']

    # Dict of Sample names (key) and paths to [tsv, [all contigs.gbk]]
    gecco_dict = walk_gecco_path(gecco_path)

    # Initiate dataframe
    gecco_out = pd.DataFrame(columns=summary_cols)

    # Add information by sample
    for i in gecco_dict.keys():
        gecco_df = pd.read_csv(gecco_dict[i][0], sep='\t').drop(unused_cols, axis=1).rename(columns=map_dict)

        # Fill columns (1 row per BGC)
        gecco_df['Sample_ID']       = i
        gecco_df['BGC_length']      = gecco_df['BGC_end']-gecco_df['BGC_start']
        gecco_df['CDS_count']       = [len(gecco_df['CDS_ID'].iloc[i].split(';')) for i in range(0, gecco_df.shape[0])] # Number of contigs in 'Annotation_ID'
        gecco_df['Prediction_tool'] = 'GECCO'

        # Add column 'InterPro_ID'
        for gbk_path in gecco_dict[i][1]:
            bgc_id = gbk_path.split('/')[-1][0:-4]
            gecco_df.loc[gecco_df['bgc_id'] == bgc_id, 'InterPro_ID'] = getInterProID(gbk_path)

        # Add empty columns with no output from GECCO
        gecco_df['BGC_complete'] = 'NA'
        gecco_df['MIBiG_ID'] = 'NA'
        gecco_out = pd.concat([gecco_out, gecco_df])

        # Fill all empty cells with NA
        for row in range(len(gecco_df['PFAM_domains'])):
            if gecco_out['PFAM_domains'].isnull().values[row]:
                gecco_out.loc[row, 'PFAM_domains'] = "NA"

    # Return data frame with ordered columns
    gecco_out = gecco_out[summary_cols]

    if verbose:
        print("Done.")

    return gecco_out


########################
# MAIN
########################

if __name__ == "__main__":

    tools = {"antiSMASH": input_antismash,
            "deepBGC": input_deepbgc,
            "GECCO": input_gecco}
    tools_provided = {}

    for tool in tools.keys():
        if tools[tool] != "":
            tools_provided[tool] = tools[tool]

    if verbose:
        print(welcome)
        print("\nYou provided directories for: " + ", ".join(tools_provided.keys()))

    # Aggregate BGC information into data frame
    summary_antismash = pd.DataFrame()
    summary_deepbgc = pd.DataFrame()
    summary_gecco = pd.DataFrame()

    for tool in tools_provided.keys():
        if tool == "antiSMASH":
            summary_antismash = antismash_workflow(input_antismash)
        elif tool == "deepBGC":
            summary_deepbgc = deepbgc_workflow(input_deepbgc)
        elif tool == "GECCO":
            summary_gecco = gecco_workflow(input_gecco)

    # Summarize and sort data frame
    summary_all = pd.concat([summary_antismash, summary_deepbgc, summary_gecco])
    summary_all.sort_values(by=["Sample_ID", "Contig_ID", "BGC_start", "BGC_length", "Prediction_tool"], axis=0, inplace=True)

    # Write results to TSV
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    summary_all.to_csv(os.path.join(outdir, 'combgc_summary.tsv'), sep='\t', index=False)
    print("Your BGC summary file is: " + os.path.join(outdir, "combgc_summary.tsv"))
