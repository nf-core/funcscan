#!/usr/bin/env python3

from multiprocessing.connection import wait
import os
import re
import sys
from Bio import SeqIO
import pandas as pd

print("Running comBGC")

# Function extract_gbk_info:
# - Iterates over list of sample folders with antiSMASH output.
# - Opens summary GBK and parses the information.
# - Extracts the knownclusterblast output from the antiSMASH folder.
# - Stores everything into a data frame which is then written to TSV file.
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

# Function: Read in the knownclusterblast txt file and extract MIBiG IDs
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

# Dictionary with all the results to be written to disk
print("\nParsing these samples:")
antismash_path = sys.argv[1]
as_results = antismash_workflow(antismash_path) # Make sure as_path ends with /

# Write the antiSMASH results into TSV (tab-separated values) file

print("\nFinished successfully.")

### TO-DO GECCO output parsing


### TO-DO DeepBGC output parsing


### TO-DO combined output TSV
# One line per BGC per tool (i.e. BGCs will be reported multiple times if found by 2 or more tools – users can filter themselves if needed)

# Example output table
# BGC_ID	Sample_ID	Prediction_tool	Contig_ID	Product_class	BGC_probability	BGC_complete	BGC_start	BGC_end	BGC_strand	BGC_length	Protein_count	Protein_ID	PFAM_domains	MIBiG_ID	InterPro_ID
# Sample_1	antiSMASH	c_001	Arylpolyene		yes	123	456	334	+	2	OGCKDNOF_00056;OGCKDNOF_00057	PF00668	BGC0001894
# Sample_1	GECCO	c_002	RiPP	0.96		123	456	334	-	1	OGCKDNOF_00056	PF00668;PF08242		IPR001031
# Sample_2	antiSMASH	c_001	NRPS		two-side-truncated	+	123	456	334	1	OGCKDNOF_00056	PF08242	BGC0001894
# Sample_2	DeepBGC	c_002	Arylpolyene	0.95		123	456	334	-	3	OGCKDNOF_00056;OGCKDNOF_00058;OGCKDNOF_00059	PF00668;PF08242;PF08243	BGC0001894
# I also have the column mappings if you need them (they are quite straightforward, though).

# Call all 3 filtering functions (antiSMASH, DeepBGC, GECCO) + collects data frames
# Join data frames (maybe sort by: 1. Sample_ID, 2. Contig_ID, 3. BGC_start) + output into TSV file
