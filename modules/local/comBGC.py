#!/usr/bin/env python3

print("Running comBGC")

from multiprocessing.connection import wait
import os
import re
import sys
from Bio import SeqIO
import pandas as pd

### TO-DO AntiSMASH parsing

# Function extract_gbk_info:
# - Iterates over list of folders with antiSMASH output.
# - Opens summary GBK and parses the information.
# - Extracts the knownclusterblast output from the antiSMASH folder .
# - Stores everything into a data frame that is then written to TSV file.
def extract_gbk_info(as_dirs, contigs):
	as_results = {}
	for as_dir in as_dirs:
		print("  - " + as_dir)

		Sample_ID = as_dir.rstrip("/") # Assuming folder name equals sample name

		gbk_path = as_dir + ".gbk"
		with open(gbk_path) as gbk:
			for record in SeqIO.parse(gbk, "genbank"):
				CDS_count = 0 # Count the open reading frames of the cluster

				for feature in record.features:

					# Use the first protocluster feature from the contig record file to extract all the infos
					if feature.type == "protocluster":

						Contig_ID = record.id
						BGC_probability = ""
						BGC_start = feature.location.start + 1 # Because zero-based start position
						BGC_end = feature.location.end
						BGC_length = feature.location.end - feature.location.start

						Product_class = feature.qualifiers["product"]
						for i in range(len(Product_class)):
							Product_class[i] = Product_class[i][0].upper() + Product_class[i][1:] # Make first letters uppercase, e.g. lassopeptide -> Lassopeptide

						if feature.qualifiers["contig_edge"] == "True":
							BGC_complete = "No"
						elif feature.qualifiers["contig_edge"] == "False":
							BGC_complete = "Yes"

					# Count functional CDSs (no pseudogenes)
					elif feature.type == "CDS" and "translation" in feature.qualifiers.keys() and feature.location.end <= BGC_end: # Make sure not to count pseudogenes (which would have no "translation tag") and within the current BGC region
						CDS_count += 1

					kcb_file = '{}knownclusterblast/{}_c1.txt'.format(as_dir, record.id) # Prbly need dir list to get file names, bc ..._c2.txt can also occur!
					MIBiG_IDs = parse_knownclusterblast(kcb_file)
					BGC_pos = record.annotations["structured_comment"]

					# Store the all the values for current GBK in a list
					as_results[Sample_ID + gbk] = [Sample_ID,
										Contig_ID,
										Product_class,
										Contig_edge,
										BGC_length,
										str(CDS_count),
										clust_IDs, # IDs and annotations of BGC BLAST matches
										clust_annotations,
										blast_num,
										blast_identity_averages,
										blast_score_cum,
										blast_score_averages,
										blast_coverage_averages,
										blast_cds_annotations]
	return as_results

# Function: Read in the antiSMASH txt file and extract MIBiG IDs
def parse_knownclusterblast(kcb_file):
    with open(kcb_file) as infile:
        hits = 0
        MIBiG_IDs = []

        for line in infile:
            if line == "Significant hits: \n" and not hits:
                hits = 1
            elif line == "\n" and hits:
                break

            if hits:
                MIBiG_ID = re.match("(BGC\d+)", line).group(1)
                MIBiG_IDs.append(MIBiG_ID)

    return MIBiG_IDs

# Dictionary with all the results to be written to disk
print("\nParsing these samples:")
as_results = extract_gbk_info(as_path, dirs, contigs)

# Write the antiSMASH results into TSV (tab-separated values) file
with open(outpath, "w", newline='') as outfile:
	fieldnames = ['Sample_ID', 'Contig_ID', 'Product_class', 'Contig_edge', 'BGC_length', 'CDS_count', 'PFAM_domain', 'Annotation_ID', 'MIBiG_ID', 'MIBiG_ID', 'InterPro_ID', 'Prediction_tool']
	outfile.write("\t".join(fieldnames) + "\n") # Write column names

	for sample in sorted(as_results.keys()): # Sort the samples not to have them in arbitrary order
		out_rows = len(as_results[sample][8]) # Determine if BLAST results were found - if yes: out_rows will be > 0
		if out_rows > 0: # If BLAST results are present: Write each result into a line in the TSV file
			for i in range(out_rows):
				for j in range(len(as_results[sample][4])): # If several products were found in the GBK file, write a line for each one
					outline1 = "\t".join(as_results[sample][:4])
					outline2 = as_results[sample][4][j]
					outline3 = "\t".join(as_results[sample][5:8])
					outline4 = as_results[sample][8][i]
					outline5 = as_results[sample][9][i]
					outline6 = as_results[sample][10][i]
					outline7 = as_results[sample][11][i]
					outline8 = as_results[sample][12][i]
					outline9 = as_results[sample][13][i]
					outline10 = ", ".join(as_results[sample][14][i])
					outfile.write("\t".join([outline1, outline2, outline3, outline4, outline5, outline6, outline7, outline8, outline9, outline10, "-", "antiSMASH"]) + "\n")
		else: # If no BLAST results are present: Write only the BGC into a line in the TSV file
			for j in range(len(as_results[sample][3])):	 # If several products were found in the GBK file, write a line for each one
				outline1 = "\t".join(as_results[sample][:4])
				outline2 = as_results[sample][4][j]
				outline3 = "\t".join(as_results[sample][5:8])
				outline4 = ""
				outline5 = ""
				outline6 = ""
				outline7 = ""
				outline8 = ""
				outline9 = ""
				outline10 = ""
				outline11 = ""
				outfile.write("\t".join([outline1, outline2, outline3, outline4, outline5, outline6, outline7, outline8, outline9, outline10, outline11]) + "\n")

print("\nFinished successfully. Your output is here: " + outpath)

### TO-DO GECCO output parsing


### TO-DO DeepBGC output parsing


### TO-DO combined output TSV
# One line per BGC per tool (i.e. BGCs will be reported multiple times if found by 2 or more tools – users can filter themselves if needed)

# Example output table
# BGC_ID	Sample_ID	Prediction_tool	Contig_ID	Product_class	BGC_probability	BGC_complete	BGC_start	BGC_end	BGC_strand	BGC_length	Protein_count	Protein_ID	PFAM_ID	MIBiG_ID	InterPro_ID
# Sample_1	antiSMASH	c_001	Arylpolyene		yes	123	456	334	+	2	OGCKDNOF_00056;OGCKDNOF_00057	PF00668	BGC0001894
# Sample_1	GECCO	c_002	RiPP	0.96		123	456	334	-	1	OGCKDNOF_00056	PF00668;PF08242		IPR001031
# Sample_2	antiSMASH	c_001	NRPS		two-side-truncated	+	123	456	334	1	OGCKDNOF_00056	PF08242	BGC0001894
# Sample_2	DeepBGC	c_002	Arylpolyene	0.95		123	456	334	-	3	OGCKDNOF_00056;OGCKDNOF_00058;OGCKDNOF_00059	PF00668;PF08242;PF08243	BGC0001894
# I also have the column mappings if you need them (they are quite straightforward, though).

# Call all 3 filtering functions (antiSMASH, DeepBGC, GECCO) + collects data frames
# Join data frames (maybe sort by: 1. Sample_ID, 2. Contig_ID, 3. BGC_start) + output into TSV file
