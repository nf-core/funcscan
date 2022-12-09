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

					infile = '{}knownclusterblast/{}_c1.txt'.format(as_dir, record.id)
					clust_IDs, clust_annotations, blast_num, blast_identity_averages, blast_score_cum, blast_score_averages, blast_coverage_averages, blast_cds_annotations = parse_cluster(infile)
					BGC_pos = record.annotations["structured_comment"]
					Contig_edge, error_id = edge_position(contigs, antismash_output_folder, record.id, BGC_pos) # Get the side of the contig edge (left/right/both/none) #### Buth this into as_results and TSV output

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

# Function: Determine on which side of the contig the BGC stops. Region of 50 base pairs on each end counts as truncated
def edge_position(contigs, sample, id, bgc_pos):
	bgc_start = int(bgc_pos["antiSMASH-Data"]["Orig. start"]) # BGC start position in contig
	bgc_end = int(bgc_pos["antiSMASH-Data"]["Orig. end"]) # BGC stop position in contig

	try:
		if bgc_start <= 50 and bgc_end >= contigs[sample][id] - 50: # Contig edge on both sides
			side = "left, right"
		elif bgc_start > 50 and bgc_end < contigs[sample][id] - 50: # No contig edge
			side = "no"
		elif bgc_start <= 50 and bgc_end < contigs[sample][id] - 50: # Contig edge on the left
			side = "left"
		elif bgc_start > 50 and bgc_end >= contigs[sample][id] - 50: # Contig edge on the right
			side = "right"
		error = ""
	except KeyError: # Return contig IDs for which no entry in the prokka file was found
		side = "NA"
		error = id

	return side, error

# Function: Read in the antiSMASH txt file and split it into the records (divided by ">>")
def file_to_list(file_path):
	hits = open(file_path).read()

	# Divide the hits into a list of lists:
	# First list level: All the records divided by ">>"
		# Second list level: Metadata at the top of each record, genes with annotations, significant BLAST hits
		# Final list looks like this: [[Metadata1, genes with annotation1, BLAST hits1], [Metadata2, genes with annotation2, BLAST hits2], [...]]
	hits = hits.split("\n\n>>\n")[1:] # Split by record
	for i in range(len(hits)): # Split each record into metadata, annotations and significant BLAST hits
		hits[i] = hits[i].split("\nTable of genes, locations, strands and annotations of subject cluster:\n")
		metadata = hits[i][0].rstrip().split("\n")

		annotation_data = hits[i][1].split("\nTable of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):\n")
		annotations = annotation_data[0].rstrip().split("\n")
		annotations = [gene.split("\t") for gene in annotations]

		blast_hits = annotation_data[1].rstrip().split("\n")
		blast_hits = [hit.split("\t") for hit in blast_hits]
		hits[i] = [metadata, annotations, blast_hits]
	return[hits]

# Function: Parse the antiSMASH KnownclusterBLAST txt files (MIBiG IDs)
def parse_cluster(file_path): # Function for calculating the average BLAST hits identity value and extracting datails from the txt files

	# Get a list of all the records in the antiSMASH txt file
	records = file_to_list(file_path)

	# Extract and calculate the results from that list
	## (Known)clusterblast ID
	clust_IDs = [re.search("\d+\. (.*)", metaline[0][0]).group(1) for metaline in records[0]]

	## (Known)clusterblast annotations
	clust_annotations = [re.search("Source: (.*)", metaline[0][1]).group(1) for metaline in records[0]]
	for j in range(len(clust_annotations)):  # Make annotations uppercase
		clust_annotations[j] = clust_annotations[j][0].upper() + clust_annotations[j][1:]

	## Cumulative BLAST score
	blast_score_cum = [re.search("Cumulative BLAST score: (\d+)", metaline[0][4]).group(1) for metaline in records[0]]

	## Average BLAST identity to (Known)clusterblast genes
	blast_num = []
	blast_identity_averages = []
	blast_score_averages = []
	blast_coverage_averages = []
	blast_cds_annotations = []

	for record in records[0]:

		# Store (Known)clusterblast IDs and corresponding annotations in dictionary (if present)
		try:
			gene_match_ids = {line[0]: line[5] for line in record[1]} # Dictionary of (Known)clusterblast gene IDs and annotations
		except IndexError:
			gene_match_ids = {}

		# Extract annotation IDs, identity values, BLAST scores and subject coverage from BLAST hits table
		blast_hit_ids = [line[1] for line in record[2]]
		blast_hit_identities = [float(line[2]) for line in record[2]]
		blast_hit_scores = [float(line[3]) for line in record[2]]
		blast_hit_coverages = [float(line[4]) for line in record[2]]

		# Get number of BLAST hits per record
		blast_num.append(str(len(record[2])))

		# Get (Known)clusterblast cds_annotations
		blast_cds_annotations.append([gene_match_ids[hit_id] for hit_id in blast_hit_ids if hit_id in gene_match_ids.keys()])

		# Calculate avarage identity score
		blast_identity_averages.append(str(round(sum(blast_hit_identities)/len(blast_hit_identities), 2)))
		blast_score_averages.append(str(round(sum(blast_hit_scores)/len(blast_hit_scores), 2)))
		blast_coverage_averages.append(str(round(sum(blast_hit_coverages)/len(blast_hit_coverages), 2)))

	return clust_IDs, clust_annotations, blast_num, blast_identity_averages, blast_score_cum, blast_score_averages, blast_coverage_averages, blast_cds_annotations

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
