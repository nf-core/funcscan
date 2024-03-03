#!/usr/bin/env python3

# Author: @darcy220606
# Date: March 2024
# Version: 0.1.0

# Required modules
import sys
import os
import pandas as pd
import numpy as np
import argparse

tool_version = "0.1.0"
#########################################
# TOP LEVEL: AMPCOMBI
#########################################
parser = argparse.ArgumentParser(prog = 'merge_taxonomy', formatter_class=argparse.RawDescriptionHelpFormatter,
                                usage='%(prog)s [options]',
                                description=('''\
    .............................................................................
                                    *merge_taxonomy*
    .............................................................................
                This script merges all three funcscan workflows with
    MMseqs2 taxonomy results. This is done in three submodules that can be
    activated seperately.
    .............................................................................'''),
                                epilog='''Thank you for running taxonomy_merge!''',
                                add_help=True)
parser.add_argument('--version', action='version', version='merge_taxonomy ' + tool_version)

#########################################
# SUBPARSERS
#########################################
subparsers = parser.add_subparsers(required=True)

#########################################
# SUBPARSERS : AMPCOMBI
#########################################
ampcombi_parser = subparsers.add_parser('ampcombi_taxa')

ampcombi_parser.add_argument("--ampcombi", dest="amp", nargs='?', help="Enter the path to the ampcombi_complete_summary.tsv' \n (default: %(default)s)",
                    type=str, default='ampcombi_complete_summary.csv')
ampcombi_parser.add_argument("--taxonomy", dest="taxa1", nargs='+', help="Enter the list of taxonomy files for all samples. ")

#########################################
# SUBPARSERS : COMBGC
#########################################
combgc_parser = subparsers.add_parser('combgc_taxa')

combgc_parser.add_argument("--combgc", dest="bgc", nargs='?', help="Enter the path to the combgc_complete_summary.tsv' \n (default: %(default)s)",
                    type=str, default='combgc_complete_summary.csv')
combgc_parser.add_argument("--taxonomy", dest="taxa2", nargs='+', help="Enter the list of taxonomy files for all samples. ")

#########################################
# SUBPARSERS : HAMRONIZATION
#########################################
hamronization_parser = subparsers.add_parser('hamronization_taxa')

hamronization_parser.add_argument("--hamronization", dest="arg", nargs='?', help="Enter the path to the hamronization_complete_summary.tsv' \n (default: %(default)s)",
                    type=str, default='hamronization_complete_summary.csv')
hamronization_parser.add_argument("--taxonomy", dest="taxa3",nargs='+', help="Enter the list of taxonomy files for all samples. ")

#########################################
# TAXONOMY
#########################################
def reformat_mmseqs_taxonomy(mmseqs_taxonomy):
    mmseqs2_df = pd.read_csv(mmseqs_taxonomy, sep='\t', header=None, names=['contig_id', 'taxid', 'rank_label', 'scientific_name', 'lineage', 'mmseqs_lineage_contig'])
    # remove the lineage column
    mmseqs2_df.drop('lineage', axis=1, inplace=True)
    mmseqs2_df['mmseqs_lineage_contig'].unique()
    # convert any classification that has Eukaryota/root to NaN as funcscan targets bacteria ONLY **
    for i, row in mmseqs2_df.iterrows():
        lineage = str(row['mmseqs_lineage_contig'])
        if 'Eukaryota' in lineage or 'root' in lineage:
            mmseqs2_df.at[i, 'mmseqs_lineage_contig'] = np.nan
            #mmseqs2_df['mmseqs_lineage_contig'].unique()
    # insert the sample name in the first column according to the file basename
    file_basename = os.path.basename(mmseqs_taxonomy)
    filename = os.path.splitext(file_basename)[0]
    mmseqs2_df.insert(0, 'sample_id', filename)
    return mmseqs2_df

#########################################
# FUNCTION : AMPCOMBI
#########################################
def ampcombi_taxa(args):
    merged_df = pd.DataFrame()

    # assign input args to variables
    ampcombi = args.amp
    taxa_list = args.taxa1

    # prepare the taxonomy files
    taxa_df = pd.DataFrame()
    # append the dfs to the taxonomy_files_combined
    for file in taxa_list: # list of taxa files ['','']
        df = reformat_mmseqs_taxonomy(file)
        taxa_df = pd.concat([taxa_df, df])

    # filter the tool df
    tool_df = pd.read_csv(ampcombi, sep=',') #current ampcombi version is comma sep. CHANGE WITH VERSION 0.2.0
    # make sure 1st and 2nd column have the same column labels
    tool_df.rename(columns={tool_df.columns[0]: 'sample_id'}, inplace=True)
    tool_df.rename(columns={tool_df.columns[1]: 'contig_id'}, inplace=True)
    # grab the real contig id in another column copy for merging
    tool_df['contig_id_merge'] = tool_df['contig_id'].str.rsplit('_', 1).str[0]

    # merge rows from taxa to ampcombi_df based on substring match in sample_id
    # grab the unique sample names from the taxonomy table
    samples_taxa = taxa_df['sample_id'].unique()
    # for every sampleID in taxadf merge the results
    for sampleID in samples_taxa:
        # subset ampcombi
        subset_tool = tool_df.loc[tool_df['sample_id'].str.contains(sampleID)]
        # subset taxa
        subset_taxa = taxa_df.loc[taxa_df['sample_id'].str.contains(sampleID)]
        # merge
        subset_df = pd.merge(subset_tool, subset_taxa, left_on = 'contig_id_merge', right_on='contig_id', how='left')
        # cleanup the table
        columnsremove = ['contig_id_merge','contig_id_y', 'sample_id_y']
        subset_df.drop(columnsremove, axis=1, inplace=True)
        subset_df.rename(columns={'contig_id_x': 'contig_id', 'sample_id_x':'sample_id'},inplace=True)
        # append in the combined_df
        merged_df = merged_df.append(subset_df, ignore_index=True)

    # write to file
    merged_df.to_csv('ampcombi_complete_summary_taxonomy.tsv', sep='\t', index=False)

#########################################
# FUNCTION : COMBGC
#########################################
def combgc_taxa(args):
    merged_df = pd.DataFrame()

    # assign input args to variables
    combgc = args.bgc
    taxa_list = args.taxa2

    # prepare the taxonomy files
    taxa_df = pd.DataFrame()
    # append the dfs to the taxonomy_files_combined
    for file in taxa_list: # list of taxa files ['','']
        df = reformat_mmseqs_taxonomy(file)
        taxa_df = pd.concat([taxa_df, df])

    # filter the tool df
    tool_df = pd.read_csv(combgc, sep='\t')
    # make sure 1st and 2nd column have the same column labels
    tool_df.rename(columns={tool_df.columns[0]: 'sample_id'}, inplace=True)
    tool_df.rename(columns={tool_df.columns[1]: 'contig_id'}, inplace=True)

    # merge rows from taxa to ampcombi_df based on substring match in sample_id
    # grab the unique sample names from the taxonomy table
    samples_taxa = taxa_df['sample_id'].unique()
    # for every sampleID in taxadf merge the results
    for sampleID in samples_taxa:
        # subset ampcombi
        subset_tool = tool_df.loc[tool_df['sample_id'].str.contains(sampleID)]
        # subset taxa
        subset_taxa = taxa_df.loc[taxa_df['sample_id'].str.contains(sampleID)]
        # merge
        subset_df = pd.merge(subset_tool, subset_taxa, left_on = 'contig_id', right_on='contig_id', how='left')
        # cleanup the table
        columnsremove = ['sample_id_y']
        subset_df.drop(columnsremove, axis=1, inplace=True)
        subset_df.rename(columns={'sample_id_x':'sample_id'},inplace=True)
        # append in the combined_df
        merged_df = merged_df.append(subset_df, ignore_index=True)

    # write to file
    merged_df.to_csv('combgc_complete_summary_taxonomy.tsv', sep='\t', index=False)

#########################################
# FUNCTION : HAMRONIZATION
#########################################
def hamronization_taxa(args):
    merged_df = pd.DataFrame()

    # assign input args to variables
    hamronization = args.arg
    taxa_list = args.taxa3

    # prepare the taxonomy files
    taxa_df = pd.DataFrame()
    # append the dfs to the taxonomy_files_combined
    for file in taxa_list: # list of taxa files ['','']
        df = reformat_mmseqs_taxonomy(file)
        taxa_df = pd.concat([taxa_df, df])

    # filter the tool df
    tool_df = pd.read_csv(hamronization, sep='\t')
    # rename the columns
    tool_df.rename(columns={'input_file_name':'sample_id', 'input_sequence_id':'contig_id'}, inplace=True)
    # reorder the columns
    new_order = ['sample_id', 'contig_id'] + [col for col in tool_df.columns if col not in ['sample_id', 'contig_id']]
    tool_df = tool_df.reindex(columns=new_order)
    # grab the real contig id in another column copy for merging
    tool_df['contig_id_merge'] = tool_df['contig_id'].str.rsplit('_', 1).str[0]

    # merge rows from taxa to ampcombi_df based on substring match in sample_id
    # grab the unique sample names from the taxonomy table
    samples_taxa = taxa_df['sample_id'].unique()
    # for every sampleID in taxadf merge the results
    for sampleID in samples_taxa:
        # subset ampcombi
        subset_tool = tool_df.loc[tool_df['sample_id'].str.contains(sampleID)]
        # subset taxa
        subset_taxa = taxa_df.loc[taxa_df['sample_id'].str.contains(sampleID)]
        # merge
        subset_df = pd.merge(subset_tool, subset_taxa, left_on = 'contig_id_merge', right_on='contig_id', how='left')
        # cleanup the table
        columnsremove = ['contig_id_merge','contig_id_y', 'sample_id_y']
        subset_df.drop(columnsremove, axis=1, inplace=True)
        subset_df.rename(columns={'contig_id_x': 'contig_id', 'sample_id_x':'sample_id'},inplace=True)
        # append in the combined_df
        merged_df = merged_df.append(subset_df, ignore_index=True)

    # write to file
    merged_df.to_csv('hamronization_complete_summary_taxonomy.tsv', sep='\t', index=False)

#########################################
# SUBPARSERS : DEFAULT
#########################################
ampcombi_parser.set_defaults(func=ampcombi_taxa)
combgc_parser.set_defaults(func=combgc_taxa)
hamronization_parser.set_defaults(func=hamronization_taxa)

if __name__ == '__main__':
    args = parser.parse_args()
    args.func(args)  # call the default function
