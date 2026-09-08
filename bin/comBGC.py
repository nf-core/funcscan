#!/usr/bin/env python3

# Written by Jasmin Frangenberg and released under the MIT license.
# See below for full license text.

from Bio import SeqIO
import pandas as pd
import argparse
import os
import re

"""
===============================================================================
MIT License
===============================================================================

Copyright (c) 2023 Jasmin Frangenberg

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

tool_version = "0.6.3"
welcome = """\
                ........................
                    * comBGC v.{version} *
                ........................
    This tool aggregates the results of BGC prediction tools:
                antiSMASH, deepBGC, and GECCO
    For detailed usage documentation please refer
    to https://nf-co.re/funcscan
    .........................................................""".format(
    version=tool_version
)

# Initialize parser
parser = argparse.ArgumentParser(
    prog="comBGC",
    formatter_class=argparse.RawTextHelpFormatter,
    description=(welcome),
    add_help=True,
)

# Input options
parser.add_argument(
    "-i",
    "--input",
    metavar="PATH(s)",
    dest="input",
    nargs="*",
    help="""path(s) to the required output file(s) of antiSMASH, DeepBGC and/or GECCO
these can be:
- antiSMASH: <sample name>.gbk and (optional) knownclusterblast/ directory
- DeepBGC:   <sample name>.bgc.tsv
- GECCO:     <sample name>.clusters.tsv
Note: Please provide files from a single sample only. If you would like to
summarize multiple samples, please see the --antismash_multiple_samples flag.""",
)
parser.add_argument(
    "-o",
    "--outdir",
    metavar="PATH",
    dest="outdir",
    nargs="?",
    help="directory for comBGC output. Default: current directory",
    type=str,
    default=".",
)
parser.add_argument(
    "-a",
    "--antismash_multiple_samples",
    metavar="PATH",
    dest="antismash_multiple_samples",
    nargs="?",
    help="""directory of antiSMASH output. Should contain subfolders (one per
sample). Can only be used if --input is not specified.""",
    type=str,
)
parser.add_argument(
    "-vv", "--verbose", help="increase output verbosity", action="store_true"
)
parser.add_argument(
    "-v", "--version", help="show version number and exit", action="store_true"
)

# Get command line arguments
args = parser.parse_args()

# Assign input arguments to variables
input = args.input
dir_antismash = args.antismash_multiple_samples
outdir = args.outdir
verbose = args.verbose
version = args.version

if version:
    exit("comBGC {version}".format(version=tool_version))

input_antismash = []
input_deepbgc = []
input_gecco = []

# Assign input files to respective tools
if input:
    for path in input:
        if path.endswith(".gbk") and not re.search("region\d\d\d\.gbk$", path): # Make sure to only fetch relevant GBK files, i.e. those containing all collective antiSMASH BGCs
            with open(path) as infile:
                for line in infile:
                    if re.search("##GECCO-Data-START##", line):
                        input_gecco.append(path)
                        break
                    elif re.search("##antiSMASH-Data-START##", line):
                        input_antismash.append(path)
                        break
        elif path.endswith("bgc.tsv"):
            input_deepbgc = path
        elif path.endswith("clusters.tsv"):
            input_gecco.append(path)
        elif path.endswith("knownclusterblast/"):
            input_antismash.append(path)

if input and dir_antismash:
    exit(
        "The flags --input and --antismash_multiple_samples are mutually exclusive.\nPlease use only one of them (or see --help for how to use)."
    )

# Make sure that at least one input argument is given
if not (input_antismash or input_gecco or input_deepbgc or dir_antismash):
    exit(
        "Please specify at least one input file (i.e. output from antismash, deepbgc, or gecco) or see --help"
    )

########################
# ANTISMASH FUNCTIONS
########################


def prepare_multisample_input_antismash(antismash_dir):
    """
    Prepare string of input paths of a given antiSMASH output folder (with sample subdirectories)
    """
    sample_paths = []
    for root, subdirs, files in os.walk(antismash_dir):
        antismash_file = "/".join([root, "index.html"])
        if os.path.exists(antismash_file):
            sample = root.split("/")[-1]
            gbk_path = "/".join([root, sample]) + ".gbk"
            kkb_path = "/".join([root, "knownclusterblast"])
            if os.path.exists(kkb_path):
                sample_paths.append([gbk_path, kkb_path])
            else:
                sample_paths.append([gbk_path])
    return sample_paths


def parse_knownclusterblast(kcb_file_path):
    """
    Extract MIBiG IDs from knownclusterblast TXT file.
    """

    with open(kcb_file_path) as kcb_file:
        hits = 0
        MIBiG_IDs = []

        for line in kcb_file:
            if line == "Significant hits: \n" and not hits:
                hits = 1  # Indicating that the following lines contain relevant information
            elif line == "\n" and hits:
                break
            elif line != "Significant hits: \n" and hits:
                MIBiG_ID = re.search("(BGC\d+)", line).group(1)
                MIBiG_IDs.append(MIBiG_ID)
    return MIBiG_IDs


def antismash_workflow(antismash_paths):
    """
    Create data frame with aggregated antiSMASH output:
    - Open summary GBK and grab relevant information.
    - Extract the knownclusterblast output from the antiSMASH folder (MIBiG annotations) if present.
    - Return data frame with aggregated info.
    """

    antismash_sum_cols = [
        "Sample_ID",
        "Prediction_tool",
        "Contig_ID",
        "Product_class",
        "BGC_probability",
        "BGC_complete",
        "BGC_start",
        "BGC_end",
        "BGC_length",
        "CDS_ID",
        "CDS_count",
        "PFAM_domains",
        "MIBiG_ID",
        "InterPro_ID",
    ]
    antismash_out = pd.DataFrame(columns=antismash_sum_cols)

    CDS_ID = []
    CDS_count = 0

    # Distinguish input files (i.e. GBK file and "knownclusterblast" folder)
    kcb_path = []
    for path in antismash_paths:
        if re.search("knownclusterblast", path):
            kcb_path = re.search(".*knownclusterblast.*", path).group()
        else:
            gbk_path = path

        kcb_files = []
        if kcb_path:
            kcb_files = [
                file
                for file in os.listdir(kcb_path)
                if file.startswith("c") and file.endswith(".txt")
            ]

        # Aggregate information
        Sample_ID = gbk_path.split("/")[-1].split(".gbk")[
            -2
        ]  # Assuming file name equals sample name
        if verbose:
            print("\nParsing antiSMASH file(s): " + Sample_ID + "\n... ", end="")

        with open(gbk_path) as gbk:
            for record in SeqIO.parse(
                gbk, "genbank"
            ):  # GBK records are contigs in this case
                # Initiate variables per contig
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
                    # Extract relevant infos from the first protocluster feature from the contig record
                    if feature.type == "protocluster":
                        if (
                            antismash_out_line
                        ):  # If there is more than 1 BGC per contig, reset the output line for new BGC. Assuming that BGCs do not overlap.
                            if not CDS_ID:
                                CDS_ID = ["NA"]
                            antismash_out_line = {  # Create dictionary of BGC info
                                "Sample_ID": Sample_ID,
                                "Prediction_tool": "antiSMASH",
                                "Contig_ID": Contig_ID,
                                "Product_class": ";".join(Product_class),
                                "BGC_probability": "NA",
                                "BGC_complete": BGC_complete,
                                "BGC_start": BGC_start,
                                "BGC_end": BGC_end,
                                "BGC_length": BGC_length,
                                "CDS_ID": ";".join(CDS_ID),
                                "CDS_count": CDS_count,
                                "PFAM_domains": ";".join(PFAM_domains),
                                "MIBiG_ID": MIBiG_ID,
                                "InterPro_ID": "NA",
                            }
                            antismash_out_line = pd.DataFrame([antismash_out_line])
                            antismash_out = pd.concat(
                                [antismash_out, antismash_out_line], ignore_index=True
                            )
                            antismash_out_line = {}

                            # Reset variables per BGC
                            CDS_ID = []
                            CDS_count = 0
                            PFAM_domains = []

                        # Extract all the BGC info
                        Product_class = feature.qualifiers["product"]
                        for i in range(len(Product_class)):
                            Product_class[i] = (
                                Product_class[i][0].upper() + Product_class[i][1:]
                            )  # Make first letters uppercase, e.g. lassopeptide -> Lassopeptide

                        if feature.qualifiers["contig_edge"] == ["True"]:
                            BGC_complete = "No"
                        elif feature.qualifiers["contig_edge"] == ["False"]:
                            BGC_complete = "Yes"

                        BGC_start = (
                            feature.location.start + 1
                        )  # +1 because zero-based start position
                        BGC_end = feature.location.end
                        BGC_length = feature.location.end - feature.location.start

                        # If there are knownclusterblast files for the BGC, get MIBiG IDs of their homologs
                        if kcb_files:
                            print(kcb_files)
                            kcb_file = "{}_c{}.txt".format(
                                record.id, str(cluster_num)
                            )  # Check if this filename is among the knownclusterblast files
                            if kcb_file in kcb_files:
                                MIBiG_IDs = ";".join(
                                    parse_knownclusterblast(
                                        os.path.join(kcb_path, kcb_file)
                                    )
                                )
                                if MIBiG_IDs != "":
                                    MIBiG_ID = MIBiG_IDs
                                cluster_num += 1

                    # Count functional CDSs (no pseudogenes) and get the PFAM annotation
                    elif (
                        feature.type == "CDS"
                        and "translation" in feature.qualifiers.keys()
                        and BGC_start != ""
                    ):  # Make sure not to count pseudogenes (which would have no "translation tag") and count no CDSs before first BGC
                        if (
                            feature.location.end <= BGC_end
                        ):  # Make sure CDS is within the current BGC region
                            if "locus_tag" in feature.qualifiers:
                                CDS_ID.append(feature.qualifiers["locus_tag"][0])
                            CDS_count += 1
                            if "sec_met_domain" in feature.qualifiers.keys():
                                for PFAM_domain in feature.qualifiers["sec_met_domain"]:
                                    PFAM_domain_name = re.search(
                                        "(.+) \(E-value", PFAM_domain
                                    ).group(1)
                                    PFAM_domains.append(PFAM_domain_name)

                # Create dictionary of BGC info
                if not CDS_ID:
                    CDS_ID = ["NA"]
                antismash_out_line = {
                    "Sample_ID": Sample_ID,
                    "Prediction_tool": "antiSMASH",
                    "Contig_ID": Contig_ID,
                    "Product_class": ";".join(Product_class),
                    "BGC_probability": "NA",
                    "BGC_complete": BGC_complete,
                    "BGC_start": BGC_start,
                    "BGC_end": BGC_end,
                    "BGC_length": BGC_length,
                    "CDS_ID": ";".join(CDS_ID),
                    "CDS_count": CDS_count,
                    "PFAM_domains": ";".join(PFAM_domains),
                    "MIBiG_ID": MIBiG_ID,
                    "InterPro_ID": "NA",
                }

                if BGC_start != "":  # Only keep records with BGCs
                    antismash_out_line = pd.DataFrame([antismash_out_line])
                    antismash_out = pd.concat(
                        [antismash_out, antismash_out_line], ignore_index=True
                    )

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


def deepbgc_workflow(deepbgc_path):
    """
    Create data frame with aggregated deepBGC output.
    """

    if verbose:
        print("\nParsing deepBGC file\n... ", end="")

    # Prepare input and output columns
    deepbgc_map_dict = {
        "sequence_id": "Contig_ID",
        "nucl_start": "BGC_start",
        "nucl_end": "BGC_end",
        "nucl_length": "BGC_length",
        "num_proteins": "CDS_count",
        "deepbgc_score": "BGC_probability",
        "product_class": "Product_class",
        "protein_ids": "CDS_ID",
        "pfam_ids": "PFAM_domains",
    }
    deepbgc_sum_cols = [
        "Sample_ID",
        "Prediction_tool",
        "Contig_ID",
        "Product_class",
        "BGC_probability",
        "BGC_complete",
        "BGC_start",
        "BGC_end",
        "BGC_length",
        "CDS_ID",
        "CDS_count",
        "PFAM_domains",
        "MIBiG_ID",
        "InterPro_ID",
    ]
    deepbgc_unused_cols = [
        "detector_version",
        "detector_label",
        "bgc_candidate_id",
        "num_domains",
        "num_bio_domains",
        "product_activity",
        "antibacterial",
        "cytotoxic",
        "inhibitor",
        "antifungal",
        "Alkaloid",
        "NRP",
        "Other",
        "Polyketide",
        "RiPP",
        "Saccharide",
        "Terpene",
        "bio_pfam_ids",
    ]

    # Grab deepBGC sample ID
    sample = os.path.basename(deepbgc_path).rsplit(".bgc", 1)[0]

    # Initiate dataframe
    deepbgc_out = pd.DataFrame(columns=deepbgc_sum_cols)

    # Add relevant deepBGC output columns per BGC
    deepbgc_df = (
        pd.read_csv(deepbgc_path, sep="\t")
        .drop(deepbgc_unused_cols, axis=1)
        .rename(columns=deepbgc_map_dict)
    )
    deepbgc_df["Sample_ID"] = sample
    deepbgc_df["Prediction_tool"] = "deepBGC"
    deepbgc_df["BGC_complete"] = "NA"
    deepbgc_df["MIBiG_ID"] = "NA"
    deepbgc_df["InterPro_ID"] = "NA"

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
    """
    Retrieve InterPro IDs from GECCO GBK file.
    """

    with open(gbk_path) as gbk:
        ip_ids = []
        id_pattern = 'InterPro\:(.*)"'

        for line in gbk:
            if line.find("InterPro:") != -1:
                new_id = re.search(id_pattern, line).group(1)
                ip_ids.append(new_id)
        ipids_str = ";".join(map(str, ip_ids))
    return ipids_str


def gecco_workflow(gecco_paths):
    """
    Create data frame with aggregated GECCO output.
    """

    if verbose:
        print("\nParsing GECCO files\n... ", end="")

    # GECCO output columns that can be mapped (comBGC:GECCO)
    map_dict = {
        "sequence_id": "Contig_ID",
        "bgc_id": "cluster_id",
        "type": "Product_class",
        "average_p": "BGC_probability",
        "start": "BGC_start",
        "end": "BGC_end",
        "domains": "PFAM_domains",
        "proteins": "CDS_ID",
    }
    summary_cols = [
        "Sample_ID",
        "Prediction_tool",
        "Contig_ID",
        "Product_class",
        "BGC_probability",
        "BGC_complete",
        "BGC_start",
        "BGC_end",
        "BGC_length",
        "CDS_ID",
        "CDS_count",
        "PFAM_domains",
        "MIBiG_ID",
        "InterPro_ID",
    ]
    unused_cols = [
        "max_p",
        "alkaloid_probability",
        "polyketide_probability",
        "ripp_probability",
        "saccharide_probability",
        "terpene_probability",
        "nrp_probability",
    ]

    tsv_path = ""
    gbk_paths = []

    for path in gecco_paths:
        if path.endswith(".tsv"):
            tsv_path = path
        else:
            gbk_paths.append(path)

    # Initiate dataframe
    gecco_out = pd.DataFrame(columns=summary_cols)

    # Add sample information
    sample = tsv_path.split("/")[-1].split(".")[0]
    gecco_df = (
        pd.read_csv(tsv_path, sep="\t")
        .drop(unused_cols, axis=1)
        .rename(columns=map_dict)
    )

    # Fill columns (1 row per BGC)
    gecco_df["Sample_ID"] = sample
    gecco_df["BGC_length"] = gecco_df["BGC_end"] - gecco_df["BGC_start"]
    gecco_df["CDS_count"] = [
        len(gecco_df["CDS_ID"].iloc[i].split(";")) for i in range(0, gecco_df.shape[0])
    ]  # Number of contigs in 'Annotation_ID'
    gecco_df["Prediction_tool"] = "GECCO"

    # Add column 'InterPro_ID'
    for gbk_path in gbk_paths:
        bgc_id = gbk_path.split("/")[-1][0:-4]
        gecco_df.loc[gecco_df["cluster_id"] == bgc_id, "InterPro_ID"] = getInterProID(
            gbk_path
        )

    # Add empty columns with no output from GECCO
    gecco_df["BGC_complete"] = "NA"
    gecco_df["MIBiG_ID"] = "NA"
    gecco_out = pd.concat([gecco_out, gecco_df])

    # Fill all empty cells with NA
    for row in range(len(gecco_df["PFAM_domains"])):
        if gecco_out["PFAM_domains"].isnull().values[row]:
            gecco_out.loc[row, "PFAM_domains"] = "NA"

    # Return data frame with ordered columns
    gecco_out = gecco_out[summary_cols]

    if verbose:
        print("Done.")

    return gecco_out


########################
# MAIN
########################

if __name__ == "__main__":
    if input_antismash:
        tools = {
            "antiSMASH": input_antismash,
            "deepBGC": input_deepbgc,
            "GECCO": input_gecco,
        }
    elif dir_antismash:
        tools = {"antiSMASH": dir_antismash}
    else:
        tools = {"deepBGC": input_deepbgc, "GECCO": input_gecco}

    tools_provided = {}

    for tool in tools.keys():
        if tools[tool]:
            tools_provided[tool] = tools[tool]

    if verbose:
        print(welcome)
        print("\nYou provided input for: " + ", ".join(tools_provided.keys()))

    # Aggregate BGC information into data frame
    summary_antismash = pd.DataFrame()
    summary_deepbgc = pd.DataFrame()
    summary_gecco = pd.DataFrame()

    for tool in tools_provided.keys():
        if tool == "antiSMASH":
            if dir_antismash:
                antismash_paths = prepare_multisample_input_antismash(dir_antismash)
                for input_antismash in antismash_paths:
                    summary_antismash_temp = antismash_workflow(input_antismash)
                    summary_antismash = pd.concat(
                        [summary_antismash, summary_antismash_temp]
                    )
            else:
                summary_antismash = antismash_workflow(input_antismash)
        elif tool == "deepBGC":
            summary_deepbgc = deepbgc_workflow(input_deepbgc)
        elif tool == "GECCO":
            summary_gecco = gecco_workflow(input_gecco)

    # Summarize and sort data frame
    summary_all = pd.concat([summary_antismash, summary_deepbgc, summary_gecco])
    summary_all.sort_values(
        by=["Sample_ID", "Contig_ID", "BGC_start", "BGC_length", "Prediction_tool"],
        axis=0,
        inplace=True,
    )

    # Rearrange and rename the columns in the summary df
    summary_all = summary_all.iloc[:, [0, 2, 1] + list(range(3, len(summary_all.columns)))]
    summary_all.rename(columns={'Sample_ID':'sample_id', 'Contig_ID':'contig_id', 'CDS_ID':'BGC_region_contig_ids'}, inplace=True)

    # Write results to TSV
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    summary_all.to_csv(
        os.path.join(outdir, "combgc_summary.tsv"), sep="\t", index=False
    )
    print("Your BGC summary file is: " + os.path.join(outdir, "combgc_summary.tsv"))
