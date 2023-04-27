#!/usr/bin/env python3

#########################################
# Authors: [Anan Ibrahim](https://github.com/brianjohnhaas), [Louisa Perelo](https://github.com/louperelo)
# File: amp_database.py
# Source: https://github.com/Darcy220606/AMPcombi/blob/main/ampcombi/amp_database.py
# Source+commit: https://github.com/Darcy220606/AMPcombi/commit/a75bc00c32ecf873a133b18cf01f172ad9cf0d2d/ampcombi/amp_database.py
# Download Date: 2023-03-08, commit: a75bc00c
# This source code is licensed under the MIT license
#########################################

# TITLE: Download the DRAMP database if input db empty AND and make database compatible for diamond

import pandas as pd
import requests
import os
from datetime import datetime
import subprocess
from Bio import SeqIO
import tempfile
import shutil


########################################
#  FUNCTION: DOWNLOAD DRAMP DATABASE AND CLEAN IT
#########################################
def download_DRAMP(db):
    ##Download the (table) file and store it in a results directory
    url = "http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.xlsx"
    r = requests.get(url, allow_redirects=True)
    with open(db + "/" + "general_amps.xlsx", "wb") as f:
        f.write(r.content)
    ##Convert excel to tab sep file and write it to a file in the DRAMP_db directly with the date its downloaded
    date = datetime.now().strftime("%Y_%m_%d")
    ref_amps = pd.read_excel(db + "/" + r"general_amps.xlsx")
    ref_amps.to_csv(db + "/" + f"general_amps_{date}.tsv", index=None, header=True, sep="\t")
    ##Download the (fasta) file and store it in a results directory
    urlfasta = (
        "http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.fasta"
    )
    z = requests.get(urlfasta)
    fasta_path = os.path.join(db + "/" + f"general_amps_{date}.fasta")
    with open(fasta_path, "wb") as f:
        f.write(z.content)
    ##Cleaning step to remove ambigous aminoacids from sequences in the database (e.g. zeros and brackets)
    new_fasta = db + "/" + f"general_amps_{date}_clean.fasta"
    seq_record = SeqIO.parse(open(fasta_path), "fasta")
    with open(new_fasta, "w") as f:
        for record in seq_record:
            id, sequence = record.id, str(record.seq)
            letters = [
                "A",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H",
                "I",
                "K",
                "L",
                "M",
                "N",
                "P",
                "Q",
                "R",
                "S",
                "T",
                "V",
                "W",
                "Y"
            ]
            new = "".join(i for i in sequence if i in letters)
            f.write(">" + id + "\n" + new + "\n")
    return os.remove(fasta_path), os.remove(db + "/" + r"general_amps.xlsx")


download_DRAMP("amp_ref_database")
