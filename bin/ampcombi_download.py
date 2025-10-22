#!/usr/bin/env python3

#########################################
# Authors: [Anan Ibrahim](https://github.com/Darcy220606/AMPcombi), [Louisa Perelo](https://github.com/louperelo)
# File: amp_database.py
# Source: https://github.com/Darcy220606/AMPcombi/blob/main/ampcombi/amp_database.py
# This source code is licensed under the MIT license
#########################################

# TITLE: Download the reference database specified by the user.

import pandas as pd
import requests
import os
import re
import subprocess
import argparse

from datetime import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

########################################
#  FUNCTION: DOWNLOAD DATABASES AND CLEAN DRAMP and APD
#########################################
def download_ref_db(database, threads):
    """
    Downloads a specified AMP (antimicrobial peptide) reference database based on the
    provided database name and saves it to the specified directory.
    This supports downloading databases only from DRAMP, APD, and UniRef100.
    Parameters:
    ----------
    db : str
        The directory path where the downloaded database should be saved.
    database : str
        The name of the database to download. Must be one of 'DRAMP', 'APD', or 'UniRef100'.
    threads : int
        Number of threads to use when downloading the UniRef100 database with `mmseqs`.
    """
    # Check which database was given
    if database == 'DRAMP':
        # Create dir
        db = 'amp_DRAMP_database'
        os.makedirs(db, exist_ok=True)
        # Download the file
        try:
            url = 'http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.txt'
            response = requests.get(url, allow_redirects=True)
            response.raise_for_status()  # Check for any download errors
            date = datetime.now().strftime("%Y_%m_%d")
            with open(db + '/' + f'general_amps_{date}.txt', 'wb') as file:
                file.write(response.content)
            print(f"File downloaded successfully and saved to {db}/general_amps_{date}.txt")
            # Create fasta version and clean it
            db_df = pd.read_csv(f'{db}/general_amps_{date}.txt', sep='\t')
            records = []
            valid_sequence_pattern = re.compile("^[ACDEFGHIKLMNPQRSTVWY]+$")
            for index, row in db_df.iterrows():
                sequence = row['Sequence']
                if valid_sequence_pattern.match(sequence):
                    record = SeqRecord(Seq(sequence), id=str(row['DRAMP_ID']), description="")
                    records.append(record)
            output_file = f'{db}/general_amps_{date}.fasta'
            SeqIO.write(records, output_file, "fasta")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download DRAMP AMP general database file: {e}")
            return

    if database == 'APD':
        # Create dir
        db = 'amp_APD_database'
        os.makedirs(db, exist_ok=True)
        # Download the file
        try:
            url = 'https://aps.unmc.edu/assets/sequences/APD_sequence_release_09142020.fasta'
            response = requests.get(url, allow_redirects=True, verify=False)  # Disable SSL verification due to site certificate issue
            response.raise_for_status()
            content = response.text
            print("APD AMP database downloaded successfully.")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download content: {e}")
            return
        # Save the content line-by-line exactly as is
        try:
            with open(db + '/' + 'APD_orig.fasta', 'w') as file:
                file.write(content)
            with open(f'{db}/APD.fasta', 'w') as output_handle:
                valid_sequence_pattern = re.compile("^[ACDEFGHIKLMNPQRSTVWY]+$")
                for record in SeqIO.parse(f'{db}/APD_orig.fasta', "fasta"):
                    sequence = str(record.seq)
                    if valid_sequence_pattern.match(sequence):
                        SeqIO.write(record, output_handle, "fasta")
            os.remove(db + '/' + 'APD_orig.fasta')
            print(f"APD AMP database saved successfully to {db}/APD.fasta")
            # Fasta to table
            headers = []
            sequences = []
            seq_ids = []
            for i, record in enumerate(SeqIO.parse(f'{db}/APD.fasta', "fasta")):
                sequence_id = record.description.split('|')[0]
                headers.append(record.description)
                sequences.append(str(record.seq))
                seq_ids.append(sequence_id)
            db_df = pd.DataFrame({
                "APD_ID": seq_ids,
                "APD_Description": headers,
                "APD_Sequence": sequences})
            db_df.to_csv(f'{db}/APD.txt', sep='\t', index=False, header=True)
            os.remove(db + '/' + 'APD.fasta')
            # Table to fasta
            records = []
            for index, row in db_df.iterrows():
                sequence = row['APD_Sequence']
                record = SeqRecord(Seq(sequence), id=str(row['APD_ID']), description="")
                records.append(record)
            output_file = f'{db}/APD.fasta'
            SeqIO.write(records, output_file, "fasta")
        except Exception as e:
            print(f"Failed to save APD AMP database: {e}")

    if database == 'UniRef100':
        # Create dir
        db = 'amp_UniRef100_database'
        os.makedirs(db, exist_ok=True)
        # Download the file
        try:
            os.makedirs(f'{db}/mmseqs2', exist_ok=True)
            command = f"mmseqs databases UniRef100 {db}/mmseqs2/ref_DB {db}/mmseqs2/tmp --remove-tmp-files true --threads {threads} -v 0"
            subprocess.run(command, shell=True, check=True)
            print(f"UniRef100 protein database downloaded successfully and saved to {db}/mmseqs2/UniRef100")
        except subprocess.CalledProcessError as e:
            print(f"Failed to download UniRef100 protein database: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Downloads a specified AMP (antimicrobial peptide) reference database based on the provided database name and saves it to the specified directory.")
    parser.add_argument("--database_id", dest="database", type=str, required=True, choices=["DRAMP", "APD", "UniRef100"],
                        help="Database ID - one of DRAMP, APD, or UniRef100. This parameter is required.")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads supplied to mmseqs databases. Only relevant in the case of 'UniRef100'. Default is 4.")

    args = parser.parse_args()
    download_ref_db(args.database, args.threads)
