#!/usr/bin/env python3

"""
This script was provided by Matthias Blum on 09-07-2023
Modified on 13-09-2023 by Francesco Costa
to find proteins containing the same domain as input proteins which could be used to improve them
"""

import json
import sqlite3
import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser()

def file_path(filepath:str) -> str:
    """Validate input path"""
    if os.path.exists(filepath):
        return filepath
    else:
        raise NotADirectoryError(filepath)

parser.add_argument(
    "--proteins_table", required=True, help="File .csv containing a two columns: domain,proteins where domain is the\
        domain intended to rescue with the pipeline and proteins contains a list of : separated protein Acc. Use pfamA_accs\
        and uniprot_acc respectively. Use one row per domain. If you input more than 10 proteins for at least one domain,\
        increment max_proteins parameter accordingly.",
 type=file_path
)

parser.add_argument(
    "--database", required=False, help="Path/to/alphafold-db ", type=str, default="/hps/nobackup/agb/interpro/mblum/alphafold/pfam-alphafold_old.sqlite"
)

parser.add_argument(
    "--max_proteins", required=False, help="Number of protein candidates to consider as possible templates (>70 plddt).",
 type=int, default=10
)


def main():
    args = parser.parse_args()

    """
    Get protein domains of interest
    """
    df = pd.read_csv(args.proteins_table)
    entries = df.domain.to_list()

    con = sqlite3.connect(args.database)
    sql = """
        SELECT protein_id, reviewed, complete, dom_score, glo_score
        FROM pfam2uniprot
        WHERE entry_id = ?
    """
    outfile = open('proteins_per_domain.csv', 'w')
    for pfam_id in entries:
        print(f"Working on {pfam_id}...")
        rows = con.execute(sql, [pfam_id]).fetchall()
        
        collected_proteins = []
        
        # use domain score
        for alphafold_id, reviewed, complete, dom_score, glo_score in rows:
            if dom_score >= 70 and glo_score >= 70:
                collected_proteins.append((alphafold_id, reviewed, 1- complete, dom_score, glo_score))
        """
        Prioritize reviewed and complete proteins, then:
            - low/med: with a low mean-pLDDT score
            - high: with a high pLDDT score
        """

        collected_proteins.sort(key=lambda x: (-x[1], x[2], -x[3], -x[4]))
        outfile.write(f'{pfam_id} {",".join([e[0] for e in collected_proteins[:args.max_proteins + 10]])}\n')
        outfile.write(f'{pfam_id} {",".join([df[df.domain == pfam_id].proteins.to_list()[0].replace(":", ",")])}\n')

    outfile.close()
    con.close()


if __name__ == "__main__":
    main()

