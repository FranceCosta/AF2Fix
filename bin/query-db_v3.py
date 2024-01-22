#!/usr/bin/env python3

"""
This script was provided by Matthias Blum on 09-07-2023
Modified on 29-08-2023 by Francesco Costa
Uses peak finding algorithm to detect suitable domains for AlphaFold submission
"""

import json
import sqlite3
from scipy.signal import find_peaks
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--database", required=False, help="Path/to/alphafold-db ", type=str, default="/hps/nobackup/agb/interpro/mblum/alphafold/pfam-alphafold.sqlite"
)

parser.add_argument(
    "--number_to_consider", required=False, help="Number of protein domains to include in the pipeline", type=int, default=50
)

parser.add_argument(
    "--proteins_per_class", required=False, help="Number of protein candidates to consider per class of AlphaFold confidence (<50, 50-70, >70 plddt). 10 more will be included and discarded by later script",
 type=int, default=10
)


def main():
    args = parser.parse_args()
    con = sqlite3.connect(args.database)

    entries = []

    """
    Load the distribution of mean domain pLDDT for each Pfam entry
    """
    sql = "SELECT id, distributions FROM entry WHERE id != 'alphafold'"
    for pfam_id, json_string in con.execute(sql):
        data = json.loads(json_string)

        """
        The distribution of pLDDT scores is splitted by superkingdom.
        There is one distribution for bacteria, eukaryotes, etc,
        so we need to merge them all to have to overall distribution.

        There are also per-superkingdom distributions where fragment proteins
        were ignored, so we need to skip these.
        """

        dist = [0] * 100
        # v is the binned distribution of plddts
        for k, v in data.items():
            if k.startswith("dom_") and not k.endswith("_nofrags"):
                for i, n in enumerate(v):
                    dist[i] += n

        total = sum(dist)
        if total < 100:
            # Ignore Pfam matching less than 100 proteins with an AF model
            continue
        # normalize
        dist = [i / total for i in dist]

        # apply peak detection algoritm
        peak_indeces , peak_heights = find_peaks(dist, height = 0.03, distance = 19)

        if len(peak_indeces) >= 2 and peak_indeces[-1] >= 70 and peak_indeces[0] <= 50:
            """
            Another way to detect distribution:
            find distributions with at least one peak above treshold with high accuracy and one with low accuracy
            """
            entries.append((pfam_id, abs(sum(dist[:50]) - sum(dist[70:]))))

    # Order entries by the absolute difference between low/high ratio (so most eqaul peaks go first)
    entries.sort(key=lambda e: e[1])
    sql = """
        SELECT protein_id, reviewed, complete, dom_score
        FROM entry2alphafold
        WHERE entry_id = ?
    """
    outfile = open('proteins_per_domain.csv', 'w')
    for pfam_id, _ in entries[:args.number_to_consider]:
        rows = con.execute(sql, [pfam_id]).fetchall()

        low = []
        med = []
        high = []

        for alphafold_id, reviewed, complete, score in rows:
            if score <= 50:
                obj = low
            elif score <= 70:
                obj = med
            else:
                obj = high

            obj.append((alphafold_id, reviewed, 1 - complete, score))
        """
        Prioritize reviewed and complete proteins, then:
            - low/med: with a low mean-pLDDT score
            - high: with a high pLDDT score
        """
        low.sort(key=lambda x: (-x[1], x[2], x[3]))
        med.sort(key=lambda x: (-x[1], x[2], x[3]))
        high.sort(key=lambda x: (-x[1], x[2], -x[3]))
        outfile.write(f'{pfam_id} {",".join([e[0] for e in low[:args.proteins_per_class + 10]])}\n')
        outfile.write(f'{pfam_id} {",".join([e[0] for e in med[:args.proteins_per_class + 10]])}\n')
        outfile.write(f'{pfam_id} {",".join([e[0] for e in high[:args.proteins_per_class + 10]])}\n')

    outfile.close()
    con.close()


if __name__ == "__main__":
    main()

