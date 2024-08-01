"""Get distribution of plddt for every pfam family with at least 100 structures in the AF-db"""

import json
import sqlite3


DATABASE = "assets/pfam-alphafold.sqlite"

def main():
    con = sqlite3.connect(DATABASE)

    entries = {}

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

        entries[pfam_id] = dist
    con.close()
    # save dictionary of distributions
    with open('data/distributions.json', 'w') as handle:
        json.dump(entries, handle)

if __name__ == "__main__":
    main()

