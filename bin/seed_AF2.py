"""

    Run multiple seed AF2

    `conda activate /nfs/research/agb/research/francesco/software/localcolabfold/colabfold-conda`

"""

import random
import pandas as pd
import os

RESULT_DIR = 'output'
OUTPUT_DIR = "AF2_seed"
WEIGHTS_DIR = "tmp/params"
MSA_DIR = "output/AF2/MSAs/"
LOGS_DIR = "logs"

def main():
    AF2_df = pd.read_csv(os.path.join(RESULT_DIR, 'summary_results', 'AF2Results.csv'))
    acc_to_pfam = AF2_df[["uniprot_acc", "pfamA_acc"]].set_index(["uniprot_acc"]).to_dict()['pfamA_acc']
    to_rescue = AF2_df[(AF2_df["plddt"]<70)&(AF2_df["condition"]=="no_template")]["uniprot_acc"].to_list()
    rescued = AF2_df[(AF2_df["uniprot_acc"].isin(to_rescue))&(AF2_df["plddt"]>70)&(AF2_df["condition"]=="best_pick_structure_plddt")]["uniprot_acc"].to_list()

    # Select 50 random rescued proteins
    random.Random(42).shuffle(rescued)
    outlist = []
    for uniprot_acc in rescued[:50]:
        pfamA_acc = acc_to_pfam[uniprot_acc]
        outlist.append([pfamA_acc, uniprot_acc])

    for item in outlist:
        for seed in range(5,9):
            pfamA_acc = item[0]
            uniprot_acc = item[1]
            msa_path = os.path.join(MSA_DIR, pfamA_acc, f"{uniprot_acc}.a3m")
            output_dir = os.path.join(OUTPUT_DIR, f"{uniprot_acc}_{seed}")
            cmd = f'sbatch -t 06:00:00 --mem 20GB --gres=gpu:a100:1 \
                    -o {os.path.join(LOGS_DIR, f"{uniprot_acc}_{seed}.log")} \
                    -e {os.path.join(LOGS_DIR, f"{uniprot_acc}_{seed}.err")} \
                    -J "{uniprot_acc}_{seed}" \
                    --wrap="colabfold_batch --data {WEIGHTS_DIR} --num-recycle 5 --use-gpu-relax \
                        --num-relax 1 --model-type auto --random-seed {seed} {msa_path} {output_dir}"'
            #print(cmd)
            os.system(cmd)

if __name__ == "__main__":
    main()