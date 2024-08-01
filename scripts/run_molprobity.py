"""

    Run molprobity on output structures and collect results

"""

import subprocess
import os
import re
import glob
import pandas as pd

AF2 = "output/AF2"
CONDITIONS = ["no_template", "template_MSA", "template_single_seq"]
CACHE = "test/molprobity"
OUT_TABLE = "output/summary_results/molprobity_results.csv"

def main():
    os.makedirs(CACHE, exist_ok=True)
    for condition in CONDITIONS:
        for pfamA_acc in os.listdir(os.path.join(AF2, condition)):
            os.makedirs(os.path.join(CACHE, condition, pfamA_acc), exist_ok=True)
            for prot in glob.glob(os.path.join(AF2, condition, pfamA_acc, "*.pdb")):
                molprob_res = runAndParse(prot)
                pdb_name = os.path.basename(prot).split(".pdb")[0]
                with open(os.path.join(CACHE, condition, pfamA_acc, pdb_name+".csv"), "wt") as fh:
                    fh.write(",".join(molprob_res.keys()) + "\n")
                    fh.write(",".join(molprob_res.values()) + "\n")
    df = collect_molprobity(CACHE)
    df.to_csv(OUT_TABLE, index=False)

def collect_molprobity(molprobity_res_dir: str) -> pd.DataFrame:
    """Collect output data for molprobity_calculation"""

    df_molprobity = pd.DataFrame(columns=['condition', 'pfamA_acc', 'uniprot_acc', 'Ramachandran outliers', 'favored', 'Rotamer outliers',
                                          'C-beta deviations', 'Clashscore', 'RMS(bonds)', 'RMS(angles)','MolProbity score'])
    for condition in CONDITIONS:
        for pfamA_acc in filter(lambda x: x.startswith('PF'), os.listdir(os.path.join(molprobity_res_dir, condition))):
            for file in filter(lambda x: x.endswith('.csv'), os.listdir(os.path.join(molprobity_res_dir, condition, pfamA_acc))):
                tmp_df = pd.read_csv(os.path.join(molprobity_res_dir, condition, pfamA_acc, file))
                uniprot_acc = file.replace(".csv", "").split("_")[0]
                tmp_df['uniprot_acc'] = uniprot_acc
                tmp_df['condition'] = condition
                tmp_df['pfamA_acc'] = pfamA_acc
                df_molprobity= pd.concat([df_molprobity, tmp_df])
    
    return df_molprobity


def runAndParse(pdb_path:str):
    """
    
        Run Molprobity and parse results
    
    """
    
    patterns = {
    "Ramachandran outliers": r"Ramachandran outliers\s*=\s*([\d\.]+) %",
    "favored": r"favored\s*=\s*([\d\.]+) %",
    "Rotamer outliers": r"Rotamer outliers\s*=\s*([\d\.]+) %",
    "C-beta deviations": r"C-beta deviations\s*=\s*([\d]+)",
    "Clashscore": r"Clashscore\s*=\s*([\d\.]+)",
    "RMS(bonds)": r"RMS\(bonds\)\s*=\s*([\d\.]+)",
    "RMS(angles)": r"RMS\(angles\)\s*=\s*([\d\.]+)",
    "MolProbity score": r"MolProbity score\s*=\s*([\d\.]+)"
    }
    cmd = f"singularity exec -e /nfs/research/agb/research/francesco/singularity_env/molprobity.sif molprobity.molprobity {pdb_path}\
                output.prefix=to_be_deleted"
    try:
        cmd_output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as exc:
        print("Status : FAIL", exc.returncode, exc.output)

    with open("to_be_deleted.out", "rt") as fh:
        out = ""
        flag = False
        for line in fh.readlines():
            if "=== Summary ===" in line:
                flag = True
            if flag == False:
                continue

            out += line
    values = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, out)
        if match:
            values[key] = match.group(1)

    for tmp_file in glob.glob("to_be_deleted*"):
        os.remove(tmp_file)
    
    return values

if __name__ == "__main__":
    main()