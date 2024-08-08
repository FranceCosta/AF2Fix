#! /usr/bin/python3

"""
    Run molprobity score
    Francesco Costa fcosta@ebi.ac.uk
    19/07/2024

"""

import argparse
import subprocess
import os
import re
import glob

def isfile(x):
    """
    'Type' for argparse - checks that file or dir exists
    """
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

parser = argparse.ArgumentParser(
    description="Provide a PDB_file"
)

parser.add_argument(
    "--PDB_file", 
    required=False, 
    help="PDB file to analyze", 
    type=isfile
)

def main():
    args = parser.parse_args()
    res_dict = runAndParse(args.PDB_file)
    pdb_name = os.path.basename(args.PDB_file).split(".pdb")[0]
    with open(pdb_name+".csv", "wt") as fh:
        fh.write(",".join(res_dict.keys()) + "\n")
        fh.write(",".join(res_dict.values()) + "\n")

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
    cmd = f"molprobity.molprobity {pdb_path}\
                output.prefix=to_be_deleted"
    try:
        _ = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as exc:
        print("Status : FAIL", exc.returncode, exc.output)
        import sys
        sys.exit()

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