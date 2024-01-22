#! /usr/bin/env python3

"""
    Convert PDB to MMCIF and add metadata needed to run localcolabfol. Inspired from https://github.com/sokrypton/ColabFold/blob/074aa2bbc0093b10e70ef632708f131d02a7aa7f/colabfold/batch.pyd
    Saves .MMCIF files naming them with a classic randomly generated PDB 4-alphanumeric code, required by colabfold. The correspondance name-PDB-code is found in .mapper.json.
    Francesco Costa
"""

import numpy as np
from colabfold.utils import CFMMCIFIO
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import Select
import os, json, string, random
from pathlib import Path
import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument(
    "--AFres_dir", 
    required=True, 
    help="Path/to/AF directory with no template results", 
    type=str
)
parser.add_argument(
    "--output_dir",
    required=False,
    help="output directory",
    type=str,
    default='./'
)


modified_mapping = {
  "MSE" : "MET", "MLY" : "LYS", "FME" : "MET", "HYP" : "PRO",
  "TPO" : "THR", "CSO" : "CYS", "SEP" : "SER", "M3L" : "LYS",
  "HSK" : "HIS", "SAC" : "SER", "PCA" : "GLU", "DAL" : "ALA",
  "CME" : "CYS", "CSD" : "CYS", "OCS" : "CYS", "DPR" : "PRO",
  "B3K" : "LYS", "ALY" : "LYS", "YCM" : "CYS", "MLZ" : "LYS",
  "4BF" : "TYR", "KCX" : "LYS", "B3E" : "GLU", "B3D" : "ASP",
  "HZP" : "PRO", "CSX" : "CYS", "BAL" : "ALA", "HIC" : "HIS",
  "DBZ" : "ALA", "DCY" : "CYS", "DVA" : "VAL", "NLE" : "LEU",
  "SMC" : "CYS", "AGM" : "ARG", "B3A" : "ALA", "DAS" : "ASP",
  "DLY" : "LYS", "DSN" : "SER", "DTH" : "THR", "GL3" : "GLY",
  "HY3" : "PRO", "LLP" : "LYS", "MGN" : "GLN", "MHS" : "HIS",
  "TRQ" : "TRP", "B3Y" : "TYR", "PHI" : "PHE", "PTR" : "TYR",
  "TYS" : "TYR", "IAS" : "ASP", "GPL" : "LYS", "KYN" : "TRP",
  "CSD" : "CYS", "SEC" : "CYS"
}

class ReplaceOrRemoveHetatmSelect(Select):
  def accept_residue(self, residue):
    hetfield, _, _ = residue.get_id()
    if hetfield != " ":
      if residue.resname in modified_mapping:
        # set unmodified resname
        residue.resname = modified_mapping[residue.resname]
        # clear hetatm flag
        residue._id = (" ", residue._id[1], " ")
        t = residue.full_id
        residue.full_id = (t[0], t[1], t[2], residue._id)
        return 1
      return 0
    else:
      return 1

def convert_pdb_to_mmcif(pdb_file, output):
    """convert existing pdb files into mmcif with the required poly_seq and revision_date"""
    i = os.path.basename(output).split('.cif')[0]
    cif_file = Path(output)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(i, pdb_file)
    cif_io = CFMMCIFIO()
    cif_io.set_structure(structure)
    cif_io.save(str(cif_file), ReplaceOrRemoveHetatmSelect())

def main():
     
    args = parser.parse_args()
    AFres_directory = args.AFres_dir
    # save output in folder named as the domain
    domain = AFres_directory.split("/")[-1]
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    #output_dir = os.path.join(args.output_dir, domain)
    #print(f'Files in directory: {os.listdir("./")}\nPath to save: {output_dir}.')
    #if os.path.isdir(output_dir) == False:
    #    os.makedirs(output_dir)
    # upload table with domain info
    # some lines contain 2 or more proteins

    # get best AF models based on domain plddt and protein plddb (both must be >= 70)
    # include proteins with repetitions of the same domain if all domains are > 70
    """
    best_models = []
    for file in filter(lambda x: x.endswith('.json') and 'rank_001' in x,  os.listdir(AFres_directory)):
        with open(os.path.join(AFres_directory, file), 'r') as fp:
            tmp_dic = json.load(fp)
            uniprot_acc = file.split('_')[0]
            domain_plddts = []
            protein_pddt = np.mean(tmp_dic['plddt'])
            for index, row in protein_info[protein_info.uniprot_acc == uniprot_acc].iterrows():
               seq_start = row['seq_start']
               seq_end = row['seq_end']
               domain_plddts.append(np.mean(tmp_dic['plddt'][seq_start:seq_end]))

        if all(domain_plddt >= 70 for domain_plddt in domain_plddts) and protein_pddt >= 70:
            best_models.append(os.path.join(AFres_directory, file.replace('.json', '.pdb').replace('scores', 'relaxed')))
    """
    best_models = []
    for file in filter(lambda x: x.endswith('.json') and 'rank_001' in x,  os.listdir(AFres_directory)):
         with open(os.path.join(AFres_directory, file), 'r') as fp:
             tmp_dic = json.load(fp)
         if np.mean(tmp_dic['plddt']) >= 70:
               best_models.append(os.path.join(AFres_directory, file.replace('.json', '.pdb').replace('scores', 'unrelaxed')))

    # get mapper of existing  files
    mapper = json.load(open(os.path.join(output_dir, '.mapper.json'), 'r'))\
                       if os.path.isfile(os.path.join(output_dir, '.mapper.json')) else {}
    
    # create new fake pdb name and assign to the old one in the mapperø
    for pdb_file in best_models:
        output_name = ''.join(random.choices(string.ascii_uppercase + string.digits, k=4)).lower()
        while output_name in mapper.values():
            output_name = ''.join(random.choices(string.ascii_uppercase + string.digits, k=4))
        old_name = pdb_file.split('/')[-1].split('.pdb')[0]
        # convert pdb to mmcif if te conversion hasn't already been done
        if old_name not in mapper.keys():
            convert_pdb_to_mmcif(pdb_file, os.path.join(output_dir, output_name+'.cif'))
            mapper[old_name] = output_name
            # save mapper only if not empty
            if mapper:
                with open(os.path.join(output_dir, '.mapper.json'), 'w') as fp:
                    json.dump(mapper, fp)

if __name__ == '__main__': 
    main()
