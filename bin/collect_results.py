#! /usr/bin/env python3

"""
    Collect AF2Fix pipeline results and tables
    Francesco Costa
"""

import argparse
from genericpath import exists
import os, json
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument(
    "--AFres_dir", 
    required=True, 
    help="Path/to/AF2 results. Should contain no_template, template_MSA, template_single_seq subfolders", 
    type=str
)

parser.add_argument(
    "--AF2Rank_res_dir", 
    required=True, 
    help="Path/to/AF2Rank results. Should contain no_template, template_MSA, template_single_seq subfolders", 
    type=str
)

parser.add_argument(
    "--database_res_dir",
    required=False,
    help="Path/to/Database calling results. Should contain target_proteins.tsv and optionally domain_lenghts.csv",
    type=str,
)

parser.add_argument(
    "--molprobity_res_dir", 
    required=True, 
    help="Path/to/molprobity results dir. Should contain no_template, template_MSA, template_single_seq subfolders", 
    type=str
)

CONDITIONS = ['no_template', 'template_MSA', 'template_single_seq']

def main():
    args = parser.parse_args()
    
    # collect db info
    protein_info = collect_db_info(args.database_res_dir)
    
    # collect AF2 results
    df = collect_AF2_results(args.AFres_dir, protein_info)

    # collect AF2Rank data
    af2Rank_results = collect_AF2_rank(args.AF2Rank_res_dir)
    
    # collect Pi-Hx interactions data
    molprobity_results = collect_molprobity(args.molprobity_res_dir)
    
    # save data
    protein_info.to_csv('protein_info_data.csv', index=False)
    df.to_csv('AF2Results.csv', index=False)
    af2Rank_results.to_csv('af2Rank_results.csv', index=False)
    molprobity_results.to_csv('molprobity_results.csv', index=False)

def collect_db_info(database_res_dir: str) -> pd.DataFrame:
    """Collect pfam and af2 databases info"""
    
    os.makedirs('tmp', exist_ok=True)
    
    # some lines contain 2 or more proteins, fix it
    with open(os.path.join(database_res_dir, 'target_proteins.tsv'), 'r') as handle:
        text_to_fix = handle.readlines()
    fixed_text = "".join(text_to_fix).replace(' PF', '\nPF')
    with open('tmp/target_proteins.tsv', 'w') as handle:
        handle.write(fixed_text)

    protein_info = pd.read_table('tmp/target_proteins.tsv', sep=' ')
    protein_info['seq_lenght'] = [len(sequence) for sequence in protein_info.sequence]
    
    # safety check
    check = len(protein_info[protein_info.seq_end > protein_info.seq_lenght])
    assert  check == 0,\
          f"{check} proteins have upper domain boundary longer than protein lenght"

    # add domain lengths info
    domain_info = pd.read_csv(os.path.join(database_res_dir, 'domain_lengths.csv'))
    protein_info = pd.merge(protein_info, domain_info.rename(columns={'domain': 'pfamA_acc'}),
     on='pfamA_acc').rename(columns={'length': 'domain_length'})
    protein_info['coverage'] = protein_info.seq_lenght / protein_info.domain_length
    
    return protein_info

def collect_AF2_results(af2_res_dir: str, protein_info: pd.DataFrame) -> pd.DataFrame:
    """Get AF2 results from conditions no_template, template_MSA, template_single_sequence"""
    # get AlphaFold data

    outlist = []
    for condition in CONDITIONS:
        for pfamA_acc in filter(lambda x: x.startswith('PF'), os.listdir(os.path.join(af2_res_dir, condition))):
            for filename in filter(lambda x:
                                x.endswith('_seed_000.json') and "_scores_rank_001_alphafold2_ptm_" in x,
                                os.listdir(os.path.join(af2_res_dir, condition, pfamA_acc))):
                uniprot_acc = filename.split('_')[0]
                with open(os.path.join(af2_res_dir, condition, pfamA_acc, filename), 'r') as handle:
                    dic = json.load(handle) 
                protein_info.sort_values(by='seq_start', ascending=True)
                if uniprot_acc in protein_info.uniprot_acc.to_list():
                    for seq_start, seq_end in zip(protein_info[protein_info.uniprot_acc == uniprot_acc].seq_start.to_list(),
                                                  protein_info[protein_info.uniprot_acc == uniprot_acc].seq_end.to_list()):
                        plddt = round(np.mean(dic['plddt']), 2)
                        domain_plddt = round(np.mean(dic['plddt'][seq_start:seq_end]), 2)
                        pae = round(np.mean(dic['pae']), 2)
                        domain_pae = round(np.mean(dic['pae'][seq_start:seq_end]), 2)
                        ptm = round(dic['ptm'], 2)

                        outlist.append([pfamA_acc, uniprot_acc, seq_start, seq_end, condition, plddt, domain_plddt,
                        pae, domain_pae, ptm])

    AF_results = pd.DataFrame(outlist, columns = ['pfamA_acc', 'uniprot_acc', 'seq_start', 'seq_end', 'condition', 
                                                'plddt', 'domain_plddt', 'pae', 'domain_pae', 'ptm'])
    
    # add best pick condition. Go on appending on outlist and overwrite AF_results
    
    # get best pick based on either domain or structure
    iter_df = AF_results[['uniprot_acc', 'seq_start', 'seq_end']].drop_duplicates()
    for uniprot_acc, seq_start, seq_end in zip(iter_df.uniprot_acc.to_list(), iter_df.seq_start.to_list(),
                                            iter_df.seq_end.to_list()):
        # consider case of no templates
        if len(AF_results[(AF_results.uniprot_acc == uniprot_acc) & (AF_results.condition == 'template_single_seq')]) > 0 and\
            len(AF_results[(AF_results.uniprot_acc == uniprot_acc) & (AF_results.condition == 'template_MSA')]) > 0:

            # get which condition is the best for domain  plddt
            
            domain_condition = 'template_single_seq'
            if AF_results[(AF_results.uniprot_acc == uniprot_acc) &\
                            (AF_results.seq_start == seq_start) &\
                            (AF_results.seq_end == seq_end) &\
                            (AF_results.condition == 'template_MSA')].domain_plddt.unique()[0] > \
                            AF_results[(AF_results.uniprot_acc == uniprot_acc) &\
                            (AF_results.seq_start == seq_start) &\
                            (AF_results.seq_end == seq_end) &\
                            (AF_results.condition == 'template_single_seq')].domain_plddt.unique()[0]:
                domain_condition = 'template_MSA'
            
            pfamA_acc, uniprot_acc, seq_start, seq_end, condition, \
                plddt, domain_plddt, pae, domain_pae, ptm = AF_results[(AF_results.uniprot_acc == uniprot_acc) &\
                            (AF_results.seq_start == seq_start) &\
                            (AF_results.seq_end == seq_end) &\
                            (AF_results.condition == domain_condition)].values[0]
            condition = "best_pick_domain_plddt"
            outlist.append([pfamA_acc, uniprot_acc, seq_start, seq_end, condition, \
                plddt, domain_plddt, pae, domain_pae, ptm])
            
            # get condition in which the structure is the best
            
            protein_condition = 'template_single_seq'
            if AF_results[(AF_results.uniprot_acc == uniprot_acc) &\
                            (AF_results.seq_start == seq_start) &\
                            (AF_results.seq_end == seq_end) &\
                            (AF_results.condition == 'template_MSA')].plddt.unique()[0] > \
                            AF_results[(AF_results.uniprot_acc == uniprot_acc) &\
                            (AF_results.seq_start == seq_start) &\
                            (AF_results.seq_end == seq_end) &\
                            (AF_results.condition == 'template_single_seq')].plddt.unique()[0]:
                
                protein_condition = 'template_MSA'
            
            pfamA_acc, uniprot_acc, seq_start, seq_end, condition, \
                plddt, domain_plddt, pae, domain_pae, ptm = AF_results[(AF_results.uniprot_acc == uniprot_acc) &\
                            (AF_results.seq_start == seq_start) &\
                            (AF_results.seq_end == seq_end) &\
                            (AF_results.condition == protein_condition)].values[0]
            condition = "best_pick_structure_plddt"
            outlist.append([pfamA_acc, uniprot_acc, seq_start, seq_end, condition, \
                plddt, domain_plddt, pae, domain_pae, ptm])

    AF_results = pd.DataFrame(outlist, columns = ['pfamA_acc', 'uniprot_acc', 'seq_start', 'seq_end', 'condition', 
                                                'plddt', 'domain_plddt', 'pae', 'domain_pae', 'ptm'])

    if 'coverage' in protein_info.columns:
        return pd.merge(AF_results, protein_info[['pfamA_acc', 'uniprot_acc', 'seq_start', 'seq_end', 'coverage',
                                                   'seq_lenght', 'sequence']], 
                                                   on=['pfamA_acc', 'uniprot_acc', 'seq_start', 'seq_end',])
    else:
        return pd.merge(AF_results, protein_info[['pfamA_acc', 'uniprot_acc', 'seq_start', 'seq_end',
                                                   'seq_lenght', 'sequence']], 
                                                   on=['pfamA_acc', 'uniprot_acc', 'seq_start', 'seq_end',])
    
def collect_AF2_rank(AF2Rank_res_dir: str) -> pd.DataFrame:
    """Collect AF2Rank data for the conditions: no_template, template_MSA, template_single_seq"""
    
    df_rank = pd.DataFrame(columns=['pfamA_acc', 'uniprot_acc', 'condition','rms_out','tm_out',
                            'plddt','ptm','composite'])

    for condition in CONDITIONS:
        for pfamA_acc in filter(lambda x: x.startswith('PF'), os.listdir(os.path.join(AF2Rank_res_dir, condition))):
            for csv_table in filter(lambda x: x.endswith('.csv'), os.listdir(os.path.join(AF2Rank_res_dir, condition, pfamA_acc))):   
                tmp_df = pd.read_csv(os.path.join(AF2Rank_res_dir, condition, pfamA_acc, csv_table), 
                )
                tmp_df = tmp_df.rename(columns={'name': 'uniprot_acc'})
                tmp_df['pfamA_acc'] = pfamA_acc
                tmp_df['condition'] = condition
                df_rank = pd.concat([df_rank, tmp_df])

    return df_rank.rename(columns={'rms_out': 'af2rank_rms_out',
                                      'tm_out': 'af2rank_tm_out', 'plddt': 'af2rank_plddt',
                                      'ptm': 'af2rank_ptm', 'composite': 'af2rank_composite'})

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

if __name__ == '__main__':
    main()