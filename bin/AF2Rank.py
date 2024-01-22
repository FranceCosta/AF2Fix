#!/usr/bin/env python3

import sys
import os
import argparse
import jax
import jax.numpy as jnp
import numpy as np
import re
import subprocess
from collections import namedtuple

parser = argparse.ArgumentParser()
parser.add_argument("name", help="name to save everything under")
parser.add_argument("--input_pdb", help="Input pdb file")
parser.add_argument("--recycles", type=int, default=1, help="Number of recycles when predicting")
parser.add_argument("--model_num", type=int, default=1, help="Which AF2 model to use")
parser.add_argument("--seed", type=int, default=0, help="RNG Seed")
parser.add_argument("--deterministic", action='store_true', help="make all data processing deterministic (no masking, etc.)")
parser.add_argument("--mask_sidechains", action='store_true', help="mask out sidechain atoms except for C-Beta")
parser.add_argument("--mask_sidechains_add_cb", action='store_true', help="mask out sidechain atoms except for C-Beta, and add C-Beta to glycines")
parser.add_argument("--seq_replacement", default='', help="Amino acid residue to fill the decoy sequence with. Default keeps target sequence")
parser.add_argument("--af2_dir", default="/nfs/research/agb/research/francesco/software/alphafold_non_docker/alphafold-2.3.2/", help="AlphaFold code and weights directory")
parser.add_argument("--output_dir", default="/nfs/research/agb/research/francesco/projects/16-05-2023_FixDomainsWithAF_v1/test", help="Output")
parser.add_argument("--tm_exec", default="/nfs/research/agb/research/francesco/software/TMscore", help="TMScore executable")

args = parser.parse_args()
sys.path.append(args.af2_dir)

from alphafold.model import model
from alphafold.model import config
from alphafold.model import data

from alphafold.data import parsers
from alphafold.data import pipeline

from alphafold.common import protein
from alphafold.common import residue_constants

"""
Read in a PDB file from a path
"""
def pdb_to_string(pdb_file):
  lines = []
  for line in open(pdb_file,"r"):
    if line[:6] == "HETATM" and line[17:20] == "MSE":
      line = "ATOM  "+line[6:17]+"MET"+line[20:]
    if line[:4] == "ATOM":
      lines.append(line)
  return "".join(lines)

"""
Compute TM Scores between two PDBs and parse outputs
pdb_pred -- The path to the predicted PDB
pdb_native -- The path to the native PDB
test_len -- run asserts that the input and output should have the same length
"""
def compute_tmscore(pdb_pred, pdb_native, test_len=True):
  tm_re = re.compile(r'TM-score[\s]*=[\s]*(\d.\d+)')
  ref_len_re = re.compile(r'Length=[\s]*(\d+)[\s]*\(by which all scores are normalized\)')
  common_re = re.compile(r'Number of residues in common=[\s]*(\d+)')
  super_re = re.compile(r'\(":" denotes the residue pairs of distance < 5\.0 Angstrom\)\\n([A-Z\-]+)\\n[" ", :]+\\n([A-Z\-]+)\\n')
  
  cmd = ([args.tm_exec, pdb_pred, pdb_native])
  tmscore_output = str(subprocess.check_output(cmd))
  try:
    tm_out = float(tm_re.search(tmscore_output).group(1))
    reflen = int(ref_len_re.search(tmscore_output).group(1))
    common = int(common_re.search(tmscore_output).group(1))
    
    seq1 = super_re.search(tmscore_output).group(1)
    seq2 = super_re.search(tmscore_output).group(1)
  except Exception as e:
    print("Failed on: " + " ".join(cmd))
    raise e

  if test_len:
    assert reflen == common, cmd
    assert seq1 == seq2, cmd
    assert len(seq1) == reflen, cmd

  return tm_out

"""
Compute aligned RMSD between two corresponding sets of poitns
true -- set of reference points. Numpy array of dimension N x 3
pred -- set of predicted points, Numpy array of dimension N x 3
"""
def jnp_rmsd(true, pred):
  def kabsch(P, Q):
    V, S, W = jnp.linalg.svd(P.T @ Q, full_matrices=False)
    flip = jax.nn.sigmoid(-10 * jnp.linalg.det(V) * jnp.linalg.det(W))
    S = flip * S.at[-1].set(-S[-1]) + (1-flip) * S
    V = flip * V.at[:,-1].set(-V[:,-1]) + (1-flip) * V
    return V@W
  p = true - true.mean(0,keepdims=True)
  q = pred - pred.mean(0,keepdims=True)
  p = p @ kabsch(p,q)
  loss = jnp.sqrt(jnp.square(p-q).sum(-1).mean() + 1e-8)
  return float(loss)

'''
Function used to add C-Beta to glycine resides
input:  3 coords (a,b,c), (L)ength, (A)ngle, and (D)ihedral
output: 4th coord
'''
def extend(a,b,c, L,A,D):
  N = lambda x: x/np.sqrt(np.square(x).sum(-1,keepdims=True) + 1e-8)
  bc = N(b-c)
  n = N(np.cross(b-a, bc))
  m = [bc,np.cross(n,bc),n]
  d = [L*np.cos(A), L*np.sin(A)*np.cos(D), -L*np.sin(A)*np.sin(D)]
  return c + sum([m*d for m,d in zip(m,d)])

"""
Create an AlphaFold model runner
name -- The name of the model to get the parameters from. Options: model_[1-5]
"""
def make_model_runner(name, recycles):
  cfg = config.model_config(name)      

  cfg.data.common.num_recycle = recycles
  cfg.model.num_recycle = recycles
  cfg.data.eval.num_ensemble = 1
  if args.deterministic:
    cfg.data.eval.masked_msa_replace_fraction = 0.0
    cfg.model.global_config.deterministic = True
  params = data.get_model_haiku_params(name, args.af2_dir + 'data/')

  return model.RunModel(cfg, params)

def empty_placeholder_template_features(num_templates, num_res):
  return {
      'template_aatype': np.zeros(
          (num_templates, num_res,
           len(residue_constants.restypes_with_x_and_gap)), dtype=np.float32),
      'template_all_atom_masks': np.zeros(
          (num_templates, num_res, residue_constants.atom_type_num),
          dtype=np.float32),
      'template_all_atom_positions': np.zeros(
          (num_templates, num_res, residue_constants.atom_type_num, 3),
          dtype=np.float32),
      'template_domain_names': np.zeros([num_templates], dtype=object),
      'template_sequence': np.zeros([num_templates], dtype=object),
      'template_sum_probs': np.zeros([num_templates], dtype=np.float32),
  }

"""
Create a feature dictionary for input to AlphaFold
runner - The model runner being invoked. Returned from `make_model_runner`
sequence - The target sequence being predicted
templates - The template features being added to the inputs
seed - The random seed being used for data processing
"""
def make_processed_feature_dict(runner, sequence, name="test", templates=None, seed=0):
  feature_dict = {}
  feature_dict.update(pipeline.make_sequence_features(sequence, name, len(sequence)))

  msa = pipeline.parsers.parse_a3m(">1\n%s" % sequence)

  feature_dict.update(pipeline.make_msa_features([msa]))

  if templates is not None:
    feature_dict.update(templates)
  else:
    feature_dict.update(empty_placeholder_template_features(num_templates=0, num_res=len(sequence)))


  processed_feature_dict = runner.process_features(feature_dict, random_seed=seed)

  return processed_feature_dict

"""
Package AlphFold's output into an easy-to-use dictionary
prediction_result - output from running AlphaFold on an input dictionary
processed_feature_dict -- The dictionary passed to AlphaFold as input. Returned by `make_processed_feature_dict`.
"""
def parse_results(prediction_result, processed_feature_dict):
  b_factors = prediction_result['plddt'][:,None] * prediction_result['structure_module']['final_atom_mask']  
  dist_bins = jax.numpy.append(0,prediction_result["distogram"]["bin_edges"])
  dist_mtx = dist_bins[prediction_result["distogram"]["logits"].argmax(-1)]
  contact_mtx = jax.nn.softmax(prediction_result["distogram"]["logits"])[:,:,dist_bins < 8].sum(-1)

  out = {"unrelaxed_protein": protein.from_prediction(processed_feature_dict, prediction_result, b_factors=b_factors),
        "plddt": prediction_result['plddt'],
        "pLDDT": prediction_result['plddt'].mean(),
        "dists": dist_mtx,
        "adj": contact_mtx}

  out.update({"pae": prediction_result['predicted_aligned_error'],
              "pTMscore": prediction_result['ptm']})
  return out

"""
Ingest a decoy protein, pass it to AlphaFold as a template, and extract the parsed output
target_seq -- the sequence to be predicted
decoy_prot -- the decoy structure to be injected as a template
model_runner -- the model runner to execute
name -- the name associated with this prediction
"""
def score_decoy(target_seq, decoy_prot, model_runner, name):
  decoy_seq_in = "".join([residue_constants.restypes[x] for x in decoy_prot.aatype]) # the sequence in the decoy PDB file

  mismatch = False
  if decoy_seq_in == target_seq:
    assert jnp.all(decoy_prot.residue_index - 1 == np.arange(len(target_seq)))
  else: # case when template is missing some residues
    if args.verbose:
      print("Sequece mismatch: {}".format(name))
    mismatch=True

    assert "".join(target_seq[i-1] for i in decoy_prot.residue_index) == decoy_seq_in 
  
  # use this to index into the template features
  template_idxs = decoy_prot.residue_index-1
  template_idx_set = set(template_idxs)

  # The sequence associated with the decoy. Always has same length as target sequence.
  decoy_seq = args.seq_replacement*len(target_seq) if len(args.seq_replacement) == 1 else target_seq

  # create empty template features
  pos = np.zeros([1,len(decoy_seq), 37, 3])
  atom_mask = np.zeros([1, len(decoy_seq), 37])

  if args.mask_sidechains_add_cb:
    pos[0, template_idxs, :5] = decoy_prot.atom_positions[:,:5]

    # residues where we have all of the key backbone atoms (N CA C)
    backbone_modelled = jnp.all(decoy_prot.atom_mask[:,[0,1,2]] == 1, axis=1)
    backbone_idx_set = set(decoy_prot.residue_index[backbone_modelled] - 1)

    projected_cb = [i-1 for i,b,m in zip(decoy_prot.residue_index, backbone_modelled, decoy_prot.atom_mask) if m[3] == 0 and b]
    projected_cb_set = set(projected_cb)
    gly_idx = [i for i,a in enumerate(target_seq) if a == "G"]
    assert all([k in projected_cb_set for k in gly_idx if k in template_idx_set and k in backbone_idx_set]) # make sure we are adding CBs to all of the glycines

    cbs = np.array([extend(c,n,ca, 1.522, 1.927, -2.143) for c, n ,ca in zip(pos[0,:,2], pos[0,:,0], pos[0,:,1])])

    pos[0, projected_cb, 3] = cbs[projected_cb]
    atom_mask[0, template_idxs, :5] = decoy_prot.atom_mask[:, :5]
    atom_mask[0, projected_cb, 3] = 1

    template = {"template_aatype":residue_constants.sequence_to_onehot(decoy_seq, residue_constants.HHBLITS_AA_TO_ID)[None],
                "template_all_atom_masks": atom_mask,
                "template_all_atom_positions":pos,
                "template_domain_names":np.asarray(["None"])}
  elif args.mask_sidechains:
    pos[0, template_idxs, :5] = decoy_prot.atom_positions[:,:5]
    atom_mask[0, template_idxs, :5] = decoy_prot.atom_mask[:,:5]

    template = {"template_aatype":residue_constants.sequence_to_onehot(decoy_seq, residue_constants.HHBLITS_AA_TO_ID)[None],
                "template_all_atom_masks": atom_mask,
                "template_all_atom_positions":pos,
                "template_domain_names":np.asarray(["None"])}
  else:
    pos[0, template_idxs] = decoy_prot.atom_positions
    atom_mask[0, template_idxs] = decoy_prot.atom_mask

    template = {"template_aatype":residue_constants.sequence_to_onehot(decoy_seq, residue_constants.HHBLITS_AA_TO_ID)[None],
                "template_all_atom_masks":decoy_prot.atom_mask[None],
                "template_all_atom_positions":decoy_prot.atom_positions[None],
                "template_domain_names":np.asarray(["None"])}

  features = make_processed_feature_dict(model_runner, target_seq, name=name, templates=template, seed=args.seed)
  result = parse_results(model_runner.predict(features, random_seed=args.seed), features)
  return result


def compute_results(af_result, prot_native):
  plddt = round(np.mean(af_result['plddt']), 3)
  ptm = float(af_result["pTMscore"])
  rms_out = jnp_rmsd(prot_native.atom_positions[:,1,:], af_result['unrelaxed_protein'].atom_positions[:,1,:])
  tm_out = compute_tmscore(os.path.join(args.output_dir, args.name+'.pdb'), args.input_pdb, test_len=True)
  return [args.name, rms_out, tm_out, plddt, ptm, tm_out*plddt*ptm]


def main():
  os.makedirs(args.output_dir, exist_ok=True)
  prot = protein.from_pdb_string(pdb_to_string(args.input_pdb))
  name = args.input_pdb.split('/')[-1].split('.pdb')[0]
  seq = "".join([residue_constants.restypes[x] for x in prot.aatype])
  # init model
  runner = make_model_runner(f'model_{args.model_num}_ptm', args.recycles)
  # run AF
  af_result = score_decoy(seq, prot, runner, name)
  # write output pdb
  pdb_lines = protein.to_pdb(af_result["unrelaxed_protein"])
  pdb_out_path = os.path.join(args.output_dir, args.name + ".pdb")
  with open(pdb_out_path, 'w') as f:
    f.write(pdb_lines)
  # get metrics
  results = compute_results(af_result, prot)
  # save results
  with open(os.path.join(args.output_dir, args.name+'.csv'), 'w') as handle:
    handle.write('name,rms_out,tm_out,plddt,ptm,composite\n')
    handle.write(",".join([str(i) for i in results])+'\n')

if __name__ == "__main__":
  main()
