#!/usr/bin/python3


"""
    Script provided by Jake Kerrison on 02/09/2023 to perform Xh-Pi analysis
    modified by Francesco Costa on 04/09/2023
    modifications to original script:
        - add Argparse to manage input arguments
        - add h_angle_cutoff as input argument to check_h_geom
        - create main
        - allow using either input PDB or directory containing PDB files
        - change .report to .csv
        - disallow pymol output by default
        - check presence of input directory or pdb file
        - control if hgen is executable
        - create output directory
        - create summary.csv
"""

import shutil
import sys  # for sys.argv[]
import math  # for arcos
import os  # for getting the current working dir - this is used to enable the .pml file to open the .pdb file
import statistics  # for the mean() function
import argparse
from io import TextIOWrapper
import subprocess

def isfile(x):
    """
    'Type' for argparse - checks that file or dir exists
    """
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

parser = argparse.ArgumentParser(
    description="Provide either PDB_file or PDB_folder"

)

parser.add_argument(
    "--PDB_file", 
    required=False, 
    help="PDB file to analyze", 
    type=isfile
)

parser.add_argument(
    "--PDB_folder", 
    required=False, 
    help="Path/to/folder containing PDB files to analyze", 
    type=isfile
)

parser.add_argument(
    "--output", 
    required=False, 
    help="Path/to/output folder",
    default='xhpi_output', 
    type=str
)

parser.add_argument(
    "--summary", 
    required=False, 
    help="If --summary and --PDB_folder PDB_FOLDER, creates a summary.csv with all .csv found in the output folder",
    action='store_true'
)

parser.add_argument(
    "--dist_cutoff", 
    required=False, 
    help="Interaction distance cutoff taken from satisfying H-bond potentials", 
    type=str,
    default=4.3
)

parser.add_argument(
    "--xtheta_cutoff", 
    required=False, 
    help="Interaction angle cutoff taken from satisfying H-bond potentials", 
    type=str,
    default=25
)

parser.add_argument(
    "--h_angle_cutoff", 
    required=False, 
    help="120deg is what is used for the cannonical H bond", 
    type=str,
    default=120
)

parser.add_argument(
    "--plddt_cutoff", 
    required=False, 
    help="if either of the species have a pLDDT confidence score less than this, it wont be reported", 
    type=str,
    default=0
)

parser.add_argument(
    "--unique_trp", 
    required=False, 
    help="T means that the script will only report on the interaction with the smallest Xdist\n\
    F means that the script will report on both rings, if an X atom is in the geometric limits to\n\
    interact with both",
    type=str,
    default="F"
)

parser.add_argument(
    "--pymol_output", 
    required=False, 
    help="If --pymol_output, a pymol output is written for each input file",
    action='store_true'
)

# defining some lists and stuffs #######################################################################################


x_atoms = ["C", "N", "O", "S"]

x_not_atom_id = ["C  ", "O  "]

aro = ["PHE", "TYR", "HIS", "TRP"]

# building the dictionary used for finding hydrogen atoms given the x atom and x residue.
# so far, included hydrogen nomenclatures are:
#   HGEN

# more nomenclatures can be added by adding hydrogen atom IDs to the lists below, as long as there is no overlap.


ala_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B", "H3B"], "m"]
             }

arg_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "CG": [["H1G", "H2G"], "a"],
              "CD": [["H1D", "H2D"], "a"],
#              "NE": [["H1E"], ""],
#              "NH1": [["H1H1", "H2H1"], ""],
#              "NH2": [["H1H2", "H2H2"], ""]
             }

asn_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "ND2": [["H1D2", "H2D2"], "n"],
             }

asp_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
#              "ND2": [["H1D2", "H2D2"], ""],
             }

cys_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "SG": [["H1G"], "s"]
             }

gln_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "CG": [["H1G", "H2G"], "a"],
              "NE2": [["H1E2", "H2E2"], "n"]
             }

glu_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "CG": [["H1G", "H2G"], "a"]
             }

gly_h_dict = {"CA": [["H1A", "H2A"], "a"],
              "N": [["H1"], "n"]
             }

his_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "ND1": [["H1D"], "r"],
              "CD2": [["H1D2"], "r"],
              "CE1": [["H1E1"], "r"],
              "NE2": [["H1E2"], "r"]
             }

ile_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B"], "a"],
              "CG1": [["H1G1", "H2G1"], "a"],
              "CG2": [["H1G2", "H2G2", "H3G2"], "m"],
              "CD1": [["H1D1", "H2D1", "H3D1"], "m"]
             }

leu_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "CG": [["H1G"], "a"],
              "CD1": [["H1D1", "H2D1", "H3D1"], "m"],
              "CD2": [["H1D2", "H2D2", "H3D2"], "m"]
             }

lys_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "CG": [["H1G", "H2G"], "a"],
              "CD": [["H1D", "H2D"], "a"],
              "CE": [["H1E", "H2E"], "a"],
#              "NZ": [["H1Z", "H2Z", "H3Z"], ""]
             }

met_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "CG": [["H1G", "H2G"], "a"],
              "CE": [["H1E", "H2E", "H3E"], "m"]
             }

phe_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "CD1": [["H1D1"], "r"],
              "CD2": [["H1D2"], "r"],
              "CE1": [["H1E1"], "r"],
              "CE2": [["H1E2"], "r"],
              "CZ": [["H1Z"], "r"]
             }

pro_h_dict = {"CA": [["H1A"], "a"],
              "CB": [["H1B", "H2B"], "a"],
              "CG": [["H1G", "H2G"], "a"],
              "CD": [["H1D", "H2D"], "a"],
             }

ser_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "OG": [["H1G"], "o"]
             }

thr_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B"], "a"],
              "OG1": [["H1G1"], "o"],
              "CG2": [["H1G2", "H2G2", "H3G2"], "m"]
             }

trp_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "CD1": [["H1D1"], "r"],
              "NE1": [["H1E1"], "r"],
              "CE3": [["H1E3"], "r"],
              "CZ2": [["H1Z2"], "r"],
              "CZ3": [["H1Z3"], "r"],
              "CH2": [["H1H2"], "r"],
             }

tyr_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B", "H2B"], "a"],
              "CD1": [["H1D1"], "r"],
              "CD2": [["H1D2"], "r"],
              "CE1": [["H1E1"], "r"],
              "CE2": [["H1E2"], "r"],
              "OH": [["H1H"], "r"]
             }

val_h_dict = {"CA": [["H1A"], "a"],
              "N": [["H1"], "n"],
              "CB": [["H1B"], "a"],
              "CG1": [["H1G1", "H2G1", "H3G1"], "m"],
              "CG2": [["H1G2", "H2G2", "H3G2"], "m"]
             }

h_dict = {"ALA": ala_h_dict,
          "ARG": arg_h_dict,
          "ASN": asn_h_dict,
          "ASP": asp_h_dict,
          "CYS": cys_h_dict,
          "GLN": gln_h_dict,
          "GLU": glu_h_dict,
          "GLY": gly_h_dict,
          "HIS": his_h_dict,
          "ILE": ile_h_dict,
          "LEU": leu_h_dict,
          "LYS": lys_h_dict,
          "MET": met_h_dict,
          "PHE": phe_h_dict,
          "PRO": pro_h_dict,
          "SER": ser_h_dict,
          "THR": thr_h_dict,
          "TRP": trp_h_dict,
          "TYR": tyr_h_dict,
          "VAL": val_h_dict,
          }


# h_match_dict {x_res : {x_atom: [h_atoms]}}

# The functions ########################################################################################################


def read_pdb(file, plddt_cutoff):  # returns pdbdict
    # {ChainID_ResNum : [Res_ID, ResNum, chainID, {Atom_ID: [AtomNum, AtomType, [x,y,z], bfactor]}]}
    pdb = open(file, encoding="ISO-8859-1").readlines()  # some .pdb files have encoding that is not utf-8 for some
                                                         # reason...
    pdb_dict = {}

    for x in pdb:
        line = x.split()

        if len(line) > 5 and x[:4] == "ATOM" and x[17:20].strip() in h_dict.keys():
            # if line contains atom information and has a res ID found in h_dict.keys() i.e. its an amino acid and not h2o

            try:
                atom_num = int(x[7:11])
                atom_id = x[12:16].strip()
                res_id = x[17:20]
                chain_id = x[21:22]
                res_num = int(x[22:26])
                x_coord = float(x[30:38])
                y_coord = float(x[38:46])
                z_coord = float(x[46:54])
                coords = [x_coord, y_coord, z_coord]
                bfactor = float(x[61:66])
                atom_type = x[77:78]

                key = chain_id + "_" + str(res_num)

                if res_id in aro:

                    if key not in pdb_dict.keys():
                        pdb_dict[key] = [res_id, res_num, chain_id, {}]

                    if atom_id not in pdb_dict[key][3].keys():
                        pdb_dict[key][3][atom_id] = [atom_num, atom_type, coords, bfactor]

                elif bfactor > plddt_cutoff:

                    if key not in pdb_dict.keys():
                        pdb_dict[key] = [res_id, res_num, chain_id, {}]

                    if atom_id not in pdb_dict[key][3].keys():
                        pdb_dict[key][3][atom_id] = [atom_num, atom_type, coords, bfactor]

            except ValueError:  # skip the line if things arent where expected
                print("Poorly formatted line, skipping......\n{}" .format(x))
                continue

    return pdb_dict


def get_resolution(file):
    # pdb_id = file[:4]
    # res_file = open("pdbtosp.txt").readlines()
    # for x in res_file:
    #     if len(x) > 4 and x[:4] == pdb_id:
    #         resolution = x[16:21]
    #         print(pdb_id,": ",resolution)
    
    pdb = open(file, encoding="ISO-8859-1").readlines()
    resolution = "XXXX"

    for x in pdb:
        line = x.split()

        if len(line) > 2:
            if line[0] == "REMARK" and line[2] == "RESOLUTION.":

                try:
                    resolution = float(line[3])
                    break

                except ValueError:
                    print("failed to find resolution. Is it under REMARK 3?\nContinuing.....")
                    continue

    print(resolution)
    if resolution == "XXXX":
        print("failed to find resolution. Is it under REMARK 3?\nContinuing.....\n\n")
    return resolution;


def get_coords(ResNum, ChainID, AtomID, Dictionary):
    key = (ChainID + "_" + str(ResNum))
    try:
        coords = Dictionary[key][3][AtomID][2]
        return (coords)
    except KeyError:
        return [9999999, 9999999, 9999999]
        print("Could not find H coords for residue number {}, atom {}\nSkipping......\n\n" .format(key, AtomID))
    # if we cant find the hydrogen atom for whatever reason, return coords that obviously arent going to be in range of
    # the aromatic ring
    # pretty botch solution, but it works.


def check_h_geom(h_coords, x_coords, rc_coords, h_angle_cutoff):
    xh_vect = get_vector(x_coords, h_coords)
    xh_mag = get_magnitude(x_coords, h_coords)
    rch_vect = get_vector(rc_coords, h_coords)
    rch_mag = get_magnitude(rc_coords, h_coords)
    xh_rch_dotprod = get_dot_product(xh_vect, rch_vect)
    xh_rch_angle = get_angle(xh_rch_dotprod, xh_mag, rch_mag)  # this is the XHangle (X--H--ringcentre)

    if xh_rch_angle > h_angle_cutoff and xh_mag < 10:
        return "T"
    else:
        return "F"


def get_x_group(x_res_id, x_atom_id):
    group = h_dict[x_res_id][x_atom_id][1]
    return group

def get_sub_dict_res(res, dict):
    newdict = {}
    for i in dict.keys():
        id = dict[i][0]
        if id == res:
            if i not in newdict.keys():
                newdict[i] = dict[i]

    return newdict

def get_vector(a_coords, b_coords):  # returns a vector between 2 sets of coordinates a->b
    x_vect = (b_coords[0] - a_coords[0])
    y_vect = (b_coords[1] - a_coords[1])
    z_vect = (b_coords[2] - a_coords[2])

    vect = [x_vect, y_vect, z_vect]
    return vect


def get_magnitude(a_coords, b_coords):  # returns the magnitude of the vecotr between 2 sets of coordinates
    # |c| =  sqrt (cx^2 + cy^2 + cz^2)

    mag = math.sqrt(((a_coords[0] - b_coords[0]) ** 2) + ((a_coords[1] - b_coords[1]) ** 2) + (
            (a_coords[2] - b_coords[2]) ** 2))

    return mag


def get_midpoint(a_coords, b_coords):  # returns coordinates correponding to the midpoint of coords a and b

    mid = [(a_coords[0] + ((b_coords[0] - a_coords[0]) / 2)),
           (a_coords[1] + ((b_coords[1] - a_coords[1]) / 2)),
           (a_coords[2] + ((b_coords[2] - a_coords[2]) / 2))]

    return mid


def get_dot_product(a_vector, b_vector):  # returns the dot product (a.b) given two vectors
    # a.b = (ax x bx) + (ay x by) + (az x bz)

    dot_prod = (a_vector[0] * b_vector[0]) + (a_vector[1] * b_vector[1]) + (a_vector[2] * b_vector[2])
    return dot_prod


def get_angle(dot_prod, a_mag, b_mag):  # returns the angle in degrees given a dot product (a.b) and the magnitude of
    # vectors a and b (|a| and |b|)

    #    a.b = |a|*|b|*cos(theta)
    #    theta = arcos(a.b / |a|*|b|)
    try:
        cos_angle = (dot_prod / (a_mag * b_mag))
        theta = math.acos(cos_angle)
        theta = theta * 180 / math.pi  # math.acos returns a radian value
        return theta  # deg = rad x 180/pi

    except ValueError:
        theta = int(999)
        return theta

def get_coord_angle(coord1, coord2, coord3):  # returns the angle in deg made by coord1 ---- coord2
    vect1 = get_vector(coord1, coord2)  # \
    vect2 = get_vector(coord3, coord2)  # \
    mag1 = get_magnitude(coord1, coord2)  # \
    mag2 = get_magnitude(coord3, coord2)  # coord3
    dot_prod = get_dot_product(vect1, vect2)
    angle = get_angle(dot_prod, mag1, mag2)

    return (angle)


def get_cross_product_vector(a_vector, b_vector):
    cross_prod = [(a_vector[1] * b_vector[2]) - (a_vector[2] * b_vector[1]),
                  (a_vector[2] * b_vector[0]) - (a_vector[0] * b_vector[2]),
                  (a_vector[0] * b_vector[1]) - (a_vector[1] * b_vector[0])]

    return cross_prod


def get_vector_magnitude(vector):
    mag = math.sqrt((vector[0] ** 2) + (vector[1] ** 2) + (vector[2] ** 2))

    return mag


# returns the vector with the length, length
def get_normalised_vector(vector, magnitude, length):
    ratio = magnitude / length

    scaled = [vector[0] / ratio,
              vector[1] / ratio,
              vector[2] / ratio]

    return scaled


def write_data(file_to_write_to, x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, x_chain, pdb_file, 
            resolution, pi_res_id, pi_res_num, pi_bfactor, pi_chain, x_dist_mag, theta_angle, planar_angle, x_plane_mag,
            rc_xplane_mag, x_pos, x_group):  
    # a function that writes the current variables (below in .format()) to the given file
    file_to_write_to.write(
        "{},{},{},{},{},{},{},{},{},{},{},{:.4},{},{:.4},{:.4},{:.4},{:.4},{:.4},{},{}\n"
            .format(x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, x_chain, pdb_file, resolution,
                    pi_res_id, pi_res_num, pi_bfactor, pi_chain, x_dist_mag, theta_angle, planar_angle, x_plane_mag,
                    rc_xplane_mag, x_pos, x_group))


def write_pymol(outfile, pi_residue_number, x_residue_number):  # a function that writes pymol code that shows the 2 residues as
    # sticks, and colors them
    # test if outfile is open
    if isinstance(outfile, TextIOWrapper):
        outfile.write("select H_acceptor_{}, resi {}\n".format(pi_residue_number, pi_residue_number))
        outfile.write(
            "show sticks, H_acceptor_{};color red, H_acceptor_{}\n".format(pi_residue_number, pi_residue_number))

        outfile.write("select H_donor_{}, resi {}\n".format(x_residue_number, x_residue_number))
        outfile.write("show sticks, H_donor_{}; color blue, H_donor_{}\n".format(x_residue_number, x_residue_number))

def add_hydrogens(pdb_file:str, output_folder: str) -> str:
    """Add hydrogens to PDB file and create filename_H.pdb file in output_folder/ directory
    input: pdb_file: path to .pdb file
    input: output: folder
    output: path to file with added hydrogens
    """
    
    # strip pdb file: consider only TER and ATOM
    infile = os.path.join(output_folder, os.path.basename(pdb_file))

    new_lines = open(pdb_file, 'r').readlines()
    with open(infile, 'w') as fh:
        for line in new_lines:
            if 'ATOM' in line or 'TER' in line:
                fh.write(line)
    
    # add hydrogens with hydride
    outfile = infile.replace('.pdb', '_H.pdb')
    
    if HGEN_CALL == 'hgen':
        subprocess.run(f"{HGEN_CALL} xyzin {infile} xyzout {outfile} << EOF HYDR APPE\nEOF",
                     shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elif HGEN_CALL == "tmp/run_hgen.sh":
        with open("tmp/run_hgen.sh", "w") as fh:
                fh.write("#! /bin/bash\n")
                fh.write("source /ccp4/bin/ccp4.setup-sh\n")
                fh.write(f"hgen xyzin {infile} xyzout {outfile} << EOF HYDR APPE\nEOF")
        subprocess.run(HGEN_CALL,
                     shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return outfile
    

########################################################################################################################
def main(input_file, args):
    """"""
    # distance and angle cut off values
    # ive taken the below values from satisfying H-bond potentials

    dist_cutoff = args.dist_cutoff
    xtheta_cutoff = args.xtheta_cutoff
    h_angle_cutoff = args.h_angle_cutoff  # 120deg is what is used for the cannonical H bond
    plddt_cutoff = args.plddt_cutoff  # if either of the species have a pLDDT confidence score less than this, it wont be reported

    # how the script should handle trp interactions
    # as there are two aromatic rings in trp there are 2 ring centres etc
    # setting unique_trp = T means that the script will only report on the interaction with the smallest Xdist
    # setting unique_trp = F means that the script will report on both rings, if an X atom is in the geometric limits to
    # interact with both
    unique_trp = args.unique_trp
    output_folder = args.output
    pdb_file = add_hydrogens(input_file, output_folder)
    pdbdict = read_pdb(pdb_file, plddt_cutoff); print(f'Running {pdb_file}')
    resolution = get_resolution(pdb_file)

    # creating the output files
    report = open(pdb_file.replace("_H.pdb", ".csv"), "w")

    # writing the .backbone file header
    report.write("x_res_id,x_res_num,x_atom_id,x_atom_type,x_atom_num,x_bfactor,x_chain,pdb,resolution,pi_res_id,pi_res_num,pi_bfactor,pi_chain,Xdist,Xtheta,planar_angle,x_height,x_width,x_pos,x_group\n")

    # pymol header to make it look a bit nicer
    pymol = 'no_pymol_output'
    if args.pymol_output == True:
        pymol = open(pdb_file.replace("_H.pdb", ".pml"), "w")
        pymol.write("load " + str(os.getcwd()) + "/{}\n".format(os.path.basename(pdb_file)))
        pymol.write("bg_color white\nset depth_cue, 0\nhide all; show cartoon\nset cartoon_fancy_helices, on\ncolour "
                    "white\n set antialias = 1\n set ribbon_radius =0.07\n set gamma=1.5\n until.ray_shadows('none')\n set "
                    "ray_trace_mode=1\n")
        pymol.write("set label_font_color, black\nset label_size, -1.2\nset stick_radius, 0.22\n")

    # making sub dictionaries for each atom we need the coordinates for
    trp_dict = get_sub_dict_res('TRP', pdbdict)
    tyr_dict = get_sub_dict_res('TYR', pdbdict)
    phe_dict = get_sub_dict_res('PHE', pdbdict)
    his_dict = get_sub_dict_res('HIS', pdbdict)

    for x_res_key in pdbdict.keys():
        # PDB file is put into a dictionary organised as follows:
        # {ChainID_ResNum : [Res_ID, ResNum, chainID, {Atom_ID: [AtomNum, AtomType, [x,y,z], bfactor]}]}
        # atoms are organised into residues

        for x_atom_id in pdbdict[x_res_key][3].keys():

            x_res_id = pdbdict[x_res_key][0]
            x_bfactor = pdbdict[x_res_key][3][x_atom_id][3]

            # print(x_res_id +  str(x_res_num))

            if x_res_id in h_dict.keys() and x_atom_id in h_dict[x_res_id].keys() and float(x_bfactor) > plddt_cutoff:

                x_atom_coords = pdbdict[x_res_key][3][x_atom_id][2]
                x_res_num = pdbdict[x_res_key][1]
                x_chain = pdbdict[x_res_key][2]
                x_atom_type = pdbdict[x_res_key][3][x_atom_id][1]
                x_atom_num = pdbdict[x_res_key][3][x_atom_id][0]

                x_group = get_x_group(x_res_id, x_atom_id)

                for pi_key in tyr_dict.keys():

                    pi_chain = tyr_dict[pi_key][2]
                    pi_res_num = tyr_dict[pi_key][1]
                    pi_res_id = tyr_dict[pi_key][0]

                    if pi_key != x_res_key:  # only continue if the x atom is not part of the aromatic residue hydrogen acceptor

                        # {ChainID_ResNum : [Res_ID, ResNum, chainID, {Atom_ID: [AtomNum, AtomType, [x,y,z], bfactor]}]}

                        try:
                            #  defining all the atom coords in the current tyr residue
                            tyr_cg = tyr_dict[pi_key][3]['CG'][2]
                            tyr_cd1 = tyr_dict[pi_key][3]['CD1'][2]
                            tyr_cd2 = tyr_dict[pi_key][3]['CD2'][2]
                            tyr_ce1 = tyr_dict[pi_key][3]['CE1'][2]
                            tyr_ce2 = tyr_dict[pi_key][3]['CE2'][2]
                            tyr_cz = tyr_dict[pi_key][3]['CZ'][2]

                        except KeyError:
                            continue

                        # defining all the ring atom bfactors to calculate the average
                        tyr_cg_bfactor = tyr_dict[pi_key][3]['CG'][3]
                        tyr_cd1_bfactor = tyr_dict[pi_key][3]['CD1'][3]
                        tyr_cd2_bfactor = tyr_dict[pi_key][3]['CD2'][3]
                        tyr_ce1_bfactor = tyr_dict[pi_key][3]['CE1'][3]
                        tyr_ce2_bfactor = tyr_dict[pi_key][3]['CE2'][3]
                        tyr_cz_bfactor = tyr_dict[pi_key][3]['CZ'][3]

                        pi_bfactor = (
                                (tyr_cg_bfactor + tyr_cd1_bfactor + tyr_cd2_bfactor + tyr_ce1_bfactor + tyr_ce2_bfactor
                                + tyr_cz_bfactor) / 6)

                        # calculating the x,y,z coords of the ring centre
                        tyr_centre = [(tyr_cg[0] + ((tyr_cz[0] - tyr_cg[0]) / 2)),
                                    (tyr_cg[1] + ((tyr_cz[1] - tyr_cg[1]) / 2)),
                                    (tyr_cg[2] + ((tyr_cz[2] - tyr_cg[2]) / 2))]

                        x_dist_mag = get_magnitude(tyr_centre, x_atom_coords)

                        if x_dist_mag < dist_cutoff and pi_bfactor > plddt_cutoff:  # only do this lot of maths if distance is within cutoff
                            x_dist_vect = get_vector(tyr_centre, x_atom_coords)  # need this for theta

                            # checking that the H atom is close to the ring centre
                            h_atom_ids = h_dict[x_res_id][x_atom_id][0]

                            h_coords = []
                            for h in h_atom_ids:
                                h_coords.append(get_coords(x_res_num, x_chain, h, pdbdict))

                            h_geom = []
                            for h_coord in h_coords:
                                h_geom.append(check_h_geom(h_coord, x_atom_coords, tyr_centre, h_angle_cutoff))

                            # print(str(x_res_num) + "    " + str(x_atom_id))
                            # print(h_geom)

                            # continue if the h atom is withing the geometric cutoffs
                            if "T" in h_geom:

                                # print(x_res_id + "    " + x_atom_id)
                                # print(h_atom_keys)

                                # generating 2 vectors on the aromatic plane (required for calculating cross product (ring
                                # normal))
                                cg_cz_vect = get_vector(tyr_cg, tyr_cz)
                                cd1_ce2_vect = get_vector(tyr_cd1, tyr_ce2)

                                # getting the magnitudes of these vectors (required for calculating cross product magnitude)
                                cg_cz_mag = get_magnitude(tyr_cg, tyr_cz)
                                cd1_ce2_mag = get_magnitude(tyr_cd1, tyr_ce2)

                                # calculating the angle between these vectors (required for calculating cross product magnitude)
                                # cgcz_cd1ce2_dot_prod = get_dot_product(cg_cz_vect, cd1_ce2_vect)
                                # cgcz_cd1ce2_angle = get_angle(cgcz_cd1ce2_dot_prod, cg_cz_mag, cd1_ce2_mag)

                                # finally calculating cross product goodness (defining ring normal vector and magnitude)
                                tyr_normal_vect = get_cross_product_vector(cg_cz_vect, cd1_ce2_vect)
                                tyr_normal_mag = get_vector_magnitude(tyr_normal_vect)

                                normal_dot_prod = get_dot_product(tyr_normal_vect, cg_cz_vect)
                                normal_angle = get_angle(normal_dot_prod, tyr_normal_mag, cg_cz_mag)
                                # print("normal angle = " + str(normal_angle))

                                # calculating the angle (theta) between x-ring centre, and ring normal
                                theta_dot_prod = get_dot_product(tyr_normal_vect, x_dist_vect)
                                theta_angle = get_angle(theta_dot_prod, x_dist_mag, tyr_normal_mag)

                                if theta_angle > 90:  # if the normal is facing away from the X atom.. (if the atom is below the plane)
                                    theta_angle = 180 - theta_angle

                                if theta_angle < xtheta_cutoff:  # carry on doing lots of fun maths, if we still are in xhpi
                                    # geometry

                                    # calculate the distance between the ring plane and X atom
                                    # the angle between the aromatic plane and X-ring center vector
                                    plane_x_angle = 90 - theta_angle

                                    # x_plane_mag is the shortest distance between the plane and the X atom
                                    x_plane_mag = x_dist_mag * math.sin((plane_x_angle * (math.pi / 180)))

                                    # Defining the posiiton on the aromatic plane that is directly below the X atom
                                    # This is done by calculating the vector between the plane and the X atom, and then
                                    # moving the X atom along this vector, by x_plane_mag (distance between the X atom and
                                    # the plane)

                                    # making a vector parallel to the normal, with a length of x_plane_mag (distance between
                                    # plane and x atom)
                                    # Can then move the X atom coords by this to get the position on the plane as described above
                                    scaled_normal = get_normalised_vector(tyr_normal_vect, tyr_normal_mag, x_plane_mag)

                                    # We dont know if the X atom is above or below the plane at this point
                                    # Thus we dont know whether to add or subtract the scaled normal from the X atom
                                    # To get around this, both are done, and then the angle between the normal vector and the
                                    # vector connecting the ring centre the the planar X position is calculated for both cases
                                    # this should be 90deg (dot product = 0) for the case where the coordinates hace been moved
                                    # onto the plane.

                                    xplane_coords_sum = [(x_atom_coords[0] + scaled_normal[0]),
                                                        (x_atom_coords[1] + scaled_normal[1]),
                                                        (x_atom_coords[2] + scaled_normal[2])]

                                    rc_xplane_vect_sum = get_vector(tyr_centre, xplane_coords_sum)

                                    xplane_coords_sub = [(x_atom_coords[0] - scaled_normal[0]),
                                                        (x_atom_coords[1] - scaled_normal[1]),
                                                        (x_atom_coords[2] - scaled_normal[2])]

                                    rc_xplane_vect_sub = get_vector(tyr_centre, xplane_coords_sub)

                                    # can find out which one is on the aromatic plane now
                                    # the dot product between plane normal and the vector from the ring centre to the x plane
                                    # coords should be 0
                                    # rounding errors mean its just a very small number close to 0, so we just take the smallest

                                    dot_prod_sum = get_dot_product(rc_xplane_vect_sum, tyr_normal_vect)
                                    dot_prod_sub = get_dot_product(rc_xplane_vect_sub, tyr_normal_vect)

                                    # print(str(dot_prod_sum) + " " + str(dot_prod_sub))

                                    if (dot_prod_sum ** 2) < (dot_prod_sub ** 2):
                                        xplane_coords = xplane_coords_sum
                                        rc_xplane_vect = rc_xplane_vect_sum
                                        x_pos = "below"

                                    else:
                                        xplane_coords = xplane_coords_sub
                                        rc_xplane_vect = rc_xplane_vect_sub
                                        x_pos = "above"

                                    rc_xplane_mag = get_magnitude(tyr_centre, xplane_coords)

                                    # calculate the angle between the vector cg-cz, and the vector connecting the ring center
                                    # to the point on the plane directly below the x atom
                                    planar_dot_prod = get_dot_product(rc_xplane_vect, cg_cz_vect)
                                    planar_angle = get_angle(planar_dot_prod, rc_xplane_mag, cg_cz_mag)

                                    # the above provides the shortest angle, and looses information about whether it should be
                                    # clockwise or anticlockwise
                                    # here the clockwise side is defined as the side with cd1 and ce1 atoms
                                    # anticlockwise is the side with cd2, and ce2

                                    # define cd1 ce1 midpoint and cd2 ce2 midpoint

                                    mid1 = get_midpoint(tyr_cd1, tyr_ce1)
                                    mid2 = get_midpoint(tyr_cd2, tyr_ce2)

                                    mid1_mid2_vect = get_vector(mid1, mid2)
                                    mid1_mid2_mag = get_magnitude(mid1, mid2)

                                    mid2_mid1_vect = get_vector(mid2, mid1)
                                    mid2_mid1_mag = get_magnitude(mid2, mid1)

                                    mid21_rcxplane_dotprod = get_dot_product(rc_xplane_vect, mid2_mid1_vect)

                                    mid12_rcxplane_dotprod = get_dot_product(rc_xplane_vect, mid1_mid2_vect)

                                    mid21_rcxplane_angle = get_angle(mid21_rcxplane_dotprod, rc_xplane_mag, mid2_mid1_mag)
                                    mid12_rcxplane_angle = get_angle(mid12_rcxplane_dotprod, rc_xplane_mag, mid1_mid2_mag)

                                    # if this is then larger than phi, it is anticlockwise
                                    # clockwise = 360 - anticlockwise
                                    if mid12_rcxplane_angle < mid21_rcxplane_angle:
                                        planar_angle = 360 - planar_angle

                                    write_data(report, x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, 
                                               x_chain, os.path.basename(pdb_file).replace("_H.pdb", ""), resolution, pi_res_id, pi_res_num, pi_bfactor, pi_chain, 
                                               x_dist_mag, theta_angle, planar_angle, x_plane_mag, rc_xplane_mag, x_pos, x_group)
                                    write_pymol(pymol, pi_res_num, x_res_num)

                # same again for phe
                # both are 6 membered rings so the geometry maths is all the same
                for pi_key in phe_dict.keys():

                    pi_chain = phe_dict[pi_key][2]
                    pi_res_num = phe_dict[pi_key][1]
                    pi_res_id = phe_dict[pi_key][0]

                    if pi_key != x_res_key:  # only continue if the x atom is not part of the aromatic residue hydrogen acceptor

                        # {ChainID_ResNum : [Res_ID, ResNum, chainID, {Atom_ID: [AtomNum, AtomType, [x,y,z], bfactor]}]}

                        try:
                            #  defining all the atom coords in the current tyr residue
                            phe_cg = phe_dict[pi_key][3]['CG'][2]
                            phe_cd1 = phe_dict[pi_key][3]['CD1'][2]
                            phe_cd2 = phe_dict[pi_key][3]['CD2'][2]
                            phe_ce1 = phe_dict[pi_key][3]['CE1'][2]
                            phe_ce2 = phe_dict[pi_key][3]['CE2'][2]
                            phe_cz = phe_dict[pi_key][3]['CZ'][2]

                        except KeyError:
                            #print("Could not find all ring atoms in phe" + str(pi_res_num) + " skipping...")
                            continue  # if we cant find all ring atoms, go to the next aromatic residue

                        # defining all the ring atom bfactors to calculate the average
                        phe_cg_bfactor = phe_dict[pi_key][3]['CG'][3]
                        phe_cd1_bfactor = phe_dict[pi_key][3]['CD1'][3]
                        phe_cd2_bfactor = phe_dict[pi_key][3]['CD2'][3]
                        phe_ce1_bfactor = phe_dict[pi_key][3]['CE1'][3]
                        phe_ce2_bfactor = phe_dict[pi_key][3]['CE2'][3]
                        phe_cz_bfactor = phe_dict[pi_key][3]['CZ'][3]

                        pi_bfactor = (
                                (phe_cg_bfactor + phe_cd1_bfactor + phe_cd2_bfactor + phe_ce1_bfactor + phe_ce2_bfactor
                                + phe_cz_bfactor) / 6)

                        phe_centre = []  # calculating the x,y,z coords of the ring centre
                        phe_centre.append(phe_cg[0] + ((phe_cz[0] - phe_cg[0]) / 2))
                        phe_centre.append(phe_cg[1] + ((phe_cz[1] - phe_cg[1]) / 2))
                        phe_centre.append(phe_cg[2] + ((phe_cz[2] - phe_cg[2]) / 2))

                        x_dist_mag = get_magnitude(phe_centre, x_atom_coords)

                        if x_dist_mag < dist_cutoff and pi_bfactor > plddt_cutoff:  # only do this lot of maths if distance is within cutoff
                            x_dist_vect = get_vector(phe_centre, x_atom_coords)  # need this for theta

                            # checking that the H atom is close to the ring centre
                            h_atom_ids = h_dict[x_res_id][x_atom_id][0]

                            h_coords = []
                            for h in h_atom_ids:
                                h_coords.append(get_coords(x_res_num, x_chain, h, pdbdict))

                            h_geom = []
                            for h_coord in h_coords:
                                h_geom.append(check_h_geom(h_coord, x_atom_coords, phe_centre, h_angle_cutoff))

                            # print(str(x_res_num) + "    " + str(x_atom_id))
                            # print(h_geom)

                            # continue if the h atom is withing the geometric cutoffs
                            if "T" in h_geom:

                                # generating 2 vectors on the aromatic plane (required for calculating cross product (ring
                                # normal))
                                cg_cz_vect = get_vector(phe_cg, phe_cz)
                                cd1_ce2_vect = get_vector(phe_cd1, phe_ce2)

                                # getting the magnitudes of these vectors (required for calculating cross product magnitude)
                                cg_cz_mag = get_magnitude(phe_cg, phe_cz)
                                cd1_ce2_mag = get_magnitude(phe_cd1, phe_ce2)

                                # calculating the angle between these vectors (required for calculating cross product magnitude)
                                # cgcz_cd1ce2_dot_prod = get_dot_product(cg_cz_vect, cd1_ce2_vect)
                                # cgcz_cd1ce2_angle = get_angle(cgcz_cd1ce2_dot_prod, cg_cz_mag, cd1_ce2_mag)

                                # finally calculating cross product goodness (defining ring normal vector and magnitude)
                                phe_normal_vect = get_cross_product_vector(cg_cz_vect, cd1_ce2_vect)
                                phe_normal_mag = get_vector_magnitude(phe_normal_vect)

                                normal_dot_prod = get_dot_product(phe_normal_vect, cg_cz_vect)
                                normal_angle = get_angle(normal_dot_prod, phe_normal_mag, cg_cz_mag)
                                # print("normal angle = " + str(normal_angle))

                                # calculating the angle (theta) between x-ring centre, and ring normal
                                theta_dot_prod = get_dot_product(phe_normal_vect, x_dist_vect)
                                theta_angle = get_angle(theta_dot_prod, x_dist_mag, phe_normal_mag)

                                if theta_angle > 90:  # if the normal is facing away from the X atom.. (if the atom is below the plane)
                                    theta_angle = 180 - theta_angle

                                if theta_angle < xtheta_cutoff:  # carry on doing lots of fun maths, if we still are in xhpi
                                    # geometry

                                    # calculate the distance between the ring plane and X atom
                                    # the angle between the aromatic plane and X-ring center vector
                                    plane_x_angle = 90 - theta_angle

                                    # x_plane_mag is the shortest distance between the plane and the X atom
                                    x_plane_mag = x_dist_mag * math.sin((plane_x_angle * (math.pi / 180)))

                                    # Defining the posiiton on the aromatic plane that is directly below the X atom
                                    # This is done by calculating the vector between the plane and the X atom, and then
                                    # moving the X atom along this vector, by x_plane_mag (distance between the X atom and
                                    # the plane)

                                    # making a vector parallel to the normal, with a length of x_plane_mag (distance between
                                    # plane and x atom)
                                    # Can then move the X atom coords by this to get the position on the plane as described above
                                    scaled_normal = get_normalised_vector(phe_normal_vect, phe_normal_mag, x_plane_mag)

                                    # We dont know if the X atom is above or below the plane at this point
                                    # Thus we dont know whether to add or subtract the scaled normal from the X atom
                                    # To get around this, both are done, and then the angle between the normal vector and the
                                    # vector connecting the ring centre the the planar X position is calculated for both cases
                                    # this should be 90deg (dot product = 0) for the case where the coordinates hace been moved
                                    # onto the plane.

                                    xplane_coords_sum = [(x_atom_coords[0] + scaled_normal[0]),
                                                        (x_atom_coords[1] + scaled_normal[1]),
                                                        (x_atom_coords[2] + scaled_normal[2])]

                                    rc_xplane_vect_sum = get_vector(phe_centre, xplane_coords_sum)

                                    xplane_coords_sub = [(x_atom_coords[0] - scaled_normal[0]),
                                                        (x_atom_coords[1] - scaled_normal[1]),
                                                        (x_atom_coords[2] - scaled_normal[2])]

                                    rc_xplane_vect_sub = get_vector(phe_centre, xplane_coords_sub)

                                    # can find out which one is on the aromatic plane now
                                    # the dot product between plane normal and the vector from the ring centre to the x plane
                                    # coords should be 0
                                    # rounding errors mean its just a very small number close to 0, so we just take the smallest

                                    dot_prod_sum = get_dot_product(rc_xplane_vect_sum, phe_normal_vect)
                                    dot_prod_sub = get_dot_product(rc_xplane_vect_sub, phe_normal_vect)

                                    # print(str(dot_prod_sum) + " " + str(dot_prod_sub))

                                    if (dot_prod_sum ** 2) < (dot_prod_sub ** 2):
                                        xplane_coords = xplane_coords_sum
                                        rc_xplane_vect = rc_xplane_vect_sum
                                        x_pos = "below"

                                    else:
                                        xplane_coords = xplane_coords_sub
                                        rc_xplane_vect = rc_xplane_vect_sub
                                        x_pos = "above"

                                    rc_xplane_mag = get_magnitude(phe_centre, xplane_coords)

                                    # calculate the angle between the vector cg-cz, and the vector connecting the ring center
                                    # to the point on the plane directly below the x atom
                                    planar_dot_prod = get_dot_product(rc_xplane_vect, cg_cz_vect)
                                    planar_angle = get_angle(planar_dot_prod, rc_xplane_mag, cg_cz_mag)

                                    # the above provides the shortest angle, and looses information about whether it should be
                                    # clockwise or anticlockwise
                                    # here the clockwise side is defined as the side with cd1 and ce1 atoms
                                    # anticlockwise is the side with cd2, and ce2

                                    # define cd1 ce1 midpoint and cd2 ce2 midpoint

                                    mid1 = get_midpoint(phe_cd1, phe_ce1)
                                    mid2 = get_midpoint(phe_cd2, phe_ce2)

                                    mid1_mid2_vect = get_vector(mid1, mid2)
                                    mid1_mid2_mag = get_magnitude(mid1, mid2)

                                    mid2_mid1_vect = get_vector(mid2, mid1)
                                    mid2_mid1_mag = get_magnitude(mid2, mid1)

                                    mid21_rcxplane_dotprod = get_dot_product(rc_xplane_vect, mid2_mid1_vect)

                                    mid12_rcxplane_dotprod = get_dot_product(rc_xplane_vect, mid1_mid2_vect)

                                    mid21_rcxplane_angle = get_angle(mid21_rcxplane_dotprod, rc_xplane_mag, mid2_mid1_mag)
                                    mid12_rcxplane_angle = get_angle(mid12_rcxplane_dotprod, rc_xplane_mag, mid1_mid2_mag)

                                    # if this is then larger than phi, it is anticlockwise
                                    # clockwise = 360 - anticlockwise
                                    if mid12_rcxplane_angle < mid21_rcxplane_angle:
                                        planar_angle = 360 - planar_angle

                                    write_data(report, x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, 
                                               x_chain, os.path.basename(pdb_file).replace("_H.pdb", ""), resolution, pi_res_id, pi_res_num, pi_bfactor, pi_chain, 
                                               x_dist_mag, theta_angle, planar_angle, x_plane_mag, rc_xplane_mag, x_pos, x_group)
                                    write_pymol(pymol, pi_res_num, x_res_num)

                for pi_key in his_dict.keys():

                    pi_chain = his_dict[pi_key][2]
                    pi_res_num = his_dict[pi_key][1]
                    pi_res_id = his_dict[pi_key][0]

                    if pi_key != x_res_key:  # only continue if the x atom is not part of the aromatic residue hydrogen acceptor

                        # {ChainID_ResNum : [Res_ID, ResNum, chainID, {Atom_ID: [AtomNum, AtomType, [x,y,z], bfactor]}]}

                        try:
                            #  defining all the atom coords in the current his residue
                            his_cg = his_dict[pi_key][3]['CG'][2]
                            his_nd1 = his_dict[pi_key][3]['ND1'][2]
                            his_cd2 = his_dict[pi_key][3]['CD2'][2]
                            his_ce1 = his_dict[pi_key][3]['CE1'][2]
                            his_ne2 = his_dict[pi_key][3]['NE2'][2]

                        except KeyError:
                            continue

                        # defining all the ring atom bfactors to calculate the average
                        his_cg_bfactor = his_dict[pi_key][3]['CG'][3]
                        his_nd1_bfactor = his_dict[pi_key][3]['ND1'][3]
                        his_cd2_bfactor = his_dict[pi_key][3]['CD2'][3]
                        his_ce1_bfactor = his_dict[pi_key][3]['CE1'][3]
                        his_ne2_bfactor = his_dict[pi_key][3]['NE2'][3]

                        pi_bfactor = ((
                                            his_cg_bfactor + his_nd1_bfactor + his_cd2_bfactor + his_ce1_bfactor + his_ne2_bfactor) / 6)

                        # calculating the ring centre of his
                        # as its a 5 memebered ring, im doing this by just averaging all the coordinates

                        his_centre = [(statistics.mean([his_cg[0], his_nd1[0], his_cd2[0], his_ce1[0], his_ne2[0]])),
                                    (statistics.mean([his_cg[1], his_nd1[1], his_cd2[1], his_ce1[1], his_ne2[1]])),
                                    (statistics.mean([his_cg[2], his_nd1[2], his_cd2[2], his_ce1[2], his_ne2[2]]))]

                        x_dist_mag = get_magnitude(his_centre, x_atom_coords)

                        if x_dist_mag < dist_cutoff and pi_bfactor > plddt_cutoff:  # only do this lot of maths if distance is within cutoff
                            x_dist_vect = get_vector(his_centre, x_atom_coords)  # need this for theta

                            # checking that the H atom is close to the ring centre
                            h_atom_ids = h_dict[x_res_id][x_atom_id][0]

                            h_coords = []
                            for h in h_atom_ids:
                                h_coords.append(get_coords(x_res_num, x_chain, h, pdbdict))

                            h_geom = []
                            for h_coord in h_coords:
                                h_geom.append(check_h_geom(h_coord, x_atom_coords, his_centre, h_angle_cutoff))

                            # print(str(x_res_num) + "    " + str(x_atom_id))
                            # print(h_geom)

                            # continue if the h atom is withing the geometric cutoffs
                            if "T" in h_geom:

                                # generating 2 vectors on the aromatic plane (required for calculating cross product (ring
                                # normal))
                                cg_ce1_vect = get_vector(his_cg, his_ce1)
                                nd1_cd2_vect = get_vector(his_nd1, his_cd2)

                                # getting the magnitudes of these vectors (required for calculating cross product magnitude)
                                cg_ce1_mag = get_magnitude(his_cg, his_ce1)
                                nd1_cd2_mag = get_magnitude(his_nd1, his_cd2)

                                # calculating the angle between these vectors (required for calculating cross product magnitude)
                                # cgce1_nd1cd2_dot_prod = get_dot_product(cg_ce1_vect, nd1_cd2_vect)
                                # cgcz_cd1ce2_angle = get_angle(cgce1_nd1cd2_dot_prod, cg_ce1_mag, nd1_cd2_mag)

                                # finally calculating cross product goodness (defining ring normal vector and magnitude)
                                his_normal_vect = get_cross_product_vector(cg_ce1_vect, nd1_cd2_vect)
                                his_normal_mag = get_vector_magnitude(his_normal_vect)

                                # calculating and printing the angle between a line on the aromatic plane and the normal
                                # making sure its actually 90deg
                                # normal_dot_prod = get_dot_product(his_normal_vect, cg_ce1_vect)
                                # normal_angle = get_angle(normal_dot_prod, his_normal_mag, cg_ce1_mag)
                                # print("normal angle = " + str(normal_angle))

                                # calculating the angle (theta) between x-ring centre, and ring normal
                                theta_dot_prod = get_dot_product(his_normal_vect, x_dist_vect)
                                theta_angle = get_angle(theta_dot_prod, x_dist_mag, his_normal_mag)

                                if theta_angle > 90:  # if the normal is facing away from the X atom.. (if the atom is below the plane)
                                    theta_angle = 180 - theta_angle

                                if theta_angle < xtheta_cutoff:  # carry on doing lots of fun maths, if we still are in xhpi
                                    # geometry

                                    # calculate the distance between the ring plane and X atom
                                    # the angle between the aromatic plane and X-ring center vector
                                    plane_x_angle = 90 - theta_angle

                                    # x_plane_mag is the shortest distance between the plane and the X atom
                                    x_plane_mag = x_dist_mag * math.sin((plane_x_angle * (math.pi / 180)))

                                    # Defining the posiiton on the aromatic plane that is directly below the X atom
                                    # This is done by calculating the vector between the plane and the X atom, and then
                                    # moving the X atom along this vector, by x_plane_mag (distance between the X atom and
                                    # the plane)

                                    # making a vector parallel to the normal, with a length of x_plane_mag (distance between
                                    # plane and x atom)
                                    # Can then move the X atom coords by this to get the position on the plane as described above
                                    scaled_normal = get_normalised_vector(his_normal_vect, his_normal_mag, x_plane_mag)

                                    # We dont know if the X atom is above or below the plane at this point
                                    # Thus we dont know whether to add or subtract the scaled normal from the X atom
                                    # To get around this, both are done, and then the angle between the normal vector and the
                                    # vector connecting the ring centre the the planar X position is calculated for both cases
                                    # this should be 90deg (dot product = 0) for the case where the coordinates hace been moved
                                    # onto the plane.

                                    xplane_coords_sum = [(x_atom_coords[0] + scaled_normal[0]),
                                                        (x_atom_coords[1] + scaled_normal[1]),
                                                        (x_atom_coords[2] + scaled_normal[2])]

                                    rc_xplane_vect_sum = get_vector(his_centre, xplane_coords_sum)

                                    xplane_coords_sub = [(x_atom_coords[0] - scaled_normal[0]),
                                                        (x_atom_coords[1] - scaled_normal[1]),
                                                        (x_atom_coords[2] - scaled_normal[2])]

                                    rc_xplane_vect_sub = get_vector(his_centre, xplane_coords_sub)

                                    # can find out which one is on the aromatic plane now
                                    # the dot product between plane normal and the vector from the ring centre to the x plane
                                    # coords should be 0
                                    # rounding errors mean its just a very small number close to 0, so we just take the smallest

                                    dot_prod_sum = get_dot_product(rc_xplane_vect_sum, his_normal_vect)
                                    dot_prod_sub = get_dot_product(rc_xplane_vect_sub, his_normal_vect)

                                    # print(str(dot_prod_sum) + " " + str(dot_prod_sub))

                                    if (dot_prod_sum ** 2) < (dot_prod_sub ** 2):
                                        xplane_coords = xplane_coords_sum
                                        rc_xplane_vect = rc_xplane_vect_sum
                                        x_pos = "below"

                                    else:
                                        xplane_coords = xplane_coords_sub
                                        rc_xplane_vect = rc_xplane_vect_sub
                                        x_pos = "above"

                                    rc_xplane_mag = get_magnitude(his_centre, xplane_coords)

                                    # calculating the planar angle for his
                                    # this is a bit different as it is not 6 memebered
                                    # thus i calculate the angle relative to the vector passing through CG, and the midpoint
                                    # between NE2, and ce1
                                    # this is where cz is closest to if it were a 6 membered ring

                                    # CE1 NE2 midpoint
                                    ce1_ne2_mid = get_midpoint(his_ce1, his_ne2)

                                    cg_mid_vect = get_vector(his_cg, ce1_ne2_mid)
                                    cg_mid_mag = get_magnitude(his_cg, ce1_ne2_mid)

                                    planar_dot_prod = get_dot_product(rc_xplane_vect, cg_mid_vect)
                                    planar_angle = get_angle(planar_dot_prod, rc_xplane_mag, cg_mid_mag)

                                    # the above provides the shortest angle, and looses information about whether it should be
                                    # clockwise or anticlockwise
                                    # here the clockwise side is defined as the side with nd1 and ce1 atoms
                                    # anticlockwise is the side with cd2, and ne2

                                    # define a ce2-ce1 vector and calculate the angle between this and the plane postion of X
                                    # define a his_rc-ce1 vector and calculate the angle between this and the plane position of X

                                    rc_ce1_vect = get_vector(his_centre, his_ce1)
                                    rc_ce1_mag = get_magnitude(his_centre, his_ce1)

                                    rc_ne2_vect = get_vector(his_centre, his_ne2)
                                    rc_ne2_mag = get_magnitude(his_centre, his_ne2)

                                    rcce1_rcxplane_dotprod = get_dot_product(rc_xplane_vect, rc_ce1_vect)
                                    rcne2_rcxplane_dotprod = get_dot_product(rc_xplane_vect, rc_ne2_vect)

                                    # print("ce2ce1_rcxplane")
                                    rcce1_rcxplane_angle = get_angle(rcce1_rcxplane_dotprod, rc_xplane_mag, rc_ce1_mag)
                                    rcne2_rcxplane_angle = get_angle(rcne2_rcxplane_dotprod, rc_xplane_mag, rc_ne2_mag)

                                    # if this is then larger than phi, it is anticlockwise
                                    # clockwise = 360 - anticlockwise
                                    if rcce1_rcxplane_angle > rcne2_rcxplane_angle:
                                        planar_angle = 360 - planar_angle

                                    write_data(report, x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, 
                                               x_chain, os.path.basename(pdb_file).replace("_H.pdb", ""), resolution, pi_res_id, pi_res_num, pi_bfactor, pi_chain, 
                                               x_dist_mag, theta_angle, planar_angle, x_plane_mag, rc_xplane_mag, x_pos, x_group)
                                    write_pymol(pymol, pi_res_num, x_res_num)

                for pi_key in trp_dict.keys():

                    trp5 = "F"
                    trp6 = "F"

                    # setting up some variables for each ring
                    # this is so i can refer to them later out of their scopes for writing the data depending on whether I
                    # want to report on interactions with both rings, or just the closest one.
                    trp5_theta_angle = float()
                    trp5_planar_angle = float()
                    trp5_x_plane_mag = float()
                    trp5_rc_xplane_mag = float()
                    trp5_x_pos = str()
                    trp6_theta_angle = float()
                    trp6_planar_angle = float()
                    trp6_x_plane_mag = float()
                    trp6_rc_xplane_mag = float()
                    trp6_x_pos = str()

                    pi_chain = trp_dict[pi_key][2]
                    pi_res_num = trp_dict[pi_key][1]
                    pi_res_id = trp_dict[pi_key][0]

                    if pi_key != x_res_key:  # only continue if the x atom is not part of the aromatic residue hydrogen acceptor

                        # {ChainID_ResNum : [Res_ID, ResNum, chainID, {Atom_ID: [AtomNum, AtomType, [x,y,z], bfactor]}]}

                        try:
                            #  defining all the atom coords in the current trp residue
                            trp_cg = trp_dict[pi_key][3]['CG'][2]
                            trp_cd1 = trp_dict[pi_key][3]['CD1'][2]
                            trp_cd2 = trp_dict[pi_key][3]['CD2'][2]
                            trp_ne1 = trp_dict[pi_key][3]['NE1'][2]
                            trp_ce2 = trp_dict[pi_key][3]['CE2'][2]
                            trp_ce3 = trp_dict[pi_key][3]['CE3'][2]
                            trp_cz2 = trp_dict[pi_key][3]['CZ2'][2]
                            trp_cz3 = trp_dict[pi_key][3]['CZ3'][2]
                            trp_ch2 = trp_dict[pi_key][3]['CH2'][2]

                        except:
                            continue

                        # defining all the ring atom bfactors to calculate the average
                        trp_cg_bfactor = trp_dict[pi_key][3]['CG'][3]
                        trp_cd1_bfactor = trp_dict[pi_key][3]['CD1'][3]
                        trp_cd2_bfactor = trp_dict[pi_key][3]['CD2'][3]
                        trp_ne1_bfactor = trp_dict[pi_key][3]['NE1'][3]
                        trp_ce2_bfactor = trp_dict[pi_key][3]['CE2'][3]
                        trp_ce3_bfactor = trp_dict[pi_key][3]['CE3'][3]
                        trp_cz2_bfactor = trp_dict[pi_key][3]['CZ2'][3]
                        trp_cz3_bfactor = trp_dict[pi_key][3]['CZ3'][3]
                        trp_ch2_bfactor = trp_dict[pi_key][3]['CH2'][3]

                        # calculating the mean bfactors for the 5 and 6 membered rings
                        trp5_b_factor = statistics.mean(
                            [trp_cg_bfactor, trp_cd1_bfactor, trp_cd2_bfactor, trp_ne1_bfactor, trp_ce2_bfactor])
                        trp6_b_factor = statistics.mean(
                            [trp_cd2_bfactor, trp_ce2_bfactor, trp_cz2_bfactor, trp_ce3_bfactor, trp_cz3_bfactor,
                            trp_ch2_bfactor])

                        # calculating the ring centre for the 5 and 6 membered rings
                        # for both cases this is done by calculating the mean of all the coordinates

                        trp5_centre = [(statistics.mean([trp_cg[0], trp_cd1[0], trp_cd2[0], trp_ne1[0], trp_ce2[0]])),
                                    (statistics.mean([trp_cg[1], trp_cd1[1], trp_cd2[1], trp_ne1[1], trp_ce2[1]])),
                                    (statistics.mean([trp_cg[2], trp_cd1[2], trp_cd2[2], trp_ne1[2], trp_ce2[2]]))]

                        trp6_centre = [
                            (statistics.mean([trp_cd2[0], trp_ce2[0], trp_cz2[0], trp_ce3[0], trp_cz3[0], trp_ch2[0]])),
                            (statistics.mean([trp_cd2[1], trp_ce2[1], trp_cz2[1], trp_ce3[1], trp_cz3[1], trp_ch2[1]])),
                            (statistics.mean([trp_cd2[2], trp_ce2[2], trp_cz2[2], trp_ce3[2], trp_cz3[2], trp_ch2[2]]))]

                        trp5_x_dist_mag = get_magnitude(trp5_centre, x_atom_coords)
                        trp6_x_dist_mag = get_magnitude(trp6_centre, x_atom_coords)

                        if trp5_x_dist_mag < dist_cutoff and trp5_b_factor > plddt_cutoff:  # only do this lot of maths if distance is within cutoff
                            trp5_x_dist_vect = get_vector(trp5_centre, x_atom_coords)  # need this for theta

                            # checking that the H atom is close to the ring centre
                            h_atom_ids = h_dict[x_res_id][x_atom_id][0]

                            h_coords = []
                            for h in h_atom_ids:
                                h_coords.append(get_coords(x_res_num, x_chain, h, pdbdict))

                            h_geom = []
                            for h_coord in h_coords:
                                h_geom.append(check_h_geom(h_coord, x_atom_coords, trp5_centre, h_angle_cutoff))

                            # print(str(x_res_num) + "    " + str(x_atom_id))
                            # print(h_geom)

                            # continue if the h atom is withing the geometric cutoffs
                            if "T" in h_geom:

                                # generating 2 vectors on the aromatic plane (required for calculating cross product (ring
                                # normal))
                                cg_ne1_vect = get_vector(trp_cg, trp_ne1)
                                cd1_cd2_vect = get_vector(trp_cd1, trp_cd2)

                                # getting the magnitudes of these vectors (required for calculating cross product magnitude)
                                cg_ne1_mag = get_magnitude(trp_cg, trp_ne1)
                                cd1_cd2_mag = get_magnitude(trp_cd1, trp_cd2)

                                # calculating the angle between these vectors (required for calculating cross product magnitude)
                                # cgne1_cd1cd2_dot_prod = get_dot_product(cg_ne1_vect, cd1_cd2_vect)
                                # cgcz_cd1ce2_angle = get_angle(cgce1_nd1cd2_dot_prod, cg_ce1_mag, nd1_cd2_mag)

                                # finally calculating cross product goodness (defining ring normal vector and magnitude)
                                trp5_normal_vect = get_cross_product_vector(cg_ne1_vect, cd1_cd2_vect)
                                trp5_normal_mag = get_vector_magnitude(trp5_normal_vect)

                                # calculating and printing the angle between a line on the aromatic plane and the normal
                                # making sure its actually 90deg
                                # normal_dot_prod = get_dot_product(his_normal_vect, cg_ce1_vect)
                                # normal_angle = get_angle(normal_dot_prod, his_normal_mag, cg_ce1_mag)
                                # print("normal angle = " + str(normal_angle))

                                # calculating the angle (theta) between x-ring centre, and ring normal
                                trp5_theta_dot_prod = get_dot_product(trp5_normal_vect, trp5_x_dist_vect)
                                trp5_theta_angle = get_angle(trp5_theta_dot_prod, trp5_x_dist_mag, trp5_normal_mag)

                                if trp5_theta_angle > 90:  # if the normal is facing away from the X atom.. (if the atom is below the plane)
                                    trp5_theta_angle = 180 - trp5_theta_angle

                                if trp5_theta_angle < xtheta_cutoff:  # carry on doing lots of fun maths, if we still are in xhpi
                                    # geometry

                                    trp5 = "T"  # trp5 ring is within the defined geometry of the X atom

                                    # calculate the distance between the ring plane and X atom
                                    # the angle between the aromatic plane and X-ring center vector
                                    trp5_plane_x_angle = 90 - trp5_theta_angle

                                    # x_plane_mag is the shortest distance between the plane and the X atom
                                    trp5_x_plane_mag = trp5_x_dist_mag * math.sin((trp5_plane_x_angle * (math.pi / 180)))

                                    # Defining the posiiton on the aromatic plane that is directly below the X atom
                                    # This is done by calculating the vector between the plane and the X atom, and then
                                    # moving the X atom along this vector, by x_plane_mag (distance between the X atom and
                                    # the plane)

                                    # making a vector parallel to the normal, with a length of x_plane_mag (distance between
                                    # plane and x atom)
                                    # Can then move the X atom coords by this to get the position on the plane as described above
                                    trp5_scaled_normal = get_normalised_vector(trp5_normal_vect, trp5_normal_mag,
                                                                            trp5_x_plane_mag)

                                    # We dont know if the X atom is above or below the plane at this point
                                    # Thus we dont know whether to add or subtract the scaled normal from the X atom
                                    # To get around this, both are done, and then the angle between the normal vector and the
                                    # vector connecting the ring centre the the planar X position is calculated for both cases
                                    # this should be 90deg (dot product = 0) for the case where the coordinates hace been moved
                                    # onto the plane.

                                    trp5_xplane_coords_sum = [(x_atom_coords[0] + trp5_scaled_normal[0]),
                                                            (x_atom_coords[1] + trp5_scaled_normal[1]),
                                                            (x_atom_coords[2] + trp5_scaled_normal[2])]

                                    trp5_rc_xplane_vect_sum = get_vector(trp5_centre, trp5_xplane_coords_sum)

                                    trp5_xplane_coords_sub = [(x_atom_coords[0] - trp5_scaled_normal[0]),
                                                            (x_atom_coords[1] - trp5_scaled_normal[1]),
                                                            (x_atom_coords[2] - trp5_scaled_normal[2])]

                                    trp5_rc_xplane_vect_sub = get_vector(trp5_centre, trp5_xplane_coords_sub)

                                    # can find out which one is on the aromatic plane now
                                    # the dot product between plane normal and the vector from the ring centre to the x plane
                                    # coords should be 0
                                    # rounding errors mean its just a very small number close to 0, so we just take the smallest

                                    trp5_dot_prod_sum = get_dot_product(trp5_rc_xplane_vect_sum, trp5_normal_vect)
                                    trp5_dot_prod_sub = get_dot_product(trp5_rc_xplane_vect_sub, trp5_normal_vect)

                                    # print(str(dot_prod_sum) + " " + str(dot_prod_sub))

                                    if (trp5_dot_prod_sum ** 2) < (trp5_dot_prod_sub ** 2):
                                        trp5_xplane_coords = trp5_xplane_coords_sum
                                        trp5_rc_xplane_vect = trp5_rc_xplane_vect_sum
                                        trp5_x_pos = "below"

                                    else:
                                        trp5_xplane_coords = trp5_xplane_coords_sub
                                        trp5_rc_xplane_vect = trp5_rc_xplane_vect_sub
                                        trp5_x_pos = "above"

                                    trp5_rc_xplane_mag = get_magnitude(trp5_centre, trp5_xplane_coords)

                                    # calculating the planar angle for his
                                    # this is a bit different as it is not 6 memebered
                                    # thus i calculate the angle relative to the vector passing through CG, and the midpoint
                                    # between NE2, and ce1
                                    # this is where cz is closest to if it were a 6 membered ring

                                    # CE1 NE2 midpoint
                                    ne1_ce2_mid = get_midpoint(trp_ne1, trp_ce2)

                                    cg_mid_vect = get_vector(trp_cg, ne1_ce2_mid)
                                    cg_mid_mag = get_magnitude(trp_cg, ne1_ce2_mid)

                                    trp5_planar_dot_prod = get_dot_product(trp5_rc_xplane_vect, cg_mid_vect)
                                    trp5_planar_angle = get_angle(trp5_planar_dot_prod, trp5_rc_xplane_mag, cg_mid_mag)

                                    # the above provides the shortest angle, and looses information about whether it should be
                                    # clockwise or anticlockwise
                                    # here the clockwise side is defined as the side with ne1 and cd1 atoms
                                    # anticlockwise is the side with cd2, and ce2

                                    # define a his_rc-ne1 vector and calculate the angle between this and the plane position
                                    # of X

                                    rc_ne1_vect = get_vector(trp5_centre, trp_ne1)
                                    rc_ne1_mag = get_magnitude(trp5_centre, trp_ne1)

                                    rc_ce2_vect = get_vector(trp5_centre, trp_ce2)
                                    rc_ce2_mag = get_magnitude(trp5_centre, trp_ce2)

                                    rcne1_rcxplane_dotprod = get_dot_product(trp5_rc_xplane_vect, rc_ne1_vect)
                                    rcce2_rcxplane_dotprod = get_dot_product(trp5_rc_xplane_vect, rc_ce2_vect)

                                    # print("ce2ce1_rcxplane")
                                    rcne1_rcxplane_angle = get_angle(rcne1_rcxplane_dotprod, trp5_rc_xplane_mag, rc_ne1_mag)
                                    rcce2_rcxplane_angle = get_angle(rcce2_rcxplane_dotprod, trp5_rc_xplane_mag, rc_ce2_mag)

                                    # if this is then larger than phi, it is anticlockwise
                                    # clockwise = 360 - anticlockwise
                                    if rcne1_rcxplane_angle > rcce2_rcxplane_angle:
                                        trp5_planar_angle = 360 - trp5_planar_angle

                                    write_pymol(pymol, pi_res_num, x_res_num)


                        if trp6_x_dist_mag < dist_cutoff and trp6_b_factor > plddt_cutoff:  # only do this lot of maths if distance is within cutoff
                            trp6_x_dist_vect = get_vector(trp6_centre, x_atom_coords)  # need this for theta

                            # checking that the H atom is close to the ring centre
                            h_atom_ids = h_dict[x_res_id][x_atom_id][0]

                            h_coords = []
                            for h in h_atom_ids:
                                h_coords.append(get_coords(x_res_num, x_chain, h, pdbdict))

                            h_geom = []
                            for h_coord in h_coords:
                                h_geom.append(check_h_geom(h_coord, x_atom_coords, trp6_centre, h_angle_cutoff))

                            # print(str(x_res_num) + "    " + str(x_atom_id))
                            # print(h_geom)

                            # continue if the h atom is withing the geometric cutoffs
                            if "T" in h_geom:

                                # generating 2 vectors on the aromatic plane (required for calculating cross product (ring
                                # normal))
                                cd2_ch2_vect = get_vector(trp_cd2, trp_ch2)
                                ce2_cz3_vect = get_vector(trp_ce2, trp_cz3)

                                # getting the magnitudes of these vectors (required for calculating cross product magnitude)
                                cd2_ch2_mag = get_magnitude(trp_cd2, trp_ch2)
                                ce2_cdz3_mag = get_magnitude(trp_ce2, trp_cz3)

                                # calculating the angle between these vectors (required for calculating cross product magnitude)
                                # cgne1_cd1cd2_dot_prod = get_dot_product(cg_ne1_vect, cd1_cd2_vect)
                                # cgcz_cd1ce2_angle = get_angle(cgce1_nd1cd2_dot_prod, cg_ce1_mag, nd1_cd2_mag)

                                # finally calculating cross product goodness (defining ring normal vector and magnitude)
                                trp6_normal_vect = get_cross_product_vector(cd2_ch2_vect, ce2_cz3_vect)
                                trp6_normal_mag = get_vector_magnitude(trp6_normal_vect)

                                # calculating and printing the angle between a line on the aromatic plane and the normal
                                # making sure its actually 90deg
                                # normal_dot_prod = get_dot_product(his_normal_vect, cg_ce1_vect)
                                # normal_angle = get_angle(normal_dot_prod, his_normal_mag, cg_ce1_mag)
                                # print("normal angle = " + str(normal_angle))

                                # calculating the angle (theta) between x-ring centre, and ring normal
                                trp6_theta_dot_prod = get_dot_product(trp6_normal_vect, trp6_x_dist_vect)
                                trp6_theta_angle = get_angle(trp6_theta_dot_prod, trp6_x_dist_mag, trp6_normal_mag)

                                # if the normal is facing away from the X atom.. (if the atom is below the plane)
                                if trp6_theta_angle > 90:
                                    trp6_theta_angle = 180 - trp6_theta_angle

                                if trp6_theta_angle < xtheta_cutoff:  # carry on doing lots of fun maths, if we still are in xhpi
                                    # geometry

                                    trp6 = "T"  # trp6 ring is within the defined geometry of the X atom

                                    # calculate the distance between the ring plane and X atom
                                    # the angle between the aromatic plane and X-ring center vector
                                    trp6_plane_x_angle = 90 - trp6_theta_angle

                                    # x_plane_mag is the shortest distance between the plane and the X atom
                                    trp6_x_plane_mag = trp6_x_dist_mag * math.sin((trp6_plane_x_angle * (math.pi / 180)))

                                    # Defining the posiiton on the aromatic plane that is directly below the X atom
                                    # This is done by calculating the vector between the plane and the X atom, and then
                                    # moving the X atom along this vector, by x_plane_mag (distance between the X atom and
                                    # the plane)

                                    # making a vector parallel to the normal, with a length of x_plane_mag (distance between
                                    # plane and x atom)
                                    # Can then move the X atom coords by this to get the position on the plane as described above
                                    trp6_scaled_normal = get_normalised_vector(trp6_normal_vect, trp6_normal_mag,
                                                                            trp6_x_plane_mag)

                                    # We dont know if the X atom is above or below the plane at this point
                                    # Thus we dont know whether to add or subtract the scaled normal from the X atom
                                    # To get around this, both are done, and then the angle between the normal vector and the
                                    # vector connecting the ring centre the the planar X position is calculated for both cases
                                    # this should be 90deg (dot product = 0) for the case where the coordinates hace been moved
                                    # onto the plane.

                                    trp6_xplane_coords_sum = [(x_atom_coords[0] + trp6_scaled_normal[0]),
                                                            (x_atom_coords[1] + trp6_scaled_normal[1]),
                                                            (x_atom_coords[2] + trp6_scaled_normal[2])]

                                    trp6_rc_xplane_vect_sum = get_vector(trp6_centre, trp6_xplane_coords_sum)

                                    trp6_xplane_coords_sub = [(x_atom_coords[0] - trp6_scaled_normal[0]),
                                                            (x_atom_coords[1] - trp6_scaled_normal[1]),
                                                            (x_atom_coords[2] - trp6_scaled_normal[2])]

                                    trp6_rc_xplane_vect_sub = get_vector(trp6_centre, trp6_xplane_coords_sub)

                                    # can find out which one is on the aromatic plane now
                                    # the dot product between plane normal and the vector from the ring centre to the x plane
                                    # coords should be 0
                                    # rounding errors mean its just a very small number close to 0, so we just take the smallest

                                    trp6_dot_prod_sum = get_dot_product(trp6_rc_xplane_vect_sum, trp6_normal_vect)
                                    trp6_dot_prod_sub = get_dot_product(trp6_rc_xplane_vect_sub, trp6_normal_vect)

                                    # print(str(dot_prod_sum) + " " + str(dot_prod_sub))

                                    if (trp6_dot_prod_sum ** 2) < (trp6_dot_prod_sub ** 2):
                                        trp6_xplane_coords = trp6_xplane_coords_sum
                                        trp6_rc_xplane_vect = trp6_rc_xplane_vect_sum
                                        trp6_x_pos = "below"

                                    else:
                                        trp6_xplane_coords = trp6_xplane_coords_sub
                                        trp6_rc_xplane_vect = trp6_rc_xplane_vect_sub
                                        trp6_x_pos = "above"

                                    trp6_rc_xplane_mag = get_magnitude(trp6_centre, trp6_xplane_coords)

                                    # calculating the planar angle for trp6
                                    # this is relative to the vector passing through cd2 and ch2

                                    trp6_planar_dot_prod = get_dot_product(trp6_rc_xplane_vect, cd2_ch2_vect)
                                    trp6_planar_angle = get_angle(trp6_planar_dot_prod, trp6_rc_xplane_mag, cd2_ch2_mag)

                                    # the above provides the shortest angle, and looses information about whether it should be
                                    # clockwise or anticlockwise
                                    # here the clockwise side is defined as the side with ce2 and cz2 atoms
                                    # anticlockwise is the side with ce3, and cz3

                                    # define cd1 ce1 midpoint and cd2 ce2 midpoint

                                    mid2 = get_midpoint(trp_ce3, trp_cz3)
                                    mid3 = get_midpoint(trp_ce2, trp_cz2)

                                    mid2_mid3_vect = get_vector(mid2, mid3)
                                    mid2_mid3_mag = get_magnitude(mid2, mid3)

                                    mid3_mid2_vect = get_vector(mid3, mid2)
                                    mid3_mid2_mag = get_magnitude(mid3, mid2)

                                    mid32_rcxplane_dotprod = get_dot_product(trp6_rc_xplane_vect, mid3_mid2_vect)

                                    mid23_rcxplane_dotprod = get_dot_product(trp6_rc_xplane_vect, mid2_mid3_vect)

                                    mid32_rcxplane_angle = get_angle(mid32_rcxplane_dotprod, trp6_rc_xplane_mag,
                                                                    mid3_mid2_mag)
                                    mid23_rcxplane_angle = get_angle(mid23_rcxplane_dotprod, trp6_rc_xplane_mag,
                                                                    mid2_mid3_mag)

                                    # if this is then larger than phi, it is anticlockwise
                                    # clockwise = 360 - anticlockwise
                                    if mid23_rcxplane_angle < mid32_rcxplane_angle:
                                        trp6_planar_angle = 360 - trp6_planar_angle

                                    write_pymol(pymol, pi_res_num, x_res_num)

                            if trp5 == "T" and trp6 == "F":
                                pi_res_id = "TRP5"
                                pi_bfactor = trp5_b_factor
                                x_dist_mag = trp5_x_dist_mag
                                theta_angle = trp5_theta_angle
                                planar_angle = trp5_planar_angle
                                x_plane_mag = trp5_x_plane_mag
                                rc_xplane_mag = trp5_rc_xplane_mag
                                x_pos = trp5_x_pos

                                write_data(report, x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, 
                                               x_chain, os.path.basename(pdb_file).replace("_H.pdb", ""), resolution, pi_res_id, pi_res_num, pi_bfactor, pi_chain, 
                                               x_dist_mag, theta_angle, planar_angle, x_plane_mag, rc_xplane_mag, x_pos, x_group)


                            if trp6 == "T" and trp5 == "F":
                                pi_res_id = "TRP6"
                                pi_bfactor = trp6_b_factor
                                x_dist_mag = trp6_x_dist_mag
                                theta_angle = trp6_theta_angle
                                planar_angle = trp6_planar_angle
                                x_plane_mag = trp6_x_plane_mag
                                rc_xplane_mag = trp6_rc_xplane_mag
                                x_pos = trp6_x_pos

                                write_data(report, x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, 
                                               x_chain, os.path.basename(pdb_file).replace("_H.pdb", ""), resolution, pi_res_id, pi_res_num, pi_bfactor, pi_chain, 
                                               x_dist_mag, theta_angle, planar_angle, x_plane_mag, rc_xplane_mag, x_pos, x_group)

                            if trp5 == "T" and trp6 == "T":

                                if unique_trp == "T":

                                    if trp5_x_dist_mag < trp6_x_dist_mag:
                                        pi_res_id = "TRP5"
                                        pi_bfactor = trp5_b_factor
                                        x_dist_mag = trp5_x_dist_mag
                                        theta_angle = trp5_theta_angle
                                        planar_angle = trp5_planar_angle
                                        x_plane_mag = trp5_x_plane_mag
                                        rc_xplane_mag = trp5_rc_xplane_mag
                                        x_pos = trp5_x_pos

                                        write_data(report, x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, 
                                               x_chain, os.path.basename(pdb_file).replace("_H.pdb", ""), resolution, pi_res_id, pi_res_num, pi_bfactor, pi_chain, 
                                               x_dist_mag, theta_angle, planar_angle, x_plane_mag, rc_xplane_mag, x_pos, x_group)

                                    if trp6_x_dist_mag < trp5_x_dist_mag:
                                        pi_res_id = "TRP6"
                                        pi_bfactor = trp6_b_factor
                                        x_dist_mag = trp6_x_dist_mag
                                        theta_angle = trp6_theta_angle
                                        planar_angle = trp6_planar_angle
                                        x_plane_mag = trp6_x_plane_mag
                                        rc_xplane_mag = trp6_rc_xplane_mag
                                        x_pos = trp6_x_pos

                                        write_data(report, x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, 
                                               x_chain, os.path.basename(pdb_file).replace("_H.pdb", ""), resolution, pi_res_id, pi_res_num, pi_bfactor, pi_chain, 
                                               x_dist_mag, theta_angle, planar_angle, x_plane_mag, rc_xplane_mag, x_pos, x_group)

                                if unique_trp == "F":
                                    pi_res_id = "TRP5"
                                    pi_bfactor = trp5_b_factor
                                    x_dist_mag = trp5_x_dist_mag
                                    theta_angle = trp5_theta_angle
                                    planar_angle = trp5_planar_angle
                                    x_plane_mag = trp5_x_plane_mag
                                    rc_xplane_mag = trp5_rc_xplane_mag
                                    x_pos = trp5_x_pos

                                    write_data(report, x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, 
                                               x_chain, os.path.basename(pdb_file).replace("_H.pdb", ""), resolution, pi_res_id, pi_res_num, pi_bfactor, pi_chain, 
                                               x_dist_mag, theta_angle, planar_angle, x_plane_mag, rc_xplane_mag, x_pos, x_group)

                                    pi_res_id = "TRP6"
                                    pi_bfactor = trp6_b_factor
                                    x_dist_mag = trp6_x_dist_mag
                                    theta_angle = trp6_theta_angle
                                    planar_angle = trp6_planar_angle
                                    x_plane_mag = trp6_x_plane_mag
                                    rc_xplane_mag = trp6_rc_xplane_mag
                                    x_pos = trp6_x_pos

                                    write_data(report, x_res_id, x_res_num, x_atom_id, x_atom_type, x_atom_num, x_bfactor, 
                                               x_chain, os.path.basename(pdb_file).replace("_H.pdb", ""), resolution, pi_res_id, pi_res_num, pi_bfactor, pi_chain, 
                                               x_dist_mag, theta_angle, planar_angle, x_plane_mag, rc_xplane_mag, x_pos, x_group)


if __name__ == "__main__":
    
    # check input args
    args = parser.parse_args()
    
    # create tmp and output directories
    os.makedirs('tmp', exist_ok=True)
    os.makedirs(args.output, exist_ok=True)

    # check presence of hgen command
    HGEN_CALL = "hgen"
    try:
        subprocess.call(HGEN_CALL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        try:
            with open("tmp/run_hgen.sh", "w") as fh:
                fh.write("#! /bin/bash\n")
                fh.write("source /ccp4/bin/ccp4.setup-sh\n")
                fh.write('hgen "@"')
            HGEN_CALL = "tmp/run_hgen.sh"
            os.system(f"chmod +x {HGEN_CALL}")
            status = subprocess.call(HGEN_CALL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if status != 1:
                raise IOError("run_hgen.sh command failed")
        except IOError:
           sys.exit(f'hgen is not in PATH or ccp4 is not installed. Terminating...')

    # check input PDB/file or folder
    if all([args.PDB_file is None, args.PDB_folder is None]) == True:   
        print('No PDB file nor PDB-containing directory provided')
        sys.exit(f"{parser.print_help()}")

    if args.PDB_file is not None:
        input_file = args.PDB_file
        main(input_file, args)
    
    elif  args.PDB_folder is not None:
        for file in filter(lambda x: x.endswith(".pdb"), os.listdir(args.PDB_folder)):
            input_file = os.path.join(args.PDB_folder, file)
            main(input_file, args)

    # create summary
    if args.PDB_folder is not None and args.summary is True:
        summary = ""
        for file in filter(lambda x: x.endswith('.csv'), os.listdir(args.output)):
            file_in = open(os.path.join(args.output, file), 'r').readlines()
            header = file_in[0]
            for line in file_in[1:]:
                summary += line

        summary = header + summary
        with open(os.path.join(args.output, 'summary.csv'), 'w') as fh:
            fh.write(summary)
    
    # delete tmp directory
    shutil.rmtree('tmp')

