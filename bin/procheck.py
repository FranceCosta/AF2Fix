#! /usr/bin/python3

"""
    Script to run PROHECK using  CCP4 suite to get ramachandran outliers
    Francesco Costa fcosta@ebi.ac.uk
    28/09/2023
"""

import argparse
import subprocess
import os, sys, shutil, re

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
    help="Path/to/output folder. Created if non existing",
    default='./', 
    type=str
)

parser.add_argument(
    "--summary", 
    required=False, 
    help="If --summary and --PDB_folder PDB_FOLDER, creates a summary.csv with all .csv found in the output folder",
    action='store_true'
)


def main(infile: str, args):

    
    outdir = os.path.abspath(args.output)
    infile = os.path.abspath(infile)
    prot_name = os.path.basename(infile).replace('.pdb', '')

    # run procheck
    os.chdir('tmp')
    if PROCHECK_CALL == 'procheck':
        subprocess.run(f"{PROCHECK_CALL} {infile} A 1",
                     shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elif PROCHECK_CALL == "./run_procheck.sh":
        with open(PROCHECK_CALL, "w") as fh:
                fh.write("#! /bin/bash\n")
                fh.write("source /ccp4/bin/ccp4.setup-sh\n")
                fh.write(f"procheck {infile} A 1")
        subprocess.run(PROCHECK_CALL,
                     shell=True)#, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    results = open(f"{prot_name}.sum", 'r').readlines()
    os.chdir('../')

    # get data
    rama_ = re.compile("[0-9]{1,3}.[0-9]{1}")

    for line in results:
        if "Ramachandran plot" in line:
            core, allow, gener, disall = list(rama_.findall(line))
    print(f'{prot_name}: {core}, {allow}, {gener}, {disall}')
    with open(os.path.join(outdir, f'{prot_name}_rama.csv'), 'w') as fh:
        fh.write(f'prot_name,core,allow,gener,disall\n{prot_name},{core},{allow},{gener},{disall}')

if __name__ == "__main__":
    
    # check input args
    args = parser.parse_args()
    
    # create tmp and output directories
    os.makedirs('tmp', exist_ok=True)
    os.makedirs(args.output, exist_ok=True)
    current_dir = os.getcwd()

    # check presence of procheck command
    os.chdir('tmp')
    PROCHECK_CALL = "procheck"
    try:
        subprocess.call(PROCHECK_CALL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        try:
            PROCHECK_CALL = "./run_procheck.sh"
            with open(PROCHECK_CALL, "w") as fh:
                fh.write("#! /bin/bash\n")
                fh.write("source /ccp4/bin/ccp4.setup-sh\n")
                fh.write("procheck --help\n")

            os.system(f"chmod +x {PROCHECK_CALL}")
            print(os.listdir())
            output = subprocess.check_output(PROCHECK_CALL, stderr=subprocess.STDOUT)
            if 'procheck  /data/pdb/p1amt.pdb  A  1.5' in str(output) == False:
                raise IOError("run_procheck.sh command failed")
        except IOError:
           sys.exit(f'procheck is not in PATH or ccp4 is not installed. Terminating...')
    os.chdir(current_dir)

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