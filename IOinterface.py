import sys
import gzip
import re
from Bio import SeqIO
import os, glob
import argparse

parser = argparse.ArgumentParser(description="This program analyzes macromolecular structures")

parser.add_argument('-i', '--input-directory',
                    dest = "input_directory",
                    action = "store",
                    required = True,
                    help = "Directory containing the input PDB files.")

parser.add_argument('-s', '--stechiometry',
                    dest = "stechiometry",
                    action = "store",
                    default = None,
                    help = "Path to a file containing the stechiometry of the complex.")

parser.add_argument('-o', '--output-directory',
                    dest = "output_directory",
                    action = "store",
                    required = True,
                    help = "Directory where the output will be saved.")

parser.add_argument('-f', '--force',
                    dest = "force",
                    action = "store_true",
                    default = False,
                    help = "Overwrite the content of the output directory.")

parser.add_argument('-v', '--verbose',
                    dest = "verbose",
                    action = "store_true",
                    default = False,
                    help = "Print the progession of the execution.")

parser.add_argument('-e', '--energy',
                    dest = "energy",
                    action = "store_true",
                    default = False,
                    help = "Calculate DOPE energy and plot the result.")


options = parser.parse_args()


def get_input_file(input_path):
    """ Handling with different kind of input: only fasta or gunzip fasta files 
    for a given path or for the current directory """
    fasp = re.compile('.pdb$|.pdb.gz$')
    path = input_path
    final_prots_files = []
    try:
        if os.path.isdir(path):
            prots_files = [f for f in os.listdir(path) if fasp.search(f) is not None]
            os.chdir(path)
    except:
        raise OSError('No directory could be read')

    for file in prots_files:
        if file.endswith('.gz'):
            with gzip.open(file) as fd:
                for line in fd:
                    chains = []
                    if line.startswith("ATOM"):
                        pdb = line.split()
                        chain = pdb[4]
                        if chain not in chains:
                            chains.append(chain)
        else:
            with open(file) as fd:
                chains=[]
                for line in fd:
                    if line.startswith("ATOM"):
                        pdb = line.split()
                        chain = pdb[4]
                        if chain not in chains:
                            chains.append(chain)

        if re.search(r'\w+_%s_%s'%('|'.join(chains),'|'.join(chains)),file) is not None:
            final_prots_files.append(file)
        else:
            raise ValueError('File name %s is not correct'%(file))

    return final_prots_files

def get_output_file(output_path, force):
    """ Handling with different kind of input: only fasta or gunzip fasta files 
    for a given path or for the current directory """
    
    if os.path.is_dir(output_path) == False:
        os.mkdir(output_path)
        os.mkdir(output_path + "/structures")
        os.mkdir(output_path + "/analysis")
    else:
        if force == True:
            out1 = output_path + "/structures"
            out2 = output_path + "/analysis"
            if os.path.is_dir() == False:
                os.mkdir(output_path + "/structures")
            if os.path.is_dir() == False:
                os.mkdir(output_path + "/analysis")
            
        else:
            print("The output directory already exists.")
            sys.exit()

if __name__ == "__main__":

    prots_files = get_input_file(options.infile)
    get_output_file(options.outfile,options.force)

    path_stech = options.stechiometry
    if path_stech is not None:
        with open(path_stech) as fd:
            stech = [f for f in fd if re.search(r'\w:\w',f)]

    verbose = options.verbose
