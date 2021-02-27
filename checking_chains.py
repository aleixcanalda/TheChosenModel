#import IOinterface.py
import sys
from classes import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Chain import Chain
parser = PDBParser(PERMISSIVE=1, QUIET=True)

def all_chains(PDB_files, verbose=False):
	""" Creating an object for every unique PDB chain. """

	unique_chain_list=[]
	id_set={"A"}
	if verbose:
		print("Processing PDB files and analyzing the chains")

	for file in PDB_files:
		for model in parser.get_structure("A",file): #Obtain the structures and iterate on all models
			chain_num = 0

			for chain in model: #Obtain all chains from the models
				
				new_instance = MyChain(chain)
				chain_num += 1

				if unique_chain_list == []:
					unique_chain_list.append(new_instance) #First instance is appended to the list.

				else:
					unique = True
					for instance in unique_chain_list:
						
						if instance.compare_sequence(new_instance) == True:
							unique = False
							break
							
					if unique == True:	#If the sequence is different from all the others of the list, it's appended.
						letter = chain_id(id_set)
						new_instance.id = letter
						id_set.add(letter)
						unique_chain_list.append(new_instance)

			if chain_num != 2:
				sys.stderr.write("All PDB files must contain two chains.")


	return unique_chain_list
							


def chain_id(id_set):
	""" Function to create unique IDs """
	
	alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
				'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'AA', 'BB',
				'CC', 'DD', 'EE', 'FF', 'GG', 'HH', 'II', 'JJ', 'KK', 'LL', 'MM', 'NN',
				'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'UU', 'VV', 'WW', 'XX', 'YY', 'ZZ']
	for letter in alphabet:
		if letter not in id_set:
			return letter

def get_interactions_dict(unique_chain_list, verbose=False):
	""" 
	Obtain the dictionary of interactions among all chains

	The dictionary will look as follows:

	[Chain1.id] = [ChainX.id that interacts with Chain1] 

	"""
	if verbose:
		sys.stderr.write("Obtaining the interactions dictionary.")

	interactions_dict = {}
	for chain1 in unique_chain_list:

		interactions_dict[chain1.id] = []
	
		for chain2 in unique_chain_list:

			if chain1.id == chain2.id:
				continue

			interact_1, interact_2 = chain1.interactions(chain2)

			if interact_1 == set() and interact_2 == set():
				continue

			else:

				interactions_dict[chain1.id].append(chain2.id)

	return interactions_dict

def start_model(interactions_dict,verbose=False):
	"""Choosing the chain with most interactions as the starting model of the macrocomplex"""
	if verbose:
		sys.stderr.write("Deciding the starting model.")

	most_inter = 0

	for chain in interactions_dict.keys():

		chain_inter = len(interactions_dict[chain])

		if chain_inter > most_inter:

			most_inter = chain_inter

			starting_chain = chain

	return starting_chain

if __name__ == "__main__":
	pdb = open("1gzx_A_B.pdb")

	files=["1gzx_A_B.pdb","1gzx_A_C.pdb","1gzx_A_D.pdb"]

	parser = PDBParser(PERMISSIVE=1, QUIET=True)

	chains = []

	struct = parser.get_structure(file="1gzx_A_B.pdb",id="A")
	for model in struct:
		for chain in model:
			ch = MyChain(chain)
			chains.append(ch)

	unique_chain_list = all_chains(files)

	int_dic = get_interactions_dict(unique_chain_list)

	print(int_dic)

	start_chain = start_model(int_dic)

	print(start_chain)