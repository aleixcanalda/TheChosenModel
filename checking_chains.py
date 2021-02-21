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
						if instance.compare_sequence(new_instance.get_sequence_chain()) == True:
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
