import sys
from classes import *
#from IOinterface.py import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.PDBIO import PDBIO
import random

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

def start_model(interactions_dict,unique_chains_list,pdbfiles,verbose=False):
	"""Choosing the chain with most interactions as the starting model of the macrocomplex"""

	model = Model("A") #Create new instance of class Model
	
	if verbose:
		sys.stderr.write("Deciding the starting model.")

	most_inter = 0

	for chain in interactions_dict.keys():

		chain_inter = len(interactions_dict[chain])

		if chain_inter > most_inter:

			most_inter = chain_inter

			starting_chain = chain

	model.add(next(x for x in unique_chains_list if x.id == starting_chain))

	second_chain = False

	for inter in interactions_dict[starting_chain]:
		if second_chain == True:
			break
		for file in pdbfiles:
			if inter in str(file) and starting_chain in str(file):
				model.add(next(x for x in unique_chains_list if x.id == inter))
				second_chain = True
				break

	return model

def equal_length_chains(chain1, chain2):
	"""Equaling the lengths of two chains."""
	
	atoms1 = sorted(chain1.get_atoms())
	atoms2 = sorted(chain2.get_atoms())

	if len(atoms1) > len(atoms2):
		return atoms1[:len(atoms2)], atoms2

	elif len(atoms1) < len(atoms2):
		return atoms1, atoms2[:len(atoms1)]
	
	else:
		return atoms1, atoms2


def clashes(chain_atoms, model):
	""" """
	
	print("Hi")
	backbone = {"CA", "C1\'"}
	chain_atoms = [a for a in chain_atoms if a.id in backbone]  # Gets only the backbone atoms
	model_atoms = [a for a in model.get_atoms() if a.id in backbone]
	Nsearch = Bio.PDB.NeighborSearch(model_atoms)
	
	clash_count = 0
	for atoms in chain_atoms:
		clash = Nsearch.search(atoms.coord, 1)
		if clash != []:
			clash_count += 1
	print(clash_count/len(chain_atoms))
	if clash_count/len(chain_atoms) >= 0.05:
		return True
	else:
		return False

	
def superimpose(unique_chains_list,pdbfiles, verbose=False):
	""" """
	print(unique_chain_list)
	int_dict = get_interactions_dict(unique_chains_list, verbose)
	model = start_model(int_dict, unique_chains_list,pdbfiles,verbose)
	chain1, chain2 = model.get_chains()
	chain_in_model = [chain1.id,chain2.id]
	print(int_dict)
	print(chain1)
	print(chain2)
	#random.shuffle(unique_chains_list)
	n = 2
	while n<len(int_dict.keys()):
		for chain in unique_chains_list:
			if chain.id in chain_in_model:
				continue
			print("looking at:%s"%(chain))	
			not_added = True
			#shuffle(chain_in_model): if we want to make more models, we should shuffle the chainsin the model so the chains are superimposed to different chains.
			for chainin in chain_in_model:
				
				if chain.id in int_dict[chainin]: #if our chain interacts with a chain inside the complex
					newChain, modelChain = equal_length_chains(chain, next(x for x in unique_chains_list if x.id == chainin))
					superimpose = Bio.PDB.Superimposer()
					superimpose.set_atoms(modelChain, newChain)
					chain_copy = chain.copy()
					chain_copy.id = chain.id
					superimpose.apply(chain_copy)
					chain_atoms = sorted(chain_copy.get_atoms())
					
					if clashes(chain_atoms, model) == True:
						continue
						
					else:
						chain_copy.parent = None
						print(chain_copy)
						model.add(chain_copy)
						not_added = False
						chain_in_model.append(chain.id)
						print(chain_in_model)
						n += 1

						if verbose:
							print("%s sucessfully added to the model" %chain.id)

						break
						
			if verbose and not_added:
				print("%s could not be added to the model" %chain.id)
			
			#n += 1
			#chain_in_model.append(chain.id)
	return model
			
def save_PDB(model, output_path, verbose=False):
	if verbose:
		print("Saving model")
	io = PDBIO()
	io.set_structure(model)
	io.save("model.pdb")
	
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

	 #int_dic = get_interactions_dict(unique_chain_list)

	#start_chain = start_model(int_dic)
	model = superimpose(unique_chain_list,files)
	save_PDB(model, "home/maria/master/Second_term/PytGA_project")
