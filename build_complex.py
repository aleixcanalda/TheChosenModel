import sys
from classes import *
#from IOinterface.py import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model as modelmodel
from Bio.PDB.PDBIO import PDBIO
import random
import gzip
import re
from Bio import SeqIO
import os, glob
import argparse
#from modeller.model import Model as modelmodel


parser = PDBParser(PERMISSIVE=1, QUIET=True)

def all_chains(PDB_files, verbose=False):
	""" Creating an object for every unique PDB chain. """

	unique_chain_list=[]
	id_set=set()
	if verbose:
		print("Processing PDB files and analyzing the chains")

	for file in PDB_files:
		for model in parser.get_structure("A",file): #Obtain the structures and iterate on all models
			chain_num = 0
			my_list = []
			for chain in model: #Obtain all chains from the models
				new_instance = MyChain(chain)
				chain_num += 1
				my_list.append(new_instance)							

			if chain_num != 2:
				sys.stderr.write("All PDB files must contain two chains.")
				
			unique_chain_list.append(my_list)


	return unique_chain_list
							

def get_interactions_dict(unique_chain_list, verbose=False):
	""" 
	Obtain the dictionary of interactions among all chains
	The dictionary will look as follows:
	'Chain1.id' = [[Chain1, Chain2 that interacts with chain 1],[Chain1,Chain3],...] 
	"""
	if verbose:
		print("Obtaining the interactions dictionary.")

	interactions_dict = {}
	
	for pair in unique_chain_list:
		chain1 = pair[0]
		chain2 = pair[1]
		
		interact_1, interact_2 = chain1.interactions(chain2)
		if interact_1 == set() and interact_2 == set():
			continue

		if interactions_dict == {}:
			interactions_dict[chain1.id]= [[chain1, chain2]]
			interactions_dict[chain2.id]= [[chain2, chain1]]
		else:
			checking1 = False
			checking2 = False
			for key in interactions_dict.keys():
				check1 = chain1.compare_sequence(interactions_dict[key][0][0]) ##Dilemes de l'Aleix 2.0 (igual coordenada, igual seqüència)
				if check1 == 2:
					interactions_dict[key].append([chain1, chain2])
					checking1 = True
				check2 = chain2.compare_sequence(interactions_dict[key][0][0])
				if check2 == 2:
					interactions_dict[key].append([chain2, chain1])
					checking2 = True
			if checking1 == False: ##WHAT WHAT Dilemes de la MAria 1.0 (per que check1 ja no es true?)
				interactions_dict[chain1.id]= [[chain1, chain2]]
			if checking2 == False:
				interactions_dict[chain2.id]= [[chain2, chain1]]
	
	return interactions_dict

def get_stech_dicts(unique_chain_list, stechiometry, verbose=False):
	""" Obtain the stechiometry from the input file and from the interacting chains"""
	if verbose:
		print("Obtaining the stechiometry dictionary.")

	stech_dict = {}
	stech_file = {}
	chains_list = []
	chain_ids = []

	for chain1, chain2 in unique_chain_list:
		if chain1.compare_sequence(chain2) != 2:

			if chain1.id not in chain_ids:
				chains_list.append(chain1)
				chain_ids.append(chain1.id)
			if chain2.id not in chain_ids:
				chains_list.append(chain2)
				chain_ids.append(chain2.id)


	for chain1 in chains_list:
		for chain2 in chains_list:
			if chain1.compare_sequence(chain2) == 1:
				if chain1.id not in stech_dict.keys():
					stech_dict[chain1.id] = [chain1.id,chain2.id]
				else:
					stech_dict[chain1.id].append(chain2.id)
				

	fd = open(stechiometry, "r")
	for line in fd:
		line = line.strip()
		stech_file[line[0]] = line[2:]
	fd.close()
	
	return stech_dict, stech_file

def start_model(interactions_dict,stechiometry=None,verbose=False):
	"""Choosing the chain with most interactions as the starting model of the macrocomplex"""

	model = modelmodel("A") #Create new instance of class Model
	
	if verbose:
		print("Deciding the starting model.")

	most_inter = 0

	for chain in interactions_dict.keys():

		chain_inter = len(interactions_dict[chain])

		if chain_inter > most_inter:

			most_inter = chain_inter

			starting_chain = chain
			
	"""
	if stechiometry != None:
		if stechiometry[starting_chain][1] == interactions_dict[starting_chain][0][1].id:
			model.add(interactions_dict[starting_chain][0][0]) #We add the chain with most interactions to the model.

			model.add(interactions_dict[starting_chain][1][1])
		else:
			model.add(interactions_dict[starting_chain][0][0]) #We add the chain with most interactions to the model.

			model.add(interactions_dict[starting_chain][0][1])
	else:
		model.add(interactions_dict[starting_chain][0][0]) #We add the chain with most interactions to the model.

		model.add(interactions_dict[starting_chain][0][1])
	"""
	model.add(interactions_dict[starting_chain][0][0]) #We add the chain with most interactions to the model.
	
	if stechiometry != None:
		for interaction in interactions_dict[starting_chain]:
			if interactions_dict[starting_chain][0][0].compare_sequence(interaction[1]) == 0:
				model.add(interaction[1])
				break
			else:
				continue
	else:
		model.add(interactions_dict[starting_chain][0][1])
	
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
	""" Calculates the clashes between two chains and returns whether the claches are in more than 5% of the atoms."""
	
	backbone = {"CA", "C1\'"}
	chain_atoms = [a for a in chain_atoms if a.id in backbone]  # Gets only the backbone atoms
	model_atoms = [a for a in model.get_atoms() if a.id in backbone]
	Nsearch = Bio.PDB.NeighborSearch(model_atoms)
	
	clash_count = 0
	for atoms in chain_atoms:
		clash = Nsearch.search(atoms.coord, 2) #Returns the atoms that clash
		if clash != []: #If there are atoms that clash, add 1 the the clash counter
			clash_count += 1
	if clash_count/len(chain_atoms) >= 0.05: #If the clashes are more than 5% of the atoms
		return True
	else:
		return False

	
def superimpose(unique_chains_list,interactions_dict, verbose=False, stechiometry=None):
	""" Main function that adds all the chains to create the final model."""
	if verbose:
		print("Starting to build the model.")
	if stechiometry != None:
		stech_dict, stech_file = get_stech_dicts(unique_chains_list, stechiometry)
		model = start_model(interactions_dict,stech_dict,verbose)
	else:
		model = start_model(interactions_dict,verbose)
	chain1, chain2 = model.get_chains()
	chain_in_model = [chain1.id,chain2.id]

	if verbose:
		print("Chains %s and %s have been selected to form the starting model." %(chain1.id, chain2.id))

	if stechiometry != None:
		#stech_dict, stech_file = get_stech_dicts(unique_chains_list, stechiometry)
		problematic_keys = {}
		for key in stech_file.keys():

			if len(stech_dict[key]) > int(stech_file[key]):
				problematic_keys[key] = 0
			elif len(stech_dict[key]) < int(stech_file[key]) and verbose:
				print("It's not possible to fulfill the stechiometry %s:%s due to an insufficient number of input chains. Only %s chains will be added to the model" %(key, stech_file[key], len(stech_dict[key])))
			
		for key in problematic_keys.keys():
			if chain1.id in stech_dict[key]:
				problematic_keys[key] += 1
				
			if chain2.id in stech_dict[key]:
				problematic_keys[key] += 1

	n = 2
	while n<len(interactions_dict.keys()):
		for chain in interactions_dict.keys():
			if chain in chain_in_model:
				continue

			if stechiometry != None and problematic_keys != {}:
				key_check = ""
				for key in problematic_keys.keys():

					if chain in stech_dict[key]:
						key_check = key
						break
				
				if key_check != "":
					if problematic_keys[key_check] == int(stech_file[key_check]): #if the number of chains is equal to the stechiometry, we don't add more chains.
						n += 1
						if verbose:
							print("Not adding chain %s due to stechiometry" %chain)
						continue
				
			not_added = True
			
			for chainin in chain_in_model:
				
				for chain_model, chain_interact in interactions_dict[chainin]: #if our chain interacts with a chain inside the complex
						
					if chain_interact.id == chain:
						chaininin = [x for x in model.get_chains() if x.id == chainin]
						newChain, modelChain = equal_length_chains(chain_model, chaininin[0])
						#Equaling the lentgh of the chains (needed for the Superimposer)
						superimpose = Bio.PDB.Superimposer() #Create a Superimposer instance
						superimpose.set_atoms(modelChain, newChain) #Create superimposition matrix
						chain_copy = chain_interact.copy() #Copy the chain so that we still keep the original version.
						chain_copy.id = chain
						superimpose.apply(chain_copy) #Apply the superposition matrix
						chain_atoms = sorted(chain_copy.get_atoms())
					
						if clashes(chain_atoms, model) == True:
							continue
						
						else:
							chain_copy.parent = None
							model.add(chain_copy)
							not_added = False
							chain_in_model.append(chain)

							if stechiometry != None and problematic_keys != {}:
								for key in problematic_keys.keys():
									if chain_copy.id in stech_dict[key]:
										problematic_keys[key] += 1
							n += 1

							if verbose:
								print("Chain %s sucessfully added to the model" %chain)

							break
						
			if verbose and not_added:
				print("%s could not be added to the model" %chain)
			
	return model
			
def save_PDB(model, output_path, verbose=False):
	""" Function that saves the model into a PDB file. """
	if verbose:
		print("Saving model")
	io = PDBIO()
	io.set_structure(model)
	filename = output_path + "/structures/model.pdb"
	io.save(filename)

	return filename
	
