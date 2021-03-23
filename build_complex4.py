import sys
from classes import *
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
from Bio.PDB.Dice import *

parser = PDBParser(PERMISSIVE=1, QUIET=True)

def all_chains(PDB_files,verbose=False):
	""" Creating an object for every unique PDB chain. """

	unique_chain_list=[]
	id_set=set()
	nomen = {}
	if verbose:
		print("Processing PDB files and analyzing the chains")

	for file in PDB_files:
		for model in parser.get_structure("A",file): #Obtain the structures and iterate on all models
			chain_num = 0
			my_list = []
			for chain in model: #Obtain all chains from the models
				new_instance = MyChain(chain)
				same_same = False
				if nomen == {}:
					nomen[file[:6]] = [new_instance]

				elif chain_num == 0:
					if file[:6] not in nomen.keys():
						nomen[file[:6]] = [new_instance]
					else:
						nomen[file[:6]].append(new_instance)

				elif chain_num == 1 and len(model) == 2: #If it's the second chain and it's a protein pair.
					if file[7:13] not in nomen.keys():
						nomen[file[7:13]] = [new_instance]
					else:
						nomen[file[7:13]].append(new_instance)

				if len(model) == 3 and new_instance.get_type() == "dna":
					chain_num += 1
					my_list.append(new_instance)
					continue

				else:
					for chains in unique_chain_list:
						for chainss in chains:
							if new_instance.compare_sequence(chainss) == 2:
								same_same = True
								break
					if same_same == False:
						id_set.add(new_instance.id)
				chain_num += 1
				my_list.append(new_instance)

			if chain_num < 2:
				sys.stderr.write("All PDB files must contain at least two chains.")

			unique_chain_list.append(my_list)

	print("unique_chain_list")
	print(unique_chain_list)
	print("nomen")
	print(nomen)
	print("interactions_dict")

	return unique_chain_list, nomen

def chain_id(id_set):
	""" Function to create unique IDs """

	alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
				'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'AA', 'BB',
				'CC', 'DD', 'EE', 'FF', 'GG', 'HH', 'II', 'JJ', 'KK', 'LL', 'MM', 'NN',
				'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'UU', 'VV', 'WW', 'XX', 'YY', 'ZZ']
	for letter in alphabet:
		if letter not in id_set:
			return letter

def get_interactions_dict(unique_chain_list, PDB_files, verbose=False, template=None):
	"""
	Obtain the dictionary of interactions among all chains
	The dictionary will look as follows:
	'Chain1.id' = [[Chain1, Chain2 that interacts with chain 1],[Chain1,Chain3],...]
	.
	.
	.
	"""
	if verbose:
		print("Obtaining the interactions dictionary.")

	interactions_dict = {}

	for pair in unique_chain_list:

		if len(pair) == 3 and template != None:
			chain1 = pair[0]
			chain2 = pair[1]

			if interactions_dict == {}:
				interactions_dict[chain1.id]= [[chain1, chain2]]
			else:
				if chain1.id in interactions_dict.keys():
					interactions_dict[chain1.id].append([chain1,chain2])
				else:
					interactions_dict[chain1.id] = [[chain1, chain2]]

		elif  len(pair) == 3 and template == None:

			chain1 = pair[0]
			chain2 = pair[1]
			chain3 = pair[2]

			interact_1, interact_2 = chain1.interactions(chain2)
			if interact_1 == set() and interact_2 == set():
				continue

			if interactions_dict == {}:
				interactions_dict[chain1.id]= [[chain1, chain2],[chain1,chain3]]
				interactions_dict[chain2.id]= [[chain2, chain1],[chain2,chain3]]
				interactions_dict[chain3.id]= [[chain3, chain1],[chain3,chain2]]

			else:

				checking1 = False
				checking2 = False
				checking3 = False
				for key in interactions_dict.keys():
					check1 = chain1.compare_sequence(interactions_dict[key][0][0]) ##Dilemes de l'Aleix 2.0 (igual coordenada, igual seqüència)
					if check1 == 2 or chain1.id == key:
						interactions_dict[key].append([chain1, chain2])
						interactions_dict[key].append([chain1, chain3])
						checking1 = True
					check2 = chain2.compare_sequence(interactions_dict[key][0][0])
					if check2 == 2 or chain2.id == key:
						interactions_dict[key].append([chain2, chain1])
						interactions_dict[key].append([chain2, chain3])
						checking2 = True
					check3 = chain3.compare_sequence(interactions_dict[key][0][0])
					if check3 == 2 or chain3.id == key:
						interactions_dict[key].append([chain3, chain1])
						interactions_dict[key].append([chain3, chain2])
						checking3 = True
				if checking1 == False: ##WHAT WHAT Dilemes de la MAria 1.0 (per que check1 ja no es true?)
					interactions_dict[chain1.id]= [[chain1, chain2],[chain1,chain3]]
				if checking2 == False:
					interactions_dict[chain2.id]= [[chain2, chain1],[chain2,chain3]]
				if checking3 == False:
					interactions_dict[chain3.id]= [[chain3, chain1],[chain3,chain2]]


		else:
			chain1 = pair[0]
			chain2 = pair[1]

			#interact_1, interact_2 = chain1.interactions(chain2)
			#if interact_1 == set() and interact_2 == set():
			#	continue

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

	print(interactions_dict)
	return interactions_dict

def parse_stech(nomen,stechiometry=None,verbose=False):
	""" Obtain the stechiometry from the input file and from the interacting chains"""
	nomen_unique={}
	for keys in nomen.keys():
		for value in nomen[keys]:
			if nomen_unique == {}:
				nomen_unique[keys] = [value.id]
			else:
				if keys not in nomen_unique.keys():
					nomen_unique[keys] = [value.id]
				elif value.id not in nomen_unique[keys]:
					nomen_unique[keys].append(value.id)
	if verbose:
		print("Obtaining the stechiometry dictionary.")
	if stechiometry:
		stech_file={}
		fd = open(stechiometry, "r")
		for line in fd:
			line = line.strip()
			stech_file[line[:6]] = line[7:]
		fd.close()
		return stech_file, nomen_unique
	else:
		return nomen_unique

def start_model(interactions_dict,verbose=False, template=None):
	"""Choosing the chain with most interactions as the starting model of the macrocomplex"""

	model = modelmodel("A") #Create new instance of class Model

	if verbose:
		print("Deciding the starting model.")

	if template:
		for hey in parser.get_structure("A",template):
			print(hey)
			for chain in hey:
				print(chain)
				instance = MyChain(chain)
				model.add(instance)

	else:

		most_inter = 0

		for chain in interactions_dict.keys():

			chain_inter = len(interactions_dict[chain])

			if chain_inter > most_inter:

				most_inter = chain_inter

				starting_chain = chain

		model.add(interactions_dict[starting_chain][0][0]) #We add the chain with most interactions to the model.

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


def superimpose(model, fixed_chain, moving_chain, apply_chain):

	superimpose = Bio.PDB.Superimposer() #Create a Superimposer instance
	superimpose.set_atoms(fixed_chain, moving_chain) #Create superimposition matrix
	chain_copy = apply_chain.copy() #Copy the chain so that we still keep the original version.
	chain_copy.id = apply_chain.id
	superimpose.apply(chain_copy) #Apply the superposition matrix
	chain_atoms = sorted(chain_copy.get_atoms())
	if clashes(chain_atoms, model) == True:
		return model, True

	else:
		chain_copy.parent = None
		model.add(chain_copy)
		return model, False


def main_loop(unique_chains_list,interactions_dict, nomen,verbose=False, stechiometry=None):
	""" Main function that adds all the chains to create the final model."""
	if verbose:
		print("Starting to build the model.")
	if stechiometry != None:
		stech_file, nomen_unique = parse_stech(nomen,stechiometry)
		model = start_model(interactions_dict,verbose)
	else:
		model = start_model(interactions_dict,verbose)
	chain1, chain2 = model.get_chains()
	chain_in_model = [chain1.id,chain2.id]
	print(chain1)
	print(chain2)
	if verbose:
		print("Chains %s and %s have been selected to form the starting model." %(chain1.id, chain2.id))

	if stechiometry != None:
		#nomen, stech_file = get_nomen(unique_chains_list, stechiometry)
		problematic_keys = {}
		for key in stech_file.keys():

			if len(nomen_unique[key]) > int(stech_file[key]):
				problematic_keys[key] = 0
			elif len(nomen_unique[key]) < int(stech_file[key]) and verbose:
				print("It's not possible to fulfill the stechiometry %s:%s due to an insufficient number of input chains. Only %s chains will be added to the model" %(key, stech_file[key], len(nomen[key])))

		for key in problematic_keys.keys():
			if chain1.id in nomen_unique[key]:
				problematic_keys[key] += 1

			if chain2.id in nomen_unique[key]:
				problematic_keys[key] += 1

	n = 2
	while n<len(interactions_dict.keys()):
		for chain in interactions_dict.keys():
			if chain in chain_in_model:
				continue

			if stechiometry != None and problematic_keys != {}:
				key_check = ""
				for key in problematic_keys.keys():

					if chain in nomen_unique[key]:
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
				if chain in chain_in_model:
					break
				for chain_model, chain_interact in interactions_dict[chainin]: #if our chain interacts with a chain inside the complex

					if chain_interact.id == chain:
						chaininin = [x for x in model.get_chains() if x.id == chainin]
						newChain, modelChain = equal_length_chains(chain_model, chaininin[0])

						model, not_added = superimpose(model, modelChain, newChain, chain_interact)
						if not_added:
							continue
						else:
							chain_in_model.append(chain)

							if stechiometry != None and problematic_keys != {}:
								for key in problematic_keys.keys():
									if chain_copy.id in nomen[key]:
										problematic_keys[key] += 1
							n += 1
							if verbose:
								print("Chain %s sucessfully added to the model" %chain)

							break

			if verbose and not_added:
				print("%s could not be added to the model" %chain)
				n +=1

	return model


def template_loop(unique_chains_list, interactions_dict, nomen, template, verbose=False, stechiometry=None):
	""" Main function that adds all the chains to create the final model."""
	if verbose:
		print("Starting to build the model with template.")
	if stechiometry != None:
		stech_file, nomen_unique = parse_stech(nomen,stechiometry)
		model = start_model(interactions_dict,verbose, template)
	else:
		nomen_unique = parse_stech(nomen)
		model = start_model(interactions_dict,verbose, template)
	chain1, chain2 = model.get_chains()
	chain_in_model = [chain1.id,chain2.id]
	template_struct = Structure("B")
	template_struct.add(model)
	print(nomen_unique)
	if verbose:
		print("Template has been selected to form the starting model.")

	if stechiometry != None:
		#nomen, stech_file = get_nomen(unique_chains_list, stechiometry)
		problematic_keys = {}
		for key in stech_file.keys():

			if len(nomen_unique[key]) > int(stech_file[key]):
				problematic_keys[key] = 0
			elif len(nomen_unique[key]) < int(stech_file[key]) and verbose:
				print("It's not possible to fulfill the stechiometry %s:%s due to an insufficient number of input chains. Only %s chains will be added to the model" %(key, stech_file[key], len(nomen[key])))

		for key in problematic_keys.keys():
			if chain1.id in nomen_unique[key]:
				problematic_keys[key] += 1

			if chain2.id in nomen_unique[key]:
				problematic_keys[key] += 1


	for key in nomen_unique.keys(): #we go through eevry P19, Q189,...
		random.shuffle(nomen_unique[key]) #we shuffle the different proteins inside each P19,..
		for ids in nomen_unique[key]: #we go through every protein in P19 to add them to the template
			random.shuffle(nomen[key])
			for v in nomen[key]: #we get the id from inside the specific P19
				if v.id == ids:
					the_chosen_one = v
					break
			for probable_child in interactions_dict[the_chosen_one.id]:
				print(the_chosen_one)
				if probable_child[0] is the_chosen_one:
					chosen_brother = probable_child[1]
					break
			chosen_template_chain,start,end = chain1.compare_dna(chosen_brother, chain2)
			print(chosen_template_chain)
			print(start)
			print(end)
			Bio.PDB.Dice.extract(template_struct,chosen_template_chain.id, start, end,"template_dna_superimp.pdb")






	for chain in interactions_dict.keys():
		if chain in chain_in_model:
			continue
		if stechiometry != None and problematic_keys != {}:
			key_check = ""
			for key in problematic_keys.keys():
				if chain in nomen_unique[key]:
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
			if chain in chain_in_model:
				break
			for chain_model, chain_interact in interactions_dict[chainin]: #if our chain interacts with a chain inside the complex

				if chain_interact.id == chain:
					chaininin = [x for x in model.get_chains() if x.id == chainin]
					newChain, modelChain = equal_length_chains(chain_model, chaininin[0])

					model, not_added = superimpose(model, modelChain, newChain, chain_interact)
					if not_added:
						continue
					else:
						chain_in_model.append(chain)

						if stechiometry != None and problematic_keys != {}:
							for key in problematic_keys.keys():
								if chain_copy.id in nomen[key]:
									problematic_keys[key] += 1
						n += 1
						if verbose:
							print("Chain %s sucessfully added to the model" %chain)

						break

		if verbose and not_added:
			print("%s could not be added to the model" %chain)
			n +=1

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
