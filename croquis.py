def superimpose(verbose=False, unique_chains_list):
	""" """
	int_dict = get_interactions_dict(unique_chains_list, verbose)
	starting_model = start_model(int_dict, verbose)
	model = Model("A") #Create new instance of class Model
	model.add(next(x for x in unique_chains_list if x.id == starting_model))
	chain_in_model = [starting_model]
	
	n = 0
	while n<len(int_dict.keys()):
		for chain in unique_chains_list:
			if chain in model:
				continue

			not_added = True
			for chainin in chain_in_model:
				
				if chain.id in int_dict[chainin]: #if our chain interacts with a chain inside the complex
					newChain, modelChain = equal_length_chains(chain, next(x for x in unique_chains_list if x.id == chainin))
					superimpose = Bio.PDB.Superimposer()
					superimpose.set_atoms(modelChain, newChain)
					chain_copy = chain.copy()
					superimpose.apply(chain_copy)
					chain_atoms = sorted(chain_copy.get_atoms())
					
					if clashes(chain_atoms, model) == True:
						continue
						
					else:
						#move.parent = None
						model.add(chain_copy)
						not_added = False

						if verbose:
							print("%s sucessfully added to the model" %chain.id)
						break
			if verbose and not_added:
				print("%s could not be added to the model" %chain.id)
			
			n += 1
			chain_in_model.append(chain.id)

			
def clashes(chain_atoms, model):
	""" """
	
	backbone = {"CA", "C1\'"}
   	chain_atoms = [a for a in move_atoms if atom.id in backbone]  # Gets only the backbone atoms
    	model_atoms = [a for a in model.get_atoms() if atom.id in backbone]
	Nsearch = Bio.PDB.NeighborSearch(model_atoms)
	
	clash_count = 0
	for atoms in chain_atoms:
		clash = Nsearch.search(atoms.coord, 2)
		if clash != []:
			clash_count += 1
	if clash_count/len(chain_atoms) >= 0.05:
		return True
	else:
		return False
	
	
	
	
	
	
