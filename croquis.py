def superimpose(verbose=False, unique_chains_list):
	int_dict = get_interactions_dict(unique_chains_list, verbose)
	starting_model = start_model(int_dict, verbose)
	model = next(x for x in unique_chains_list if x.id == starting_model)
	chain_in_model = [starting_model]

	while n<len(unique_chains_list):
		for chain in unique_chains_list:
			if chain in model:
				continue
			for chain_interaction in chain_in_model:
				for chain_interaction2 in int_dict[chain.id]:
					if chain_interaction == chain_interaction2:
						







			n += 1
			chain_in_model.append(chain.id)