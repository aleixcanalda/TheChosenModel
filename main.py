from IOinterface import *
from classes import *
from build_complex import *
from energies import *

prots_files = get_input_file(options.input_directory)

get_output_file(options.output_directory, options.force)

chain_list, nomen = all_chains(prots_files, options.verbose)

int_dict = get_interactions_dict(chain_list, prots_files, options.verbose, options.template)

model = main_loop(chain_list, int_dict, options.verbose, options. stechiometry, options.template)

filename = save_PDB(model, options.output_directory, options.verbose)

if options.energy:
	energy = DOPE_Energy(filename, options.output_directory, options.verbose)
	get_profile(energy, options.output_directory)
