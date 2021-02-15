import argparse
import sys
import gzip
import re
from pathlib import Path
import os

parser = argparse.ArgumentParser(description="This program does BLA BLA BLA")

parser.add_argument('-i', '--input-directory',
					dest = "input_directory",
					action = "store",
					required = True,
					help = "Directory containing the input PDB files.")

parser.add_argument('-s', '--stechiometry',
					dest = "stechiometry",
					action = "store",
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


options = parser.parse_args()

class DirectoryError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return "%s is not a valid directory" %(self.value)

ind = Path(options.input_directory)

if ind.is_dir() == True:
	for file in os.listdir(options.input_directory):
		if re.match(r"\w+_\w_\w\.pdb[\.gz]?", file) == None:
			raise DirectoryError(options.input_directory)
else:
	raise DirectoryError(options.input_directory)

outd = Path(options.output_directory)

if outd.is_dir() == False:
	os.mkdir(options.output_directory)
	os.mkdir(options.output_directory + "/structures")
	os.mkdir(options.output_directory + "/analysis")
else:
	if options.force == True:
		out1 = Path(options.output_directory + "/structures")
		out2 = Path(options.output_directory + "/analysis")
		if out1.is_dir() == False:
			os.mkdir(options.output_directory + "/structures")
		if out2.is_dir() == False:
			os.mkdir(options.output_directory + "/analysis")
		
	else:
		print("The output directory already exists.")
		sys.exit()

style="max-width: 540px;"