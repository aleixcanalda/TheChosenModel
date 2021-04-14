# TheChosenModel
## Software Requirements

To be able to run the program smoothly and correctly, there are certain requirements that need to be met:

 * Python version 3.6

The following python modules:

 * Argparse

 * Biopython

 * re

 * random

 * gzip

 * os

 * sys

In order to run the -e --energy argument:

 * Modeller version 10.1

 * pylab

## Download and Installation

You can download the package to start its use with the following code:

```{.sh}

git clone https://github.com/aleixcanalda/TheChosenModel.git
cd TheChosenModel

```
There should be a script called setup.py with which we will do the installation. Nevertheless, before doing so, we should make sure that we have the requirements described above, otherwise the installation won't work.

```{.sh}

sudo python3 setup.py install

```
## Command-line arguments

These are all the arguments that can be introduced to our program:

```{.sh}

  -h, --help            show this help message and exit

  -i INPUT_DIRECTORY, --input-directory INPUT_DIRECTORY
                        Directory containing the input PDB files

  -s STECHIOMETRY, --stechiometry STECHIOMETRY
                        Path to a file containing the stechiometry of the complex.

  -o OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                        Directory where the output will be saved.

  -f, --force           Overwrite the content of the output directory.

  -v, --verbose         Print the progession of the execution.
  
  -e, --energy          Calculate DOPE energy and plot the result.
  
  -t --template-DNA     DNA template for models with only protein-DNA interactions.
   
  -m --models           Number of resulting models.
```

### Input (mandatory argument)

As we can see above, the input has to be a directory containing only PDB files. These files have to include two interacting chains from the model and the name of each file has to follow a set structure. For protein-protein interactions the structure is: protein1_protein2.name_chain1_chain2.pdb(.gz) where protein1 is an alphanumerical string that will be used to identify the stechiometry of the structure, the name is an alphanumerical string (name of the complex) and the chains must coincide with the IDs of the chain IDs inside the file. As it can be seen, the PDB files can either be compressed or not. Additionally, there is the option of introducing input files with DNA interactions following this structure (when there is a combination of protei-protein and protein-DNA interactions). Then, instead of protein2 being any string, it will  put "DNA" to indicate that it's a DNA strand and chain2 will still be the chain ID inside the PDB file.

Example of input files for the protein macrocomplex 1GZX:

P69905_P68871.1gzx_A_B.pdb

P69905_P69905.1gzx_A_C.pdb

P69905_P68871.1gzx_A_D.pdb

For the cases with DNA-protein interactions exclusively (there are no protein-protein interactions) the structure is: protein.DNA.name_chain1_SenseAntisense.pdb(.gz). The protein is also an alphanumerical string that will be used to identify the stechiometry of the structure, then DNA to identify that it's a protein-DNA interaction and the name of the structure (name). Chain1 is the ID of the protein chain and it must coincide with the chain ID inside the file. SenseAntisense are the IDs of the two strands of DNA and they must coincide with the ID inside the file.

Example of input files for protein-DNA interactions (from the macro-complex 2O61):

P05412.DNA.1t2k_C_EF.pdb

Q04206.DNA.5u01_A_EF.pdb

### Stechiometry

The stechiometry is a .txt file that must contain a line for each protein using the same nomenclature used in the input files and it will indicate how many times the protein has to appear in the final model. Example of an stechiometry file:

P69905:2

P68871:2

### Output directory (mandatory argument)

If the directory doesn't exist, it will be created automatically but if the user introduces an already existing directory, it will raise an error unless he specifies that he wants to overwrite the content (using the -force argument).

### Energy
It calculates the DOPE score for each model and it builds a DOPE profile used to create a plot of the DOPE score along the sequence. To use this argument is necessary to have the modules of modeller and pylab installed.

### Template DNA
It is a mandatory argument when the model consists only of protein-DNA interactions that form a long strand of DNA where the proteins interact. The input is the path to a file containing the both strands of a DNA sequence as a PDB file.

## Examples

Our first example will be with the protein macro-complex of the oxy T state of haemoglobin, with the PDB ID **1GZX**. This is a structure made with X-ray diffraction with a resolution of 2.10 Å. The structure was deposited by the authors in 2002 from their paper about the T state haemoglobin (Paoli et al., 1996).

To create the complex running TheChosenModel in the terminal, the following code has to be executed:

```{.sh}

TheChosenModel.py -i home/examples/example1/1gzx/ -o home/examples/results/1gzx -f -v -e

```
With this line what we are doing is giving as input (-i) the path to a folder containing the pairs of chain interactions that will form the macro-complex 1GZX (correct naming structure explained in the input argument section), the path to the output folder where we want to store the results (-o) and, if this folder already exists, we'll overwrite it (-f). Also, we will show on the terminal the different steps that the program is going through, to keep track of the process, with the verbose argument (-v). And lastly, we will calculate the energies of the resulting model as well as its energy profile.

The other example we will show is the macro-complex with a PDB ID **2O61**. This complex is the crystal structure, obtained with X-ray diffraction, of NFkB, IRF7 and IRF3 bound to the interferon-b enhancer. The resolution of the complex is 2.80 Å and it was deposited into the database in 2006 by the authors of a paper about the activation of the interferon beta (IFN-beta) gene. Its activation requires an assembly of an enhanceosome containing ATF-2/c-Jun, IRF-3/IRF-7, and NFkappaB. (Panne et al., 2007)

This macro-complex is formed solely by DNA-protein interactions so, in order for the program to work, it is necessary to add a DNA template that will help with the superimposition process. The commands used to run the example are:

```{.sh}

TheChosenModel.py -i home/examples/example2/2O61/ -o home/examples/results/2O61 -t home/examples/example2/template_file.pdb -m 4 -e -s home/examples/example2/stechiometry_file.txt

```
With this line what we are doing is first giving as input (-i) the folder containing the chain interactions that will form the macro-complex 2O61 and the output folder where we want to store the results (-o). Also, we will send the template ("template_file") for the DNA (-t), we will receive four different models (-m 4) and calculate the DOPE energy for each of them. Lastly, the path to an stechiometry file is provided (-s, required syntax for the stechiometry explained in the stechiometry argument section) and the program will try to satisfy it the best it can.
