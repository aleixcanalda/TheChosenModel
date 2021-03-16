---
### BEGIN-Of-YAML-Block ###
#
## ######################################################################################
##
##   README.md
##
##     A LaTeX-extended MarkDown template for MScBI-ALG.
##
## ######################################################################################
##
##                 CopyLeft 2020 (CC:BY-NC-SA) --- Josep F Abril
##
##   This file should be considered under the Creative Commons BY-NC-SA License
##   (Attribution-Noncommercial-ShareAlike). The material is provided "AS IS", 
##   mainly for teaching purposes, and is distributed in the hope that it will
##   be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
##   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
## ######################################################################################
#
# title-meta is title string without any LaTeX format to be used as pdftitle, part of emails subject...
title-meta: Macro-complex structure modeller
#
# title is the big title for the cover page, fully LaTeX formated to fit into a shortstack command...
title: |
  \textsc{Macro-complex structure modeller}
#subtitle:
#
# runtitle is the running header or footer, used i.e. by fancyheadings...
runtitle: |
  Macro-complex structure modeller
#
# author-meta sets the pdfauthor variable...
author-meta: !!str 'Aleix Canalda, Maria Díaz'
#
# authors to appear on the title page...
author:
- name: Aleix Canalda Baltrons
author:
- name: Maria Díaz Ros
#
# authorshort defines a brief list of authors, i.e. for headings
authorshort: YOURSURNAME, YOURNAMEshort
#
### end-Of-YAML-Block ###
---

# Introduction

We have developed a program to model the macro-complex structure of biomolecules, including proteins and DNA, using as an input the sequences of the pair interactions of the complex. This manual will include a tutorial explaining how to use the program using some examples, a theoretical explanation of the algorithm and the biological background behind it and, lastly, an analysis of the results of some examples, to assess the quality of the solutions our program offers.

# Tutorial

## Installation
...

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
```

### Input (mandatory argument)

As we can see above, the input has to be a directory containing only PDB files. These files have to include two interacting chains from the model and the name of each file has to follow a set structure: name_chain1_chain2.pdb(.gz) where the name is an alphanumerical string and the chains must coincide with the IDs of the chains inside the file. As it can be seen, the PDB files can either be compressed or not.

### Stechiometry

The stechiometry file must contain a line for each chain and its associated number, such as:

A:2

B:2

### Output directory (mandatory argument)

If the directory doesn't exist, it will be created automatically but if the user introduces an already existing directory, it will raise an error unless he specifies that he wants to overwrite the content (using the -force argument).

## Examples


# Algorithm

Our program obtains at least the input and output in order to work correctly and with that creates a list with all the chains that are interacting, nevertheless, the chains are saved as a new type of class which has the same functions as the PDB chain class, but with some extra ones that will help us out throughout our program. After we have all the chains in a list, we want to know which ones interact among each other and create a dictionary with all the interacting chains. After that, we want to create the starting final model, and we will do so by adding two chains that we know are interacting, from the chain that presents the most interactions. That way, we choose as a starting model a chain that is related to various other chains and will make the process of adding all the chains and searching for interacting chains onto the structure much more efficient. Moreover, by doing so we will avoid choosing a starting chain that doesn't present interactions at all. Once a starting model has been chosen, we loop through all the interactions and see if we can superimpose the rest of chains inside the interaction dictionary with the model in an iterative fashion. 

The program will superimpose two identical subunits (in other words, chains that have at least 95% similarity) and along with them, their respective interacting chain, applying to the former the transformation matrix applied to the identical chains' superimposition. This way, when superimposing a chain, our model checks whether its respective interacting chain presents clashes with the rest of the model, in order to keep the chain or not add it at all. We will end up then with all the chains that interact and don't clash with the previous model as the final model and we'll be ready to check how well the model turns out.

# Biological background

Obtaining the structure of a protein has been done for some years now through experimental procedures such as X-ray crystallography, NMR spectroscopy and electron microscopy. This has enabled to create a Protein Database (PDB) where many protein structures are stored. Thanks to this we can advance in the field of Structural Bioinformatics and obtain for example information on pairs of interacting chains of a macrocomplex which would allow us to create macrocomplex structures, which is exactly what this program is about. These macrocomplex structures are also known as the quaternary structure, where tertiary structures interact with each other to form a larger structure.

It should also be noted that it isn't as easy to obtain the structures for all proteins. For example, transmembrane proteins are harder than soluble proteins, which means that in the PDB there is a certain bias on the different proteins that we can find there. Furthermore, it is also very expensive and time-consuming to obtain PPIs.That is why being able to predict certain interactions with a program like ours can be important to help fill this gap of knowledge.

Examples of interacting forces would be hydrogen bonds, disulfide bonds but also Van der Waals forces or electrostatic forces. Examples of macrocomplexes that present these interactions would be hemoglobin, enhanceosome DNA-Protein interactions, receptors located at the membrane, etc.

Here, we will talk about the biological background that the program needed in order to be developed. 

## Macrocomplexes 

When talking about proteins, we know that most don't interact on their own, they form macrocomplexes. The chains usually interact in a way that they keep the hydrophobic residues together and expose the rest of residues to the solvent in order to obtain a more favourable structure. When taking this into account, we know that chains aren't in space on their own, they interact with each other in order to form the whole structure. For this reason, in our program we determined whether two chains were interacting or not, without taking into account if they had a joint pdb file, in order to detect all interacting chains for the whole complex, without the need to have many pdb files. We determined that two chains were interacting if they formed hydrogen bonds, which would mean that their atoms are at a length of 3.5 A or less (Narayan et al., 2000).Therefore, by knowing the pairs of interactions among them in space, we can reconstruct the whole stucture. Knowing the structures of these macrocomplexes will further our knowledge of protein-protein interactions (PPIs) and the biochemical functions they partake. 

## Superimposition

We assume that an interaction (A-B) will interact with another pair of interacting chains (for example, A-C) which we will then be able to superimpose (A-A) in order to obtain the structure with the 3 chains as shown in Figure X. 

![Figure X. Process for a superimposition of two interacting pairs of chains.](./images/superimp.png "Figure X. Process for a superimposition of two interacting pairs of chains.")

However, when doing a superimposition we have to watch out that the chain that we add does not interfere with the already created model/structure. This would be called a clash and we can see if the added chain presents a clash by seeing if the backbone atoms are in contact with each other at a distance of less than 2 A and if these clashes happen in more than 5% of the structure (Batsanov, 2001).

## Strengths

## Weaknesses

# Analysis

With our program we are able to analyze the DOPE energy of each model that we create with the optional argument "-e, --energy" and we obtain an energy profile for each model that uses MODELLER to obtain it on a window of 13 residues. The DOPE energy can be seen on the terminal screen. An energy that is very negative overall means that a good model was produced, whereas the less negative the energy, the worse the model is, which makes it easy to compare and find the model that is more suitable. 

Another way to analyze the model is, in those cases where the complete structure is available for download in the PDB, by superimposing with a program such as Chimera the model obtained with the complete known structure. By doing so, we can obtain the RMSD from this superimposition, where the smaller it is, being 0.00 the smallest value, the better the model. In this case, an RMSD of 0.00 would mean that the model is identical to the complete structure, so the reconstruction of the macrocomplex has been done perfectly.

With our examples, we decided to analyze the models both ways:

\begin{itemize}
  \item \textbf{1gzx}. This macrocomplex is also known as hemoglobin when oxygen is bound to all four chains. The main stechiometry is A2B2, however the chains in the example are presented as A, B, C, D. After running it through our program we obtained a perfect reconstruction, with a good energy profile and an RMSD of 0.00, as seen in Figures X and Y
  ![Figure X. Energy profile of 1gzx model.](./images/DOPE_energy_1gzx.png "Figure X. Energy profile of 1gzx model.")
  ![Figure Y. Chimera superimposition of 1gzx model with full structure.](./images/1gzxsuperimp.png "Figure Y. Chimera superimposition of 1gzx model with full structure.")
  \item \textbf{6gmh}. This macrocomplex is the structure of an activated transcription complex Polymerase II. Its main stechiometry is A1B1C1D1E1F1G1H1I1J1K1L1M1N1O1P1Q1R1S1T1, which are the chains that were available for our program. In this example we also obtained a good profile and an RMSD of 0.00, as seen in Figures X and Y.
  ![Figure X. Energy profile of 6gmh model.](./images/DOPE_energy_6gmh.png "Figure X. Energy profile of 6gmh model.")
  ![Figure Y. Chimera superimposition of 6gmh model with full structure.](./images/6gmhsuperimp.png "Figure Y. Chimera superimposition of 6gmh model with full structure.")



\end{itemize}

# References

Narayanan Eswar, C. Ramakrishnan, Deterministic features of side-chain main-chain hydrogen bonds in globular protein structures, Protein Engineering, Design and Selection, Volume 13, Issue 4, April 2000, Pages 227–238.

Batsanov S.S.; Van der Waals Raddi of Elements, Inorganic Materials, 2001.

```{.sh .numberlines startFrom="100"}
# 
pandoc -f $PDOCFLGS                          \
       --template=./template/readme.tex      \
       -t latex --natbib --listings          \
       --number-sections                     \
       --variable papersize:a4paper          \
       --variable toc=true                   \
       --variable lof=true                   \
       --variable lot=true                   \
       --variable geometry:margin=1.5cm      \
       --variable fontsize=10pt              \
       -o $RF.tex $RF.md;                    \
pdflatex $RF.tex; bibtex $RF; pdflatex $RF.tex;  pdflatex $RF.tex;
#
# --highlight-style pygments
# --highlight-style tango
# ...
```
