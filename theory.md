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



# Biological background

Obtaining the structure of a protein has been done for some years now through experimental procedures such as X-ray crystallography, NMR spectroscopy and electron microscopy. This has enabled to create a Protein Database (PDB) where many protein structures are stored. Thanks to this we can advance in the field of Structural Bioinformatics and obtain for example information on pairs of interacting chains of a macrocomplex which would allow us to create macrocomplex structures, which is exactly what this program is about. These macrocomplex structures are also known as the quaternary structure, where tertiary structures interact with each other to form a larger structure.

It should also be noted that it isn't as easy to obtain the structures for all proteins. For example, transmembrane proteins are harder than soluble proteins, which means that in the PDB there is a certain bias on the different proteins that we can find there. Furthermore, it is also very expensive and time-consuming to obtain PPIs.That is why being able to predict certain interactions with a program like ours can be important to help fill this gap of knowledge.

Examples of interacting forces would be hydrogen bonds, disulfide bonds but also Van der Waals forces or electrostatic forces. Examples of macrocomplexes that present these interactions would be hemoglobin, enhanceosome DNA-Protein interactions, receptors located at the membrane, etc.

Here, we will talk about the biological background that the program needed in order to be developed. 

## Macrocomplexes 

When talking about proteins, we know that most don't interact on their own, they form macrocomplexes. The chains usually interact in a way that they keep the hydrophobic residues together and expose the rest of residues to the solvent in order to obtain a more favourable structure. When taking this into account, we know that chains aren't in space on their own, they interact with each other in order to form the whole structure. For this reason, in our program we determined whether two chains were interacting or not, without taking into account if they had a joint pdb file, in order to detect all interacting chains for the whole complex, without the need to have many pdb files. We determined that two chains were interacting if they formed hydrogen bonds, which would mean that their atoms are at a length of 3.5 A or less (Narayan et al., 2000).Therefore, by knowing the pairs of interactions among them in space, we can reconstruct the whole stucture. Knowing the structures of these macrocomplexes will further our knowledge of protein-protein interactions (PPIs) and the biochemical functions they partake. 

## Superimposition

We assume that an interaction (A-B) will interact with another pair of interacting chains (for example, A-C) which we will then be able to superimpose (A-A) in order to obtain the structure with the 3 chains as shown in Figure X. 

![Figure X. Process for a superimposition of two interacting pairs of chains.](./superimp.png "Figure X. Process for a superimposition of two interacting pairs of chains.")

However, when doing a superimposition we have to watch out that the chain that we add does not interfere with the already created model/structure. This would be called a clash and we can see if the added chain presents a clash by seeing if the backbone atoms are in contact with each other at a distance of less than 2 A and if these clashes happen in more than 5% of the structure (Batsanov, 2001).

## Strengths

## Weaknesses

# Analysis

...

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
