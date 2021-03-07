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
```

### Input (mandatory argument)

As we can see above, the input has to be a directory containing only PDB files. These files have to include two interacting chains from the model and the name of each file has to follow a set structure: name_chain1_chain2.pdb(.gz) where the name is an alphanumerical string and the chains must coincide with the IDs of the chains inside the file. As it can be seen, the PDB files can either be compressed or not.

### Stechiometry

The stechiometry file must contain a line for each chain and its associated number, such as:

A:2

B:2

### Output directory (mandatory argument)

If the directory doesn't exist, it will be created automatically but if the user introduces an already existing directory, it will raise an error unless he specifies that he want to overwrite the content (using the -force argument).



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
