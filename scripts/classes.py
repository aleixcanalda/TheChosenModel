#!/usr/bin/env python
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
import Bio.PDB.NeighborSearch
import sys 
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio import pairwise2

class MyChain(Chain):
	""" A customized chain PDB class """
	protein = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N','GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W','ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', "UNK": "X"}
	dna = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T'}
	rna = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}

	def __init__(self, chain):
		self.child_list = chain.child_list
		self._id = chain._id
		self.parent = chain.parent
		self.xtra = chain.xtra
		self.level = chain.level # The other attibutes are needed for biopython to be able to work.

	def get_type(self):
		first_res = self.child_list[0].resname.strip()  # Get the first residue name to see what kind of sequence it is
		if first_res in self.protein: #if the first residue is in the protein dictionary, we have a protein sequence
			return "protein"
		
		elif first_res in self.dna: #if the first residue is in the dna dictionary, we have a dna sequence
			return "dna"

		else:
			return "rna"

	def get_sequence_chain(self):
		first_res = self.child_list[0].resname.strip()  # Get the first residue name to see what kind of sequence it is
		seq=""
		if first_res in self.protein: #if the first residue is in the protein dictionary, we have a protein sequence
			for res in self:
				if res.id[0] == " ": #to avoid getting the HETATOM
					seq += self.protein[res.resname.strip()]
		
		elif first_res in self.dna: #if the first residue is in the dna dictionary, we have a dna sequence
			for res in self:
				if res.id[0] == " ": #to avoid getting the HETATOM
					seq += self.dna[res.resname.strip()]

		else:
			for res in self:
				if res.id[0] == " ":
					seq += self.rna[res.resname.strip()]
		return seq

	def compare_sequence(self,other_chain):
		chain_seq = self.get_sequence_chain()
		other_chain_seq = other_chain.get_sequence_chain()
		alignments = pairwise2.align.globalxx(chain_seq,other_chain_seq)[0] #we'll align the two sequences, and get the score, but only the best one
		align1 = alignments[0]
		align2 = alignments[1]
		ident = sum(base1 == base2 for base1, base2 in zip(align1, align2))  # Calculate number of identities
		if ident/len(align1) >= 0.95:
			#if the chains are homolog (have the same chain in stechiometry) but are in a different position in space, we'll still keep them
			if (self.child_list[0].child_list[0].get_vector()[0],self.child_list[0].child_list[0].get_vector()[1]) == (other_chain.child_list[0].child_list[0].get_vector()[0],self.child_list[0].child_list[0].get_vector()[1]):
				return 2
			else:
				return 1
		else: 
			return 0
	
	def compare_dna(self,other_chain, template_chain):
		""" Compares a template DNA with a fragment DNA and gives out the score and the specific position """
		chain_seq = self.get_sequence_chain()
		other_chain_seq = other_chain.get_sequence_chain()
		template_chain_seq = template_chain.get_sequence_chain()
		alignments = pairwise2.align.localms(chain_seq,other_chain_seq, 10, -10, -7, -5)[0] 
		alignments2 = pairwise2.align.localms(template_chain_seq,other_chain_seq, 10, -10, -7, -5)[0] #we will use these scores and a local alignment because our sequence aligns with a part of the template
		if alignments[2] > alignments2[2]: #we will compare the alignment score to see which strand we're looking at, plus send the start and end positions
			return self,alignments[3], alignments[4]
		else:
			return template_chain,alignments2[3], alignments2[4]


	def interactions(self, other_chain):
		""" Searches all the atoms interacting between two chains"""
		atoms_self = set()
		atoms_other = set()
		atom_list = [a for a in self.get_atoms()] 
		Nsearch = Bio.PDB.NeighborSearch(atom_list) #Creates an instance of NeighborSearch with the atoms of one chain as query
		
		for atom in other_chain.get_atoms():
			interaction = Nsearch.search(atom.coord, 4) #For all atoms of the other chain, assess if they're close to the other chain
			if interaction != []:
				atoms_other.add(atom)
				atoms_self.update(interaction)

		return (atoms_self, atoms_other)

