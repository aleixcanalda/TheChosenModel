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

	def __init__(self,chain):
		self.child_list = chain.child_list
		self._id = chain._id
		self.parent = chain.parent
		self.interactions = []  # Here there will be stored the chain's known interactions as tuple of residues
		self.xtra = chain.xtra
		self.level = chain.level # The other attibutes are needed for biopython to be able to work.

	def get_sequence_chain(self):
		first_res = self.child_list[0].resname.strip()  # Get the first residue name to see what kind of sequence it is
		seq=""
		if first_res in self.protein: #if the first residue is in the protein dictionary, we have a protein sequence
			ppb=PPBuilder()
			for pp in ppb.build_peptides(self):
				seq = pp.get_sequence()
		
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
		alignments = pairwise2.align.globalxx(chain_seq,other_chain)[0] #we'll align the two sequences, and get the score, but only the best one
		align1 = alignments[0]
		align2 = alignments[1]
		ident = sum(base1 == base2 for base1, base2 in zip(align1, align2))  # Calculate number of identities
		if ident/len(align1) >= 0.95:
			return True
		else: 
			return False



if __name__ == "__main__":
	pdb = open("1gzx_A_B.pdb")

	parser = PDBParser(PERMISSIVE=1, QUIET=True)

	struct = parser.get_structure(file="1gzx_A_B.pdb",id="A")
	for model in struct:
		for chain in model:
			ch = MyChain(chain)
			print(ch.compare_sequence("AGTGCTGATGCTGTGCTAGTCGTA"))
