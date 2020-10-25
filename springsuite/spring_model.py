#! /usr/bin/env python3
import argparse
import os

class Molecule:
	def __init__(self, fileName):
		self.calpha = []
		with open(fileName) as file:
			for index, line in enumerate(file):
				if "CA" == line[13:15]:
					x = self.toFloat(line[30:38])
					y = self.toFloat(line[38:46])
					z = self.toFloat(line[46:54])
					occupancy = self.toFloat(line[54:60])
					temperature = self.toFloat(line[54:60])
					residue = line[17:20]
					self.calpha.append(dict(x=x, y=y, z=z,
											residue=residue,
											occupancy=occupancy,
											temperature=temperature))

	def toFloat(self, x):
		try:
			return float(x)
		except:
			return 0.0

class Alignment:
	def __init__(self, fileName):
		self.queryName = None
		self.queryAlignment = []
		self.queryAlignmentStart = []
		self.templateName = None
		self.templateAlignment = []
		self.templateAlignmentStart = []
		self.readFile(fileName)

	def readFile(self, fileName):
		with open(args.query) as file:
			for index, line in enumerate(file):
				cols = line.split()
				if len(cols) > 1 and cols[0] == "Query":
					self.queryName = cols[1]
				if len(cols) > 1 and cols[0].startswith(">"):
					self.templateName = cols[0][1:]
				if self.queryName and self.templateName:
					if len(cols) > 3:
						if cols[0] == "Q" and cols[1] == self.queryName:
							self.queryAlignmentStart.append(cols[2])
							self.queryAlignment.append(cols[3])
						if cols[0] == "T" and cols[1] == self.templateName:
							self.templateAlignmentStart.append(cols[2])
							self.templateAlignment.append(cols[3])
				if len(cols) > 1 and cols[0] == "No" and cols[1] == "2":
					break

def main(args):
	# load template molecule
	templateMolecule = Molecule(args.template)
	print (templateMolecule.calpha)
	alignment = Alignment(args.query)
	print ('\n'.join(alignment.templateAlignment))

def toThreeAmino (seq):
	code = dict(G="GLY", A="ALA", V="VAL", L="LEU", I="ILE", M="MET", F="PHE", P="PRO", Y="TYR", W="TRP",
            	K="LYS", S="SER", C="CYS", N="ASN", Q="GLN", H="HIS", T="THR", E="GLU", D="ASP", R="ARG")
	return code[seq] if seq in code else "XXX"

def toSingleAmino (seq):
	code = dict(GLY="G", ALA="A", VAL="V", LEU="L", ILE="I", MET="M", PHE="F", PRO="P", TYR="Y", TRP="W",
            	LYS="K", SER="S", CYS="C", ASN="N", GLN="Q", HIS="H", THR="T", GLU="E", ASP="D", ARG="R")
	return code[seq] if seq in code else "X"

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Create a 3D model from HH-search results.')
	parser.add_argument('-q', '--query', help='HHR target file result', required=True)
	parser.add_argument('-t', '--template', help='Structure template', required=True)
	args = parser.parse_args()
	main(args)