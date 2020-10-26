#! /usr/bin/env python3
import argparse
import os

class Molecule:
	def __init__(self, fileName):
		self.calpha = {}
		self.hetatm = {}
		self.biomol = {}
		self.rotmat = None
		with open(fileName) as file:
			for index, line in enumerate(file):
				key = line[0:6].strip()
				if key in ["ATOM", "HETATM"]:
					atom = line[13:15].strip()
					atomNumber = line[6:11]
					chainName = line[21:22]
					if chainName not in self.calpha:
						self.calpha[chainName] = {}
						self.hetatm[chainName] = {}
					x = self.toFloat(line[30:38])
					y = self.toFloat(line[38:46])
					z = self.toFloat(line[46:54])
					occupancy = self.toFloat(line[54:60])
					temperature = self.toFloat(line[54:60])
					residue = line[17:20]
					residueNumber = self.toInt(line[22:26])
					atomDict = dict(x=x, y=y, z=z,
                     				residue=residue, occupancy=occupancy, temperature=temperature)
					if key == "ATOM":
						if "CA" == atom:
							self.calpha[chainName][residueNumber] = atomDict
					elif atom != "HOH":
						self.hetatm[chainName][atomNumber] = atomDict
				biokey = "REMARK 350 APPLY THE FOLLOWING TO CHAINS:"
				if self.rotmat is None and line[0:len(biokey)] == biokey:
					chains = line[len(biokey):].split(",")
					chains = list(map(lambda x: x.strip(), chains))
					biomolMat1 = self.getFloats(file)
					biomolMat2 = self.getFloats(file)
					biomolMat3 = self.getFloats(file)
					matrix = [biomolMat1, biomolMat2, biomolMat3]
					self.rotmat = dict(chains=chains, matrix=matrix)
		if not self.calpha:
			raise Exception("Molecule has no atoms.")

	def getFloats(self, file):
		nextLine = next(file)
		matLine = nextLine[23:].split()
		matLine = list(map(lambda x: self.toFloat(x), matLine))
		return matLine

	def createUnit(self):
		count = 0
		for chain in self.rotmat["chains"]:
			chainCopy = self.calpha[chain].copy()
			for atomNumber in chainCopy:
				atom = chainCopy[atomNumber]
				rotmat = self.rotmat["matrix"]
				newx = atom["x"] * rotmat[0][0] + atom["y"] * rotmat[0][1] + atom["z"] * rotmat[0][2] + rotmat[0][3]
				newy = atom["x"] * rotmat[1][0] + atom["y"] * rotmat[1][1] + atom["z"] * rotmat[1][2] + rotmat[1][3]
				newz = atom["x"] * rotmat[2][0] + atom["y"] * rotmat[2][1] + atom["z"] * rotmat[2][2] + rotmat[2][3]
				atom["x"] = newx
				atom["y"] = newy
				atom["z"] = newz
			self.calpha[str(count)] = chainCopy
			count = count + 1

	def toFloat(self, x):
		try:
			return float(x)
		except:
			return 0.0

	def toInt(self, x):
		try:
			return int(x)
		except:
			return 0

	def save(self, chainName, outputName):
		print ("Writing PDB file to %s." % outputName)
		f = open(outputName, "w")
		for residueNumber in sorted(self.calpha[chainName].keys()):  
			ca = self.calpha[chainName][residueNumber]
			if ca["residue"] is not None:
				atomString = "ATOM  %5d  CA  %s A %3d    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n" % (residueNumber, ca["residue"], residueNumber, ca["x"], ca["y"], ca["z"], ca["occupancy"], ca["temperature"], "C")
				f.write(atomString)
		f.close()

class Alignment:
	def __init__(self, fileName):
		self.queryName = None
		self.queryAlignment = []
		self.queryStart = []
		self.templateName = None
		self.templateAlignment = []
		self.templateStart = []
		self.readFile(fileName)

	def readFile(self, fileName):
		with open(fileName) as file:
			for index, line in enumerate(file):
				cols = line.split()
				if len(cols) > 1 and cols[0] == "Query":
					self.queryName = cols[1]
				if len(cols) > 1 and cols[0].startswith(">"):
					self.templateName = cols[0][1:]
				if self.queryName and self.templateName:
					if len(cols) > 3:
						if cols[0] == "Q" and cols[1] == self.queryName:
							self.queryStart.append(self.toInt(cols[2]))
							self.queryAlignment.append(cols[3])
						if cols[0] == "T" and cols[1] == self.templateName:
							self.templateStart.append(self.toInt(cols[2]))
							self.templateAlignment.append(cols[3])
				if len(cols) > 1 and cols[0] == "No" and cols[1] == "2":
					break

	def createModel(self, templateChain):
		for residueNumber in templateChain:
			templateChain[residueNumber]["residue"] = None
		for i in range(len(self.templateStart)):
			start = self.templateStart[i]
			querySequence = self.queryAlignment[i]
			templateSequence = self.templateAlignment[i]
			n = len(querySequence)
			tcount = 0
			for j in range(n):
				qs = querySequence[j]
				ts = templateSequence[j]
				if qs != "-" and ts != "-":
					residueNumber = start + tcount
					if residueNumber in templateChain:
						templateChain[residueNumber]["residue"] = toThreeAmino(qs)
					else:
						print ("Warning: Skipping invalid residue identifier [%s]." % residueNumber)
				if ts != "-":
					tcount += 1

	def toInt(self, x):
		try:
			return int(x)
		except:
			return 0

def main(args):
	templateMolecule = Molecule(args.a_template)
	if args.a_chain not in templateMolecule.calpha:
		raise Exception("Chain not found in template [%s]" % args.a_chain)
	templateMolecule.createUnit()
	os.system("mkdir -p temp")
	for key in templateMolecule.calpha.keys():
		templateMolecule.save(key, "temp/%s.pdb" % key)

	queryChain = templateMolecule.calpha[args.a_chain]
	queryAlignment = Alignment(args.a_result)
	queryAlignment.createModel(queryChain)
	outputName = "%s.pdb" % queryAlignment.queryName
	templateMolecule.save(args.a_chain, "temp/query.pdb")

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
	parser.add_argument('-ar', '--a_result', help='First HHR target file result', required=True)
	parser.add_argument('-ac', '--a_chain', help='First template chain name', required=False)
	parser.add_argument('-at', '--a_template', help='First template PDB', required=True)
	parser.add_argument('-br', '--b_result', help='Second HHR target file result', required=True)
	parser.add_argument('-bc', '--b_chain', help='Second structure chain name', required=False, default="A")
	parser.add_argument('-bt', '--b_template', help='Second template PDB', required=False, default="A")
	parser.add_argument('-t', '--template', help='Structure template', required=True)
	parser.add_argument('-o', '--output', help='Output PDB file', required=False)
	args = parser.parse_args()
	main(args)
 