#! /usr/bin/env python3
import argparse
import os

class Molecule:
	def __init__(self, fileName = None):
		self.calpha = {}
		self.hetatm = {}
		self.biomol = {}
		self.rotmat = []
		self.atoms = []
		if fileName is not None:
			self.fromFile(fileName)

	def fromFile(self, fileName):
		biomolFound = False
		with open(fileName) as file:
			for index, line in enumerate(file):
				key = line[0:6].strip()
				if key in ["ATOM", "HETATM"]:
					atom = line[12:16]
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
					atomNumber = self.toInt(line[6:11])
					atomName = line[12:16]
					atomDict = dict(x=x, y=y, z=z,
									residue=residue,
									occupancy=occupancy,
									temperature=temperature,
									atomNumber=atomNumber,
									atomName=atomName,
									residueNumber=residueNumber)
					if key == "ATOM":
						if atom.strip() == "CA":
							self.calpha[chainName][residueNumber] = atomDict
					elif atom.strip() != "HOH":
						self.hetatm[chainName][atomNumber] = atomDict
					self.atoms.append(atomDict)
				biokey = "REMARK 350 BIOMOLECULE:"
				if line[0:len(biokey)] == biokey:
					biomolFound = self.toInt(line[len(biokey):]) == 1
				biokey = "REMARK 350 APPLY THE FOLLOWING TO CHAINS:"
				if biomolFound and line[0:len(biokey)] == biokey:
					chains = line[len(biokey):].split(",")
					chains = list(map(lambda x: x.strip(), chains))
					biomolMatId, biomolMat1 = self.getFloats(file)
					if biomolMatId == 1:
						biomolMatId, biomolMat2 = self.getFloats(file)
						biomolMatId, biomolMat3 = self.getFloats(file)
						matrix = [biomolMat1, biomolMat2, biomolMat3]
						self.rotmat.append(dict(chains=chains, matrix=matrix))
		if not self.calpha:
			raise Exception("Molecule has no atoms.")

	def getFloats(self, file):
		nextLine = next(file)
		matId = self.toInt(nextLine[20:23])
		matLine = nextLine[23:].split()
		matLine = list(map(lambda x: self.toFloat(x), matLine))
		return matId, matLine

	def createUnit(self):
		molecule = Molecule()
		for matrixDict in self.rotmat:
			for chain in matrixDict["chains"]:
				chainCopy = self.calpha[chain].copy()
				for atomNumber in chainCopy:
					atom = chainCopy[atomNumber]
					rotmat = matrixDict["matrix"]
					self.applyMatrix(atom, rotmat)
				molecule.calpha[chain] = chainCopy
		return molecule

	def applyMatrix(self, atom, rotmat):
		newx = atom["x"] * rotmat[0][0] + atom["y"] * rotmat[0][1] + atom["z"] * rotmat[0][2] + rotmat[0][3]
		newy = atom["x"] * rotmat[1][0] + atom["y"] * rotmat[1][1] + atom["z"] * rotmat[1][2] + rotmat[1][3]
		newz = atom["x"] * rotmat[2][0] + atom["y"] * rotmat[2][1] + atom["z"] * rotmat[2][2] + rotmat[2][3]
		atom["x"] = newx
		atom["y"] = newy
		atom["z"] = newz
		return atom

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

	def saveChain(self, chainName, outputName):
		print ("Writing PDB file to %s." % outputName)
		f = open(outputName, "w")
		for residueNumber in sorted(self.calpha[chainName].keys()):  
			ca = self.calpha[chainName][residueNumber]
			if ca["residue"] is not None:
				f.write(self.atomString(ca))
		f.close()

	def save(self, outputName):
		print ("Writing atoms to PDB file to %s." % outputName)
		f = open(outputName, "w")
		for atom in self.atoms:
			f.write(self.atomString(atom))
		f.close()

	def atomString(self, atom):
		return "ATOM  %5d %s %s A %3d    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (atom["atomNumber"], atom["atomName"], atom["residue"], atom["residueNumber"], atom["x"], atom["y"], atom["z"], atom["occupancy"], atom["temperature"])

class Alignment:
	def __init__(self, fileName):
		self.queryName = None
		self.alignment = []
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
							self.alignment.append(cols[3])
						if cols[0] == "T" and cols[1] == self.templateName:
							self.templateStart.append(self.toInt(cols[2]))
							self.templateAlignment.append(cols[3])
				if len(cols) > 1 and cols[0] == "No" and cols[1] == "2":
					break

	def createModel(self, templateChain):
		for residueNumber in templateChain:
			templateChain[residueNumber]["residue"] = None
		for i in range(len(self.templateStart)):
			templateStart = self.templateStart[i]
			templateSequence = self.templateAlignment[i]
			queryStart = self.queryStart[i]
			querySequence = self.alignment[i]
			n = len(querySequence)
			tcount = 0
			qcount = 0
			for j in range(n):
				qs = querySequence[j]
				ts = templateSequence[j]
				if qs != "-" and ts != "-":
					residueNumber = templateStart + tcount
					if residueNumber in templateChain:
						templateChain[residueNumber]["residue"] = self.toThreeAmino(qs)
						templateChain[residueNumber]["residueNumber"] = queryStart + qcount
					else:
						print ("Warning: Skipping invalid residue identifier [%s]." % residueNumber)
				if qs != "-":
					qcount = qcount + 1
				if ts != "-":
					tcount = tcount + 1

	def toInt(self, x):
		try:
			return int(x)
		except:
			return 0

	def toThreeAmino (self, seq):
		code = dict(G="GLY", A="ALA", V="VAL", L="LEU", I="ILE", M="MET", F="PHE", P="PRO", Y="TYR", W="TRP",
					K="LYS", S="SER", C="CYS", N="ASN", Q="GLN", H="HIS", T="THR", E="GLU", D="ASP", R="ARG")
		return code[seq] if seq in code else "XXX"

	def toSingleAmino (self, seq):
		code = dict(GLY="G", ALA="A", VAL="V", LEU="L", ILE="I", MET="M", PHE="F", PRO="P", TYR="Y", TRP="W",
					LYS="K", SER="S", CYS="C", ASN="N", GLN="Q", HIS="H", THR="T", GLU="E", ASP="D", ARG="R")
		return code[seq] if seq in code else "X"

def TMalign(fileA, fileB, outputName):
	os.system("tmalign/TMalign %s %s -m temp/tmalign.mat > temp/tmalign.out" % (fileA, fileB))
	rotmat = []
	with open("temp/tmalign.mat") as file:
		line = next(file)
		line = next(file)
		for i in range(3):
			line = next(file)
			rotmatLine = line[1:].split()
			rotmatLine = list(map(lambda x: float(x), rotmatLine))
			rotmatLine = [rotmatLine[1], rotmatLine[2], rotmatLine[3], rotmatLine[0]]
			rotmat.append(rotmatLine)
	molecule = Molecule(fileA)
	for atom in molecule.atoms:
		molecule.applyMatrix(atom, rotmat)
	molecule.save(outputName)
	return molecule

def buildModel(resultFile, templateFile, chainName, outputName):
	template = Molecule(templateFile)
	if chainName not in template.calpha:
		raise Exception("Chain not found in template [%s]" % chainName)
	chain = template.calpha[chainName]
	alignment = Alignment(resultFile)
	alignment.createModel(chain)
	template.saveChain(chainName, outputName)
	os.system("./pulchra/pulchra %s" % outputName)

def main(args):
	os.system("mkdir -p temp")
	templateMolecule = Molecule(args.template)
	bioMolecule = templateMolecule.createUnit()
	for key in bioMolecule.calpha.keys():
		bioMolecule.saveChain(key, "temp/template%s.pdb" % key)
	buildModel(args.a_result, args.a_template, args.a_chain, "temp/modelA.pdb")
	coreMolecule = TMalign("temp/modelA.rebuilt.pdb", "temp/template%s.pdb" % args.template_core, "temp/tmalignCore.pdb")
	buildModel(args.b_result, args.b_template, args.b_chain, "temp/modelB.pdb")
	partners = []
	for chainName in bioMolecule.calpha.keys():
		if chainName != args.template_core:
			partnerMolecule = TMalign("temp/modelB.rebuilt.pdb", "temp/template%s.pdb" % chainName, "temp/tmalign%s.pdb" % chainName)
			partners.append(partnerMolecule)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Create a 3D model from HH-search results.')
	parser.add_argument('-ar', '--a_result', help='First HHR target file result', required=True)
	parser.add_argument('-ac', '--a_chain', help='First template chain name', required=False)
	parser.add_argument('-at', '--a_template', help='First template PDB', required=True)
	parser.add_argument('-br', '--b_result', help='Second HHR target file result', required=True)
	parser.add_argument('-bc', '--b_chain', help='Second structure chain name', required=False, default="A")
	parser.add_argument('-bt', '--b_template', help='Second template PDB', required=False, default="A")
	parser.add_argument('-ts', '--template', help='Structure template', required=True)
	parser.add_argument('-tc', '--template_core', help='Core template chain name', required=True)
	parser.add_argument('-o', '--output', help='Output PDB file', required=False)
	args = parser.parse_args()
	main(args)
 