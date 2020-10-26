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
