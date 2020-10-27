class Alignment:
	def __init__(self, fileName):
		self.queryName = None
		self.queryStart = list()
		self.queryAlignment = list()
		self.templateName = None
		self.templateStart = list()
		self.templateAlignment = list()
		self.templateDssp = list()
		self.readFile(fileName)

	def readFile(self, fileName):
		with open(fileName) as file:
			for index, line in enumerate(file):
				cols = line.split()
				if len(cols) > 1 and cols[0] == "Query":
					self.queryName = cols[1].split()[0][0:14]
				if len(cols) > 1 and cols[0].startswith(">"):
					self.templateName = cols[0][1:]
				if self.queryName and self.templateName:
					if len(cols) > 2:
						if cols[0] == "Q" and cols[1] == self.queryName:
							self.queryStart.append(self.getStart(cols[2]))
							self.queryAlignment.append(cols[3])
						if cols[0] == "T" and cols[1] == self.templateName:
							self.templateStart.append(self.getStart(cols[2]))
							self.templateAlignment.append(cols[3])
						if cols[0] == "T" and cols[1] == "ss_dssp":
							self.templateDssp.append(cols[2])
				if len(cols) > 1 and cols[0] == "No" and cols[1] == "2":
					break

	def createModel(self, templateChain):
		previousResidue = dict()
		for residueNumber in templateChain:
			previousResidue[residueNumber] = templateChain[residueNumber]["residue"]
			templateChain[residueNumber]["residue"] = None
		tcount = 0
		for i in range(len(self.templateStart)):
			templateStart = self.templateStart[i]
			templateSequence = self.templateAlignment[i]
			templateDssp = self.templateDssp[i]
			queryStart = self.queryStart[i]
			querySequence = self.queryAlignment[i]
			n = len(querySequence)
			qcount = 0
			for j in range(n):
				qs = querySequence[j]
				ts = templateSequence[j]
				dssp = templateDssp[j]
				if dssp != "-" and qs != "-" and ts != "-":
					residueNumber = tcount
					if residueNumber in templateChain:
						pr = previousResidue[residueNumber]
						if pr != self.toThreeAmino(ts):
							print ("Warning: Ignoring mismatching residue [%s != %s]." % (pr, self.toThreeAmino(ts)))
						templateChain[residueNumber]["residue"] = self.toThreeAmino(qs)
						templateChain[residueNumber]["residueNumber"] = queryStart + qcount
					else:
						print ("Warning: Skipping missing residue [%s]." % residueNumber)
				if qs != "-":
					qcount = qcount + 1
				if dssp != "-":
					tcount = tcount + 1

	def getStart(self, x):
		try:
			return int(x)
		except:
			raise Exception("Invalid start index in alignment [%s]." % x)

	def toThreeAmino (self, seq):
		code = dict(G="GLY", A="ALA", V="VAL", L="LEU", I="ILE", M="MET", F="PHE", P="PRO", Y="TYR", W="TRP",
					K="LYS", S="SER", C="CYS", N="ASN", Q="GLN", H="HIS", T="THR", E="GLU", D="ASP", R="ARG")
		return code[seq] if seq in code else "XXX"

	def toSingleAmino (self, seq):
		code = dict(GLY="G", ALA="A", VAL="V", LEU="L", ILE="I", MET="M", PHE="F", PRO="P", TYR="Y", TRP="W",
					LYS="K", SER="S", CYS="C", ASN="N", GLN="Q", HIS="H", THR="T", GLU="E", ASP="D", ARG="R")
		return code[seq] if seq in code else "X"
