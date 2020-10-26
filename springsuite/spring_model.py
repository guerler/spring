#! /usr/bin/env python3
import argparse
import os
from spring_package.Molecule import Molecule
from spring_package.Alignment import Alignment

def buildModel(resultFile, templateFile, chainName, outputName):
	template = Molecule(templateFile)
	if chainName not in template.calpha:
		raise Exception("Chain not found in template [%s]" % chainName)
	chain = template.calpha[chainName]
	alignment = Alignment(resultFile)
	alignment.createModel(chain)
	template.saveChain(chainName, outputName)
	os.system("./build/pulchra %s" % outputName)

def TMalign(fileA, fileB):
	os.system("build/TMalign %s %s -m temp/tmalign.mat > temp/tmalign.out" % (fileA, fileB))
	rotmat = list()
	with open("temp/tmalign.mat") as file:
		line = next(file)
		line = next(file)
		for i in range(3):
			line = next(file)
			rotmatLine = line[1:].split()
			rotmatLine = list(map(lambda x: float(x), rotmatLine))
			rotmatLine = [rotmatLine[1], rotmatLine[2], rotmatLine[3], rotmatLine[0]]
			rotmat.append(rotmatLine)
	with open("temp/tmalign.out") as file:
		for i in range(14):
			line = next(file)
		try:
			tmscore = float(line[9:17])
		except:
			raise Exception("TMalign::Failed to retrieve TMscore.")
	molecule = Molecule(fileA)
	for atom in molecule.atoms:
		molecule.applyMatrix(atom, rotmat)
	return tmscore, molecule

def main(args):
	os.system("mkdir -p temp")
	buildModel(args.a_result, args.a_template, args.a_chain, "temp/modelA.pdb")
	buildModel(args.b_result, args.b_template, args.b_chain, "temp/modelB.pdb")
	templateMolecule = Molecule(args.template)
	bioMolecule = templateMolecule.createUnit()
	for key in bioMolecule.calpha.keys():
		bioMolecule.saveChain(key, "temp/template%s.pdb" % key)
	coreTMscore, coreMolecule = TMalign("temp/modelA.rebuilt.pdb", "temp/template%s.pdb" % args.template_core)
	for chainName in bioMolecule.calpha.keys():
		if chainName != args.template_core:
			print("Evaluating chain %s..." % chainName) 
			partnerTMscore, partnerMolecule = TMalign("temp/modelB.rebuilt.pdb", "temp/template%s.pdb" % chainName)
			minTMscore = min(coreTMscore, partnerTMscore)
			print("min-TMscore: %2.5f" % minTMscore)
			partnerMolecule.save("temp/tmalign%s.pdb" % chainName)

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
 