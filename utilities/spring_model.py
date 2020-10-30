#! /usr/bin/env python3
import argparse
import os

from spring_package.Alignment import Alignment
from spring_package.Energy import Energy
from spring_package.Molecule import Molecule

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
    baseA = os.path.basename(fileA)
    baseB = os.path.basename(fileB)
    baseA = os.path.splitext(baseA)[0]
    baseB = os.path.splitext(baseB)[0]
    tmName = "temp/tmalign.%s.%s" % (baseA, baseB)
    os.system("build/TMalign %s %s -m %s.mat > %s.out" % (fileA, fileB, tmName, tmName))
    rotmat = list()
    with open("%s.mat" % tmName) as file:
        line = next(file)
        line = next(file)
        for i in range(3):
            line = next(file)
            rotmatLine = line[1:].split()
            rotmatLine = list(map(lambda x: float(x), rotmatLine))
            rotmatLine = [rotmatLine[1], rotmatLine[2], rotmatLine[3], rotmatLine[0]]
            rotmat.append(rotmatLine)
    with open("%s.out" % tmName) as file:
        for i in range(14):
            line = next(file)
        try:
            tmscore = float(line[9:17])
            line = next(file)
            tmscore = max(tmscore, float(line[9:17]))
        except:
            raise Exception("TMalign::Failed to retrieve TMscore.")
    molecule = Molecule(fileA)
    for atom in molecule.atoms:
        molecule.applyMatrix(atom, rotmat)
    return tmscore, molecule

def main(args):
    print("SPRING Model")
    print("Sequence A: %s" % args.a_result)
    print("Sequence B: %s" % args.b_result)
    print("Template C: %s" % args.c_template)
    os.system("mkdir -p temp")
    os.system("rm -f temp/*.*")
    buildModel(args.a_result, args.a_template, args.a_chain, "temp/monomerA.pdb")
    buildModel(args.b_result, args.b_template, args.b_chain, "temp/monomerB.pdb")
    interfaceEnergy = Energy()
    templateMolecule = Molecule(args.c_template)
    maxScore = -9999
    maxMolecule = None
    maxTemplate = None
    for biomolNumber in range(len(templateMolecule.rotmat.keys())):
        os.system("rm -f temp/template*.pdb")
        if biomolNumber == 0:
            bioMolecule = templateMolecule
        else:
            bioMolecule = templateMolecule.createUnit(biomolNumber)
        if len(bioMolecule.calpha.keys()) > 1 and args.c_chain in bioMolecule.calpha:
            for chainName in bioMolecule.calpha.keys():
                bioMolecule.saveChain(chainName, "temp/template%s.pdb" % chainName)
            coreTMscore, coreMolecule = TMalign("temp/monomerA.rebuilt.pdb", "temp/template%s.pdb" % args.c_chain)
            for chainName in bioMolecule.calpha.keys():
                if chainName != args.c_chain and len(bioMolecule.calpha[chainName]) > 0:
                    print("Evaluating chain %s..." % chainName)
                    try:
                        partnerTMscore, partnerMolecule = TMalign("temp/monomerB.rebuilt.pdb", "temp/template%s.pdb" % chainName)
                    except:
                        print("Warning: Failed TMalign [%s]." % chainName)
                        continue
                    TMscore = min(coreTMscore, partnerTMscore)
                    print("min-TMscore: %5.5f" % TMscore)
                    energy = -interfaceEnergy.get(coreMolecule, partnerMolecule)
                    print("Interaction: %5.5f" % energy)
                    springScore = TMscore * args.wtm + energy * args.wenergy
                    print("SpringScore: %5.5f" % springScore)
                    if springScore > maxScore:
                        maxScore = springScore
                        maxMolecule = partnerMolecule
                        maxTemplate = bioMolecule
            outputName = "%s_%1.2f.pdb" % (args.output, maxScore)
            coreMolecule.save(outputName, chainName="0")
            maxMolecule.save(outputName, chainName="1", append=True)
            maxTemplate.save(outputName, append=True)
    if maxMolecule is None:
        print ("Warning: Failed to determine model.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a 3D model from HH-search results.')
    parser.add_argument('-ar', '--a_result', help='First HHR target file result', required=True)
    parser.add_argument('-ac', '--a_chain', help='First template chain name', required=True)
    parser.add_argument('-at', '--a_template', help='First template PDB', required=True)
    parser.add_argument('-br', '--b_result', help='Second HHR target file result', required=True)
    parser.add_argument('-bc', '--b_chain', help='Second structure chain name', required=True)
    parser.add_argument('-bt', '--b_template', help='Second template PDB', required=True)
    parser.add_argument('-ct', '--c_template', help='Structure template', required=True)
    parser.add_argument('-cc', '--c_chain', help='Core template chain name', required=True)
    parser.add_argument('-o', '--output', help='Output model file', required=True)
    parser.add_argument('-wt', '--wtm', help='Weight TM-score', type=float, default=12.4, required=False)
    parser.add_argument('-we', '--wenergy', help='Weight Energy term', type=float, default=0.2, required=False)
    parser.add_argument('-tp', '--temp', help='Temporary directory', required=False, default="temp/")
    parser.add_argument('-st', '--show_template', help='Add template to model structure', type=bool, required=False, default=True)
    args = parser.parse_args()
    main(args)
 