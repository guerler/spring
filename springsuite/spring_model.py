#! /usr/bin/env python3
import argparse
from os import system
from os.path import basename, splitext

from spring_package.Alignment import Alignment, getTemplates
from spring_package.DBKit import DBKit
from spring_package.Energy import Energy
from spring_package.Molecule import Molecule
from spring_package.Utilities import getChain, getCrossReference, getName


def createPDB(identifier, pdbDatabase, outputName):
    pdb = getName(identifier)
    pdbDatabaseId = "%s.pdb" % pdb
    pdbDatabase.createFile(pdbDatabaseId, outputName)


def buildModel(resultFile, identifier, pdbDatabase, outputName):
    createPDB(identifier, pdbDatabase, outputName)
    template = Molecule(outputName)
    pdbChain = getChain(identifier)
    if pdbChain not in template.calpha:
        raise Exception("Chain not found in template [%s]" % pdbChain)
    chain = template.calpha[pdbChain]
    alignment = Alignment(resultFile)
    alignment.createModel(chain)
    template.saveChain(pdbChain, outputName)
    system("./build/pulchra %s" % outputName)


def TMalign(fileA, fileB):
    baseA = basename(fileA)
    baseB = basename(fileB)
    baseA = splitext(baseA)[0]
    baseB = splitext(baseB)[0]
    tmName = "temp/tmalign.%s.%s" % (baseA, baseB)
    system("build/TMalign %s %s -m %s.mat > %s.out" % (fileA, fileB, tmName, tmName))
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
        except Exception:
            raise Exception("TMalign::Failed to retrieve TMscore.")
    molecule = Molecule(fileA)
    for atom in molecule.atoms:
        molecule.applyMatrix(atom, rotmat)
    return tmscore, molecule


def getFrameworks(aTemplates, bTemplates, crossReference, minScore=10, maxTries=5):
    for aTemplate in aTemplates:
        if aTemplate in crossReference:
            partners = crossReference[aTemplate]
            for pTemplate in partners:
                if pTemplate in bTemplates:
                    minZ = min(aTemplates[aTemplate], bTemplates[pTemplate])
                    if minZ > minScore and maxTries > 0:
                        maxTries = maxTries - 1
                        yield aTemplate


def main(args):
    print("SPRING Model")
    print("Sequence A: %s" % args.a_hhr)
    print("Sequence B: %s" % args.b_hhr)
    aTop, aTemplates = getTemplates(args.a_hhr)
    bTop, bTemplates = getTemplates(args.b_hhr)
    system("mkdir -p temp")
    system("rm -f temp/*.*")
    pdbDatabase = DBKit(args.index, args.database)
    crossReference = getCrossReference(args.cross)
    interfaceEnergy = Energy()
    buildModel(args.a_hhr, aTop, pdbDatabase, "temp/monomerA.pdb")
    buildModel(args.b_hhr, bTop, pdbDatabase, "temp/monomerB.pdb")
    maxScore = -9999
    maxMolecule = None
    maxTemplate = None
    for aTemplate in getFrameworks(aTemplates, bTemplates, crossReference):
        templateFile = "temp/template.pdb"
        templateChain = getChain(aTemplate)
        createPDB(aTemplate, pdbDatabase, templateFile)
        templateMolecule = Molecule(templateFile)
        for biomolNumber in range(len(templateMolecule.biomol.keys())):
            system("rm -f temp/template*.pdb")
            if biomolNumber == 0:
                bioMolecule = templateMolecule
            else:
                bioMolecule = templateMolecule.createUnit(biomolNumber)
            if len(bioMolecule.calpha.keys()) > 1 and templateChain in bioMolecule.calpha:
                for chainName in bioMolecule.calpha.keys():
                    bioMolecule.saveChain(chainName, "temp/template%s.pdb" % chainName)
                coreTMscore, coreMolecule = TMalign("temp/monomerA.rebuilt.pdb", "temp/template%s.pdb" % templateChain)
                for chainName in bioMolecule.calpha.keys():
                    if chainName != templateChain and len(bioMolecule.calpha[chainName]) > 0:
                        print("Evaluating chain %s..." % chainName)
                        try:
                            partnerTMscore, partnerMolecule = TMalign("temp/monomerB.rebuilt.pdb", "temp/template%s.pdb" % chainName)
                        except Exception:
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
        print("Warning: Failed to determine model.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a 3D model from HH-search results.')
    parser.add_argument('-a', '--a_hhr', help='First HHR target file result', required=True)
    parser.add_argument('-b', '--b_hhr', help='Second HHR target file result', required=True)
    parser.add_argument('-i', '--index', help='PDB Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='PDB Database files (ffdata)', required=True)
    parser.add_argument('-c', '--cross', help='PDB Cross Reference', required=True)
    parser.add_argument('-o', '--output', help='Output model file', required=True)
    parser.add_argument('-t', '--temp', help='Temporary directory', required=False, default="temp")
    parser.add_argument('-wt', '--wtm', help='Weight TM-score', type=float, default=1.0, required=False)
    parser.add_argument('-we', '--wenergy', help='Weight Energy term', type=float, default=0.0, required=False)
    parser.add_argument('-sr', '--show_reference', help='Add reference template to model structure', type=bool, required=False, default=True)
    args = parser.parse_args()
    main(args)
