#! /usr/bin/env python3
import argparse
from os import system
from os.path import basename, isfile, splitext

from spring_package.Alignment import Alignment
from spring_package.DBKit import DBKit
from spring_package.Energy import Energy
from spring_package.Molecule import Molecule
from spring_package.Utilities import getChain, getCrossReference, getName, getTemplates


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
    tmName = "temp/tmalign"
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


def getFrameworks(aTemplates, bTemplates, crossReference, minScore, maxTries):
    templateHits = list()
    for aTemplate in aTemplates:
        if aTemplate in crossReference:
            partners = crossReference[aTemplate]["partners"]
            templates = crossReference[aTemplate]["templates"]
            for pIndex in range(len(partners)):
                pTemplate = partners[pIndex]
                templatePair = templates[pIndex]
                if pTemplate in bTemplates:
                    minZ = min(aTemplates[aTemplate], bTemplates[pTemplate])
                    templateHits.append(dict(templatePair=templatePair, score=minZ))
    templateList = sorted(templateHits, key=lambda item: item["score"], reverse=True)
    print("Found %d templates." % len(templateList))
    for templateHit in templateList:
        if templateHit["score"] < minScore or maxTries == 0:
            break
        maxTries = maxTries - 1
        yield templateHit["templatePair"]


def main(args):
    print("SPRING - Complex Model Creation")
    print("Sequence A: %s" % args.a_hhr)
    print("Sequence B: %s" % args.b_hhr)
    aTop, aTemplates = getTemplates(args.a_hhr)
    bTop, bTemplates = getTemplates(args.b_hhr)
    system("mkdir -p temp")
    system("rm -f temp/*.*")
    outputName = args.output
    pdbDatabase = DBKit(args.index, args.database)
    crossReference = getCrossReference(args.cross)
    interfaceEnergy = Energy()
    buildModel(args.a_hhr, aTop, pdbDatabase, "temp/monomerA.pdb")
    buildModel(args.b_hhr, bTop, pdbDatabase, "temp/monomerB.pdb")
    maxScore = -9999
    maxInfo = None
    minScore = float(args.minscore)
    maxTries = int(args.maxtries)
    for [aTemplate, bTemplate] in getFrameworks(aTemplates, bTemplates, crossReference, minScore=minScore, maxTries=maxTries):
        print("Evaluating Complex Template: %s." % aTemplate)
        templateFile = "temp/template.pdb"
        createPDB(aTemplate, pdbDatabase, templateFile)
        templateMolecule = Molecule(templateFile)
        aTemplateChain = getChain(aTemplate)
        bTemplateChain = getChain(bTemplate)
        print("Evaluating chain %s and %s..." % (aTemplate, bTemplate))
        for biomolNumber in range(len(templateMolecule.biomol.keys())):
            system("rm -f temp/template_*.pdb")
            if biomolNumber == 0:
                bioMolecule = templateMolecule
            else:
                bioMolecule = templateMolecule.createUnit(biomolNumber)
            if len(bioMolecule.calpha.keys()) > 1 and aTemplateChain in bioMolecule.calpha and bTemplateChain in bioMolecule.calpha:
                print("Evaluating biomolecule %i..." % biomolNumber)
                bioMolecule.saveChain(aTemplateChain, "temp/template_0.pdb")
                bioMolecule.saveChain(bTemplateChain, "temp/template_1.pdb")
                try:
                    coreTMscore, coreMolecule = TMalign("temp/monomerA.rebuilt.pdb", "temp/template_0.pdb")
                    partnerTMscore, partnerMolecule = TMalign("temp/monomerB.rebuilt.pdb", "temp/template_1.pdb")
                except Exception:
                    print("Warning: Failed TMalign [%s]." % bTemplateChain)
                    continue
                TMscore = min(coreTMscore, partnerTMscore)
                print("  min-TMscore: %5.5f" % TMscore)
                energy = -interfaceEnergy.get(coreMolecule, partnerMolecule)
                print("  Interaction: %5.5f" % energy)
                springScore = TMscore * args.wtm + energy * args.wenergy
                print("  SpringScore: %5.5f" % springScore)
                if springScore > maxScore:
                    maxScore = springScore
                    maxInfo = "%s\t %5.5f\t %5.5f\n" % (outputName, TMscore, energy)
                    coreMolecule.save(outputName, chainName="0")
                    partnerMolecule.save(outputName, chainName="1", append=True)
                    if args.showtemplate == "true":
                        bioMolecule.save(outputName, append=True)
    if maxInfo is not None:
        print("Completed.")
        print("SpringScore: %5.5f" % maxScore)
        print("Result stored to %s" % outputName)
        logExists = isfile(args.log)
        logFile = open(args.log, "a+")
        if not logExists:
            logFile.write("# Description: Name, TMscore, Energy\n")
        logFile.write(maxInfo)
        logFile.close()
    else:
        print("Warning: Failed to determine model.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a 3D model from HH-search results.')
    parser.add_argument('-a', '--a_hhr', help='First HHR target file result', required=True)
    parser.add_argument('-b', '--b_hhr', help='Second HHR target file result', required=True)
    parser.add_argument('-i', '--index', help='PDB Database Index file (ffindex)', required=True)
    parser.add_argument('-d', '--database', help='PDB Database files (ffdata)', required=True)
    parser.add_argument('-c', '--cross', help='PDB Cross Reference', required=True)
    parser.add_argument('-o', '--output', help='Output model file', required=True)
    parser.add_argument('-g', '--log', help='Log file', required=True)
    parser.add_argument('-wt', '--wtm', help='Weight TM-score', type=float, default=1.0, required=False)
    parser.add_argument('-we', '--wenergy', help='Weight Energy term', type=float, default=0.0, required=False)
    parser.add_argument('-ms', '--minscore', help='Minimum min-Z score threshold', type=float, default=10.0, required=False)
    parser.add_argument('-mt', '--maxtries', help='Maximum number of templates', type=int, default=50, required=False)
    parser.add_argument('-sr', '--showtemplate', help='Add reference template to model structure', required=False, default="true")
    args = parser.parse_args()
    main(args)
