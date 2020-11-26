from os.path import isfile


def validateIdentifier(identifier):
    if len(identifier) != 6 or identifier[4:5] != "_":
        raise Exception("Invalid list entry (`PDB_CHAIN`): %s." % identifier)


def getId(identifier):
    identifier = identifier.strip()
    validateIdentifier(identifier)
    return identifier[:4].upper() + identifier[4:6]


def getChain(identifier):
    validateIdentifier(identifier)
    pdbChain = identifier[5:6]
    return pdbChain


def getName(identifier):
    validateIdentifier(identifier)
    pdb = identifier[:4].lower()
    return pdb


def getCrossReference(crossReferenceFile):
    crossReference = dict()
    with open(crossReferenceFile) as file:
        for line in file:
            columns = line.split()
            if len(columns) < 2:
                raise Exception("Invalid Cross Reference Entry %s." % line)
            core = columns[0]
            partner = columns[-1]
            if core not in crossReference:
                crossReference[core] = []
            crossReference[core].append(partner)
    return crossReference


def getTemplates(hhrFile, minScore=10):
    result = dict()
    topTemplate = None
    if isfile(hhrFile):
        with open(hhrFile) as file:
            for index, line in enumerate(file):
                if index > 8:
                    if not line.strip():
                        break
                    templateId = line[4:10]
                    templateScore = float(line[57:63])
                    if templateScore > minScore:
                        if topTemplate is None:
                            topTemplate = templateId
                        result[templateId] = templateScore
    return topTemplate, result
