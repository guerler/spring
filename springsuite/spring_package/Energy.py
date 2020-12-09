import math

NTYPE = 21
NDIST = 20
NSCALE = 2.0


class Energy:
    def __init__(self):
        self.dfire = list()
        with open("spring_package/dfire/dfire.txt") as file:
            for line in file:
                self.dfire.append(float(line))

    def get(self, moleculeA, moleculeB):
        result = 0
        chainA = list(moleculeA.calpha.keys())[0]
        chainB = list(moleculeB.calpha.keys())[0]
        for i in moleculeA.calpha[chainA]:
            atomA = moleculeA.calpha[chainA][i]
            indexA = self.toResCode(atomA["residue"])
            for j in moleculeB.calpha[chainB]:
                atomB = moleculeB.calpha[chainB][j]
                indexB = self.toResCode(atomB["residue"])
                dist2 = ((atomA["x"] - atomB["x"]) ** 2 +
                         (atomA["y"] - atomB["y"]) ** 2 +
                         (atomA["z"] - atomB["z"]) ** 2)
                dist = int((math.sqrt(dist2) * NSCALE))
                if dist < NDIST:
                    index = indexA * NTYPE * NDIST + indexB * NDIST + dist
                    result = result + self.dfire[index]
        return result

    def toResCode(self, seq):
        code = dict(ALA=0, CYS=1, ASP=2, GLU=3, PHE=4, GLY=5, HIS=6, ILE=7, LYS=8, LEU=9, MET=10,
                    ASN=11, PRO=12, GLN=13, ARG=14, SER=15, THR=16, VAL=17, TRP=18, TYR=19)
        return code[seq] if seq in code else 20
