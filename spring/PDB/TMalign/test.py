import sys
import os
currentPath = os.path.dirname(os.path.abspath(__file__))
pdbPath = os.path.dirname(currentPath)
print(currentPath)
print(pdbPath)

sys.path.append(pdbPath)

from pdbchains import PDBchains
from tmalign import _tmalign

pdb1 = PDBchains()
pdb1.Read('./PDB1.pdb', ca_only = True)
seq1 = pdb1.sequence[0]

pdb2 = PDBchains()
pdb2.Read("./PDB1.pdb", ca_only = True)
seq2 = pdb2.sequence[0]
x = _tmalign(seq1,pdb1.coord[0],seq2,pdb2.coord[0])

print(x)
