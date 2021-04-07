import pandas as pd
import glob
from bromine_generation import bromines_from_com
from com_files_from_smiles import VEHICLe_string_to_com
from xyzfile_generation import xyz_from_com
from xyz2mol import read_xyz_file
from xyz2mol import xyz2mol
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from anion_generation import anion_from_com
from rdkit.Chem import AllChem
from com_files_from_smiles import write_indexfile

VEHICLe=pd.read_csv('VEHICLe.csv')

VEHICLe_1=VEHICLe[:31]

write_indexfile(VEHICLe_1)

#mol=Chem.MolFromSmiles('c1cc[nH]c1')
#molH=Chem.AddHs(mol)
#molc=AllChem.Compute2DCoords(molH)
#a=molH.GetAtoms()

#for c in molH.GetConformers():
 #   index=[]
  #  atomicnum = []
   # coords = list(c.GetPositions())

    #for j, atoms in enumerate(a):
     #   b=atoms.GetAtomicNum()
      #  d = atoms.GetIdx()
       # index.append(d)
        #atomicnum.append(b)

    #zip=zip(atomicnum, index, coords)

    #for u, t, i in zip:
     #   print(u, t, i)


        #if b == 1:
         #   x=atoms.GetNeighbors()
          #  for i, k in enumerate(x):
           #     print(k, k.GetIdx())
