import glob
import sys
import cirpy
from datetime import date
import rdkit
import pandas as pd
import glob
from VEHICLe_read import read_csv_to_dataframe
from VEHICLe_read import read_list_of_smiles
from VEHICLe_read import read_mol_from_smiles
import numpy as np
import rdkit
from com_files_from_smiles import VEHICLe_string_to_com
#from mol_translator.structure.structure_write import write_mol_toxyz
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
from anion_generation import anion_from_com
from xyzfile_generation import xyz_from_com
from bromine_generation import bromines_from_com
from xyz2mol import xyz2mol
from xyz2mol import read_xyz_file
from xyz2mol import AC2mol
from xyz2mol import xyz2AC
#from xyz2mol import x
import pytest


files=glob.glob('xyz_files\\*')

ffn=glob.glob('neutral_comfiles/*')
ffb=glob.glob('bromine_comfiles/*')

#for i in ffb:
 #   xyz_from_com(i)

for x in files:
    print(x)
    split = x.split('\\')
    name = split[1].split('_')
    xyz_read = read_xyz_file(x)

    mols = xyz2mol(xyz_read[0], xyz_read[2], xyz_read[1], use_huckel=True)
    for i in mols:

        try:
            rdkit.Chem.Draw.MolToFile(i, name[0] + name[2] + '.png', kekulize=True)
        except 'RDKit ERROR' :
            print(file, 'has incorrect bromine')
#for k in ffn:
#    bromines_from_com(k)

#for j in ffb:
 #   xyz_from_com(j)



#AC=xyz2AC(xyz_read[0], xyz_read[2], xyz_read[1])
#mol=AC2mol(AC[1], AC[0], xyz_read[0], xyz_read[1])
#mols=xyz2mol(xyz_read[0], xyz_read[2], xyz_read[1])



