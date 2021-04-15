import pandas as pd
import glob
from datetime import date
from bromine_generation import bromines_from_com
from com_files_from_smiles import VEHICLe_string_to_com
from xyzfile_generation import xyz_from_com
from xyz2mol import read_xyz_file
from xyz2mol import xyz2mol
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from anion_generation import anions_from_smiles
from rdkit.Chem import AllChem
from com_files_from_smiles import write_indexfile

VEHICLe=pd.read_csv('VEHICLe.csv')

VEHICLe_1=VEHICLe[:2]

anions_from_smiles(VEHICLe_1)




