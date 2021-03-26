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
from bromine_generation import bromine_comfiles_from_xyz
from xyzfile_generation import xyz_from_com
from xyz2mol import xyz2mol
from xyz2mol import read_xyz_file
import pytest

#files=glob.glob('neutral_comfiles/*')
#for i in files:
#    xyz_from_com(i)

xyz_from_com('neutral_comfiles\\S1_2021-02-26_wb97xd_631gd_opt.com')
#print(read_xyz_file('xyz_files/S1_2021-03-26.xyz'))

b#romine_comfiles_from_xyz('xyz_files\\S1_2021-03-26_2021-02-26.xyz')