import pandas as pd
import glob
from datetime import date
from bromine_generation import bromines_from_smiles
from xyzfile_generation import xyz_from_com
from xyz2mol import read_xyz_file
from xyz2mol import xyz2mol
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from anion_generation import anions_from_smiles
from rdkit.Chem import AllChem
from VEHICLe_filters import nitrogen_only_filter
from VEHICLe_filters import NCN_filter
from plots import plot_activation_map
import numpy as np
from rdkit.Chem.FragmentMatcher import FragmentMatcher
from rdkit.Chem import rdRGroupDecomposition
from print_file import print_smiles_to_gridimage
from plots import plot_3Dactivation_map
from plots import plot_activation_map

VEHICLe=pd.read_csv('VEHICLe.csv')
VEHICLe_first30=VEHICLe[1:31]
smniles = VEHICLe['Smiles']
interesting_structures=pd.read_csv('interesting_structures.csv')
nitrogen_only=pd.read_csv('Nitrogen_only.csv')
#NCN=pd.read_csv('NCN.csv')
NOnly_data=pd.read_csv('VEHICLe_NOnly_fulldata.csv')

plot_activation_map(NOnly_data, '2D_activation_map')
#plot_2Dactivation_map(NOnly_data, '2D_activation_eaoverAcidity', 'n')

print_smiles_to_gridimage(VEHICLe_first30, 'VEHICLe_first30.png', 30)