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
from workflow import comfile_workflow

VEHICLe=pd.read_csv('VEHICLe.csv')

nitrogen_only=pd.read_csv('Nitrogen_only.csv')

S1773=nitrogen_only[nitrogen_only.Regid == 'S1773']

S123_S13281=nitrogen_only[38:40]

S13282=nitrogen_only[40:41]

bromines_from_smiles(S13282)



