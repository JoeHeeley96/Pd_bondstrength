import pandas as pd
import glob
from datetime import date
from bromine_generation import bromines_from_smiles
from xyzfile_generation import xyz_from_com
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from anion_generation import anions_from_smiles
from rdkit.Chem import AllChem
from VEHICLe_filters import NCN_filter
from plots import plot_activation_map
import numpy as np
from rdkit.Chem.FragmentMatcher import FragmentMatcher
from rdkit.Chem import rdRGroupDecomposition
from print_file import print_smiles_to_gridimage
from plots import plot_3Dactivation_map
from plots import plot_activation_map
from VEHICLe_filters import fulldata_filter
import glob
import copy
import os
import re
from sklearn import metrics
import sys
from impression_input import write_imp_input
import math
from mol_translator.imp_converter import dataframe_write as dfw
from tqdm import tqdm
from impression_input import logfile_to_aemol
from impression_input import xyzfile_to_aemol
from mol_translator.properties.structure.bond_angle import get_bond_angles
from mol_translator.properties.structure.dihedral_angle import get_dihedral_angle
from workflow import property_calculate_workflow
import pickle
from mol_translator.aemol import aemol
import pybel as pyb
from logfile_read import energy_readdata
from mol_translator import aemol
from property_calculate import calculate_relative_properties
from workflow import impression_input_workflow
from VEHICLe_filters import dj2_filter


VEHICLe = pd.read_csv('VEHICLe.csv')
dj1=pd.read_csv('dj1.csv')
dj2=pd.read_csv('dj2.csv')
dj1_fulldata=pd.read_csv('dj1_fulldata.csv')
comfiles = glob.glob('neutral_comfiles/*') + glob.glob('anion_comfiles/*') + glob.glob('bromine_comfiles/*')
logfiles= glob.glob('logfiles/*')
xyzfiles = glob.glob('xyzfiles/*')
anion_logfiles = glob.glob('logfiles/*anion*')
bromine_logfiles = glob.glob('logfiles/*bromine*')
neutral_logfiles = list(set(logfiles) - set(anion_logfiles) - set(bromine_logfiles))

#training_set_first200 = neutral_logfiles[0:200]
#test_set_last55 = neutral_logfiles[200:]

#training_set_every2nd = sorted(neutral_logfiles, key=str)[0::2]
#test_set_every2nd = sorted(neutral_logfiles, key=str)[1::2]

#impression_input_workflow(training_set_every2nd, dj1_fulldata, 'dj1_every2nd_train')
#impression_input_workflow(test_set_every2nd, dj1_fulldata, 'dj1_every2nd_test')

print(len(dj2))