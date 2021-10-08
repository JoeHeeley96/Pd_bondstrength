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
from VEHICLe_filters import nitrogen_only_filter
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

import pickle
from mol_translator.aemol import aemol
import pybel as pyb
from mol_translator import aemol

VEHICLe_Nonly=pd.read_csv('Nitrogen_only_VEHICLe.csv')
nonly_fulldata=pd.read_csv('Nonly_fulldata.csv')
comfiles = glob.glob('neutral_comfiles/*') + glob.glob('anion_comfiles/*') + glob.glob('bromine_comfiles/*')
logfiles= glob.glob('logfiles/*')
xyzfiles = glob.glob('xyzfiles/*')
anion_logfiles = glob.glob('logfiles/*anion*')
bromine_logfiles = glob.glob('logfiles/*bromine*')
neutral_logfiles = list(set(logfiles) - set(anion_logfiles) - set(bromine_logfiles))
sdf_files = glob.glob('aemols/*')
training_set_first200 = neutral_logfiles[0:200]
test_set_last55 = neutral_logfiles[200:]

S1035_log = ['logfiles\\S1035_2anion_2021-04-21_wb97xd-631gd_opt.log']
S1035_xyz = ['xyzfiles\\S1035_2anion_2021-09-30_wb97xd-631gd_opt.xyz']

#aemols = logfile_to_aemol(training_set)
#inputs = write_imp_input(aemols, nonly_fulldata)

def output_structure_check(xyzfiles, logfiles):
    '''
    requires inputs in alphabetical order
    '''

    for i,j in zip(xyzfiles, logfiles):
        in_file = i.split('\\')
        in_split = in_file[1].split('-')

        out_file = j.split('\\')
        out_split = out_file[1].split('-')
        if in_split[0] != out_split[0]:
            raise NameError('Files dont match')

        input_aemol = aemol(in_file[1].split('_')[0])
        input_aemol.from_file(i, ftype = 'xyz')

        output_aemol = aemol(in_file[1].split('_')[0])
        output_aemol.from_file(j, ftype='log')

        input_structure = output_aemol.structure['conn']
        output_structure = input_aemol.structure['conn']
        #print(in_split[0])
        print(input_structure)
        print(output_structure)

        for x,y in zip(input_structure, output_structure):
            print(np.count_nonzero(x == y))


output_structure_check(xyzfiles[:1], logfiles[:1])
