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
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from mol_translator.imp_converter import dataframe_write as dfw
from tqdm import tqdm
from impression_input import logfile_to_aemol
from impression_input import xyzfile_to_aemol
from mol_translator.properties.structure.bond_angle import get_bond_angles
from mol_translator.properties.structure.dihedral_angle import get_dihedral_angle
from workflow import logfile_analysis_workflow
import pickle
from property_calculate import calculate_properties
from mol_translator.aemol import aemol
import pybel as pyb
from logfile_read import energy_readdata
from mol_translator import aemol
import sklearn
from property_calculate import calculate_relative_properties
from workflow import impression_input_workflow
from VEHICLe_filters import dj2_filter
from sampling import get_tm_df
from calculation_check import output_structure_check
from calculation_check import opt_check
from random import randint
from rdkit import DataStructs
from sampling import get_fps_ids
from workflow import comfile_generation_workflow
from workflow import logfile_fingerprint_sampling_workflow
from workflow import vehicle_fingerprint_sampling_workflow
from calculation_check import comfile_check
from neutral_generation import basehet_from_smiles
from sampling import random_sample
from VEHICLe_read import read_mol_from_smiles
from plots import plot_binary_activation_map

VEHICLe = pd.read_csv('VEHICLe.csv')

dj1 = pd.read_csv('dj1/dj1.csv')
dj1_clean = pd.read_csv('dj1/dj1_clean.csv')
dj1_clean_fulldata = pd.read_csv('dj1/data/dj1_clean_fulldata.csv')
dj1logfiles = glob.glob('dj1/logfiles/*')
dj1anion_logfiles = glob.glob('dj1/logfiles/*anion*')
dj1bromine_logfiles = glob.glob('dj1/logfiles/*bromine*')
dj1neutral_logfiles = list(set(dj1logfiles) - set(dj1anion_logfiles) - set(dj1bromine_logfiles))
independent_testset = pd.read_csv('imp_inputs/dj1clean_random_testset_for_sampling_atom_df.csv')
trainingset_for_sampling = pd.read_csv('dj1/dj1_training_for_sampling.csv')

dj2 = pd.read_csv('dj2/dj2.csv')
dj2_comfiles = glob.glob('dj2/anion_comfiles/*') + glob.glob('dj2/bromine_comfiles/*') + glob.glob('dj2/neutral_comfiles/*')
dj2_fpsampled1000 = pd.read_csv('dj2/dj2_fpsampled1000.csv')
dj2_randsample1000 = pd.read_csv('dj2/dj2_randsample1000.csv')
dj2_logfiles = glob.glob('dj2/logfiles/*')

dj3 = pd.read_csv('dj3/dj3.csv')
dj3_comfiles = glob.glob('dj3/neutral_comfiles/*') + glob.glob('dj3/anion_comfiles/*') + glob.glob('dj3/bromine_comfiles/*')
dj3_logfiles = glob.glob('dj3/lofiles/*')

xyzfiles = glob.glob('xyzfiles/*')


'''structures= logfile_analysis_workflow('dj2/logfiles', outname=None, calc_check=True, calc_data=False, plot_map=False, write=False)

with open('dj2_checked_files.txt', 'a+') as f:
    for i in structures[0] + structures[1]:
        for line in f:
            if line.endswith(i):
                break
            else:
                print(date.today(), i, file=f)
'''

'''anion_diff = []
bromine_diff = []

regid = dj1_clean_fulldata['Regid']

for i in tqdm(regid):
    df = dj1_fulldata[dj1_fulldata['Regid'] == i]
    df_nan = df.fillna(0)
    anions = []
    bromines = []

    for j in df_nan.columns:
        if 'anion' in j:
            anions.append([x for x in df_nan[j].values if df_nan[j].values != 0])

        elif 'bromine' in j:
            bromines.append([x for x in df_nan[j].values if df_nan[j].values != 0])


    anion_diff = anion_diff + (list(np.array(sorted(anions)[-1]) - np.array(sorted(anions)[-2])))
    bromine_diff = bromine_diff + (list(np.array(sorted(bromines, reverse=True)[0]) - np.array(sorted(bromines, reverse=True)[1])))

print(anion_diff)
print(bromine_diff)
print('average lowest anions:', sum(anion_diff)/len(anion_diff))
print('average highest bromines:', sum(bromine_diff)/len(bromine_diff))'''

successful_activations = ['S1080', 'S1875', 'S1863', 'S5857', 'S6018', 'S5854', 'S5703', 'S1869']

fulldata_filter(VEHICLe_dataframe=VEHICLe, calculate_properties_dataframe=dj1_clean_fulldata, outname='dj1/images/dj1_FDfilter_ActMap_A0.91-2.5_Ea0.91-1.2_TopRight',
                a_upper=2.5, a_lower=0.91, ea_upper=1.2, ea_lower=0.91)
