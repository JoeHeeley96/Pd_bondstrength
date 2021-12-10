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

VEHICLe = pd.read_csv('VEHICLe.csv')
dj1 = pd.read_csv('dj1/dj1.csv')
dj1_clean = pd.read_csv('dj1/dj1_clean.csv')
dj2 = pd.read_csv('dj2/dj2.csv')
dj1_fulldata = pd.read_csv('dj1/data/dj1_clean_fulldata.csv')
comfiles = glob.glob('neutral_comfiles/*') + glob.glob('anion_comfiles/*') + glob.glob('bromine_comfiles/*')
dj1logfiles = glob.glob('dj1/logfiles/*')
xyzfiles = glob.glob('xyzfiles/*')
dj1anion_logfiles = glob.glob('dj1/logfiles/*anion*')
dj1bromine_logfiles = glob.glob('dj1/logfiles/*bromine*')
aemols = glob.glob('aemols/*')
dj1neutral_logfiles = list(set(dj1logfiles) - set(dj1anion_logfiles) - set(dj1bromine_logfiles))
anion_xyzfiles = glob.glob('xyzfiles/*anion*')
dj1_clean_fulldata = pd.read_csv('dj1/data/dj1_clean_fulldata.csv')

independent_testset = pd.read_csv('imp_inputs/dj1clean_random_testset_for_sampling_atom_df.csv')
#trainingset_for_sampling = pd.read_csv('dj1_training_for_sampling.csv')

dj2_fpsampled1000 = pd.read_csv('dj2/dj2_fpsampled1000.csv')
dj2_randsample1000 = pd.read_csv('dj2/dj2_randsample1000.csv')

dj1_fullatomdf = pd.read_csv('dj1clean_full_atom_df.csv')
dj1_fullpairdf = pd.read_csv('dj1clean_full_pair_df.csv')


#vehicle_fingerprint_sampling_workflow(dj2, 'dj2_fpsampled', num=1000, write=True)

dj1_training = pd.DataFrame()



#regid = dj2_fpsampled1000['Regid']
#logfile_analysis_workflow('logfiles', None, calc_check=True, calc_data=False,
#                          plot_map=False, write=False)



#comfile_generation_workflow(dj2_fpsampled1000)
#comfile_check(comfiles=comfiles, sampled_structures=regid)

'''structures= logfile_analysis_workflow('logfiles', outname=None, calc_check=True, calc_data=False, plot_map=False, write=False)

with open('dj2_checked_files.txt', 'r+') as f:
    for i in structures[1]:
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


sample_size = [10, 20, 40, 60, 80, 100, 150, 200]

'''for i in sample_size:

    sampled_atom = pd.DataFrame()
    sampled_pair = pd.DataFrame()

    fpsamp_fulldata = pd.DataFrame()
    fpsamp = vehicle_fingerprint_sampling_workflow(vehicle_dataframe=trainingset_for_sampling, outname='dj1_clean_fpsampled',
                                                   num=i, write=False)

    for j in fpsamp['Regid']:
        sampled_atom = sampled_atom.append(dj1_fullatomdf[dj1_fullatomdf['molecule_name'] == j])
        sampled_pair = sampled_pair.append(dj1_fullpairdf[dj1_fullpairdf['molecule_name'] == j])

    with open('dj1clean_fpsamp' + str(i) + '_atom_df.csv', 'w') as f:
        print(sampled_atom.to_csv(sep=','), file=f)

    with open('dj1clean_fpsamp' + str(i) + '_pair_df.csv', 'w') as f:
        print(sampled_pair.to_csv(sep=','), file=f)'''

comfile_generation_workflow(dj2_randsample1000)

