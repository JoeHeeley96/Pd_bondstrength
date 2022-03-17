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
from mol_translator.structure.rdkit_converter import aemol_to_rdmol
import sklearn
from property_calculate import calculate_relative_properties
from workflow import impression_input_workflow
from VEHICLe_filters import dj2_filter
from sampling import get_tm_df
from sampling import every2nd_sampling
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
from VEHICLe_filters import selectivity_filter
from property_calculate import find_average_diff
from add_substituent import add_methyls

VEHICLe = pd.read_csv('VEHICLe.csv')

dj1 = pd.read_csv('dj1/dj1.csv')
dj1_fulldata = pd.read_csv('dj1/dj1_clean_fulldata.csv')
dj1_clean = pd.read_csv('dj1/dj1_clean.csv')
dj1_clean_fulldata = pd.read_csv('dj1/data/dj1_clean_fulldata.csv', index_col=0)
dj1_logfiles = glob.glob('dj1/logfiles/*')
dj1anion_logfiles = glob.glob('dj1/logfiles/*anion*')
dj1bromine_logfiles = glob.glob('dj1/logfiles/*bromine*')
dj1neutral_logfiles = list(set(dj1_logfiles) - set(dj1anion_logfiles) - set(dj1bromine_logfiles))
independent_testset = pd.read_csv('dj1/imp_inputs/dj1clean_random_testset_for_sampling_atom_df.csv')
trainingset_for_sampling = pd.read_csv('dj1/dj1_training_for_sampling.csv')

dj2 = pd.read_csv('dj2/dj2.csv')
dj2_comfiles = glob.glob('dj2/comfiles/*')
dj2_fpsampled1000 = pd.read_csv('dj2/dj2_fpsampled1000.csv')
dj2_randsample1000 = pd.read_csv('dj2/dj2_randsample1000.csv')
dj2_logfiles = glob.glob('dj2/logfiles/*')
dj2anion_logfiles = glob.glob('dj2/logfiles/*anion*')
dj2bromine_logfiles = glob.glob('dj2/logfiles/*bromine*')
dj2neutral_logfiles = list(set(dj2_logfiles) - set(dj2anion_logfiles) - set(dj2bromine_logfiles))
dj2_fulldata = pd.read_csv('dj2/dj2_fulldata.csv', index_col=0)

dj3 = pd.read_csv('dj3/dj3.csv')
#dj3_for_dj1 = pd.read_csv('dj3/dj3_for_dj1.csv')
dj3_comfiles = glob.glob('dj3/comfiles/*')
dj3_logfiles = glob.glob('dj3/logfiles/*')
dj3_fulldata = pd.read_csv('dj3/dj3_fulldata.csv')
dj3anion_logfiles = glob.glob('dj3/logfiles/*anion*')
dj3bromine_logfiles = glob.glob('dj3/logfiles/*bromine*')
dj3neutral_logfiles = list(set(dj3_logfiles) - set(dj3anion_logfiles) - set(dj3bromine_logfiles))

dj3_full_atom = pd.read_csv('dj3/dj3_full_atom_df.csv')
dj3_full_pair = pd.read_csv('dj3/dj3_full_pair_df.csv')
dj3_train = pd.read_csv('dj3/dj3_train.csv_atom_df.csv')
dj3_test_atom = pd.read_csv('dj3/dj3_test.csv_pair_df.csv')


xyzfiles = glob.glob('xyzfiles/*')
dj1_dj2_fulldata = pd.concat([dj1_clean_fulldata, dj2_fulldata], axis=0)
dj1_dj2_neutral_logfiles = dj1neutral_logfiles + dj2neutral_logfiles

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
selectivity_targets = ['S2741', 'S1869', 'S5703', 'S1875', 'S5855', 'S2748', 'S5854', 'S6018', 'S5857']
exploration_targets = ['S7959', 'S5853', 'S9705', 'S54', 'S1072', 'S1873', 'S5704', 'S5707', 'S5706']

dj3_overlapped_structures = ['DR001', 'DR002', 'DR003', 'DR005', 'DR006', 'DR007', 'DR010', 'DR011', 'DR013', 'DR015', 'DR020', 'DR021',
                  'DR028', 'DR032', 'DR034', 'DR040', 'DR041', 'DR045', 'DR046', 'DR056', 'DR057', 'DR058', 'DR059',
                  'DR063', 'DR064']

dj3_for_dj1_structures = [x for x in dj3.Regid.unique() if x not in dj3_overlapped_structures]


#fulldata_filter(VEHICLe_dataframe=VEHICLe, calculate_properties_dataframe=dj1_clean_fulldata, outname='dj1/images/dj1_FDfilter_ActMap_A-0.6-1.5_Ea-0.6-0.8',
 #               a_upper=1.5, a_lower=0.0, ea_upper=0.8, ea_lower=0.6)

#selectivity_filter(VEHICLe_dataframe=VEHICLe, calculate_properties_dataframe=dj1_clean_fulldata,
#                   outname='dj1/images/dj1_SFilter_2.0-2.5', s_lower=2.0, s_upper=2.5, write=True)

#logfile_analysis_workflow(logfilelocation='dj2/logfiles', xyzfilelocation='dj2/xyzfiles', outname='dj2/dj2', calc_check=False,
#                          calc_data=True, calc_rel_data=True, plot_map=False, write=True)

#plot_activation_map(calculate_properties_dataframe=dj3_fulldata, outname='dj3/dj3_activationmap_unbound')

#fulldata_filter(VEHICLe_dataframe=dj3, calculate_properties_dataframe=dj3_fulldata, outname='dj3/images/check',
#                a_upper=100, a_lower=3.0, ea_upper=2.00, ea_lower=0.0)

#impression_input_workflow(logfiles=dj3neutral_logfiles, fulldata_df=dj3_fulldata, outname='dj3/dj3_full')

#samp = every2nd_sampling(dj1neutral_logfiles)

#impression_input_workflow(logfiles=samp[0], fulldata_df=dj1_fulldata, outname='dj1/dj1_train', write=True)
#impression_input_workflow(logfiles=samp[1], fulldata_df=dj1_fulldata, outname='dj1/dj1_test', write=True)

#plot_binary_activation_map(calculate_properties_dataframe=dj1_clean_fulldata, successful_activations=successful_activations,
#                           selectivity_targets=[], exploration_targets=[],
#                           outname='dj1/dj1_binary_activation_map', limit=1.2, write=True)

#impression_input_workflow(logfiles=dj1_dj2_neutral_logfiles, fulldata_df=dj1_dj2_fulldata, outname='dj1_dj2_full')


#s_upper = [0.9, 1.1, 1.3, 1.5, 1.7, 2.0, 2.5, 3.0]
#s_lower = [0.5, 0.9, 1.1, 1.3, 1.5, 1.7, 2.0, 2.5]

#for j, k in zip(s_lower, s_upper):

 #   outname = 'dj2/images/dj2_selectivity_filter_' + str(j) + '-' + str(k)

  #  selectivity_filter(VEHICLe_dataframe=VEHICLe, calculate_properties_dataframe=dj2_fulldata, outname=outname,
   #                        s_upper=k, s_lower=j, write=True)


def classifier_input(vehicle_dataframe, calculate_properties_dataframe, outname, write=True):

    regid = vehicle_dataframe['Regid']

    class_dataframe = pd.DataFrame(columns=['Regid', 'Smiles', 'Activation'])

    for i in regid:
        activation = ''
        props = calculate_properties_dataframe[calculate_properties_dataframe['Regid'] == i]
        prop_nan = props.fillna(0)

        anions = [j for j in prop_nan.columns if 'anion' in j]
        bromines = [k for k in prop_nan.columns if 'bromine' in k]

        if round(float(max(anions), 2))/round(float(max(bromines), 2)) < 1.20:
            activation = 'Acidic'

        else:
            activation = 'Nucleophilic'

        data = {'Regid': i, 'Smiles': vehicle_dataframe[vehicle_dataframe['Regid'] == i].Smiles.value,
                'Activation': activation}

        class_dataframe = class_dataframe.append(data)

        print(data)


#classifier_input(vehicle_dataframe=dj1_clean, calculate_properties_dataframe=dj1_clean_fulldata, outname=None)

aemols = logfile_to_aemol(glob.glob('dj2\\logfiles\\S447_2021-12-10_wb97xd_631gd_opt.log'))

for i in aemols:
    add_methyls(i)

