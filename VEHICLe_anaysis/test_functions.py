import pandas as pd
import glob
from datetime import date
from bromine_generation import bromines_from_smiles
from xyzfile_generation import xyz_from_com
from VEHICLe_filters import djr_filter
from VEHICLe_filters import atom_filter
from anion_generation import anions_from_smiles
from plots import plot_activation_map
import numpy as np
from VEHICLe_filters import synthesis_filter
from rdkit.Chem.FragmentMatcher import FragmentMatcher
from plots import plot_3Dactivation_map
from plots import plot_activation_map
from sklearn import metrics
from impression_input import write_imp_input
from tqdm import tqdm
import rdkit
from impression_input import exclude_structures
from impression_input import file_to_aemol
from impression_input import xyzfile_to_aemol
from workflow import logfile_analysis_workflow
from mol_translator.aemol import aemol
from workflow import impression_input_workflow
from sampling import get_tm_df
from sampling import every2nd_sampling
from workflow import comfile_generation_workflow
from workflow import logfile_fingerprint_sampling_workflow
from workflow import vehicle_fingerprint_sampling_workflow
from VEHICLe_filters import atom_filter

VEHICLe = pd.read_csv('VEHICLe.csv')

dj1 = pd.read_csv('dj1/dj1.csv')
dj1_fulldata = pd.read_csv('dj1/data/dj1_clean_fulldata.csv')
dj1_clean = pd.read_csv('dj1/dj1_clean.csv')
dj1_clean_fulldata = pd.read_csv('dj1/data/dj1_clean_fulldata.csv', index_col=0)
dj1_clean_fulldata_relative = pd.read_csv('dj1/data/dj1_fulldata_relative.csv')
dj1_logfiles = glob.glob('dj1/logfiles/*')
dj1anion_logfiles = glob.glob('dj1/logfiles/*anion*')
dj1bromine_logfiles = glob.glob('dj1/logfiles/*bromine*')
dj1neutral_logfiles = list(set(dj1_logfiles) - set(dj1anion_logfiles) - set(dj1bromine_logfiles))
dj1_sub_fulldata = pd.read_csv('dj1/dj1.1/DJ1.1_FINAL_fulldata.csv')
independent_testset = pd.read_csv('dj1/imp_inputs/dj1clean_random_testset_for_sampling_atom_df.csv')

dj2 = pd.read_csv('dj2/dj2.csv')
dj2_comfiles = glob.glob('dj2/comfiles/*')
dj2_fpsampled1000 = pd.read_csv('dj2/dj2_fpsampled1000.csv')
dj2_randsample1000 = pd.read_csv('dj2/dj2_randsample1000.csv')
dj2_logfiles = glob.glob('dj2/logfiles/*')
dj2anion_logfiles = glob.glob('dj2/logfiles/*anion*')
dj2bromine_logfiles = glob.glob('dj2/logfiles/*bromine*')
dj2neutral_logfiles = list(set(dj2_logfiles) - set(dj2anion_logfiles) - set(dj2bromine_logfiles))
dj2methyl_logfiles = glob.glob('dj2/logfiles/*Methylated*')
dj2fluorine_logfiles = glob.glob('dj2/logfiles/*Fluorinated*')
dj2_fulldata = pd.read_csv('dj2/data/dj2_fulldata.csv', index_col=0)
dj2_neutrallogfiles_unsub = list(set(dj2neutral_logfiles) - set(dj2methyl_logfiles) - set(dj2fluorine_logfiles))
dj2_sub_fulldata = pd.read_csv('dj2/dj2.1/dj2_sub_fulldata.csv')

dj3 = pd.read_csv('dj3/dj3.csv')
#dj3_for_dj1 = pd.read_csv('dj3/dj3_for_dj1.csv')
dj3_comfiles = glob.glob('dj3/comfiles/*')
dj3_logfiles = glob.glob('dj3/logfiles/*')
dj3_fulldata = pd.read_csv('dj3/dj3_fulldata.csv')
dj3anion_logfiles = glob.glob('dj3/logfiles/*anion*')
dj3bromine_logfiles = glob.glob('dj3/logfiles/*bromine*')
dj3neutral_logfiles = list(set(dj3_logfiles) - set(dj3anion_logfiles) - set(dj3bromine_logfiles))

#dj3_full_atom = pd.read_csv('dj3/dj3_full_atom_df.csv')
#dj3_full_pair = pd.read_csv('dj3/dj3_full_pair_df.csv')
#dj3_train = pd.read_csv('dj3/dj3_train.csv_atom_df.csv')
#dj3_test_atom = pd.read_csv('dj3/dj3_test.csv_pair_df.csv')


xyzfiles = glob.glob('xyzfiles/*')
dj1_dj2_fulldata = pd.concat([dj1_clean_fulldata, dj2_fulldata], axis=0)
dj1_dj2_neutral_logfiles = dj1neutral_logfiles + dj2_neutrallogfiles_unsub
test = np.array(dj1neutral_logfiles + dj2_neutrallogfiles_unsub)

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
dj1_sub_and_dj2_fulldata = dj1_sub_fulldata.append(dj2_fulldata)

dj1sub_atom = pd.read_csv('dj1/dj1.1/dj1_ext_atom_df.csv')
dj2_atom = pd.read_csv('dj2/imp_inputs/dj2_full_atom_df.csv')

dj1_s = dj1sub_atom['molecule_name'].unique()
dj2_s = dj2_atom['molecule_name'].unique()

dj2_real = [x for x in dj2_s if x not in dj1_s]

dj2_logfiles_for_dj1 = []
k = [x for x in dj2_s if x in dj1_s]

for i in dj2_neutrallogfiles_unsub:
    for j in dj2_real:
        if j + '_' in i:

            dj2_logfiles_for_dj1.append(i) if i not in dj2_logfiles_for_dj1 else dj2_logfiles_for_dj1


dj1_sub_and_dj2_unsub = dj1neutral_logfiles + dj2_logfiles_for_dj1

dj2_sub_atom = pd.read_csv('dj2/dj2.1/dj2_sub_atom_df.csv')
dj2_sub_pair = pd.read_csv('dj2/dj2.1/dj2_sub_pair_df.csv')

dj1_atom = pd.read_csv('dj1/imp_inputs/dj1_full_atom_df.csv')
dj1_pair = pd.read_csv('dj1/imp_inputs/dj1_full_pair_df.csv')

dj2_atom = pd.read_csv('dj2/imp_inputs/dj2_full_atom_df.csv')

exc = dj2_atom['molecule_name'].unique()

#exclude_structures(exc, dj1_atom, dj1_pair, outname='dj2/dj2.1/dj1_unsub_test_for_dj2_sub')

syn = synthesis_filter(dj2, export=False)

djr_filter(syn, dj2_fulldata, outname='dj2/dj2_synthesis_filter_split_0.8', djr_upper=1000, djr_lower=0, write=True)

#impression_input_workflow(neutral_logfiles=dj2neutral_logfiles, fulldata_df=dj2_sub_fulldata,
 #                         outname='dj2/dj2.1/dj2_sub', write=True)

#logfile_analysis_workflow('dj2/logfiles', 'dj2/xyzfiles', 'dj2/dj2_sub', calc_check=False, plot_map=True)
