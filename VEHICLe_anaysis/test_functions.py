import pandas as pd
import glob
import numpy as np
from plots import plot_activation_map

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
dj3_comfiles = glob.glob('dj3/comfiles/*')
dj3_logfiles = glob.glob('dj3/logfiles/*')
dj3_fulldata = pd.read_csv('dj3/dj3_fulldata.csv')
dj3anion_logfiles = glob.glob('dj3/logfiles/*anion*')
dj3bromine_logfiles = glob.glob('dj3/logfiles/*bromine*')
dj3neutral_logfiles = list(set(dj3_logfiles) - set(dj3anion_logfiles) - set(dj3bromine_logfiles))

dj1_dj2_fulldata = pd.concat([dj1_clean_fulldata, dj2_fulldata], axis=0)
dj1_dj2_neutral_logfiles = dj1neutral_logfiles + dj2_neutrallogfiles_unsub
test = np.array(dj1neutral_logfiles + dj2_neutrallogfiles_unsub)

dj1_rel = pd.read_csv('dj1/data/dj1_fulldata_relative.csv')
dj2_rel = pd.read_csv('dj2/data/dj2_fulldata_relative.csv')
dj3_rel = pd.read_csv('dj3/dj3_fulldata_relative.csv')

dj3_overlapped_structures = ['DR001', 'DR002', 'DR003', 'DR005', 'DR006', 'DR007', 'DR010', 'DR011', 'DR013', 'DR015', 'DR020', 'DR021',
                  'DR028', 'DR032', 'DR034', 'DR040', 'DR041', 'DR045', 'DR046', 'DR056', 'DR057', 'DR058', 'DR059',
                  'DR063', 'DR064']

dj3_for_dj1_structures = [x for x in dj3.Regid.unique() if x not in dj3_overlapped_structures]
dj1_sub_and_dj2_fulldata = dj1_sub_fulldata.append(dj2_fulldata)


succ_act = ['S1080', 'S1875', 'S1863', 'S5857', 'S6018', 'S5854', 'S5703', 'S1869', 'S5853', 'S1072']

unsucc_act = ['S5704', 'S2741']

targeted = ['S291', 'S5706', 'S5707', 'S2741', 'S7959', 'S54', 'S8063', 'S7959', 'S1873',
            'S5706', 'S5855', 'S2746', 'S5704', 'S49', 'S50']

paper_targets = ['S5704', 'S2741', 'S7959', 'S54', 'S1873']

plot_activation_map(dj1_clean_fulldata, outname='DJ1_bare_activation_map', plot_title='DJ1 Activation Map', highlight=False)

dj3_exc = ['DR014', 'DR042', 'DR043', 'DR044', 'DR079', 'DR080', 'DR081', 'DR082', 'DR083', 'DR084', 'DR085',
           'DR086', 'DR087', 'DR088', 'DR001', 'DR002', 'DR003', 'DR005', 'DR006', 'DR007', 'DR010', 'DR011', 'DR013',
           'DR015', 'DR020', 'DR021', 'DR028', 'DR032', 'DR034', 'DR040', 'DR041', 'DR045', 'DR046', 'DR056', 'DR057',
           'DR058', 'DR059', 'DR063', 'DR064', 'DR009', ]

'''from impression_input import classifier_input
from add_substituent import add_substituent

rdmol = rdkit.Chem.MolFromSmiles('CF')
smart = rdkit.Chem.MolToSmarts(rdmol)

print('molsmart to add:', smart)

samp = random_sample(dj1_clean['Regid'], length=50)

add_substituent(dj1_clean, sample_pool=samp, outname='dj1_CFH2.csv',
                frag_smarts=smart, write=True, generate_com=True)'''


