import sys
import glob
from mol_translator import aemol
from mol_translator.structure.find_paths import pybmol_find_paths, pybmol_find_all_paths
import pandas as pd
import numpy as np
from tqdm import tqdm
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
from random import randint
import os
import shutil
import copy

def first200_sampling(neutral_logfiles):

    training_set = neutral_logfiles[0:200]
    test_set = []

    for i in neutral_logfiles:
        if i not in training_set:
            test_set.append(i)

    for j in test_set:
        if j in training_set:
            raise ValueError('Same molecules in the test and training set')

    print('training set contains: ', len(training_set), 'molecules')
    print('test set contains: ', len(test_set), 'molecules')

    return training_set, test_set


def every2nd_sampling(neutral_logfiles):

    training_set = sorted(neutral_logfiles, key=str)[0::2]
    test_set = []

    for i in neutral_logfiles:
        if i not in training_set:
            test_set.append(i)

    print('training set contains: ', len(training_set), 'molecules')
    print('test set contains: ', len(test_set), 'molecules')

    return training_set, test_set



def random_sample(sample_pool, length):
    '''
    :param sample_pool: a list of .log files you want to sample from
    :param length: how many structures you want to sample
    :return: a list of .log files that have been randomly sampled
    '''

    sampled_structures = []

    for j in range(400):
        value = randint(0, len(sample_pool))

        if len(sampled_structures) < length:
            file = sample_pool[value]
            sampled_structures.append(file) if file not in sampled_structures else sampled_structures

    print('Sampled', len(sampled_structures), 'structures')

    return sampled_structures

def get_tm_df(training_structures, vehicle_df):
    '''
    Should return a dataframe of rdkit fingerprints and an array of fingerprint similarity
    :param training_structures: should be a list of regids you want to sample from
    :param VEHICLe_subset: should contain the SMILES strings for all your regids
                (you can use either the subset you're working from or just the full VEHICLe database)

    :return: tm_df (dataframe of molecular fingerprints)
    :return: tm_array (fingerprint similarity array)
    '''

    names = []
    ecfps = []
    for j in tqdm(training_structures):
        smiles = vehicle_df.loc[vehicle_df['Regid'] == j].Smiles.item()
        mol = Chem.MolFromSmiles(smiles)
        molH = Chem.AddHs(mol)
        fp = AllChem.GetMorganFingerprintAsBitVect(molH, radius=2, nBits=2048)
        names.append(j)
        ecfps.append(fp)

    tm_array = np.zeros(shape=(len(names), len(names)), dtype='float32')
    for i1, (id1, fp1) in tqdm(enumerate(zip(names, ecfps))):
        for i2, (id2, fp2) in enumerate(zip(names, ecfps)):
            tm_array[i1][i2] = DataStructs.FingerprintSimilarity(fp1, fp2)

    tm_df = pd.DataFrame(list(zip(names, ecfps)), columns=['molecule_name', 'ecfp'])

    return tm_df, tm_array

def get_fps_ids(tm_df, tm_array, outname, num=100, write=False):
    '''
    :param tm_df: output [0] from get_tm_df
    :param tm_array: output [1] from get_tm_df
    :param num: number of sampled you want to sample

    :return: list of len(num) of structures to sample,
                also writes the list as a .txt file if write=True
    '''
    pbar = tqdm(total=num)
    selected = list(np.random.choice(tm_df.molecule_name.unique(), 1))
    for molid in selected:
        pbar.update(1)

    while len(selected) < num:
        selected_ids = [x in selected for x in tm_df.molecule_name.values]
        big_d = 0.0
        to_pick = ''
        tm_array = tm_array
        for molid in tm_df.molecule_name.unique():
            if molid in selected:
                continue
            id1 = tm_df.loc[(tm_df.molecule_name == molid)].index.values[0]
            distance = np.sum(tm_array[id1][selected_ids])
            if distance > big_d:
                big_d = distance
                to_pick = molid
        pbar.update(1)
        selected.append(to_pick)

    print('selected', len(selected), 'structures')

    if write:
        with open(outname + str(num) + '.txt', 'w') as f:
            for j, molid in enumerate(selected):
                string = "{0:<10d}\t,\t{1:<10s}".format(j, molid)
                print(molid, file=f)

    return selected