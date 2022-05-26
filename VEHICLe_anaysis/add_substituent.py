from random import randint
import pandas as pd
from rdkit import Chem as chm
from workflow import comfile_generation_workflow

def frag_type_array():

    frag_type_array = {'C': 'Methylated', 'F': 'Fluorinated'}

    return frag_type_array

def add_substituent(dataframe, outname, frag_smarts='C', write=True, generate_com=True):

    'needs documenting, defaults to methylating all substrates randomly'

    smiles = dataframe['Smiles']

    substituted_dataframe = pd.DataFrame()

    for i in smiles:
        regid = dataframe[dataframe['Smiles'] == i].Regid.values
        mol = chm.MolFromSmiles(i)
        molh = chm.AddHs(mol)
        h_frag = chm.MolFromSmarts('[H]')
        replace_frag = chm.MolFromSmarts(frag_smarts)

        mod_mol = chm.ReplaceSubstructs(molh, h_frag, replace_frag)

        structure_select = randint(0, len(mod_mol) - 1)

        structureh = mod_mol[structure_select]
        structure = chm.RemoveHs(structureh)

        sub_array = frag_type_array()
        sub_type = sub_array.get(frag_smarts)

        substituted_dataframe = substituted_dataframe.append({'Regid': regid[0] + '-' + sub_type, 'Smiles': chm.MolToSmiles(structure)}, ignore_index=True)

    if write:
        with open(outname, 'w') as f:
            print(substituted_dataframe.to_csv(), file=f)

    if generate_com:
        comfile_generation_workflow(substituted_dataframe)




