from random import randint
import pandas as pd
from rdkit import Chem as chm
from workflow import comfile_generation_workflow

def frag_type_array():

    frag_type_array = {'C': 'Methylated',
                       'F': 'Fluorinated',
                       '[Pd]': 'Palladated',
                       '[#6](-[#6])=[#8]': 'Acetylated',
                       '[#6]-[#6]': 'Et',
                       '[#6]-[#6]-[#6]': 'nPr',
                       '[#6]-[#6]-[#6]-[#6]': 'nBu',
                       '[#6](-[#6])(-[#6])-[#6]': 'tBu',
                       '[#6](-[#6])-[#6]': 'iPr',
                       '[#6]-[#9]': 'CFH2'}

    return frag_type_array

def add_substituent(dataframe, sample_pool, outname, frag_smarts='C', write=True, generate_com=True):

    'needs documenting, defaults to methylating all substrates randomly'

    sub_array = frag_type_array()

    substituted_dataframe = pd.DataFrame(columns=['Regid', 'Smiles'])

    for i in sample_pool:
        smiles = dataframe[dataframe['Regid'] == i].Smiles.values
        mol = chm.MolFromSmiles(smiles[0])
        molh = chm.AddHs(mol)
        h_frag = chm.MolFromSmarts('[H]')
        replace_frag = chm.MolFromSmarts(frag_smarts)

        mod_mol = chm.ReplaceSubstructs(molh, h_frag, replace_frag)

        structure_select = randint(0, len(mod_mol) - 1)

        structureh = mod_mol[structure_select]
        structure = chm.RemoveHs(structureh)

        sub_type = sub_array.get(frag_smarts)

        print(i + '-' + sub_type, chm.MolToSmiles(structure))

        sub_row = pd.DataFrame({'Regid': i + '-' + sub_type, 'Smiles': chm.MolToSmiles(structure)}, index=[0])

        substituted_dataframe = pd.concat([substituted_dataframe, sub_row])

    if write:
        with open(outname, 'w') as f:
            print(substituted_dataframe.to_csv(), file=f)

    if generate_com:
        comfile_generation_workflow(substituted_dataframe)
