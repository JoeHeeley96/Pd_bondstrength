import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import cirpy
import urllib
import pymongo
from pymongo import MongoClient
from pymongo import MongoClient


def mol2_to_mongo(mol2file):
    # Make an empty dictionary and make it accessible outside the function
    global d
    d = {}
    # Read in mol2 files as rdkit objects and returns error if it fails
    file = Chem.MolFromMol2File(mol2file, removeHs=False)
    if file == None:
        print('File input error')
    else:
        # Define properties you want from the redkit.Chem. as a dictionary (make sure to add file extension if needed!)
        rdkitStrings = {'SMILES': 'MolToSmiles',
                        'XYZcoords': 'MolToXYZBlock',
                        'Molecular Formula': 'rdMolDescriptors.CalcMolFormula',
                        'Molecular Mass': 'rdMolDescriptors.CalcExactMolWt'}

        # Define properties you want from cirpy as a dictionary
        cirpyStrings = {'Name': 'names',
                        'CAS': 'cas'}

        # Use the rdkitStrings to give us info from rdkit and add to dictionary d
        for key in rdkitStrings.keys():
            rdstring = rdkitStrings[key]
            a = 'rdkit.Chem.' + rdstring + '(file)'
            b = eval(a)

            if key not in d.values():
                d[key] = b

        # Get the info from cirpy using the smiles string we got above and add them to dictionary d
        for key in cirpyStrings.keys():
            smiles = d['SMILES']
            cstring = cirpyStrings[key]
            x = cirpy.resolve(smiles, cstring)

            if key not in d.values():
                d[key] = x
        # Lets us know where the gaps are in cirpy
        for x, y in d.items():
            if y == None:
                d[x] = 'Not Retrievable'

        print(d)

    prompt=input('Enter Data into Heterocycle Database?')
    if prompt == 'no':
        print('Data not Entered')
