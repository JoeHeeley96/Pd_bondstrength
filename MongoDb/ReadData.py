from typing import Any
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import cirpy
from cirpy import resolve
from pymongo import MongoClient
import time

def ReadData(mol2file):
	# Make an empty dictionary and make it accessible outside the function
	global data
	data = {}
	# Read in mol2 files as rdkit objects and returns error if it fails
	file = Chem.MolFromMol2File(mol2file, removeHs=False)
	if file == None:
		print('File input error')
	else:
		# Define properties you want from the rdkit.Chem. as a dictionary (make sure to add file extension if needed!)
		rdkitStrings = {'SMILES': 'MolToSmiles',
						'XYZcoords': 'MolToXYZBlock',
						'Molecular Formula': 'rdMolDescriptors.CalcMolFormula',
						'Molecular Mass': 'rdMolDescriptors.CalcExactMolWt'}

		# Define properties you want from cirpy as a dictionary
		cirpyStrings = {'Name': 'names',
						'CAS': 'cas'}

		# Get the info from cirpy using the smiles string we got above and add them to dictionary d
		for key in cirpyStrings.keys():
			smiles = rdkit.Chem.MolToSmiles(file)
			cstring = cirpyStrings[key]
			x = cirpy.resolve(smiles, cstring)

			if key not in data.values():
				data[key] = x

		# Use the rdkitStrings to give us info from rdkit and add to dictionary d
		for key in rdkitStrings.keys():
			rdstring = rdkitStrings[key]
			a = 'rdkit.Chem.' + rdstring + '(file)'
			b = eval(a)

			if key not in data.values():
				data[key] = b

		# Lets us know where the gaps are in cirpy
		for x, y in data.items():
			if y == None:
				data[x] = 'Not Retrievable'

	#Print data so we can have a look at it
	print('Printing data:')
	time.sleep(1)
	print(data)
	time.sleep(1)

	for key in data.keys():
		if data[key]=='Not Retrievable':
			print('WARNING!', key, 'Data not found')


	#Asks user whether to submit data to Heterocycle Database or not
	prompt: Any=input('Submit data to Heterocycle Database?')
	if prompt =='no':
		print('Data not submitted')

	#Submits to Hdb (this is a test, so actually goes to testdb, CHANGE ME!)
	elif prompt =='yes':
		print('Submitting data to Heterocycle Database...')
		cluster=MongoClient('mongodb+srv://yv19151:#Googleplex1213@cluster1.lwhvo.gcp.mongodb.net/HeterocycleDatabase?retryWrites=true&w=majority')
		database=cluster['HeterocycleDatabase']
		collection=database['Heterocycles']

		query=collection.find_one({'SMILES': data['SMILES']})
		if query == None:
			collection.insert_one(data)
			print('Data submitted!')
		else:
			print('This data has already been submitted!')


FileNames=['pyridazine.mol2', 'pyridine.mol2', 'pyrido_32d_pyrimidine.mol2',
		   'pyrimidine.mol2', 'quinazoline.mol2', 'quinoline.mol2', 'quinoxaline.mol2', '5-azaindole.mol2',
		   '7-azaindazole.mol2']

ReadData('data/15-napthyridine.mol2')

'''for file in FileNames:
	string='data/' + file
	ReadData(string)'''

#To do:
#Popuate HDb chemicals collection
#Create another collection for the results for calculations
