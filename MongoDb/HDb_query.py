import pymongo
from pymongo import MongoClient
from pprint import pprint
import rdkit
from rdkit import Chem

cluster=MongoClient('mongodb+srv://yv19151:#Googleplex1213@cluster1.lwhvo.gcp.mongodb.net/HeterocycleDatabase?retryWrites=true&w=majority')
database=cluster['HeterocycleDatabase']
collection=database['Heterocycles']

p=collection.find({'Name': 'indole'})
