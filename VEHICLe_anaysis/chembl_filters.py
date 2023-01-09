import pandas as pd
from rdkit import Chem
import numpy as np
from tqdm import tqdm
import rdkit.Chem

def isRingAromatic(mol, bondRing):
    for id in bondRing:
        if not mol.GetBondWithIdx(id).GetIsAromatic():
            return False
    return True

def dj1_point2_filter(chembl_data, outname):

    regid = []
    smiles = []
    amino_acid = Chem.MolFromSmiles('NCC(NCC(NCC(N)=O)=O)=O')

    with open(chembl_data, 'r') as f:
        for line in f:
            goodmol = True
            ch = False
            linsplit = line.split()

            if 'canonical_smiles' not in linsplit[1]:
                rdmol = Chem.MolFromSmiles(linsplit[1])

                try:
                    for molatom in rdmol.GetAtoms():
                        if molatom.GetAtomicNum() not in range(6, 10):
                            goodmol = False

                except AttributeError:
                    goodmol = False
                    print('Error:', linsplit[0], linsplit[1], 'Not Recognized')

                if goodmol:
                    if rdmol.HasSubstructMatch(amino_acid):
                        pass

                    elif len(rdmol.GetAtoms()) > 50:
                        pass

                    else:
                        hmol = Chem.AddHs(rdmol)
                        ri = rdmol.GetRingInfo()

                        for atoms, bonds in zip(ri.AtomRings(), ri.BondRings()):
                            if isRingAromatic(rdmol, bonds):
                                aromatic_type_array = []
                                for atom in atoms:
                                    atomtype = rdmol.GetAtomWithIdx(atom).GetAtomicNum()
                                    if atomtype == 6:
                                        for n in hmol.GetAtomWithIdx(atom).GetNeighbors():
                                            if n.GetAtomicNum() == 1:
                                                ch = True

                                    aromatic_type_array.append(atomtype)

                                if 16 not in aromatic_type_array:
                                    if 8 not in aromatic_type_array:
                                        if ch:
                                            regid.append(linsplit[0])
                                            smiles.append(linsplit[1])

    data = pd.DataFrame({'Regid': regid, 'Smiles': smiles})
    print(len(data['Regid'].unique()))

    with open(outname + '.csv', 'w') as p:
        print(data.to_csv(), file=p)
