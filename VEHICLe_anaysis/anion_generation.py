import glob
from rdkit import Chem
from rdkit.Chem import AllChem
from datetime import date

def anions_from_smiles(dataframe):

    regid = dataframe['Regid']

    for l in regid:
        comfiles = glob.glob('neutral_comfiles/' + str(l) + '_*')
        index_file = glob.glob('txt_files/' + str(l) + '_indexing.txt')
        H_index_list = []
        C_index_list = []
        C_numbering = []
        row = dataframe[dataframe.Regid == l]
        smiles = row.Smiles.item()

        mol = Chem.MolFromSmiles(smiles)
        molH = Chem.AddHs(mol)
        molc = AllChem.Compute2DCoords(molH)
        a = molH.GetAtoms()

        for j, atoms in enumerate(a):
            b = atoms.GetAtomicNum()
            d = atoms.GetIdx()
            if b == 1:
                e = atoms.GetNeighbors()
                for j, k in enumerate(e):
                    if k.GetAtomicNum() == 6:
                        if str(k.GetHybridization()) == 'SP2':
                            H_index_list.append(d)
                            C_index = k.GetIdx()
                            C_index_list.append(C_index)

        for x in range(len(C_index_list)):
            C_numbering.append(x + 1)

        match_Cnum_Cindex = zip(H_index_list, C_numbering)

        for file in comfiles:
            todays_date = str(date.today())
            split = file.split('\\')
            name = split[1].split('_')
            for m in index_file:
                for n in match_Cnum_Cindex:
                    with open(m) as p:
                        for line in p:
                            comindex = line[2] + line[3]
                            if float(n[0]) == float(comindex):
                                xyz = line[4:]
                                anion_comfilename = 'anion_comfiles/' + name[0] + '_' + str(n[1]) + 'anion_' + todays_date + '_' + name[2] + '-' + name[3] + '_opt.com'
                                anion_chkfilename = name[0] + '_' + str(n[1]) + 'anion_' + todays_date + '_' + name[ 2] + '-' + name[3] + '_opt.chk'
                                with open(anion_comfilename, 'w') as q:
                                    with open(file) as f:
                                        for comline in f:
                                            com_xyz=comline[2:]
                                            if comline.startswith('%chk='):
                                                print('%chk=' + anion_chkfilename, file=q)
                                            elif '0 1\n' in comline:
                                                print('-1 1', file=q)
                                            elif str(xyz).lstrip() != com_xyz:
                                                print(comline.strip('\n'), file=q)
