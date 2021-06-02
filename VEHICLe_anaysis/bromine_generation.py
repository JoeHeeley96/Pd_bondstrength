from xyzfile_generation import xyz_from_com
from xyz2mol import xyz2mol
from datetime import date
import glob
from rdkit import Chem
from rdkit.Chem import AllChem


def bromines_from_smiles(dataframe):
    regid = dataframe['Regid']

    for l in regid:
        C_index_list = []
        C_numbering = []
        comfiles = glob.glob('neutral_comfiles/' + str(l) + '_*')
        index_file = glob.glob('txt_files/' + str(l) + '_indexing.txt')
        row = dataframe[dataframe.Regid == l]
        smiles = row.Smiles.item()

        mol = Chem.MolFromSmiles(smiles)
        molH = Chem.AddHs(mol)
        molc = AllChem.Compute2DCoords(molH)
        a = molH.GetAtoms()

        for j, atoms in enumerate(a):
            b = atoms.GetAtomicNum()
            if b == 1:
                e = atoms.GetNeighbors()
                for h, k in enumerate(e):
                    if k.GetAtomicNum() == 6:
                        if str(k.GetHybridization()) == 'SP2':
                            C_index = k.GetIdx()
                            C_index_list.append(C_index)

        for x in range(len(C_index_list)):
            C_numbering.append(x + 1)

        match_Cnum_Cindex = zip(C_index_list, C_numbering)

        for file in comfiles:
            todays_date = str(date.today())
            split = file.split('\\')
            name = split[1].split('_')
            for m in index_file:
                for n, o in match_Cnum_Cindex:
                    with open(m) as p:
                        for line in p:
                            comindex = line[2] + line[3]
                            if float(n) == float(comindex):
                                bromine_comfilename = 'bromine_comfiles/' + name[0] + '_' + str(o) + 'bromine_' + todays_date + '_' + name[2] + '-' + name[3] + '_opt.com'
                                bromine_chkfilename = name[0] + '_' + str(o) + 'bromine_' + todays_date + '_' + name[
                                    2] + '-' + name[3] + '_opt.chk'
                                split = line.split(' ')
                                Br_xcoord = float(split[2])
                                Br_ycoord = float(split[3])
                                Br_zcoord = float(split[4])
                                Br_coords = (str(Br_xcoord) + ' ' + str(Br_ycoord) + ' ' + str(Br_zcoord))

                        with open(bromine_comfilename, 'w') as q:
                            with open(file) as f:
                                for comline in f:
                                    split=comline.split()

                                    if comline.startswith('%chk='):
                                        print('%chk=' + bromine_chkfilename, file=q)

                                    elif comline == '0 1\n':
                                        print('1 1', file=q)

                                    elif Br_coords in comline:
                                        comsplit=comline.split()
                                        if Br_xcoord == float(comsplit[1]):
                                            if Br_ycoord == float(comsplit[2]):
                                                if Br_zcoord == float(comsplit[3]):
                                                    print(comline.strip('\n'), file=q)
                                                    print('Br', Br_xcoord, Br_ycoord, (Br_zcoord + 2), file=q)

                                    else:
                                        print(comline.strip('\n'), file=q)