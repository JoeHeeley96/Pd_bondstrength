import pandas as pd
import glob

def g16_scfread(file):
    energy = -404
    with open(file, 'r') as f:
        for line in f:
            if 'SCF Done' in line:
                items = line.split()
                energy = float(items[4])

    return energy

def energy_readdata(file_locations, outname):

        energies = []
        date=[]
        structure=[]
        RegId = []

        list_of_files=glob.glob(file_locations + '/*')
        for file in list_of_files:
            energy = g16_scfread(file)
            if energy != -404:
                energies.append(energy)
                x = file.split('\\')
                z = x[1].split('_')
                RegId.append(z[0])
                print('WARNING: CHECK THIS LOOP IT HASNT BEEN TESTED')
                if z[1].str.contains('*anion'):
                    structure.append(z[1])
                    date.append(z[2])
                elif z[1].str.contains('*bromine'):
                    structure.append(z[2])
                    date.append(z[2])

                else:
                    structure.append('Base_het')
                    date.append(z[1])

            else:
                print(file, 'Energy not found')



        dict = {'Regid': RegId, 'Structure': structure, 'Energy': energies, 'Date': date}
        data = pd.DataFrame(dict)

        with open(outname, 'w') as y:
            print(data.to_csv(sep=','), file=y)


def energy_readreferences(reference_locations):
    Ref_E = []
    ExpNo = []
    functional = []
    basisset = []
    structure = []
    calctype = []

    for ref in reference_locations:
        ref_file = glob.glob(ref + '/*')
        for file in ref_file:
            refs = g16_scfread(file)
            if refs != -404:
                Ref_E.append(refs)

                x = file.split('\\')
                z = x[1].split('_')
                ExpNo.append(z[1])
                functional.append(z[-3])
                basisset.append(z[-2])
                structure.append(z[-4])
                calctype.append(z[-1])
            else:
                print(file, 'Energy not found')

    ref_dict = {'Energy': Ref_E, 'ExpNo': ExpNo, 'Functional': functional, 'Basis Set': basisset,
                'Calculation': calctype,
                'Structure': structure}

    global Reference_df
    Reference_df = pd.DataFrame(ref_dict)

    with open('reference_calcs.csv', 'w') as y:
        print(Reference_df.to_csv(sep=','), file=y)