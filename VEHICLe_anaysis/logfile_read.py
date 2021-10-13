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

def energy_readdata(file_locations, outname, write=True):

        energies = []
        date=[]
        structure=[]
        RegId = []
        theory=[]
        calc=[]

        list_of_files=glob.glob(file_locations + '/*')
        for file in list_of_files:
            energy = g16_scfread(file)
            if energy != -404:
                energies.append(energy)
                x = file.split('\\')
                z = x[1].split('_')
                RegId.append(z[0])
                theory.append(z[3])
                calc.append(z[-1])

                if 'anion' in z[1]:
                    structure.append(z[1])
                    date.append(z[2])
                elif 'bromine' in z[1]:
                    structure.append(z[1])
                    date.append(z[2])

                else:
                    structure.append('Base_het')
                    date.append(z[1])

            else:
                print(file, 'Energy not found')

        dict = {'Regid': RegId, 'Structure': structure, 'Energy': energies, 'Date': date, 'Functional/basisset': theory, 'Calculation': calc}
        data = pd.DataFrame(dict)

        if write:
            with open(outname, 'w') as f:
                print(data.to_csv(sep=','), file=f)

        return data

