from datetime import date

def xyz_from_type_array(filename, zipped_type_array):
    print('This function needs a fix: needs to include charge and number of atoms in the writing of the xyz file')
    Atoms = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 16: 'S'}
    with open(filename, 'w') as f:
        for i in zipped_type_array:
            for atomicnum, atom in Atoms.items():
                if i[0] == atomicnum:
                    print(atom, *i[1], sep=' ', file=f)

def xyz_from_com(comfilename):
    todays_date = date.today()
    split = comfilename.split('\\')
    name = split[1].split('_')

    with open(comfilename) as p:
        xyz=[]
        read = p.readlines()[5:]
        remove_command=read[1:] ### Why did I do this? ###
        charge_multip= read[0]
        charge=list(charge_multip)[0]


        for line in remove_command:
            if not line.strip():
                continue
            else:
                xyz.append(line.strip('\n'))
        numatoms= len(xyz)

        with open('xyz_files/' + name[0] +'_' + str(todays_date) + '_' + name[1] + '.xyz', 'w') as g:
            print(numatoms,'\n', 'charge=' + charge + '=', file=g, sep='')

            for x in xyz:
                line=str(x)
                print(line.strip('\n'), file=g)
