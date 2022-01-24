from datetime import date


def xyzcoords_from_type_array(filename, zipped_type_array):
    Atoms = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 16: 'S'}
    with open(filename, 'w') as f:
        for i in zipped_type_array:
            for atomicnum, atom in Atoms.items():
                if i[0] == atomicnum:
                    print(atom, i[1], *i[2], file=f)

def xyz_from_com(comfilename, save_loc):
    split = comfilename.split('\\')
    name = split[1].split('.')
    atoms=['C', 'H', 'N', 'O', 'S', 'Br']
    with open(comfilename) as p:
        xyz=[]
        read = p.readlines()[5:]
        remove_command=read[1:] ### Why did I do this? ###
        charge_multip= read[0]
        charge = int(charge_multip[0] + charge_multip[1])

        for line in remove_command:
            if not line.strip():
                continue
            else:
                for i in atoms:
                    if i in line:
                        xyz.append(line.strip('\n'))
        numatoms= len(xyz)

        with open(save_loc + '/' + name[0] + '.xyz', 'w') as g:
            print(numatoms,'\n', 'charge=' + str(charge) + '=', file=g, sep='')

            for x in xyz:
                line=str(x)
                print(line.strip('\n'), file=g)


