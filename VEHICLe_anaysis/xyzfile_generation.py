from datetime import date
def xyz_from_com(comfilename):
    todays_date = date.today()
    split = comfilename.split('\\')
    name = split[1].split('_')
    with open(comfilename) as p:
        xyz=[]
        read = p.readlines()[5:]
        remove_command=read[1:]
        charge_multip= read[0]
        charge=charge_multip[0]

        for line in remove_command:
            if not line.strip():
                continue
            else:
                xyz.append(line.strip('\n'))
        numatoms= len(xyz)

        with open('xyz_files/' + name[0] +'_' + str(todays_date) + '_' + name[1] + '.xyz', 'w') as g:
            print(numatoms, '\n',charge, file=g)
            for x in xyz:
                line=str(x)
                print(line.strip('\n'), file=g)
