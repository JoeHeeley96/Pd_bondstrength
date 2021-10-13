def output_structure_check(xyzfiles, logfiles):
    '''
    requires inputs in alphabetical order
    '''

    for i,j in zip(xyzfiles, logfiles):
        in_file = i.split('\\')
        in_split = in_file[1].split('-')

        out_file = j.split('\\')
        out_split = out_file[1].split('-')
        if in_split[0] != out_split[0]:
            raise NameError('Files dont match')

        input_aemol = aemol(in_file[1].split('_')[0])
        input_aemol.from_file(i, ftype = 'xyz')

        output_aemol = aemol(in_file[1].split('_')[0])
        output_aemol.from_file(j, ftype='log')

        input_structure = output_aemol.structure['conn']
        output_structure = input_aemol.structure['conn']
        #print(in_split[0])
        print(input_structure)
        print(output_structure)

        for x,y in zip(input_structure, output_structure):
            print(np.count_nonzero(x == y))