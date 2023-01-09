import pandas as pd


def calculate_properties(dataset, outname, write=True):
    regId = list(set(dataset.Regid))
    property_data=pd.DataFrame(columns=['Regid'])

    for id in regId:
        regid_dict= {'Regid': id}
        regid_data = dataset[dataset['Regid'].str.fullmatch(id)]
        alist = list(regid_data.itertuples(index=False))
        basehet = (regid_data[regid_data['Structure'].str.contains('Base_het')]).Energy

        for j in basehet:
            for i in alist:
                if 'anion' in i.Structure:
                    acidity_au =(i[3] + -229.010318493) - (j + -228.403959536)
                    acidity_kcalpermol = acidity_au * 627.5
                    regid_dict[i[2] + '_acidity_kcalpermol'] = acidity_kcalpermol

                elif 'bromine' in i.Structure:
                    elec_affin_au = (j + -2571.19136375) - i[3]
                    elec_affin_kcalpermol = elec_affin_au * 627.5
                    regid_dict[i[2] + '_electrophile_affinity_kcalpermol'] = elec_affin_kcalpermol

        property_data=property_data.append(regid_dict, ignore_index=True)

    if write:
        with open(outname, 'w') as f:
            print(property_data.to_csv(sep=','), file=f)
    return property_data

def calculate_deltaG(deltaG_rawdata, outname, write=True):

    regId = list(set(deltaG_rawdata.Regid))
    property_data = pd.DataFrame(columns=['Regid'])

    for id in regId:
        regid_dict = {'Regid': id}
        regid_data = deltaG_rawdata[deltaG_rawdata['Regid'].str.fullmatch(id)]
        alist = list(regid_data.itertuples(index=False))
        basehet = (regid_data[regid_data['Structure'].str.contains('Base_het')]).Energy

        for j in basehet:
            for i in alist:
                if 'anion' in i.Structure:
                    acidity_au = (i[-1] + -229.046860464999) - (j + -228.477217341)
                    acidity_kcalpermol = acidity_au * 627.5
                    regid_dict[i[2] + '_acidity_kcalpermol'] = acidity_kcalpermol

                elif 'bromine' in i.Structure:
                    elec_affin_au = (j + -2571.16789474) - i[-1]
                    elec_affin_kcalpermol = elec_affin_au * 627.5
                    regid_dict[i[2] + '_electrophile_affinity_kcalpermol'] = elec_affin_kcalpermol

        property_data = pd.concat([property_data, pd.DataFrame(regid_dict, index=[0])], ignore_index=True)

    if write:
        with open(outname + '.csv', 'w') as f:
            print(property_data.to_csv(sep=','), file=f)
    return property_data

def calculate_relative_properties(calculate_properties_dataframe, outname, write=True):
    regid=list(set(calculate_properties_dataframe.Regid))
    rel_props_df = pd.DataFrame([])

    for i in regid:
        rel_props = {}
        rel_props['Regid'] = i
        df = calculate_properties_dataframe.loc[calculate_properties_dataframe['Regid'] == i]

        for j in df.columns:
            if 'anion' in j:
                acidity = float(df[j].values)
                rel_props['rel_' + j] = (24.044391262487466/acidity)

            elif 'bromine' in j:
                elec_aff = float(df[j].values)
                rel_props['rel_' + j] = (elec_aff/183.647244829844)
        rel_props_df = rel_props_df.append(rel_props, ignore_index=True)

    fill_nan = rel_props_df.fillna(0)

    if write:
        with open(outname, 'w') as f:
            print(fill_nan.to_csv(sep=','), file=f)

    return fill_nan


def find_average_diff(calculate_properties_dataframe):

    regids = calculate_properties_dataframe['Regid']
    diff_acidity = []
    diff_elec = []

    for i in regids:
        acidity = []
        elec_aff = []
        props = calculate_properties_dataframe[calculate_properties_dataframe['Regid'] == i]
        prop_nan = props.fillna(0)

        anions = [j for j in prop_nan.columns if 'anion' in j]
        bromines = [k for k in prop_nan.columns if 'bromine' in k]

        for l in anions:
            acidity = acidity + prop_nan[l].values.tolist()

        for m in bromines:
            elec_aff = elec_aff + prop_nan[m].values.tolist()

        a_max = max(acidity)
        a2_max = max([x for x in acidity if x != a_max])

        b_max = max(elec_aff)
        b2_max = max([y for y in elec_aff if y != b_max])

        if a2_max != 0:
            diff_acidity.append(a_max - a2_max)

        if b2_max != 0:
            diff_elec.append(b_max - b2_max)

    print('Average delta acidity: ', sum(diff_acidity) / len(diff_acidity))
    print('Average delta Electrophile Affinity: ', sum(diff_elec) / len(diff_elec))