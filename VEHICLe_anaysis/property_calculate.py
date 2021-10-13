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
                    acidity_au =(i[3] + -229.010318493) - (j + -228.425492536)
                    acidity_kcalpermol = acidity_au * 627.5
                    regid_dict[i[2] + '_acidity_kcalpermol'] = acidity_kcalpermol

                elif 'bromine' in i.Structure:
                    elec_affin_au = (j + -2571.17518775) - i[3]
                    elec_affin_kcalpermol = elec_affin_au * 627.5
                    regid_dict[i[2] + '_electrophile_affinity_kcalpermol'] = elec_affin_kcalpermol

        property_data=property_data.append(regid_dict, ignore_index=True)

    if write:
        with open(outname, 'w') as f:
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
                rel_props['rel_' + j] = (24.038503/acidity)

            elif 'bromine' in j:
                elec_aff = float(df[j].values)
                rel_props['rel_' + j] = (elec_aff/183.651714)
        rel_props_df = rel_props_df.append(rel_props, ignore_index=True)

    fill_nan = rel_props_df.fillna(0)

    if write:
        with open(outname, 'w') as f:
            print(fill_nan.to_csv(sep=','), file=f)

    return fill_nan