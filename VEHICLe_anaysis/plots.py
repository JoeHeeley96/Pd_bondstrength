import itertools

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import math
from impression_input import exclude_structures_prop_dfs
from mpl_toolkits import mplot3d
from matplotlib import animation
import numpy as np
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import matplotlib.lines as mlines
from sklearn import linear_model
from sklearn.metrics import r2_score
#from rotation import rotanimate
from mpl_toolkits.mplot3d import Axes3D

def plot_animation_3D(xdata, ydata, zdata, regids, outname):
    fig = plt.figure()
    ax1 = plt.subplot(111, projection='3d')
    num_colours = 250
    cm = plt.get_cmap('gist_rainbow')
    cNorm = colors.Normalize(vmin=0, vmax=num_colours - 3)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    ax1.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(num_colours)])

    def init():
        for x,y,z,regid in zip(xdata,ydata,zdata,regids):
            if x < 1.2:
                if y > 0:
                    plot = ax1.scatter(x, y, z, label=regid, marker='x', zorder=2)
        plt.xlabel('Relative Electrophile Affinity')
        plt.ylabel('Relative Acidity')
        ax1.set_zlabel('Relative Differentiation')
        plt.title('Heterocycle properties relative to Pyrazolo[1,5-a]pyrimidine')
        plt.autoscale(enable=True, axis='both',tight=None)
        ax1.legend(bbox_to_anchor=(1.05, 1, 2, 1))
        return fig,

    def animate(i):
        ax1.view_init(elev=10, azim=i * 4)
        return fig,

    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=90, interval=50, blit=True)
    fn = outname
    ani.save(fn + '.gif', writer='pillow', fps=100)


def plot_activation_map(calculate_properties_dataframe, outname, plot_title='Activation_map', highlight=True):

    plt.clf()
    fig1, ax1=plt.subplots()
    num_colours = 250

    cm = plt.get_cmap('gist_rainbow')
    cNorm = colors.Normalize(vmin=0, vmax=num_colours - 3)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    ax1.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(num_colours)])
    regid=list(set(calculate_properties_dataframe.Regid))

    for i in regid:
        acidities= {}
        elec_affs= {}
        relative_acidities= []
        relative_electro_affs=[]

        df = calculate_properties_dataframe.loc[calculate_properties_dataframe['Regid'] == i]
        for j in df.columns:
            if 'anion' in j:
                acidities[j] = float(df[j].values)
            elif 'bromine' in j:
                elec_affs[j] = float(df[j].values)

        for m, c in acidities.items():
            if c == (min(acidities.values())):
                relative_acidities.append(24.044391262487466/c) if 24.044391262487466/c not in relative_acidities else relative_acidities
                position=list(m.split('_')[0])

                if len(relative_acidities) == 1:
                    break

        for f,l in elec_affs.items():
            if l == max(elec_affs.values()):
                relative_electro_affs.append(l/183.647244829844) if l/183.647244829844 not in relative_electro_affs else relative_electro_affs
                position = list(f.split('_')[0])


        if len(relative_acidities) != len(relative_electro_affs):
            print('There is a problem with: ', i, 'please check!')

        plot=ax1.scatter(relative_electro_affs, relative_acidities, label=i, marker='x', zorder=2)

    plt.xlabel('Relative Nucleophilicity')
    plt.ylabel('Relative Acidity')
    plt.title(plot_title)
    plt.xlim(0.4, 1.3)
    plt.ylim(0.2, 2.5)
    plt.savefig('dj1/dj1_nothighlighted.png', dpi=(1200))

    plt.gca().add_patch(Rectangle((0, 1), 2, 11, color='salmon'))
    plt.gca().add_patch(Rectangle((1, 0), 1, 5, color='lightblue'))
    plt.gca().add_patch(Rectangle((1, 1), 1, 11, color='mediumpurple'))
    #ax1.legend(bbox_to_anchor=(1.05,1,2,1))


    plt.savefig(outname, dpi=(1200))


def plot_3Dactivation_map(calculate_properties_dataframe, outname, animate):
    fig = plt.figure()
    ax1 = plt.subplot(111, projection='3d')
    num_colours = 250

    cm = plt.get_cmap('gist_rainbow')
    cNorm = colors.Normalize(vmin=0, vmax=num_colours - 3)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    ax1.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(num_colours)])
    regid=list(set(calculate_properties_dataframe.Regid))
    acidities = {}
    elec_affs = {}
    relative_acidities = []
    relative_electro_affs = []
    relative_selectivity = []
    regids=[]

    for i in regid:
        A = []
        Ea = []
        regids.append(i)
        df = calculate_properties_dataframe.loc[calculate_properties_dataframe['Regid'] == i]
        for j in df.columns:
            if 'anion' in j:
                acidities[j] = float(df[j].values)
            elif 'bromine' in j:
                elec_affs[j] = float(df[j].values)

        for m, c in acidities.items():
            if c == (min(acidities.values())):
                    A.append(c)
                    relative_acidities.append(24.044391262487466/c) if 24.044391262487466/c not in relative_acidities else relative_acidities
                    position=list(m.split('_')[0])
                    #for a,b in elec_affs.items():
                     #   if position[0] in a:
                      #      relative_electro_affs.append(b/183.651714) if b/183.651714 not in relative_electro_affs else relative_electro_affs


        for f,l in elec_affs.items():
            if l == max(elec_affs.values()):
                Ea.append(l)
                relative_electro_affs.append(l/183.647244829844) if l/183.647244829844not in relative_electro_affs else relative_electro_affs
                position = list(f.split('_')[0])

                #if len(relative_electro_affs) == 2:
                 #   for w, e in acidities.items():
                  #      if position[0] in w:
                   #         relative_acidities.append(e / 24.038503) if e/24.038503 not in relative_acidities else relative_acidities


        if len(relative_acidities) != len(relative_electro_affs):
            print('There is a problem with: ', i, 'please check!')

        relative_selectivity.append((Ea[0]/A[0])/(183.651714/24.038503))
    print(relative_selectivity)

    if animate == 'yes':
        plot_animation_3D(relative_electro_affs, relative_acidities, relative_selectivity, regids, outname)
    else:
        for i, j in zip(relative_electro_affs, relative_acidities):
            if i < 1.2:
                if j > 0:
                    for w,x,y,z in zip(relative_electro_affs, relative_acidities, relative_selectivity, regids):
                        plot=ax1.scatter(w, x, y, label=z, marker='x', zorder=2)
        plt.xlabel('Relative Electrophile Affinity')
        plt.ylabel('Relative Acidity')
        ax1.set_zlabel('Relative Differentiation')
        plt.title('Heterocycle properties relative to Pyrazolo[1,5-a]pyrimidine')
        ax1.legend(bbox_to_anchor=(1.05,1,2,1))
        #ax1.set_xlim3d(0.0, 1.5)
        #ax1.set_ylim3d(0.7, 1.5)
        #ax1.set_zlim3d(0,3)
        plt.autoscale(enable=True, axis='both',tight=None)

        plt.savefig(outname, dpi=(600))


def plot_static_3D_plots(calculate_properties_dataframe, outname):
    print('fill me in to make a static plot')

def plot_moving_3D_plots(calculate_properties_dataframe, outname):
    xdata=[]
    ydata=[]
    zdata=[]
    regids=[]

    fig = plt.figure()
    ax1 = plt.subplot(111, projection='3d')
    num_colours = 250
    cm = plt.get_cmap('gist_rainbow')
    cNorm = colors.Normalize(vmin=0, vmax=num_colours - 3)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    ax1.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(num_colours)])

    regid = list(set(calculate_properties_dataframe.Regid))

    for i in regid:
        regids.append(i)
        acidities= {}
        elec_affs= {}

        df = calculate_properties_dataframe.loc[calculate_properties_dataframe['Regid'] == i]
        for j in df.columns:
            if 'anion' in j:
                acidities[j] = float(df[j].values)
            elif 'bromine' in j:
                elec_affs[j] = float(df[j].values)

        sort_acidity=sorted(x for x in acidities.values() if math.isnan(x) == False)
        sort_electrophile_affinity=sorted(y for y in elec_affs.values() if math.isnan(y) == False)
        if len(sort_acidity) > 1:
            delta_acidities=sort_acidity[1] - sort_acidity[0]
            delta_electro_aff=sort_electrophile_affinity[-1] - sort_electrophile_affinity[-2]
            elec_minus_aciditiy= sort_electrophile_affinity[-1]-sort_acidity[0]

            xdata.append(delta_acidities)
            ydata.append(delta_electro_aff)
            zdata.append(elec_minus_aciditiy)


    plot_animation_3D(xdata,ydata,zdata,regids, outname)


def plot_binary_activation_map(calculate_properties_dataframe, successful_activations, selectivity_targets,
                               exploration_targets, outname, limit=1.21, write=True):
    '''
    This function plots the same activation map as in the plots.plot_activation_map function, but also draw the boundary
    between acidic and nucleophilic actiations depending on the value of limit (this should be determined experimentally)

    :param calculate_properties_dataframe: Output of the property_calculate.calculate_properties function or
                                        workflow.logfile_analysis_workflow function. Labelled as "outname"_fulldata.csv
    :param successful_activations: A list of Regids that have been tested in the lab and agree with predictions
    :param selectivity_targets: A list of Regids that have been chosen to investigate the selectivity metric
    :param exploration_targets: A list of Regids that have been chosen to fill out unexplored areas of the map
    :param outname: the name the graph will be saved under.
    :param limit: Where the boundary should be between acidic and nucleophilic activation.
    :param write: Bool, saves graph as outname.png if True.
    :return:
    '''

    fig1, ax1=plt.subplots()
    regid=list(set(calculate_properties_dataframe.Regid))

    for i in regid:
        acidities= {}
        elec_affs= {}
        relative_acidities= []
        relative_electro_affs=[]

        df = calculate_properties_dataframe.loc[calculate_properties_dataframe['Regid'] == i].round(decimals=15)
        for j in df.columns:
            if 'anion' in j:
                acidities[j] = float(df[j].values)
            elif 'bromine' in j:
                elec_affs[j] = float(df[j].values)

        for m, c in acidities.items():
            if c == (min(acidities.values())):
                relative_acidities.append(24.044391262487466/c) if 24.044391262487466/c not in relative_acidities else relative_acidities
                position=list(m.split('_')[0])

                if len(relative_acidities) == 1:
                    break

        for f, l in elec_affs.items():
            if l == max(elec_affs.values()):
                relative_electro_affs.append(l/183.64724482984454) if l/183.64724482984454 not in relative_electro_affs else relative_electro_affs
                position = list(f.split('_')[0])

        for j in successful_activations:
            if i == j:
                ax1.scatter(relative_electro_affs, relative_acidities,
                            edgecolors='green', facecolors='none', zorder=4, s=140)

            else:

                if relative_acidities[0] and relative_electro_affs[0] == 1:
                    ax1.scatter(relative_electro_affs[0], relative_acidities[0],
                                marker='x', color='green', zorder=4)

        if len(relative_acidities) == 0:
            print('no ch bond', i)

        elif relative_acidities[0]/relative_electro_affs[0] <= limit:
            ax1.scatter(relative_electro_affs, relative_acidities,
                        marker='x', color='blue', zorder=2)

        elif relative_acidities[0]/relative_electro_affs[0] > limit:
            ax1.scatter(relative_electro_affs, relative_acidities,
                        marker='x', color='red', zorder=2)

        for k in selectivity_targets:
            if i == k:
                ax1.scatter(relative_electro_affs, relative_acidities,
                            edgecolors='black', facecolors='none', zorder=3, s=140)

        for l in exploration_targets:
            if i == l:
                ax1.scatter(relative_electro_affs, relative_acidities,
                            edgecolors='black', facecolors='none', zorder=3, s=140)


    x = np.linspace(0, 100.0, 500)
    y = limit*x

    boundary_leg = mlines.Line2D([], [], color='black', linestyle='--',
                          label='Boundary')

    Ea_cross = mlines.Line2D([], [], color='blue', marker='x', linestyle='None', markersize=7,
                          label='Predicted Electrophilic Activation')

    A_cross = mlines.Line2D([], [], color='red', marker='x', linestyle='None', markersize=7,
                          label='Predicted Acidic Activation')

    G_cross = mlines.Line2D([], [], color='green', marker='x', linestyle='None', markersize=7,
                          label='Heterocyclic Standard')

    G_ring = mlines.Line2D([], [], markeredgecolor='green', markerfacecolor='None', marker='o', linestyle='None',
                           markersize=7, label='Experimentally Observed')

    O_ring = mlines.Line2D([], [], markeredgecolor='black', markerfacecolor='None', marker='o', linestyle='None',
                           markersize=7, label='Exploration Targets')

    #P_ring = mlines.Line2D([], [], markeredgecolor='white', markerfacecolor='None', marker='o', linestyle='None',
     #                      markersize=4, label='Exploration Targets')

    plt.plot(x, y, label='Boundary', color='black', linestyle='--', zorder=3)

    plt.xlabel('Relative Electrophile Affinity')
    plt.ylabel('Relative Acidity')
    plt.title('DJ1 Activation Map')
    plt.gca().add_patch(Polygon([[0.1, (0.1*limit)], [3.0, (3.0*limit)], [0.9, 11]], facecolor='salmon', edgecolor='salmon', zorder=1, closed=True))
    plt.gca().add_patch(Polygon([[0.1, (0.1*limit)], [3.0, (3.0*limit)], [0.9, -11]], color='lightblue', zorder=1))
    plt.xlim(0.4, 1.3)
    plt.ylim(0.2, 2.5)
    #plt.legend(handles=[Ea_cross, A_cross, G_cross, G_ring, O_ring], prop={'size': 7.8}, bbox_to_anchor=(0.997, 0.997), loc=1, borderaxespad=0)

    if write:
        plt.savefig(outname, dpi=(1200))

def plot_all_maps(dj1_fulldata_rel, dj2_full_data_rel, dj3_fulldata_rel, outname):
    '''
    plots all datasets onto one activation map
    :param dj1_fulldata_rel:
    :param dj2_full_data_rel:
    :param dj3_fulldata_rel:
    :param outname:
    :return:
    '''
    plt.clf()
    fig1, ax1 = plt.subplots()

    dj1_fulldata_rel['set'] = 'DJ1'
    dj2_full_data_rel['set'] = 'DJ2'
    dj3_fulldata_rel['set'] = 'DJ3'

    color_dict = {'DJ1': '#00B050', 'DJ2': '#FF00FF', 'DJ3': '#ED7D31'}

    exc = exclude_structures_prop_dfs(dj1_fulldata_rel['Regid'], dj2_full_data_rel, dj3_fulldata_rel, write=False)

    mega_df = pd.concat([dj1_fulldata_rel, exc[0], exc[1]], axis=0)

    for i in mega_df['Regid']:
        zorder = 2
        acidities = []
        elec_affs = []

        df = mega_df.loc[mega_df['Regid'] == i]

        for j in df.columns:

            if 'anion' in j:
                acidities.append(float(df[j].values))

            elif 'bromine' in j:
                elec_affs.append(float(df[j].values))

        if df.iloc[0]['set'] == 'DJ1':
            zorder=3

        elif df.iloc[0]['set'] == 'DJ3':
            zorder = 3

        if max(acidities) == 0:
            pass
        else:
            plt.plot(max(elec_affs), max(acidities),  marker='x', color=color_dict.get(df.iloc[0]['set']), label=df.iloc[0]['set'], zorder=zorder)


    plt.xlabel('Relative Electrophile Affinity')
    plt.ylabel('Relative Acidity')
    plt.title('Highlighted Heteroaromatic Activation Map')
    #plt.gca().add_patch(Rectangle((0, 1), 2, 11, color='darksalmon'))
    #plt.gca().add_patch(Rectangle((1, 0), 1, 5, color='lightblue'))
    #plt.gca().add_patch(Rectangle((1, 1), 1, 11, color='mediumpurple'))

    DJ1 = mlines.Line2D([], [], color='#00B050', marker='x', linestyle='None', markersize=7,
                             label='DJ1')

    DJ2 = mlines.Line2D([], [], color='#FF00FF', marker='x', linestyle='None', markersize=7,
                            label='DJ2')

    DJ3 = mlines.Line2D([], [], color='#ED7D31', marker='x', linestyle='None', markersize=7,
                            label='DJ3')
    plt.legend(handles=[DJ1, DJ2, DJ3], prop={'size': 7.8}, bbox_to_anchor=(1,1))
    plt.xlim(0.4, 1.3)
    plt.ylim(0.1, 10)

    plt.savefig(outname + '.png', dpi=(600))


def plot_acidity_vs_nucleophilicity(dj1, dj2, dj3, outname, write=True):
    plt.clf()
    fig1, ax1 = plt.subplots()

    dj1['set'] = 'DJ1'
    dj2['set'] = 'DJ2'
    dj3['set'] = 'DJ3'

    color_dict = {'DJ1': '#00B050', 'DJ2': '#FF00FF', 'DJ3': '#ED7D31'}

    exc = exclude_structures_prop_dfs(dj1['Regid'], dj2, dj3, write=False)

    mega_df = pd.concat([dj1, exc[0], exc[1]], axis=0)

    for i in mega_df['Regid']:
        zorder = 2
        acidities = []
        elec_affs = []

        df = mega_df.loc[mega_df['Regid'] == i]

        for j in df.columns:

            if 'anion' in j:
                acidities.append(float(df[j].values))

            elif 'bromine' in j:
                elec_affs.append(float(df[j].values))

        if df.iloc[0]['set'] == 'DJ1':
            zorder = 3

        elif df.iloc[0]['set'] == 'DJ3':
            zorder = 3

        plt.scatter(elec_affs, acidities, marker='x', color=color_dict.get(df.iloc[0]['set']),
                 label=df.iloc[0]['set'], zorder=zorder)


    plt.xlabel('Electrophile Affinity [Kcal $mol^{-1}]$')
    plt.ylabel('Acidity [Kcal $mol^{-1}]$')
    plt.title('Acidity vs Electrophile Affinity')
    plt.ylim(-35, 75)
    plt.xlim(35, 270)
    #plt.gca().add_patch(Rectangle((0, 1), 2, 11, color='darksalmon'))
    #plt.gca().add_patch(Rectangle((1, 0), 1, 5, color='lightblue'))
    #plt.gca().add_patch(Rectangle((1, 1), 1, 11, color='mediumpurple'))

    DJ1 = mlines.Line2D([], [], color='#00B050', marker='x', linestyle='None', markersize=7,
                             label='DJ1')

    DJ2 = mlines.Line2D([], [], color='#FF00FF', marker='x', linestyle='None', markersize=7,
                            label='DJ2')

    DJ3 = mlines.Line2D([], [], color='#ED7D31', marker='x', linestyle='None', markersize=7,
                            label='DJ3')
    plt.legend(handles=[DJ1, DJ2, DJ3], prop={'size': 7.8})

    plt.savefig(outname + '.png', dpi=(600))


def plot_quadrant_activation_map(relative_properties_dataframe, successful_activations, unsuccessful_activations,
                                 targeted_activations, outname, acidity_limit=0.68, nucleophilicity_limit=0.84,
                                 djr_limit=1.22, write=True):
    plt.clf()
    fig1, ax1 = plt.subplots()

    for i in relative_properties_dataframe['Regid']:
        acidities = []
        elec_affs = []

        df = relative_properties_dataframe.loc[relative_properties_dataframe['Regid'] == i]

        for j in df.columns:

            if 'anion' in j:
                acidities.append(float(df[j].values))

            elif 'bromine' in j:
                elec_affs.append(float(df[j].values))

        max_e = round(max(elec_affs), 2)
        max_a = round(max(acidities), 2)

        if max_e < nucleophilicity_limit and max_a < acidity_limit:
            plt.plot(max_e, max_a, marker='x', color='darkgrey')

        elif max_e >= nucleophilicity_limit and max_a < acidity_limit:
            plt.plot(max_e, max_a, marker='x', color='blue')

        elif max_e < nucleophilicity_limit and max_a > acidity_limit:
            plt.plot(max_e, max_a, marker='x', color='red')

        elif max_e >= nucleophilicity_limit and max_a > acidity_limit:

            if max_a/max_e < djr_limit:
                plt.plot(max_e, max_a, marker='x', color='blue')

            else:
                plt.plot(max_e, max_a, marker='x', color='red')

        for h in successful_activations:
            if i == h:
                ax1.scatter(max_e, max_a,
                            edgecolors='green', facecolors='none', zorder=3, s=140)

        for k in unsuccessful_activations:
            if i == k:
                ax1.scatter(max_e, max_a,
                            edgecolors='goldenrod', facecolors='none', zorder=3, s=140)

        for l in targeted_activations:
            if i == l:
                ax1.scatter(max_e, max_a,
                            edgecolors='black', facecolors='none', zorder=3, s=140)

    x = np.linspace(nucleophilicity_limit, 100.0, 500)
    y = djr_limit*x
    plt.plot(x, y, color='black', linestyle='--')

    boundary_line = mlines.Line2D([], [], color='black', linestyle='--',
                                 label='Selectivity Metric Boundary')

    Ea_cross = mlines.Line2D([], [], color='blue', marker='x', linestyle='None', markersize=7,
                             label='Predicted Nucleophilic Activation')

    A_cross = mlines.Line2D([], [], color='red', marker='x', linestyle='None', markersize=7,
                            label='Predicted Acidic Activation')

    G_cross = mlines.Line2D([], [], color='darkgrey', marker='x', linestyle='None', markersize=7,
                            label='Predicted No Activation')

    G_ring = mlines.Line2D([], [], markeredgecolor='green', markerfacecolor='None', marker='o', linestyle='None',
                           markersize=10, label='Experiment Matches Prediction')

    B_ring = mlines.Line2D([], [], markeredgecolor='black', markerfacecolor='None', marker='o', linestyle='None',
                           markersize=10, label='Targeted Substrates')

    gold_ring = mlines.Line2D([], [], markeredgecolor='goldenrod', markerfacecolor='None', marker='o', linestyle='None',
                           markersize=10, label='Experiment Contradicts Prediction')



    plt.xlabel('Relative Nucleophilicity')
    plt.ylabel('Relative Acidity')
    plt.title('DJ1 Activation Map')
    plt.gca().add_patch(Rectangle((0, 0), nucleophilicity_limit, acidity_limit, color='lightgrey'))
    plt.gca().add_patch(Rectangle((0, acidity_limit), nucleophilicity_limit, 11, color='salmon'))
    plt.gca().add_patch(Rectangle((nucleophilicity_limit, 0), 11, acidity_limit, color='lightblue'))
    plt.gca().add_patch(Rectangle((nucleophilicity_limit, acidity_limit), 11, 11, color='mediumpurple'))

    plt.xlim(0.4, 1.3)
    plt.ylim(0.2, 2.5)
    plt.legend(handles=[Ea_cross, A_cross, G_cross, G_ring, #gold_ring,
                        B_ring, boundary_line], prop={'size': 7.8}, bbox_to_anchor=(0.455, 0.997),
               loc=1, borderaxespad=0)

    plt.savefig(outname + '.png', dpi=(1200))


def dataset_count(dict_of_datasets, outname):
    '''

    :param dict_of_datasets: should be a dictionary of the datasets you want ot count as atom_dfs
    :return: writed bar chart of counts for acidities and nucleophilicities
    '''
    plt.clf()
    acidities = []
    elec_affs = []

    ind = np.arange(3)
    width = 0.17

    for label, data in dict_of_datasets.items():
        h_df = data[data['typestr'] == 'H']
        a = h_df[h_df['shift'] !=0]
        acidities.append(len(a))

        c_df = data[data['typestr'] == 'C']
        ea = c_df[c_df['shift'] != 0]
        elec_affs.append(len(ea))

        print(label, 'contains: ', len(a), 'positions')

    ax, fig = plt.subplots()
    bar1 = plt.bar(ind, acidities, width, color='red', align='center', label='Acidities')
    bar2 = plt.bar(ind+width, elec_affs, width, color='blue', align='center', label='Electrophile Affinities')

    plt.ylabel('Count')
    plt.xticks(ind + width, labels=('DJ1', 'DJ2', 'DJ3'))
    plt.xlabel('Data Set')
    plt.title('Count of Heterocyclic Properties')
    plt.legend(title='Property', bbox_to_anchor=(0.25, 1), prop={'size': 8})
    plt.savefig(outname + '.png')


def plot_djr_barchart(list_of_structures, rel_prop_df, outname):

    NMR_yield_A = [60, -62, 51, 57, 69, -62, -0.5, -6]
    NMR_yield_B = [48, 45, 57, 0.5, 21, -0.5, -0.5, -50]
    NMR_yield_C = [0.5, 42, 68, 48, 63, -90, -43, -87]
    x_ticks = []
    exp_select = []

    for i in list_of_structures:

        df = rel_prop_df[rel_prop_df.Regid == i]

        acidity = -100000
        nuc = -100000

        for j in df.columns:
            if 'anion' in j:
                df_acidity = float(df[j].values)
                if df_acidity > acidity:
                    acidity = round(df_acidity, 2)

            elif 'bromine' in j:
                df_nuc = float(df[j].values)
                if df_nuc > nuc:
                    nuc = round(df_nuc, 2)

        comp = round(acidity/nuc, 2)
        x_ticks.append(comp)

        if comp > 1.22:
            exp_select.append(-1)

        else:
            exp_select.append(+1)
        print(i, comp)
    #print(comp_select)
    #print(exp_select)

    exp_select.sort(reverse=True)
    x_ticks.sort()

    print(exp_select)
    print(x_ticks)

    ind = np.arange(8)
    width = 0.17

    bar1 = plt.bar(ind, NMR_yield_A, width=width, align='center', label='Conditions A', color=np.where(np.array(NMR_yield_A)>0, '#1AC2EA', '#FE2712'))
    plt.bar(ind+width, NMR_yield_B, width=width, align='center', label='Conditions B', color=np.where(np.array(NMR_yield_B)>0, '#1121F3', '#DB1717'))
    plt.bar(ind+width*2, NMR_yield_C, width=width, align='center', label='Conditions C', color=np.where(np.array(NMR_yield_C)>0, '#050557', '#950F0F'))
    plt.axhline(0, color='grey', lw=0.5)
    plt.axvline(4.7, color='black', linestyle='--')
    plt.ylim(-90, 90)
    plt.xticks(ind+width, labels=x_ticks)
    plt.xlabel('Selectivity Metric')
    plt.ylabel('% NMR Yield')
    plt.title('Experimental Selectivity vs Selectivity Metric')
    yticks = plt.yticks()
    plt.yticks(yticks[0], labels=[str(abs(x)) for x in yticks[0]])

    legA1 = plt.bar(0, 0, 0, color='#1AC2EA')
    legA2 = plt.bar(0, 0, 0, color='#FE2712')

    legB1 = plt.bar(0, 0, 0, color='#1121F3')
    legB2 = plt.bar(0, 0, 0, color='#DB1717')

    legC1 = plt.bar(0, 0, 0, color='#050557')
    legC2 = plt.bar(0, 0, 0, color='#950F0F')

    plt.legend([(legA1, legA2), (legB1, legB2), (legC1, legC2)], ['Conditions A', 'Conditions B', 'Conditions C'],
               numpoints=3, handler_map={tuple: HandlerTuple(ndivide=None)})
    plt.savefig(outname + '.png', dpi=(600))


def deltaG_deltaE_compare(deltaE_fulldata, deltaG_fulldata, outname):

    plt.clf()
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    num_colours = 29

    cm = plt.get_cmap('hsv')
    cNorm = colors.Normalize(vmin=0, vmax=num_colours - 3)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    ax1.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(num_colours)])
    ax2.set_prop_cycle(color=[scalarMap.to_rgba(i) for i in range(num_colours)])

    regids = deltaG_fulldata['Regid']

    for i in regids:

        deltaE_data = deltaE_fulldata[deltaE_fulldata['Regid'] == i]
        deltaG_data = deltaG_fulldata[deltaG_fulldata['Regid'] == i]

        deltaE_ac = []
        deltaG_ac = []

        deltaE_ea = []
        deltaG_ea = []

        for column in deltaG_data.columns:
            if 'anion' in column:
                deltaE_ac.extend(list(deltaE_data[column].values))
                deltaG_ac.extend(list(deltaG_data[column].values))

            elif 'bromine' in column:
                deltaE_ea.extend(list(deltaE_data[column].values))
                deltaG_ea.extend(list(deltaG_data[column].values))


        ax1.plot(deltaG_ac, deltaE_ac, label=i, marker='x')
        ax2.plot(deltaG_ea, deltaE_ea, label=i, marker='x')

    ax1.set_xlabel('Delta G')
    ax1.set_ylabel('Delta E')
    ax1.set_title('Delta E vs Delta G for DJ1 Acidities')
    fig1.savefig(outname + '_ACIDITIES.png', dpi=(600))

    ax2.set_xlabel('Delta G')
    ax2.set_ylabel('Delta E')
    ax2.set_title('Delta E vs Delta G for DJ1 Nucleophilicities')

    ax1.legend(ncol=2, loc='upper left')
    fig1.savefig(outname + '_ACIDITIES__With_LEG.png', dpi=(600))
    fig2.savefig(outname + '_NUCLEOPHILICITIES.png', dpi=(600))
