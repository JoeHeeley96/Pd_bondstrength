import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import math
from mpl_toolkits import mplot3d
from matplotlib import animation
import numpy as np
import matplotlib.lines as mlines
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
    ani.save(fn + '.gif', writer='pillow', fps=1000 / 50)


def plot_activation_map(calculate_properties_dataframe, outname):

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
                #for a,b in elec_affs.items():
                 #   if position[0] in a:
                  #      relative_electro_affs.append(b/183.651714) if b/183.651714 not in relative_electro_affs else relative_electro_affs

                if len(relative_acidities) == 1:
                    break

        for f,l in elec_affs.items():
            if l == max(elec_affs.values()):
                relative_electro_affs.append(l/183.647244829844) if l/183.647244829844 not in relative_electro_affs else relative_electro_affs
                position = list(f.split('_')[0])
                #if len(relative_electro_affs) == 2:
                 #   for w, e in acidities.items():
                  #      if position[0] in w:
                   #         relative_acidities.append(e / 24.038503) if e/24.038503 not in relative_acidities else relative_acidities


        if len(relative_acidities) != len(relative_electro_affs):
            print('There is a problem with: ', i, 'please check!')
        if relative_acidities[0] >=1:
            if relative_electro_affs[0] >=1:
                print(i)

        plot=ax1.scatter(relative_electro_affs, relative_acidities, label=i, marker='x', zorder=2)
    plt.xlabel('Relative Electrophile Affinity')
    plt.ylabel('Relative Acidity')
    plt.title('Heterocycle properties relative to Pyrazolo[1,5-a]pyrimidine')
    plt.gca().add_patch(Rectangle((0, 1), 2, 4, color='salmon'))
    plt.gca().add_patch(Rectangle((1, 0), 1, 5, color='lightblue'))
    plt.gca().add_patch(Rectangle((1, 1), 1, 4, color='mediumpurple'))
    ax1.legend(bbox_to_anchor=(1.05,1,2,1))
    plt.xlim(0.6, 1.25)
    plt.ylim(0, 2.6)

    plt.savefig(outname, dpi=(600))

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

def plot_binary_activation_map(calculate_properties_dataframe, successful_activations, outname,
                               limit=1.2, write=True):
    '''
    This function plots the same activation map as in the plots.plot_activation_map function, but also draw the boundary
    between acidic and nucleophilic actiations depending on the value of limit (this should be determined experimentally)

    :param calculate_properties_dataframe: Output of the property_calculate.calculate_properties function or
                                        workflow.logfile_analysis_workflow function. Labelled as "outname"_fulldata.csv
    :param successful_activations: A list of Regids that have been tested in the lab and agree with predictions
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
                            edgecolors='green', facecolors='none', zorder=3, s=140)
                if i =='S1875':
                    print(relative_acidities[0]/relative_electro_affs[0])
            else:

                if relative_acidities[0] and relative_electro_affs[0] == 1:
                    ax1.scatter(relative_electro_affs[0], relative_acidities[0],
                                marker='x', color='green', zorder=3)

                elif relative_acidities[0]/relative_electro_affs[0] <= limit:
                    ax1.scatter(relative_electro_affs, relative_acidities,
                                marker='x', color='blue', zorder=2)

                elif relative_acidities[0]/relative_electro_affs[0] > limit:
                    ax1.scatter(relative_electro_affs, relative_acidities,
                                marker='x', color='red', zorder=2)


    x = np.linspace(0, 100.0, 500)
    y = limit*x

    Ea_cross = mlines.Line2D([], [], color='blue', marker='x', linestyle='None', markersize=7,
                          label='Predicted Electrophilic Activation')

    A_cross = mlines.Line2D([], [], color='red', marker='x', linestyle='None', markersize=7,
                          label='Predicted Acidic Activation')

    G_cross = mlines.Line2D([], [], color='green', marker='x', linestyle='None', markersize=7,
                          label='Heterocyclic Standard')

    G_ring = mlines.Line2D([], [], markeredgecolor='green', markerfacecolor='None', marker='o', linestyle='None',
                           markersize=7, label='Experimentally Observed')

    plt.plot(x, y, label='Boundary', color='black', linestyle='--', zorder=3)

    plt.xlabel('Relative Electrophile Affinity')
    plt.ylabel('Relative Acidity')
    plt.title('DJ1 Binary Activation "Heat" Map')
    plt.gca().add_patch(Polygon([[0.1, (0.1*limit)], [3.0, (3.0*limit)], [0.9, 11]], facecolor='salmon',edgecolor='salmon', zorder=1, closed=True))
    plt.gca().add_patch(Polygon([[0.1, (0.1*limit)], [3.0, (3.0*limit)], [0.9, -11]], color='lightblue', zorder=1))
    plt.xlim(0.6, 1.25)
    plt.ylim(0, 2.6)
    plt.legend(handles=[Ea_cross, A_cross, G_cross, G_ring], prop={'size': 7.8}, bbox_to_anchor=(0.997, 0.997), loc=1, borderaxespad=0)

    if write:
        plt.savefig(outname, dpi=(600))