import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import animation

def plot_3D_map_with_exchange_E(relative_directed_activation_data, rel_prop_data, outname, animate=False):

    print('HARD CODED LIMIT FOR REL Gexch on line 337 and 355, SORT THIS OUT!!!')
    structures = []
    rel_acidity = []
    rel_nucleophilicity = []
    rel_Gexch = []
    red = ['S5854', 'S5857', 'S6018']
    blue = ['S1863', 'S1869', 'S5853', 'S1873', 'S5703', 'S1875', 'S5704']
    purple = ['S1080', 'S1072']
    grey = ['S7959', 'S54', 'S8063', 'S49', 'S50']
    color = []
    alpha = []


    for i in relative_directed_activation_data['Structure']:

        structures.append(i.split('-')[1]) if i.split('-')[1] not in structures else structures

    for regid in structures:

        if regid in red:
            color.append('red')
            alpha.append(1)

        elif regid in blue:
            color.append('blue')
            alpha.append(1)

        elif regid in purple:
            color.append('purple')
            alpha.append(1)

        elif regid in grey:
            color.append('black')
            alpha.append(1)

        else:
            color.append('grey')
            if animate:
                alpha.append(0.2)
            else:
                alpha.append(0.4)

        id_G = relative_directed_activation_data[relative_directed_activation_data['Structure'].str.endswith(regid)]

        id_rel_props = rel_prop_data[rel_prop_data['Regid'].str.fullmatch(regid)]

        het_anions = []
        het_bromines = []
        het_exchG = id_G['ExchG_relative'].values

        for j in id_rel_props.columns:
            if 'acidity' in j:
                het_anions.append(float(id_rel_props[j].values))

            if 'electrophile_affinity' in j:
                het_bromines.append(float(id_rel_props[j].values))


        rel_acidity.append(max(het_anions))
        rel_nucleophilicity.append(max(het_bromines))
        rel_Gexch.append(max(het_exchG))

        #print(regid, max(het_anions), max(het_bromines), max(het_exchG))

    data = pd.DataFrame({'Structures': structures, 'Rel_acidity': rel_acidity,
                         'Rel_nucleophilicity': rel_nucleophilicity, 'Rel_donicity': rel_Gexch})
    fig = plt.figure()
    ax1 = plt.subplot(111, projection='3d')

    if animate:
        plot_animation_3D(rel_acidity, rel_nucleophilicity, rel_Gexch, color, alpha, structures, outname)

    else:
        for a, b, c, d, e in zip(rel_acidity, rel_nucleophilicity, rel_Gexch, color, alpha):

            if c < 15:
                ax1.scatter(a, b, c, color=d, marker='x', zorder=2, alpha=e)
        plt.xlabel('Relative Acidity')
        plt.ylabel('Relative Nucleophilicity')
        ax1.set_zlabel('Relative Gexch')
        plt.title('3D Activation Map - Exchange Energies')
        plt.autoscale(enable=True, axis='both', tight=None)

        plt.savefig(outname + '.png', dpi=(1200))

    with open('acidity_nucleophilicity_donicity_3D_data.csv', 'w') as f:
        print(data.to_csv(), file=f)



def plot_animation_3D(xdata, ydata, zdata, colors, alpha, structures, outname):

    fig = plt.figure()
    ax1 = plt.subplot(111, projection='3d')

    def init():
        for x, y, z, color, alph, id in zip(xdata, ydata, zdata, colors, alpha, structures):
            if z < 15:
                plot = ax1.scatter(x, y, z, color=color, marker='x', zorder=2, alpha=alph)

            #if color == 'blue':
             #   print(id, x, y, z)
              #  ax1.text(x,y,z, id, color=color)

        plt.xlabel('Relative Acidity')
        plt.ylabel('Relative Nucleophilicity')
        ax1.set_zlabel('$Relative E_{exch}$')
        plt.title('3D Activation Map - Exchange Energies')
        plt.autoscale(enable=True, axis='both', tight=None)
        return fig,

    def animate(i):
        ax1.view_init(elev=10, azim=i * 4)
        return fig,

    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=90, interval=75, blit=True)
    fn = outname
    ani.save(fn + '.gif', writer='pillow', fps=100)