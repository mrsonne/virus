import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

def plot(times, sick, hospitalized,
         ventilator, recovered, dead,
         vimpgrid_xmpgrid_ympgrid_zntilator_capacity=None,
         show_recovered=False,
         title=''):

    max_hospitalized = np.max(hospitalized)
    max_ventilator = np.max(ventilator)
    max_recovered = np.max(recovered)
    max_sick = np.max(sick)
    total_deaths = dead[-1]

    fig, ax = plt.subplots(figsize=(15, 8))

    ax.plot(times, sick, 'b-', linewidth=2,
            label='Sick (max={:5.1e})'.format(max_sick))
    ax.plot(times, hospitalized, 'b--', linewidth=2,
            label='Hospitalized (max={:5.1e})'.format(max_hospitalized))
    ax.plot(times, ventilator, 'b:', linewidth=2,
            label='Ventilator (max={:5.1e})'.format(max_ventilator))
    ax.plot(times, dead, 'r-', linewidth=2,
            label='Dead (total={:5.1e})'.format(total_deaths))

    if show_recovered:
        alpha = 1
        ylim = None
    else:
        ylim = ax.get_ylim()
        alpha = 0.2

    ax.plot(times, recovered, 'g-', linewidth=2, alpha=alpha, label='Recovered (max={:5.1e})'.format(max_recovered))
    if ylim: ax.set_ylim(ylim)


    if ventilator_capacity is not None:
        ax.hlines(ventilator_capacity, xmin=min(t_eval), xmax=max(t_eval), color='black', linewidth=2,
                label='Ventilator capacity')


    # TODO: show parameters in plot

#     ax.legend(loc='upper left')
    font_size = 16
    ax.set_xlabel('Time (day)', fontsize=font_size)
    ax.set_ylabel('Individuals', fontsize=font_size)
    ax.legend(fontsize=font_size) 
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    plt.title(title, fontsize=font_size*1.25)
    plt.show()


def cplot(x, y, z, title=''):
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    grid_x, grid_y = np.meshgrid(xi, yi)
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

    fig, ax = plt.subplots(figsize=(8, 8))
    cs = ax.contour(grid_x, grid_y, grid_z, colors='black')
    font_size = 16
    ax.clabel(cs, inline=1, fontsize=font_size, fmt='%1.0f')
    ax.set_xlabel('$p_d$ (%)', fontsize=font_size)
    ax.set_ylabel('τ (days)', fontsize=font_size)
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    plt.title(title, fontsize=font_size*1.25)
    plt.show()
