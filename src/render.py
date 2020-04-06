import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

def plot(times, sick, hospitalized,
         ventilator, recovered, dead,
         ventilator_capacity=None,
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


def cplot(x, y, z, p_d_nom, tau_nom, title=''):
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    grid_x, grid_y = np.meshgrid(xi, yi)
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

    fig, ax = plt.subplots(figsize=(8, 8))
    cs = ax.contour(grid_x, grid_y, grid_z, colors='black')
    ax.scatter(p_d_nom, tau_nom, s=50, c='none', edgecolors='black')
    font_size = 16
    ax.clabel(cs, inline=1, fontsize=font_size, fmt='%1.0f')
    ax.set_xlabel('$p_d$ (%)', fontsize=font_size)
    ax.set_ylabel('τ (days)', fontsize=font_size)
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    plt.title(title, fontsize=font_size*1.25)
    plt.show()


def ua_plot(xvals, yvals, ypercentiles, xnames, yname, percentiel_name, n_bins=25, title=''):
    font_size = 16
    nattrs = xvals.shape[0]
    fig = plt.figure(figsize=(14, 7), constrained_layout=True)
    gs = fig.add_gridspec(ncols=3, nrows=nattrs)
    axs = []
    for row in range(2):
        axs.append(fig.add_subplot(gs[row, 0]))

    ax_big = fig.add_subplot(gs[:, 1:])

    for ax, _xvals, name in zip(axs, xvals, xnames):
       ax.hist(_xvals, bins=n_bins, density=True, align='mid')
       ax.set_xlabel(name, fontsize=font_size)
       ax.set_ylabel('Density', fontsize=font_size)
       ax.tick_params(axis='x', labelsize=font_size)
       ax.tick_params(axis='y', labelsize=font_size)


    ax_big.axvspan(ypercentiles[0], ypercentiles[-1], alpha=0.3, color='seagreen',
                   label='{} % CI ({:5.1e} to {:5.1e})'.format(percentiel_name, ypercentiles[0], ypercentiles[-1]))
    ax_big.hist(yvals, bins=n_bins, density=True, align='mid')
    ax_big.axvline(ypercentiles[1], color='black',
                   label='Median ({:5.1e})'.format(ypercentiles[1]))
    ax_big.set_xlabel(yname, fontsize=font_size)
    ax_big.set_ylabel('Density', fontsize=font_size)
    ax_big.legend(fontsize=font_size) 
    ax_big.tick_params(axis='x', labelsize=font_size)
    ax_big.tick_params(axis='y', labelsize=font_size)

    plt.show()