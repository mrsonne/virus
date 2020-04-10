import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from .solver import get_kIplus, get_kIminus

def table_str(header, rows, title):
    fstrs = '{:<10} {:^10} {:^10} {:^5}'
    fstrn = '{:<10} {:^10} {:^10} {:^5}'

    lines = [fstrn.format(*row) for row in rows]
    lines.insert(0, fstrs.format(*header))


    sepchr = '-'
    length = max([len(line) for line in lines])
    sep = length*sepchr
    lines.insert(0, sep)
    lines.insert(2, sep)
    lines.append(sep)
    lines.insert(0, title)
    return '\n'.join(lines)


def par_table(population, n_infected_init, infections_at_tau, pars):

    kIplus = get_kIplus(pars['E'], pars['p_t'])
    kIminus = get_kIminus(infections_at_tau, pars['tau'])


    fstr = '{:20} {:7} {}'
    ostrs = []
    ostrs.append(fstr.format('Encounters',  pars['E'], '/day'))
    ostrs.append(fstr.format('Population', population, ''))
    ostrs.append(fstr.format('Ventilators', pars['ventilator_capacity'], ''))
    ostrs.append(fstr.format('Infected at day 0', n_infected_init, ''))
    ostrs.append(fstr.format('Infection time τ', pars['tau'], 'day'))

    fstr = '{:20} {:7.2f} {}'
    ostrs.append(fstr.format('k_I+', kIplus, '/day'))
    ostrs.append(fstr.format('k_I-', kIminus, '/day'))
    ostrs.append(fstr.format('-k_I+ / k_I-', -kIplus / kIminus, ''))

    fstr = '{:20} {:7.1f} %'
    ostrs.append(fstr.format('Infections at τ (%)', infections_at_tau*100))
    ostrs.append(fstr.format('p_t', pars['p_t']*100))
    ostrs.append(fstr.format('p_h', pars['p_h']*100))
    ostrs.append(fstr.format('p_d', pars['p_d']*100))
    ostrs.append(fstr.format('p_v', pars['p_v']*100))
    ostrs.append( fstr.format('p_d,nv', pars['p_dnv']*100))

    ostrs.insert(0, 'Parameters')
    max_length = max([len(ostr) for ostr in ostrs])
    ostrs.insert(0, '-'*max_length) 
    ostrs.insert(2, '-'*max_length) 
    ostrs.append('-'*max_length) 

    return '\n'.join(ostrs)


def plot(times, infected, hospitalized,
         ventilator, recovered, dead,
         ventilator_capacity=None,
         show_recovered=False,
         title=''):

    max_hospitalized = np.max(hospitalized)
    max_ventilator = np.max(ventilator)
    max_recovered = np.max(recovered)
    max_infected = np.max(infected)
    total_deaths = dead[-1]

    fig, ax = plt.subplots(figsize=(15, 8))

    ax.plot(times, infected, 'b-', linewidth=2,
            label='Infected (max={:5.1e})'.format(max_infected))
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


def cplot(x, y, z, x_nom, y_nom, xlabel, ylabel, title=''):
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    grid_x, grid_y = np.meshgrid(xi, yi)
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

    fig, ax = plt.subplots(figsize=(8, 8))
    cs = ax.contour(grid_x, grid_y, grid_z, colors='black')
    ax.scatter(x_nom, y_nom, s=50, c='none', edgecolors='black')
    font_size = 16
    ax.clabel(cs, inline=1, fontsize=font_size, fmt='%1.0f')
    ax.set_xlabel(xlabel, fontsize=font_size)
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    plt.title(title, fontsize=font_size*1.25)
    plt.show()


def ua_plot(xvals, yvals, ypercentiles, xnames, yname, ventilator_capacity, pct_above_maxvent, n_bins=25, title=''):
    font_size = 16
    nattrs = xvals.shape[0]
    fig = plt.figure(figsize=(14, 7), constrained_layout=True)
    gs = fig.add_gridspec(ncols=3, nrows=nattrs)
    axs = []
    for row in range(nattrs):
        axs.append(fig.add_subplot(gs[row, 0]))

    ax_big = fig.add_subplot(gs[:, 1:])

    for ax, _xvals, name in zip(axs, xvals, xnames):
       ax.hist(_xvals, bins=n_bins, density=True, align='mid')
       ax.set_xlabel(name, fontsize=font_size)
       ax.set_ylabel('Density', fontsize=font_size)
       ax.tick_params(axis='x', labelsize=font_size)
       ax.tick_params(axis='y', labelsize=font_size)


    y_avg = np.average(yvals)
    y_p50 = np.percentile(yvals, 50)
    ax_big.hist(yvals, bins=n_bins, density=True, align='mid')
    xlim = ax_big.get_xlim()
    ax_big.axvspan(ventilator_capacity, xlim[1], alpha=0.3, color='tomato',
                   label='Capacity exceeded (p={:<4.1f} %)'.format(pct_above_maxvent))
    ax_big.axvline(y_avg, color='black',
                   label='Average: {:<4.0f}'.format(y_avg))
    # ax_big.axvline(y_p50, color='black',
    #                label='Median: {:<4.0f}'.format(y_p50))

    ax_big.set_xlim(xlim)
    ax_big.set_xlabel(yname, fontsize=font_size)
    ax_big.set_ylabel('Density', fontsize=font_size)
    ax_big.legend(fontsize=font_size) 
    ax_big.tick_params(axis='x', labelsize=font_size)
    ax_big.tick_params(axis='y', labelsize=font_size)
    plt.show()


def ua_timeseries(times, values):
    font_size = 16
    fig, ax = plt.subplots(figsize=(15, 8))
    avg = np.average(values, axis=1)
    avg_max = max(avg)
    p50 = np.percentile(values, 50, axis=1)
    p95 = np.percentile(values, 95, axis=1)
    p05 = np.percentile(values, 5, axis=1)
    ax.plot(times, values, color="black", alpha=0.1, linewidth=2)
    ax.plot(times, values[:, 0], color="seagreen", alpha=0.2, linewidth=2, label='Samples')
    ax.fill_between(times, p05, p95, color="orange", alpha=0.5, linewidth=1, label='CI 90 %', zorder=9999)
    # ax.plot(times, p50, color="black", alpha=1, linewidth=2, label='Median', zorder=10000)
    ax.plot(times, avg, color="black", alpha=1, linewidth=2,
            label='Average (max={:5.1e})'.format(avg_max), zorder=10000)
    ax.legend(fontsize=font_size, loc='upper left')
    ax.set_xlabel('Time (day)', fontsize=font_size)
    ax.set_ylabel('Ventilators required', fontsize=font_size)
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    plt.show()
