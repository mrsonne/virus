import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from scipy.interpolate import griddata
from scipy.stats import percentileofscore

from .solver import get_kIplus, get_rate_Iminus

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


def par_table(population, n_infected_init, infections_at_tau, k, pars):

    kIplus = get_kIplus(pars['E'], pars['p_t'])
    kIminus = get_rate_Iminus(pars['tau'], infections_at_tau, k)


    fstr = '{:18} {:7} {}'
    ostrs = []
    ostrs.append(fstr.format('Encounters',  pars['E'], '/day'))
    ostrs.append(fstr.format('Population', population, ''))
    ostrs.append(fstr.format('Vent. capacity', pars['ventilator_capacity'], ''))
    ostrs.append(fstr.format('Infected at day 0', n_infected_init, ''))
    ostrs.append(fstr.format('Infctn time (τ)', pars['tau'], 'day'))
    ostrs.append(fstr.format('Infctns at τ', infections_at_tau*100, '%'))
    ostrs.append(fstr.format('Infctn stages (k)', k, ''))

    fstr = '{:18} {:8.3f} {}'
    ostrs.append(fstr.format('Mean infctn time', k/kIminus, 'day'))
    ostrs.append(fstr.format('r_I+', kIplus, '/day'))
    ostrs.append(fstr.format('r_I-', kIminus, '/day'))
    ostrs.append(fstr.format('r_I+ / r_I-', kIplus / kIminus, ''))

    fstr = '{:18} {:8.1f} %'
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
         ventilated, recovered, dead,
         ventilator_capacity=None,
         show_recovered=False,
         title=''):

    max_hospitalized = np.max(hospitalized)
    max_ventilated = np.max(ventilated)
    max_recovered = np.max(recovered)
    max_infected = np.max(infected)
    total_deaths = dead[-1]

    fig, ax = plt.subplots(figsize=(15, 8))

    ax.plot(times, infected, 'b-', linewidth=2,
            label='Infected (max={:5.1e})'.format(max_infected))
    ax.plot(times, hospitalized, 'b--', linewidth=2,
            label='Hospitalized (max={:5.1e})'.format(max_hospitalized))
    ax.plot(times, ventilated, 'b:', linewidth=2,
            label='Ventilated (max={:5.1e})'.format(max_ventilated))
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


def ua_plot(xvals, yvals, xnames, yname,
            response_ts=None, pars=None, # specify these together
            n_bins=25, title='', threshold=None):
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
    n_zero_to_threshold = 20

    if threshold is None:
        _threshold = np.max(yvals) 
    else:
        _threshold = threshold
        pct_above_threshold = 100. - percentileofscore(yvals, threshold)
        label = 'Capacity exceeded (p={:<4.1f} %)'.format(pct_above_threshold)

    width = float(_threshold)/n_zero_to_threshold
    bins = np.arange(0, np.max(yvals), step=width)
    hist_values, bin_edges = np.histogram(yvals, bins=bins, density=True)
    width = bin_edges[1] - bin_edges[0]
    ax_big.bar(bin_edges[:-1][:n_zero_to_threshold], hist_values[:n_zero_to_threshold],
               width=width, align='edge', label='Capacity OK')

    if threshold is not None:
        ax_big.bar(bin_edges[:-1][n_zero_to_threshold:], hist_values[n_zero_to_threshold:],
                color='tomato', width=width, align='edge', label=label)

    ax_big.axvline(y_p50, color='black', linestyle='-',
                   label='Median: {:<4.0f}'.format(y_p50))
    ax_big.axvline(y_avg, color='black', linestyle=':',
                   label='Average: {:<4.0f}'.format(y_avg))
    
    ax_big.set_xlabel(yname, fontsize=font_size)
    ax_big.set_ylabel('Density', fontsize=font_size)
    ax_big.legend(fontsize=font_size) 
    ax_big.tick_params(axis='x', labelsize=font_size)
    ax_big.tick_params(axis='y', labelsize=font_size)
    plt.show()


def ua_timeseries(times, values):

    # This dictionary defines the colormap
    cdict = {'red':  ((0.0, 0.0, 0.0),   # no red at 0
                    (0.5, 1.0, 1.0),   # orange at 0.5
                    (1.0, 1.0, 1.0)),  # set to 0.8 so its not too bright at 1

            'green': ((0.0, 0.8, 0.8),   # set to 0.8 so its not too bright at 0
                    (0.5, 0.4, 0.4),   # orange at 0.5
                    (1.0, 0.0, 0.0)),  # no green at 1

            'blue':  ((0.0, 0.0, 0.0),   # no blue at 0
                    (0.5, 0.0, 0.0),   # orange at 0.5
                    (1.0, 0.0, 0.0))   # no blue at 1
        }

    cmap = colors.LinearSegmentedColormap('my_colormap', cdict, 100)

    font_size = 16
    fig, ax = plt.subplots(figsize=(15, 8))
    avg = np.average(values, axis=1)
    avg_max = max(avg)
    pvals = [0, 20, 40, 60, 80, 95]
    percentiles = np.percentile(values, pvals, axis=1)
    p50 = np.percentile(values, 50, axis=1)
    p50_max = max(p50)
    for lower, upper, p_low, p_up in zip(percentiles[:-1], percentiles[1:], pvals[:-1], pvals[1:]):
        p_mean = 0.5*(p_low + p_up)
        ax.fill_between(times, lower, upper, color=cmap(p_mean/100),
                        alpha=0.5, linewidth=2, zorder=10000,
                        label='Percentile {:.0f} % - {:.0f} %'.format(p_low, p_up))

    # Samples
    ax.plot(times, values, color="gray", alpha=0.1, linewidth=1)
    # make a single label
    ax.plot(times, values[:, 0], color="gray", alpha=0.2, linewidth=1, label='Samples') 

    ax.plot(times, p50, color="black", alpha=1, linewidth=2,
            label='Median (max={:5.1e})'.format(p50_max), zorder=10000)
    ax.plot(times, avg, color="black", linestyle=':', alpha=1, linewidth=2,
            label='Average (max={:5.1e})'.format(avg_max), zorder=10000)
    ax.legend(fontsize=font_size, loc='upper right')
    ylim_max = np.max(percentiles[-1])
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.set_ylim(0, ylim_max*1.05)
    ax.set_xlabel('Time (day)', fontsize=font_size)
    ax.set_ylabel('Ventilators required', fontsize=font_size)
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    plt.show()
