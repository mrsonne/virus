import itertools
import numpy as np
from scipy.stats import percentileofscore
from .data import COUNTRYFUN_FOR_COUNTRYID
from .solver import solve, get_kIplus
from . import render

def frc_to_pct(val):
    return 100*val


def dummy_transform(val):
    return val


def get_last(time_series):
    return time_series[-1]


def get_pars(virus_id, country_id, encounters_per_day, tspan):


    try:
        pars, population = COUNTRYFUN_FOR_COUNTRYID[country_id](virus_id)
    except KeyError:
        err_str = 'Unknown country ID "{}". Available IDs: {}'
        raise KeyError(err_str.format(country_id,
                                      ', '.join(list(COUNTRYFUN_FOR_COUNTRYID.keys()))))


    if encounters_per_day is None:
        pars['E'] = -kIminus/get_kIplus(1., p_t) # gives constant infected count
    else:
        pars['E'] = encounters_per_day

    if tspan:
        pars['tspan'] = tspan

    return pars, population


def virus(virus_id, country_id, encounters_per_day=None,
                show_recovered=False, tspan=None):
    """
    Virus simulation
    """

    # 80 % infections disappeared after infection time
    infections_at_tau = 0.2

    pars, population = get_pars(virus_id, country_id,
                                encounters_per_day,
                                tspan)


    n_infected_init = 5
    
    print(render.par_table(population, n_infected_init, infections_at_tau, pars))

    # TODO: subtract n_infected_init
    y0 = [population, n_infected_init, 0, 0]

    tss = solve(y0, infections_at_tau, **pars)

    title = '{} ($E$={}/day)'.format(virus_id, pars['E'])
    render.plot(tss["times"], 
                tss["infected"],
                tss["hospitalized"],
                tss["ventilator"],
                tss["recovered"],
                tss["dead"],
                show_recovered=show_recovered, title=title)



def contour(virus_id, country_id, par1, par2, response,
            encounters_per_day=None, tspan=None,
            nsteps=25):
    """Grid
    """

    # 80 % infections disappeared after infection time
    infections_at_tau = 0.2

    pars, population = get_pars(virus_id, country_id,
                                encounters_per_day,
                                tspan)

    parstr1 = par1['name']
    parstr2 = par2['name']
    try:
        ftrans1 = par1["transform"]
    except KeyError:
        ftrans1 = dummy_transform

    try:
        ftrans2 = par2["transform"]
    except KeyError:
        ftrans2 = dummy_transform

    # Save nominel values
    par1_nom = pars[parstr1]
    par2_nom = pars[parstr2]

    n_infected_init = 5
    y0 = [population, n_infected_init, 0, 0]

    
    pars1 = np.linspace(par1['min'], par1['max'], nsteps)
    pars2 = np.linspace(par2['min'], par2['max'], nsteps)
    pars1_grid, pars2_grid = [], []
    resp_grid = []
    for i, (p1, p2) in enumerate(itertools.product(pars1, pars2)):
        pars[parstr1] = p1
        pars[parstr2] = p2
        time_series = solve(y0, infections_at_tau, **pars)[response["name"]]
        resp_grid.append(response["transform"](time_series))
        pars1_grid.append(p1)
        pars2_grid.append(p2)

    title = '{} ({}, $E$={}/day)'.format(response["title"], virus_id, pars['E'])
    render.cplot(ftrans1(np.array(pars1_grid)), 
                 ftrans2(np.array(pars2_grid)),
                 resp_grid,
                 ftrans1(par1_nom), 
                 ftrans2(par2_nom), 
                 par1["axlabel"], par2["axlabel"], 
                 title=title)


def ua(virus_id, country_id,
       encounters_per_day,
       nsamples=1000,
       tspan=None,
       infections_at_tau=0.2):
    """Uncertainty analysis
    """
    # Seed so I dont have to change the text for each run
    np.random.seed(0)

    # 80 % infections disappeared after infection time
    infections_at_tau = 0.2

    pars, population = get_pars(virus_id, country_id,
                                encounters_per_day,
                                tspan)


    # Reduce onset time
    n_infected_init = 250
    y0 = [population, n_infected_init, 0, 0]

    xnames = r'$p_{\rm{h}}$ (%)', r'$\tau$ (days)', '$E$ (day\u207B\u00B9)'
    mean = [pars["p_h"] , pars["tau"], pars["E"]]
    cov = [[0.00001, 0. , 0], [0., 4., 0], [0., 0, 16]]
    xvals = np.random.multivariate_normal(mean, cov, nsamples).T

    ventilators = []
    ventilator_series = []
    for p_h, tau, E in xvals.T:
        pars['p_h'] = p_h
        pars['tau'] = tau
        pars['E'] = E

        tss = solve(y0, infections_at_tau, **pars)
        ventilators.append(max(tss["ventilators_required"]))
        ventilator_series.append(tss["ventilators_required"])

    times = tss["times"]

    yname = 'Ventilators required'
    # Convert p_d to percent for plotting
    xvals[0, :] *= 100
    percentiles = np.percentile(ventilators, q=(5, 50, 95))
    pct_above_maxvent = 100. - percentileofscore(ventilators, pars["ventilator_capacity"])
    title = '{}, $E$={}'.format(virus_id, E)
    render.ua_plot(xvals, ventilators, percentiles, xnames, yname, pars["ventilator_capacity"], pct_above_maxvent,
                   title=title)
    return times, np.array(ventilator_series)
