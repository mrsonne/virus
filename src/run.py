import itertools
import numpy as np
from scipy.stats import percentileofscore
from .data import COUNTRYFUN_FOR_COUNTRYID
from .solver import solve, get_kIplus
from . import render

def get_pars(virus_id, country_id, encounters_per_day, tspan):


    try:
        pars, population = COUNTRYFUN_FOR_COUNTRYID[country_id](virus_id)
    except KeyError:
        err_str = 'Unknown country ID "{}". Available IDs: {}'
        raise KeyError(err_str.format(country_id,
                                      ', '.join(list(COUNTRYFUN_FOR_COUNTRYID.keys()))))


    if encounters_per_day is None:
        pars['E'] = -kIminus/get_kIplus(1., p_t) # constant sick count
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


    n_sick_init = 5
    
    print(render.par_table(population, n_sick_init, infections_at_tau, pars))

    # TODO: subtract n_sick_init
    y0 = [population, n_sick_init, 0, 0]

    (times, 
    sick,
    hospitalized,
    ventilator,
    recovered,
    dead,
    _) = solve(y0, infections_at_tau, **pars)

    title = '{} ($E$={}/day)'.format(virus_id, pars['E'])
    render.plot(times, sick, hospitalized, ventilator, recovered, dead,
                show_recovered=show_recovered, title=title)




def contour(virus_id, country_id, encounters_per_day=None,
                tspan=None):
    """Grid
    """

    # 80 % infections disappeared after infection time
    infections_at_tau = 0.2

    pars, population = get_pars(virus_id, country_id,
                                encounters_per_day,
                                tspan)

    # Save nominel values
    p_d_nom = pars['p_d'] 
    tau_nom = pars['tau']

    n_sick_init = 5
    y0 = [population, n_sick_init, 0, 0]

    deads = []
    nsteps = 25
    _ps_d = np.linspace(0.001, 0.01, nsteps)
    _taus = np.linspace(6, 21, nsteps)
    taus = []
    ps_d = []
    for i, (p_d, tau) in enumerate(itertools.product(_ps_d, _taus)):
        pars['p_d'] = p_d
        pars['tau'] = tau
        # Use only last element of "deads" time series
        dead_count = solve(y0, infections_at_tau, **pars)[-2][-1]
        deads.append(dead_count)
        taus.append(tau)
        ps_d.append(p_d)

    title = 'Deads ({}, $E$={}/day)'.format(virus_id, pars['E'])
    render.cplot(np.array(ps_d)*100, np.array(taus), deads, p_d_nom*100, tau_nom, title=title)


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
    n_sick_init = 250
    y0 = [population, n_sick_init, 0, 0]

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

        (times,
        _,
        _,
        _,
        _,
        _,
        ventilators_required) = solve(y0, infections_at_tau, **pars)
        ventilators.append(max(ventilators_required))
        ventilator_series.append(ventilators_required)

    yname = 'Ventilators required'
    # Convert p_d to percent for plotting
    xvals[0, :] *= 100
    percentiles = np.percentile(ventilators, q=(5, 50, 95))
    pct_above_maxvent = 100. - percentileofscore(ventilators, pars["ventilator_capacity"])
    title = '{}, $E$={}'.format(virus_id, E)
    render.ua_plot(xvals, ventilators, percentiles, xnames, yname, pars["ventilator_capacity"], pct_above_maxvent,
                   title=title)
    return times, np.array(ventilator_series)
