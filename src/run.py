import itertools
import numpy as np
from scipy.stats import percentileofscore
from .data import COUNTRYFUN_FOR_COUNTRYID
from .solver import solve, get_kIplus
from . import render

def get_pars(virus_id, country_id, encounters_per_day, infections_at_tau, tspan):


    try:
        (population,
        ventilator_capacity,
        _tspan,
        p_t,
        p_v,
        p_h,
        p_d,
        p_dnv,
        tau) = COUNTRYFUN_FOR_COUNTRYID[country_id](virus_id)
    except KeyError:
        err_str = 'Unknown country ID "{}". Available IDs: {}'
        raise KeyError(err_str.format(country_id,
                                      ', '.join(list(COUNTRYFUN_FOR_COUNTRYID.keys()))))

    kIminus = np.log(infections_at_tau)/tau

    if encounters_per_day is None:
        E = -kIminus/get_kIplus(1., p_t) # constant sick count
    else:
        E = encounters_per_day

    if not tspan:
        tspan = _tspan

    return (E, 
            population,
            ventilator_capacity,
            _tspan,
            p_t,
            p_v,
            p_h,
            p_d,
            p_dnv,
            tau, kIminus, tspan)


def virus(virus_id, country_id, encounters_per_day=None,
                show_recovered=False, tspan=None):
    """
    Virus simulation
    """

    # 80 % infections disappeared after infection time
    infections_at_tau = 0.2

    (E, 
    population,
    ventilator_capacity,
    _tspan,
    p_t,
    p_v,
    p_h,
    p_d,
    p_dnv,
    tau,
    kIminus,
    tspan) = get_pars(virus_id, country_id,
                      encounters_per_day,
                      infections_at_tau, tspan)


    n_sick_init = 5
    
    kIplus = get_kIplus(E, p_t)
    print(render.par_table(E, population, ventilator_capacity,
                           n_sick_init, tau, kIplus, kIminus, infections_at_tau,
                           p_t, p_h, p_d, p_v, p_dnv))

    y0 = [population, n_sick_init, 0, 0]

    (times, 
    sick,
    hospitalized,
    ventilator,
    recovered,
    dead,
    _) = solve(kIplus,
               kIminus,
               p_d,
               p_dnv,
               p_h,
               p_v,
               tspan, 
               y0,
               ventilator_capacity)

    title = '{} ($E$={}/day)'.format(virus_id, E)
    render.plot(times, sick, hospitalized, ventilator, recovered, dead,
                show_recovered=show_recovered, title=title)




def contour(virus_id, country_id, encounters_per_day=None,
                tspan=None):
    """Grid
    """

    # 80 % infections disappeared after infection time
    infections_at_tau = 0.2

    (E, 
    population,
    ventilator_capacity,
    _tspan,
    p_t,
    p_v,
    p_h,
    p_d_nom,
    p_dnv,
    tau_nom,
    kIminus,
    tspan) = get_pars(virus_id, country_id,
                      encounters_per_day,
                      infections_at_tau, tspan)


    n_sick_init = 5
    y0 = [population, n_sick_init, 0, 0]

    kIplus = get_kIplus(E, p_t)

    deads = []
    nsteps = 25
    _ps_d = np.linspace(0.001, 0.01, nsteps)
    _taus = np.linspace(6, 21, nsteps)
    taus = []
    ps_d = []
    for i, (p_d, tau) in enumerate(itertools.product(_ps_d, _taus)):
        kIminus = np.log(infections_at_tau)/tau
        (_, 
        _,
        _,
        _,
        _,
        dead,
        _) = solve(kIplus,
                   kIminus,
                   p_d,
                   p_dnv,
                   p_h,
                   p_v,
                   tspan, 
                   y0,
                   ventilator_capacity)
        deads.append(dead[-1])
        taus.append(tau)
        ps_d.append(p_d)

    title = 'Deads ({}, $E$={}/day)'.format(virus_id, E)
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

    (E, 
    population,
    ventilator_capacity,
    _tspan,
    p_t,
    p_v,
    p_h_nom,
    p_d,
    p_dnv,
    tau_nom,
    kIminus,
    tspan) = get_pars(virus_id, country_id,
                      encounters_per_day,
                      infections_at_tau, tspan)

    # Reduce onset time
    n_sick_init = 250
    y0 = [population, n_sick_init, 0, 0]

    xnames = r'$p_{\rm{h}}$ (%)', r'$\tau$ (days)', '$E$ (day\u207B\u00B9)'
    mean = [p_h_nom , tau_nom, E]
    cov = [[0.00001, 0. , 0], [0., 4., 0], [0., 0, 16]]
    xvals = np.random.multivariate_normal(mean, cov, nsamples).T

    ventilators = []
    ventilator_series = []
    for p_h, tau, E in xvals.T:
        kIplus = get_kIplus(E, p_t)
        kIminus = np.log(infections_at_tau)/tau
        (times,
        _,
        _,
        _,
        _,
        _,
        ventilators_required) = solve(kIplus,
                                      kIminus,
                                      p_d,
                                      p_dnv,
                                      p_h,
                                      p_v,
                                      tspan, 
                                      y0, ventilator_capacity)
        ventilators.append(max(ventilators_required))
        ventilator_series.append(ventilators_required)

    yname = 'Ventilators required'
    # Convert p_d to percent for plotting
    xvals[0, :] *= 100
    percentiles = np.percentile(ventilators, q=(5, 50, 95))
    pct_above_maxvent = 100. - percentileofscore(ventilators, ventilator_capacity)
    title = '{}, $E$={}'.format(virus_id, E)
    render.ua_plot(xvals, ventilators, percentiles, xnames, yname, ventilator_capacity, pct_above_maxvent,
                   title=title)
    return times, np.array(ventilator_series)
