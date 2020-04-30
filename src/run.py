import itertools
import numpy as np
from scipy.stats import percentileofscore
from .data import COUNTRYFUN_FOR_COUNTRYID
from .solver import solve, get_rate_Iplus, get_y0, get_rate_Iminus
from . import render

def frc_to_pct(val):
    return 100*val


def dummy_transform(val):
    return val


def get_last(time_series):
    return time_series[-1]


def get_first(vals):
    return vals[0]


def get_max(time_series):
    return np.max(time_series)


def get_pars(virus_id, country_id, encounters_per_day,
             tspan, survival_at_tau, k, p_t):


    try:
        pars, population = COUNTRYFUN_FOR_COUNTRYID[country_id](virus_id, p_t)
    except KeyError:
        err_str = 'Unknown country ID "{}". Available IDs: {}'
        raise KeyError(err_str.format(country_id,
                                      ', '.join(list(COUNTRYFUN_FOR_COUNTRYID.keys()))))

    if encounters_per_day is None:
        pars['Es'] = get_rate_Iminus(pars['tau'], survival_at_tau, k)/get_rate_Iplus(1., pars['p_t']) # gives constant infected count
    else:
        pars['Es'] = encounters_per_day

    if tspan:
        pars['tspan'] = tspan

    return pars, population


def virus(virus_id, country_id, encounters_per_day=None,
                show_recovered=False, tspan=None,
                survival_at_tau=0.2, k=1, p_t=None):
    """
    Virus simulation
    """

    _encounters_per_day = np.atleast_1d(encounters_per_day)
    pars, population = get_pars(virus_id, country_id,
                                _encounters_per_day,
                                tspan, survival_at_tau, k,
                                p_t)


    n_infected_init = 5

    # DEBUG: see pure survival function
    # n_infected_init = population

    if len(_encounters_per_day) == 1:
        print(render.par_table(population, n_infected_init, survival_at_tau, k, pars))

    y0 = get_y0(population, n_infected_init, k)

    tss = solve(y0, survival_at_tau, **pars)

    title = '{} ($E$={}/day)'.format(virus_id, ','.join([str(E) for E in pars['Es']]))
    render.plot(tss["times"], 
                tss["infected"],
                tss["hospitalized"],
                tss["ventilator"],
                tss["recovered"],
                tss["dead"],
                show_recovered=show_recovered, title=title)



def contour(virus_id, country_id, par1, par2, response,
            encounters_per_day=None, tspan=None,
            survival_at_tau=0.2, k=1,
            nsteps=25):
    """Grid
    """
    _encounters_per_day = np.atleast_1d(encounters_per_day)
    p_t = None
    pars, population = get_pars(virus_id, country_id,
                                _encounters_per_day,
                                tspan, survival_at_tau, k,
                                p_t)

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
    y0 = get_y0(population, n_infected_init, k)

    
    pars1 = np.linspace(par1['min'], par1['max'], nsteps)
    pars2 = np.linspace(par2['min'], par2['max'], nsteps)
    pars1_grid, pars2_grid = [], []
    resp_grid = []
    for i, (p1, p2) in enumerate(itertools.product(pars1, pars2)):
        pars[parstr1] = p1
        pars[parstr2] = p2
        time_series = solve(y0, survival_at_tau, **pars)[response["name"]]
        resp_grid.append(response["transform"](time_series))
        pars1_grid.append(p1)
        pars2_grid.append(p2)

    title = '{} ({}, $E$={}/day)'.format(response["title"], virus_id, ','.join([str(E) for E in pars['Es']]))
    render.cplot(ftrans1(np.array(pars1_grid)), 
                 ftrans2(np.array(pars2_grid)),
                 resp_grid,
                 ftrans1(par1_nom), 
                 ftrans2(par2_nom), 
                 par1["axlabel"], par2["axlabel"], 
                 title=title)


def ua(virus_id, country_id,
       encounters_per_day,
       parobjs, response,
       nsamples=1000,
       tspan=None,
       survival_at_tau=0.2, k=1):
    """Uncertainty analysis
    """
    _encounters_per_day = np.atleast_1d(encounters_per_day)
    p_t = None
    pars, population = get_pars(virus_id, country_id,
                                _encounters_per_day,
                                tspan, survival_at_tau, k,
                                p_t)


    response_ftrans = response["transform"]

    # Reduce onset time by increasing initially infected individuals
    n_infected_init = 5
    y0 = get_y0(population, n_infected_init, k)

    # Prepare sampling
    xvals_for_parname = {}
    for parobj in parobjs:
        if "mean" in parobj:
            mean = parobj["mean"]
        else: 
            mean = pars[parobj["name"]]

        try:
            cov = np.diag([parobj['std']**2 for _ in mean])
            xvals_for_parname[parobj["name"]] = np.random.multivariate_normal(mean, cov, nsamples)
        except (ValueError, TypeError):
           xvals_for_parname[parobj["name"]] = np.random.normal(mean, parobj['std'], nsamples).T

    # Initialize transformed response and full time series
    response_trns = []
    response_ts = []

    n_time_eval = 1000

    # Loop over parameter sets and solve
    for idxsmpl in range(nsamples):

        # Tranfer parameters
        for parobj in parobjs:
            pars[parobj["name"]] = xvals_for_parname[parobj["name"]][idxsmpl]
                
        tss = solve(y0, survival_at_tau, **pars, n_time_eval=n_time_eval)
        if len(tss[response["name"]]) < n_time_eval: continue
        response_trns.append(response_ftrans(tss[response["name"]]))
        response_ts.append(tss[response["name"]])

    times = tss["times"]
    yname = response["title"]

    # Transform
    for parobj in parobjs:
        fun = parobj["transform"] if "transform" in parobj else dummy_transform
        xvals_for_parname[parobj["name"]] = fun(xvals_for_parname[parobj["name"]])

    # Plotting
    title = '{}'.format(virus_id)
    render.ua_plot(xvals_for_parname, response_trns, parobjs, yname,
                   response_ts=response_ts,
                   pars=pars,
                   title=title,
                   threshold=pars["ventilator_capacity"])
    return times, np.stack(response_ts, axis=0)
