import numpy as np
from scipy.integrate import solve_ivp
from scipy import optimize
from scipy.stats import erlang


# Indices in state array
idx_for_state = {'healthy' : 0,
                 'recovered': 1,
                 'dead': 2,
                 'infected_1': 3, # holds the total number of infected individuals
                }
NSTATES = len(idx_for_state)
N_INFECTED_STAGES = 1


def get_y0(population, n_infected_init, infected_stages=1):
    """
    TODO: subtract n_infected_init
    """
    global NSTATES
    NSTATES += infected_stages - 1
    y0 = [0]*NSTATES
    y0[idx_for_state['healthy']] = max(population - n_infected_init, 0)

    # Add stages to index mapper and state array
    global N_INFECTED_STAGES
    N_INFECTED_STAGES = infected_stages
    if infected_stages > 0:
        idx_stage1 = idx_for_state['infected_1']
        for istage in range(1, infected_stages):
            idx_for_state['infected_{}'.format(istage + 1)] = idx_stage1 + istage
    else:
        raise BaseException('Number of infection stages musr be 1 or greater')

    y0[idx_for_state['infected_1']] = min(n_infected_init, population)

    return y0


def get_I_tot(y):
    """
    Sum infection stages to give total number of infected individuals
    """
    idx_1 = idx_for_state['infected_1']
    return sum(y[idx_1 : idx_1 + N_INFECTED_STAGES])


def get_last_iidx():
    return idx_for_state['infected_1'] + N_INFECTED_STAGES - 1


def get_first_iidx():
    return idx_for_state['infected_1']


def get_I_last(y):
    """
    Get last infection stage
    """
    return y[get_last_iidx()]


def get_ventilators_required(yinfected, ytot, p_h, p_v):
    hospitalized = yinfected*ytot*p_h
    return hospitalized*p_v


def dydt(t, y, kIminus, kIplus, p_r, p_d, p_dnv,
         ventilator_capacity, p_h, p_v, ytot):
    """
    Model 1: no collateral effect of using all ventilators. Keep ICU capacity
    """
    dydt = [0]*NSTATES
    I_tot = get_I_tot(y)
    I_last = get_I_last(y)
    # Ventilator can be required in all stages since the data that we have are pooled for all hospitalized
    ventilators_required = get_ventilators_required(I_tot,
                                                    ytot,
                                                    p_h,
                                                    p_v)

    ventilators_missing = ventilators_required - ventilator_capacity
    ventilators_missing = max(ventilators_missing, 0) / ytot

    # All stages can transmit virus
    dIplus_dt = kIplus*I_tot*y[idx_for_state['healthy']]
    dydt[idx_for_state['healthy']] = -dIplus_dt

    # You only die or recover from the last stage... A bit odd when ventilators_missing is calculated based on I_tot
    # dIminus_dt is zero if more ventilators are required than there are individuals in the last stage... 
    dIminus_dt = kIminus*( max(I_last - ventilators_missing, 0 ))
    dydt[idx_for_state['dead']] = p_d * dIminus_dt
    # If I_last < ventilators_missing only I_last will die
    dydt[idx_for_state['dead']] += p_dnv * kIminus * min(ventilators_missing, I_last)
    # dydt[idx_for_state['dead']] += p_dnv * kIminus * ventilators_missing
    dydt[idx_for_state['recovered']] = dIminus_dt - dydt[idx_for_state['dead']]
    
    dydt[get_first_iidx()] += dIplus_dt
    # dydt[get_last_iidx()] -=  dIminus_dt
    dydt[get_last_iidx()] -=  dydt[idx_for_state['dead']] +  dydt[idx_for_state['recovered']]

    # print(list(range(get_first_iidx(), get_last_iidx())))
    # print(list(range(get_first_iidx() + 1, get_last_iidx() + 1)))

    # Exclude last since it's already handled
    for iidx in range(get_first_iidx(), get_last_iidx()):
        dydt[iidx] -= kIminus*y[iidx]

    # Exclude first since it's already handled
    for iidx in range(get_first_iidx() + 1, get_last_iidx() + 1):
        dydt[iidx] += kIminus*y[iidx - 1]

    return dydt



def exceed_ventilator_capacity(t, y, kIminus,
                               kIplus, p_r, p_d, p_dnv,
                               ventilator_capacity, p_h, p_v, ytot):
    ventilators_required = get_ventilators_required(get_I_tot(y),
                                                    ytot,
                                                    p_h,
                                                    p_v)
    return ventilator_capacity - ventilators_required


def get_rate_Iplus(encounters_per_day, p_t):
    """
    Equivalent for infected to meet healthy and vice versa
    p_S = 0.5 (p_H = 0.5)
        S         H
    S   0.25      0.25
    H   0.25      0.25
    """
    return 2*encounters_per_day*p_t


def obj(rate, tau, sf_at_tau, k):
    return erlang.sf(tau, a=k, scale=1./rate) - sf_at_tau


def get_rate_Iminus(tau, sf_at_tau, k=1):
    """
    Determine the rate based on the survival function at tau
    """
    
    if sf_at_tau == 'tau_is_mean':
        return k/tau
    else:
        r_exp = -np.log(sf_at_tau) / tau
        if k == 1:
            return r_exp
        else:
            r_est = r_exp*k
            sol = optimize.root_scalar(obj, bracket=[1e-9, 5*r_est],
                                    x0=r_est,
                                    args=(tau, sf_at_tau, k),
                                    method='brentq')
            return sol.root



def extract_time_series(sol, ytot, p_h, p_v, ventilator_capacity):
    infected = get_I_tot(sol.y)*ytot
    recovered = sol.y[idx_for_state['recovered']]*ytot
    dead = sol.y[idx_for_state['dead']]*ytot
    hospitalized = infected*p_h
    ventilator = np.minimum(hospitalized*p_v, ventilator_capacity)
    ventilators_required = hospitalized*p_v
    return sol.t, infected, recovered, dead, hospitalized, ventilator, ventilators_required


def solve(y0,
          survival_at_tau,
          E,
          p_t,
          tau,
          p_d,
          p_dnv,
          p_h,
          p_v,
          tspan,
          ventilator_capacity,
          n_time_eval=1000):

    kIplus = get_rate_Iplus(E, p_t)
    kIminus = get_rate_Iminus(tau, survival_at_tau, N_INFECTED_STAGES)

    ytot = np.sum( y0 )

    # normalize to get probabilities
    y0 /= ytot
    exceed_ventilator_capacity.terminal = False
    exceed_ventilator_capacity.direction = -1

    p_r = 1. - p_d
    times = np.linspace(*tspan, n_time_eval)

    # Legacy from time where solve_ivp could be called multiple times with changing parameters
    ts, infected, recovered, dead, hospitalized, ventilator, ventilators_required = [], [], [], [], [], [], []

    sol = solve_ivp(dydt, tspan, y0,
                    args=(kIminus, kIplus, p_r,
                        p_d,
                        p_dnv,
                        ventilator_capacity,
                        p_h,
                        p_v, ytot),
                    t_eval=times,
                    method='Radau',
                    rtol=1e-6, # tighten tolerance to avoid negative population count
                    events=exceed_ventilator_capacity
                )

    # Check closure
    # print('Closure OK', np.allclose(np.sum(sol.y, axis=0) - 1, 0))

    (_ts, _infected, _recovered,
    _dead, _hospitalized, _ventilator,
    _ventilators_required) = extract_time_series(sol,
                                                 ytot,
                                                 p_h,
                                                 p_v,
                                                 ventilator_capacity)

    ts = np.append(ts, _ts)
    infected = np.append(infected, _infected)
    recovered = np.append(recovered, _recovered)
    dead = np.append(dead, _dead)
    hospitalized = np.append(hospitalized, _hospitalized)
    ventilator = np.append(ventilator, _ventilator)
    ventilators_required = np.append(ventilators_required, _ventilators_required)

    return {'times': ts,
            'infected': infected,
            'hospitalized': hospitalized,
            'ventilator': ventilator,
            'recovered': recovered,
            'dead': dead,
            'ventilators_required': ventilators_required}
