import numpy as np
from scipy.integrate import solve_ivp

# Indices in state array
idx_for_state = {'healthy' : 0,
                 'recovered': 1,
                 'dead': 2,
                 'infected': 3, # holds the total number of infected individuals
                }
nstates = len(idx_for_state)
N_INFECTED_STAGES = 1


def get_y0(population, n_infected_init, infected_stages=1):
    """
    TODO: subtract n_infected_init
    """
    y0 = [0]*nstates
    y0[idx_for_state['healthy']] = population

    # Add stages to index mapper and state array
    global N_INFECTED_STAGES
    N_INFECTED_STAGES = infected_stages
    if infected_stages > 0:
        for istage in range(infected_stages):
            idx_for_state['infected_{}'.format(istage + 1)] = idx_for_state['infected'] + istage

        del idx_for_state['infected']

    else:
        raise BaseException('Number of infection stages musr be 1 or greater')

    y0[idx_for_state['infected_1']] = n_infected_init

    return y0


def get_Itot(y):
    idx_1 = idx_for_state['infected_1']
    return sum(y[idx_1 : idx_1 + N_INFECTED_STAGES])



def get_ventilators_required(yinfected, ytot, p_h, p_v):
    hospitalized = yinfected*ytot*p_h
    return hospitalized*p_v


def dydt(t, y, kIminus, kIplus, p_r, p_d, p_dnv,
         ventilator_capacity, p_h, p_v, ytot):
    """
    Model 1: no collateral effect of using all ventilators. Keep ICU capacity
    """
    dydt = [0]*nstates
    ventilators_required = get_ventilators_required(y[idx_for_state['infected_1']],
                                                    ytot,
                                                    p_h,
                                                    p_v)

    ventilators_missing = ventilators_required - ventilator_capacity
    ventilators_missing = max(ventilators_missing, 0) / ytot

    # All stages can transmit virus
    I_tot = get_Itot(y)
    dIplus_dt = kIplus*I_tot*y[idx_for_state['healthy']]
    dydt[idx_for_state['healthy']] = -dIplus_dt

    dIminus_dt = kIminus*( y[idx_for_state['infected_1']] - ventilators_missing )
    dydt[idx_for_state['dead']] = p_d * dIminus_dt
    dydt[idx_for_state['dead']] += p_dnv * kIminus * ventilators_missing
    dydt[idx_for_state['recovered']] = dIminus_dt - dydt[idx_for_state['dead']]
    
    dydt[idx_for_state['infected_1']] = dIplus_dt - dIminus_dt
    return dydt



def exceed_ventilator_capacity(t, y, kIminus,
                               kIplus, p_r, p_d, p_dnv,
                               ventilator_capacity, p_h, p_v, ytot):
    ventilators_required = get_ventilators_required(y[idx_for_state['infected_1']],
                                                    ytot,
                                                    p_h,
                                                    p_v)
    return ventilator_capacity - ventilators_required


def get_kIplus(encounters_per_day, p_t):
    """
    Equivalent for infected to meet healthy and vice versa
    p_S = 0.5 (p_H = 0.5)
        S         H
    S   0.25      0.25
    H   0.25      0.25
    """
    return 2*encounters_per_day*p_t



def get_kIminus(infections_at_tau, tau):
    return -np.log(infections_at_tau)/tau


def extract_time_series(sol, ytot, p_h, p_v, ventilator_capacity):
    infected = sol.y[idx_for_state['infected_1']]*ytot
    recovered = sol.y[idx_for_state['recovered']]*ytot
    dead = sol.y[idx_for_state['dead']]*ytot
    hospitalized = sol.y[idx_for_state['infected_1']]*ytot*p_h
    ventilator = np.minimum(hospitalized*p_v, ventilator_capacity)
    ventilators_required = hospitalized*p_v
    return sol.t, infected, recovered, dead, hospitalized, ventilator, ventilators_required


def solve(y0,
          infections_at_tau,
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

    kIplus = get_kIplus(E, p_t)
    kIminus = get_kIminus(infections_at_tau, tau)

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
                    events=exceed_ventilator_capacity
                )

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
