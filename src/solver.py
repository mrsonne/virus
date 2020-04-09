import numpy as np
from scipy.integrate import solve_ivp

# Indices in state array
healthy_idx = 0
sick_idx = 1
recovered_idx = 2
dead_idx = 3


def get_ventilators_required(ysick, ytot, p_h, p_v):
    hospitalized = ysick*ytot*p_h
    return hospitalized*p_v


def dydt(t, y, kIminus, kIplus, p_r, p_d, p_dnv,
         ventilator_capacity, p_h, p_v, ytot):
    """
    Model 1: no collateral effect of using all ventilators. Keep ICU capacity
    """
    ventilators_required = get_ventilators_required(y[sick_idx],
                                                    ytot,
                                                    p_h,
                                                    p_v)

    ventilators_missing = ventilators_required - ventilator_capacity
    ventilators_missing = max(ventilators_missing, 0) / ytot

    dIplus_dt = kIplus*y[sick_idx]*y[healthy_idx]
    dhealthy_dt = -dIplus_dt

    dIminus_dt = kIminus*( y[sick_idx] - ventilators_missing )
    ddead_dt = -p_d * dIminus_dt
    ddead_dt -= p_dnv * kIminus * ventilators_missing
    drecovered_dt = -dIminus_dt - ddead_dt
    
    dsick_dt = dIplus_dt + dIminus_dt
    return [dhealthy_dt, dsick_dt, drecovered_dt, ddead_dt]



def exceed_ventilator_capacity(t, y, kIminus,
                               kIplus, p_r, p_d, p_dnv,
                               ventilator_capacity, p_h, p_v, ytot):
    ventilators_required = get_ventilators_required(y[sick_idx],
                                                    ytot,
                                                    p_h,
                                                    p_v)
    return ventilator_capacity - ventilators_required


def get_kIplus(encounters_per_day, p_t):
    """
    Equivalent for sick to meet healthy and vice versa
    p_S = 0.5 (p_H = 0.5)
        S         H
    S   0.25      0.25
    H   0.25      0.25
    """
    return 2*encounters_per_day*p_t



def get_kIminus(infections_at_tau, tau):
    return np.log(infections_at_tau)/tau


def extract_time_series(sol, ytot, p_h, p_v, ventilator_capacity):
    sick = sol.y[sick_idx]*ytot
    recovered = sol.y[recovered_idx]*ytot
    dead = sol.y[dead_idx]*ytot
    hospitalized = sol.y[sick_idx]*ytot*p_h
    ventilator = np.minimum(hospitalized*p_v, ventilator_capacity)
    ventilators_required = hospitalized*p_v
    return sol.t, sick, recovered, dead, hospitalized, ventilator, ventilators_required


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
    ts, sick, recovered, dead, hospitalized, ventilator, ventilators_required = [], [], [], [], [], [], []

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

    (_ts, _sick, _recovered,
    _dead, _hospitalized, _ventilator,
    _ventilators_required) = extract_time_series(sol,
                                                 ytot,
                                                 p_h,
                                                 p_v,
                                                 ventilator_capacity)

    ts = np.append(ts, _ts)
    sick = np.append(sick, _sick)
    recovered = np.append(recovered, _recovered)
    dead = np.append(dead, _dead)
    hospitalized = np.append(hospitalized, _hospitalized)
    ventilator = np.append(ventilator, _ventilator)
    ventilators_required = np.append(ventilators_required, _ventilators_required)

    return {'times': ts,
            'infected': sick,
            'hospitalized': hospitalized,
            'ventilator': ventilator,
            'recovered': recovered,
            'dead': dead,
            'ventilators_required': ventilators_required}
