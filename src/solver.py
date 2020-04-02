import numpy as np
from scipy.integrate import solve_ivp

# Indices in state array
healthy_idx = 0
sick_idx = 1
recovered_idx = 2
dead_idx = 3

def get_ventilators_required(ysick, ytot, p_sick_to_hospitalized, p_hospitalized_to_ventilator):
    hospitalized = ysick*ytot*p_sick_to_hospitalized
    return hospitalized*p_hospitalized_to_ventilator


def dydt(t, y, kIminus, kIplus, p_sick_to_recovered, p_sick_to_dead, p_sick_to_dead_no_ventilator,
         ventilator_capacity, p_sick_to_hospitalized, p_hospitalized_to_ventilator, ytot):
    """
    Model 1: no collateral effect of using all ventilators. Keep ICU capacity
    """
    ventilators_required = get_ventilators_required(y[sick_idx],
                                                    ytot,
                                                    p_sick_to_hospitalized,
                                                    p_hospitalized_to_ventilator)

    ventilators_missing = ventilators_required - ventilator_capacity
    ventilators_missing = max(ventilators_missing, 0) / ytot

    dIplus_dt = kIplus*y[sick_idx]*y[healthy_idx]
    dhealthy_dt = -dIplus_dt

    dIminus_dt = kIminus*( y[sick_idx] - ventilators_missing )
    ddead_dt = p_sick_to_dead * dIminus_dt
    ddead_dt += p_sick_to_dead_no_ventilator * kIminus * ventilators_missing
    drecovered_dt = dIminus_dt - ddead_dt
    
    dsick_dt = dIplus_dt - dIminus_dt
    return [dhealthy_dt, dsick_dt, drecovered_dt, ddead_dt]



def exceed_ventilator_capacity(t, y, kIminus,
                               kIplus, p_sick_to_recovered, p_sick_to_dead, p_sick_to_dead_no_ventilator,
                               ventilator_capacity, p_sick_to_hospitalized, p_hospitalized_to_ventilator, ytot):
    ventilators_required = get_ventilators_required(y[sick_idx],
                                                    ytot,
                                                    p_sick_to_hospitalized,
                                                    p_hospitalized_to_ventilator)
    return ventilator_capacity - ventilators_required


def get_kIplus(encounters_per_day, p_transmision, tau):
    """
    Equivalent for sick to meet healthy and vice versa
    p_S = 0.5 (p_H = 0.5)
        S         H
    S   0.25      0.25
    H   0.25      0.25
    """
    return 2*encounters_per_day*p_transmision*tau


def extract_time_series(sol, ytot, p_sick_to_hospitalized, p_hospitalized_to_ventilator, ventilator_capacity):
    sick = sol.y[sick_idx]*ytot
    recovered = sol.y[recovered_idx]*ytot
    dead = sol.y[dead_idx]*ytot
    hospitalized = sol.y[sick_idx]*ytot*p_sick_to_hospitalized
    ventilator = np.minimum(hospitalized*p_hospitalized_to_ventilator, ventilator_capacity)
    return sol.t, sick, recovered, dead, hospitalized, ventilator


def solve(encounters_per_day,
          p_transmision,
          tau,
          kIminus,
          p_sick_to_dead,
          p_sick_to_dead_no_ventilator,
          p_sick_to_hospitalized,
          p_hospitalized_to_ventilator,
          tspan, y0, ventilator_capacity,
          n_time_eval=1000):

    kIplus = get_kIplus(encounters_per_day, p_transmision, tau)

    ytot = np.sum( y0 )

    # normalize to get probabilities
    y0 /= ytot
    exceed_ventilator_capacity.terminal = False
    exceed_ventilator_capacity.direction = -1

    p_sick_to_recovered = 1. - p_sick_to_dead
    times = np.linspace(*tspan, n_time_eval)

    run = True
    ts, sick, recovered, dead, hospitalized, ventilator = [], [], [], [], [], []

    sol = solve_ivp(dydt, tspan, y0,
                    args=(kIminus, kIplus, p_sick_to_recovered,
                        p_sick_to_dead,
                        p_sick_to_dead_no_ventilator,
                        ventilator_capacity,
                        p_sick_to_hospitalized,
                        p_hospitalized_to_ventilator, ytot),
                    t_eval=times,
                    method='Radau',
                    events=exceed_ventilator_capacity
                )

    _ts, _sick, _recovered, _dead, _hospitalized, _ventilator = extract_time_series(sol,
                                                                            ytot,
                                                                            p_sick_to_hospitalized,
                                                                            p_hospitalized_to_ventilator,
                                                                            ventilator_capacity)

    ts = np.append(ts, _ts)
    sick = np.append(sick, _sick)
    recovered = np.append(recovered, _recovered)
    dead = np.append(dead, _dead)
    hospitalized = np.append(hospitalized, _hospitalized)
    ventilator = np.append(ventilator, _ventilator)


    return ts, sick, hospitalized, ventilator, recovered, dead
