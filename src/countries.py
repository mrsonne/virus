import numpy as np
from .data import PARSFUN_FOR_VIRUSID
from .solver import solve, get_kIplus
from . import render

def _get_denmark():
    population = int(5e6)
    ventilator_capacity = 1000
    return population, ventilator_capacity


COUNTRYFUN_FOR_COUNTRYID = {'denmark': _get_denmark}


def run_country(virus_id, country_id, encounters_per_day=None,
                show_recovered=False, tspan=None):
    """
    Danish data in Table 2 (~1000 deads per year)
    https://www.ssi.dk/sygdomme-beredskab-og-forskning/sygdomsovervaagning/i/influenzasaesonen---opgoerelse-over-sygdomsforekomst-2018-19
    """


    try:
        (_tspan,
        p_transmision,
        p_hospitalized_to_ventilator,
        p_sick_to_hospitalized,
        p_sick_to_dead,
        p_sick_to_dead_no_ventilator,
        tau
        ) = PARSFUN_FOR_VIRUSID[virus_id]()
    except KeyError:
        err_str = 'Unknown virus ID "{}". Available IDs: {}'
        raise KeyError(err_str.format(virus_id,
                                      ', '.join(list(PARSFUN_FOR_VIRUSID.keys()))))

    if not tspan:
        tspan = _tspan
    # Always 80 % loose infection after infection duration
    infections_at_tau = 0.2
    kIminus = -np.log(infections_at_tau)/tau

    if encounters_per_day is None:
        _encounters_per_day = kIminus/get_kIplus(1., p_transmision) # constant sick count
    else:
        _encounters_per_day = encounters_per_day

    try:
        population, ventilator_capacity = COUNTRYFUN_FOR_COUNTRYID[country_id]()
    except KeyError:
        err_str = 'Unknown country ID "{}". Available IDs: {}'
        raise KeyError(err_str.format(country_id,
                                      ', '.join(list(COUNTRYFUN_FOR_COUNTRYID.keys()))))


    n_sick_init = 5
    
    kIplus = get_kIplus(_encounters_per_day, p_transmision)

    fstr = '{:20} {:7} {}'
    ostrs = []
    ostrs.append(fstr.format('Encounters per day',  _encounters_per_day, ''))
    ostrs.append(fstr.format('Population', population, ''))
    ostrs.append(fstr.format('Ventilators', ventilator_capacity, ''))
    ostrs.append(fstr.format('Sick at day 0', n_sick_init, ''))
    ostrs.append(fstr.format('Infection time τ', tau, 'day'))

    fstr = '{:20} {:7.2f} {}'
    ostrs.append(fstr.format('k_I+', kIplus, '/day'))
    ostrs.append(fstr.format('k_I-', kIminus, '/day'))
    ostrs.append(fstr.format('k_I+ / k_I-', kIplus / kIminus, ''))

    fstr = '{:20} {:7.1f} %'
    ostrs.append(fstr.format('Infections at τ (%)', infections_at_tau*100))
    ostrs.append(fstr.format('p_t', p_transmision*100))
    ostrs.append(fstr.format('p_h', p_sick_to_hospitalized*100))
    ostrs.append(fstr.format('p_d', p_sick_to_dead*100))
    ostrs.append(fstr.format('p_v', p_hospitalized_to_ventilator*100))
    ostrs.append( fstr.format('p_d,nv', p_sick_to_dead_no_ventilator*100))

    ostrs.insert(0, 'Parameters')
    max_length = max([len(ostr) for ostr in ostrs])
    ostrs.insert(0, '-'*max_length) 
    ostrs.insert(2, '-'*max_length) 
    ostrs.append('-'*max_length) 

    print('\n'.join(ostrs))

    y0 = [population, n_sick_init, 0, 0]

    (times, 
    sick,
    hospitalized,
    ventilator,
    recovered,
    dead) = solve(_encounters_per_day,
                         p_transmision,
                         tau,
                         kIminus,
                         p_sick_to_dead,
                         p_sick_to_dead_no_ventilator,
                         p_sick_to_hospitalized,
                         p_hospitalized_to_ventilator,
                         tspan, 
                         y0, ventilator_capacity)

    title = '{} ($E$={})'.format(virus_id, encounters_per_day)
    render.plot(times, sick, hospitalized, ventilator, recovered, dead,
                show_recovered=show_recovered, title=title)
