def get_covid19_parameters():
    # Assumes same hospitalization rate as flu. Data is not reliable for COVID-19 yet. Parameter is not critical
    _, p_transmision, _, p_sick_to_hospitalized, p_sick_to_dead, _, _ = get_flu_parameters()
    tspan = [0, 200]
    p_sick_to_dead = .005 
    # p_sick_to_dead_no_ventilator =  p_sick_to_dead # For testing: no effect of exceeding ventilator capacity

    # https://medium.com/@tomaspueyo/coronavirus-the-hammer-and-the-dance-be9337092b56
    p_sick_to_dead_no_ventilator = .9 # 
    p_hospitalized_to_ventilator = 0.2 # of hospitalized
    tau = 14 #days
    return tspan, p_transmision, p_hospitalized_to_ventilator, p_sick_to_hospitalized, p_sick_to_dead, p_sick_to_dead_no_ventilator, tau


def get_flu_parameters():
    # Get flu data and use as baseline
    p_sick_to_hospitalized, p_sick_to_dead = get_us_flu_data() 
    tspan = [0, 700]
    p_sick_to_dead_no_ventilator = 0.
    p_transmision = 0.0026 # adjusted to fit flu data
    p_hospitalized_to_ventilator = 0. # unknown
    tau = 7 #days
    return tspan, p_transmision, p_hospitalized_to_ventilator, p_sick_to_hospitalized, p_sick_to_dead, p_sick_to_dead_no_ventilator, tau


def get_us_flu_data():
    # US flu data
    # https://www.cdc.gov/flu/about/burden/preliminary-in-season-estimates.htm
    # population = 250.e6
    cases = 45000000. # per year
    p_sick_to_hospitalized = 600000. / cases # ~1.3 %
    p_sick_to_dead = 40000. / cases # ~0.1 %
    return p_sick_to_hospitalized, p_sick_to_dead

PARSFUN_FOR_VIRUSID = dict(flu=get_flu_parameters,
                           covid19=get_covid19_parameters)


def _get_denmark(virus_id):
    """
    Danish flu data in Table 2 (~1000 deads per year)
    https://www.ssi.dk/sygdomme-beredskab-og-forskning/sygdomsovervaagning/i/influenzasaesonen---opgoerelse-over-sygdomsforekomst-2018-19
    """
    population = int(5e6)
    ventilator_capacity = 1000

    try:
        (tspan,
        p_t,
        p_v,
        p_h,
        p_d,
        p_dnv,
        tau
        ) = PARSFUN_FOR_VIRUSID[virus_id]()
    except KeyError:
        err_str = 'Unknown virus ID "{}". Available IDs: {}'
        raise KeyError(err_str.format(virus_id,
                                      ', '.join(list(PARSFUN_FOR_VIRUSID.keys()))))

    return {'ventilator_capacity': ventilator_capacity,
           'tspan': tspan,
           'p_t': p_t,
           'p_v': p_v,
           'p_h': p_h,
           'p_d': p_d,
           'p_dnv': p_dnv,
           'tau': tau}, population


COUNTRYFUN_FOR_COUNTRYID = {'denmark': _get_denmark}
