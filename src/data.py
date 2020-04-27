def get_covid19_parameters():
    # Assumes same hospitalization rate as flu. Data is not reliable for COVID-19 yet. Parameter is not critical
    _, p_t, _, p_h, p_d, _, _ = get_flu_parameters()
    tspan = [0, 200]
    p_d = .005 
    # p_dnv =  p_d # For testing: no effect of exceeding ventilator capacity

    # https://medium.com/@tomaspueyo/coronavirus-the-hammer-and-the-dance-be9337092b56
    p_dnv = .9 # 
    p_v = 0.2 # of hospitalized
    tau = 14 #days
    return tspan, p_t, p_v, p_h, p_d, p_dnv, tau


def get_flu_parameters():
    # Get flu data and use as baseline
    p_h, p_d = get_us_flu_data() 
    tspan = [0, 700]
    p_dnv = 0.
    p_t = 0.005026 # adjusted to fit flu data with k=1
    p_v = 0. # unknown
    tau = 7 #days
    return tspan, p_t, p_v, p_h, p_d, p_dnv, tau


def get_us_flu_data():
    # US flu data
    # https://www.cdc.gov/flu/about/burden/preliminary-in-season-estimates.htm
    cases = 45000000. # per year
    p_h = 650000. / cases # ~1.3 %
    p_d = 47000. / cases # ~0.1 %
    return p_h, p_d

PARSFUN_FOR_VIRUSID = dict(flu=get_flu_parameters,
                           covid19=get_covid19_parameters)


def _get_denmark(virus_id, p_t_spec):
    """
    Danish flu data in Table 2 (~1000 deads per year)
    https://www.ssi.dk/sygdomme-beredskab-og-forskning/sygdomsovervaagning/i/influenzasaesonen---opgoerelse-over-sygdomsforekomst-2018-19
    """
    population = int(5.5e6)
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

    _p_t = p_t_spec or p_t 

    return {'ventilator_capacity': ventilator_capacity,
           'tspan': tspan,
           'p_t': _p_t,
           'p_v': p_v,
           'p_h': p_h,
           'p_d': p_d,
           'p_dnv': p_dnv,
           'tau': tau}, population



def _get_usa(virus_id, p_t_spec):
    """
    """
    population = int(327e6)
    ventilator_capacity = population # unknown

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

    _p_t = p_t_spec or p_t 
    return {'ventilator_capacity': ventilator_capacity,
           'tspan': tspan,
           'p_t': _p_t,
           'p_v': p_v,
           'p_h': p_h,
           'p_d': p_d,
           'p_dnv': p_dnv,
           'tau': tau}, population


COUNTRYFUN_FOR_COUNTRYID = {'denmark': _get_denmark,
                            'usa': _get_usa}
