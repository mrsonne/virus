def get_covid19_parameters():
    # Assumes same hospitalization rate as flu. Data is not reliable for COVID-19 yet. Parameter is not critical
    _, p_transmision, _, p_sick_to_hospitalized, p_sick_to_dead, _, _ = get_flu_parameters()
    tspan = [0, 20]
    p_sick_to_dead *= 6. # relative to flu
    # p_sick_to_dead_no_ventilator =  p_sick_to_dead # For testing: no effect of exceeding ventilator capacity

    # https://medium.com/@tomaspueyo/coronavirus-the-hammer-and-the-dance-be9337092b56
    # p_sick_to_dead_no_ventilator =  0.04 # TODO: Of all patients including collateral effects. Source??
    p_sick_to_dead_no_ventilator = .9 # 
    p_hospitalized_to_ventilator = 0.2 # of hospitalized
    tau = 14 #days
    return tspan, p_transmision, p_hospitalized_to_ventilator, p_sick_to_hospitalized, p_sick_to_dead, p_sick_to_dead_no_ventilator, tau


def get_flu_parameters():
    # Get flu data and use as baseline
    p_sick_to_hospitalized, p_sick_to_dead = get_us_flu_data() 
    tspan = [0, 100]
    p_sick_to_dead_no_ventilator = 0. # TODO: of the patients requirering ventilator?
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
