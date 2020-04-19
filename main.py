from src import run, render, multistage, solver
import numpy as np
# Seed so I dont have to change the text for each run
np.random.seed(0)

# TODO 
# Notebook
# virus_id as cmd line input
# Compare multiple scenarios
# Make two population groups with different parameters


# How many people do you meet per day
encounters_per_day = 50 # adjusted to fit flu data (adjusted with transmission probability)
# encounters_per_day = 15 # damped infected count
# encounters_per_day = 35 # flattened

# encounters_per_day = None # constant infected count

# run.virus('flu', 'denmark', encounters_per_day, k=1)
# run.virus('flu', 'usa', encounters_per_day)
# run.virus('covid19', 'denmark', None)
# run.virus('covid19', 'denmark', encounters_per_day, k=1)
# run.virus('covid19', 'denmark', encounters_per_day, show_recovered=True, tspan=[0,320])

# par1 = dict(axlabel=r'$p_{\rm{d}}$ (%)',
#             name='p_d',
#             min=0.001,
#             max=0.01,
#             transform=run.frc_to_pct
#             )
# par2 = dict(axlabel='τ (days)',
#             name='tau',
#             min=6,
#             max=21,
#             )

# response = dict(name='dead',
#                 transform=run.get_last,
#                 title="Dead",
#                )

# run.contour('covid19', 'denmark', par1, par2, response,
#             encounters_per_day, tspan=[0, 400])
# run.contour('flu', 'denmark', encounters_per_day, tspan=[0, 400])

    # xnames = r'$p_{\rm{h}}$ (%)', r'$\tau$ (days)', '$E$ (day\u207B\u00B9)'
    # mean = [pars["p_h"] , pars["tau"], pars["E"]]
    # cov = [[0.00001, 0. , 0], [0., 4., 0], [0., 0, 16]]

par1 = dict(axlabel=r'$p_{\rm{h}}$ (%)',
            name='p_h',
            std=0.003,
            transform=run.frc_to_pct
            )

par2 = dict(axlabel=r'$\tau$ (days)',
            name='tau',
            std=1.5,
            )

par3 = dict(axlabel='$E$ (day\u207B\u00B9)',
            name='E',
            std=4.,
            )

response = dict(name='ventilators_required',
                transform=run.get_max,
                title="Ventilator required",
               )

pars = par1, par2, par3
times, ventilator_tseries = run.ua('covid19', 'denmark', 30,
                                    pars, response,
                                    nsamples=100, tspan=[0, 400])
render.ua_timeseries(times, ventilator_tseries.T)


# def multistage_ex():
#     # multistage example
#     k = 5
#     scale = 1.5
#     rate = 1./scale
#     # Adjust exp rate to same mean as Erlang
#     rate_exp = rate/k 
#     multistage.example(rate, rate_exp, k)

#     # Adjust rate to get same survival function value (0.2) after two weeks
#     rate_exp = solver.get_rate_Iminus(14, 0.2)
#     rate = solver.get_rate_Iminus(14, 0.2, k=k)
#     multistage.example(rate, rate_exp, k)


# multistage_ex()