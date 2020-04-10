from src import run, render

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

# run.virus('flu', 'denmark', encounters_per_day)
# run.virus('covid19', 'denmark', None)
run.virus('covid19', 'denmark', encounters_per_day)
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

times, ventilator_tseries = run.ua('covid19', 'denmark', 30,
                                    nsamples=500, tspan=[0, 400])

render.ua_timeseries(times, ventilator_tseries.T)
