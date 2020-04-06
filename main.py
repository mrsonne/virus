from src import countries

# TODO 
# Notebook
# virus_id as cmd line input
# Compare multiple scenarios
# Make two population groups with different parameters


# How many people do you meet per day
encounters_per_day = 50 # adjusted to fit flu data (adjusted with transmission probability)
# encounters_per_day = 15 # damped sick count
# encounters_per_day = 35 # flattened

# encounters_per_day = None # constant sick count

# countries.run_country('flu', 'denmark', encounters_per_day)
# countries.run_country('covid19', 'denmark', encounters_per_day)
# countries.run_country('covid19', 'denmark', encounters_per_day, show_recovered=True, tspan=[0,320])


# countries.map_country('covid19', 'denmark', encounters_per_day, tspan=[0, 400])
# countries.map_country('flu', 'denmark', encounters_per_day, tspan=[0, 400])

countries.ua('covid19', 'denmark', encounters_per_day, 
              nsamples=100, tspan=[0, 400])