import matplotlib.pyplot as plt
import numpy as np

def plot(times, sick, hospitalized,
         ventilator, recovered, dead,
         ventilator_capacity=None):

    max_hospitalized = int(np.max(hospitalized))
    max_ventilator = int(np.max(ventilator))
    total_deaths = int(dead[-1])

    fig, ax = plt.subplots(figsize=(15, 8))

    ax.plot(times, sick, 'b-', linewidth=2,
            label='Sick')
    ax.plot(times, hospitalized, 'b--', linewidth=2,
            label='Hospitalized (max={:})'.format(max_hospitalized))
    ax.plot(times, ventilator, 'b:', linewidth=2,
            label='Ventilator (max={})'.format(max_ventilator))
    ax.plot(times, dead, 'r-', linewidth=2,
            label='Dead (total={})'.format(total_deaths))
#     ax.plot(times, recovered, 'g-', linewidth=2,
#             label='Recovered')
    if ventilator_capacity is not None:
        ax.hlines(ventilator_capacity, xmin=min(t_eval), xmax=max(t_eval), color='black', linewidth=2,
                label='Ventilator capacity')


    # TODO: show parameters in plot

    ax.set_xlabel('Time (τ)')
    ax.set_ylabel('Individuals')
    ax.legend(loc='upper left')
    plt.show()
