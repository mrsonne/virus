import matplotlib.pyplot as plt
import numpy as np

def plot(times, sick, hospitalized,
         ventilator, recovered, dead,
         ventilator_capacity=None,
         show_recovered=False,
         title=''):

    max_hospitalized = np.max(hospitalized)
    max_ventilator = np.max(ventilator)
    max_recovered = np.max(recovered)
    total_deaths = dead[-1]

    fig, ax = plt.subplots(figsize=(15, 8))

    ax.plot(times, sick, 'b-', linewidth=2,
            label='Sick')
    ax.plot(times, hospitalized, 'b--', linewidth=2,
            label='Hospitalized (max={:5.1e})'.format(max_hospitalized))
    ax.plot(times, ventilator, 'b:', linewidth=2,
            label='Ventilator (max={:5.1e})'.format(max_ventilator))
    ax.plot(times, dead, 'r-', linewidth=2,
            label='Dead (total={:5.1e})'.format(total_deaths))

    if show_recovered:
        alpha = 1
        ylim = None
    else:
        ylim = ax.get_ylim()
        alpha = 0.2

    ax.plot(times, recovered, 'g-', linewidth=2, alpha=alpha, label='Recovered (max={:5.1e})'.format(max_recovered))
    if ylim: ax.set_ylim(ylim)


    if ventilator_capacity is not None:
        ax.hlines(ventilator_capacity, xmin=min(t_eval), xmax=max(t_eval), color='black', linewidth=2,
                label='Ventilator capacity')


    # TODO: show parameters in plot

#     ax.legend(loc='upper left')
    font_size = 16
    ax.set_xlabel('Time (τ)', fontsize=font_size)
    ax.set_ylabel('Individuals', fontsize=font_size)
    ax.legend(fontsize=font_size) 
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    plt.title(title, fontsize=font_size*1.25)
    plt.show()
