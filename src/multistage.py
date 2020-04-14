import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import expon, erlang
import matplotlib.pyplot as plt


# def dydt_exp(t, y, r, y0):
#     """ Two equivalent ways
#         Probability of decaying in time interval delta_t
#           int_t^(t + delta_t) p(t) dt = p(t)delta_t for delta_t -> 0
#         Thus the change in population (remaining survivors)
#         dS(t) = -p(t)*dt = -r*S(t)dt
#
#     """
#     # val = -r*y
#     val = [-y0*p_erlang(t, r, 1)] # assumes whole population have the same "age"
#     return val


# def dydt_erlang(t, y, r, k, y0):
#   """"
#   dS(t) = -p(t)*dt cannot be manipulated into depending on the current state only 
#   """"
#     return [-y0*p_erlang(t, r, k)] # assumes whole population have the same "age"


# def p_erlang(t, r, k):
#     """pdf erlang
#     """
#     return r**k * t**(k-1) * np.exp(-r*t) / np.math.factorial(k - 1)


# def s_erlang(t, r, k):
#     """Survival function
#     """
#     return sum([(r*t)**n * np.exp(-r*t) / np.math.factorial(n) for n in range(k)])


def dydt_lct(t, y, r, k):
    """linear chain trick
    """
    dydt = -r*y
    dydt[1:] -= dydt[:-1]
    return dydt


def example():
    tspan = [0, 12.]
    t_eval = np.linspace(*tspan, 35)
    y0 = [1,]

    k = 5
    # scale to same mean
    scale = 1.5
    r = 1./scale

    # Exp
    r_exp = r/k
    r_exp
    pdf_exp = expon.pdf(t_eval, scale=1./r_exp)
    survival_exp = y0[0]*expon.sf(t_eval, scale=1./r_exp)

    y0_lct = np.zeros(k)
    y0_lct[0] = y0[0]
    sol_lct = solve_ivp(dydt_lct, tspan, y0_lct, t_eval=t_eval, args=(r, k), method='Radau')

    # print('mean exp: ', 1/r_exp)
    # print('rate erlang: ', r)
    # print('scale erlang: ', 1./r)
    # print('mean erlang: ', k/r)

    pdf_erlang = erlang.pdf(t_eval, a=k, scale=1./r)
    survival_erlang = y0[0]*erlang.sf(t_eval, a=k, scale=1./r)

    fig, axs = plt.subplots(2, 1, figsize=(15, 15), sharex=True)

    pdf_idx = 0
    survival_idx = 1

    exp_label = 'Exp(x; $\lambda$={:5.3f})'.format(r_exp)
    erlang_label = 'Erlang(x; $\lambda$={:5.3f},$k$={})'.format(r, k)

    axs[pdf_idx].plot(t_eval, pdf_exp, 'k-', label='PDF {}'.format(exp_label))
    axs[pdf_idx].plot(t_eval, pdf_erlang, 'k--', label='PDF {}'.format(erlang_label))

    axs[survival_idx].plot(t_eval, survival_exp, 'k-', label='{}'.format(exp_label))
    axs[survival_idx].plot(t_eval, survival_erlang, 'k--', label='{}'.format(erlang_label) )

    axs[survival_idx].plot(t_eval, np.sum(sol_lct.y, axis=0), marker='o', linewidth=0,
                           markeredgecolor='k', markerfacecolor='none', 
                           label=r'$I = I_1 + I_2 + I_3 + I_4 + I_5$')

    for istage in range(k):
        axs[survival_idx].plot(t_eval, sol_lct.y[istage], ':', label='$I_{}$'.format(istage + 1))


    # Styling
    font_size = 16
    for ax in axs:
        ax.legend(fontsize=font_size) 

        ax.tick_params(axis='x', labelsize=font_size)
        ax.tick_params(axis='y', labelsize=font_size)


    axs[-1].set_xlabel('x', fontsize=font_size)
    plt.show()


if __name__ == '__main__':
    example()