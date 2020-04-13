import numpy as np

def dydt_lct(t, y, r, k):
    """linear chain trick
    """
    dydt = -r*y
    dydt[1:] -= dydt[:-1]
    return dydt


def dydt_erlang(t, y, r, k, y0):
    return [-y0*p_erlang(t, r, k)] # experimental and broken


def dydt_exp(t, y, r, y0):
    """ Two equivalent ways
    """
    # val = -r*y
    val = [-y0*p_erlang(t, r, 1)]
    return val


def p_erlang(t, r, k):
    """pdf erlang
    """
    return r**k * t**(k-1) * np.exp(-r*t) / np.math.factorial(k - 1)


def s_erlang(t, r, k):
    """Survival function
    """
    return sum([(r*t)**n * np.exp(-r*t) / np.math.factorial(n) for n in range(k)])


def get_ynew(dt, y, r, k):
    return y*s_erlang(dt, r, k)


def example():
    from scipy.integrate import solve_ivp
    import matplotlib.pyplot as plt

    tspan = [0, 12.]
    t_eval = np.linspace(*tspan, 50)
    y0 = [18,]

    k = 5
    # scale to same mean
    scale = 1.5
    r = 1./scale

    # Exp
    r_exp = r/k
    print('mean exp: ', 1/r_exp)
    sol_exp = solve_ivp(dydt_exp, tspan, y0, t_eval=t_eval, args=(r_exp, y0[0]), method='Radau')
    pdf_exp = p_erlang(t_eval, r_exp, 1)
    survival_exp = y0[0]*s_erlang(t_eval, r_exp, 1)

    sol = solve_ivp(dydt_erlang, tspan, y0, t_eval=t_eval, args=(r, k, y0[0]), method='Radau')

    y0_lct = np.zeros(k)
    y0_lct[0] = y0[0]
    sol_lct = solve_ivp(dydt_lct, tspan, y0_lct, t_eval=t_eval, args=(r, k), method='Radau')

    yold = y0[0]
    dts = np.diff(t_eval)
    ys_fd_exp = y0

    print('rate erlang: ', r)
    print('scale erlang: ', 1./r)
    print('mean erlang: ', k/r)
    for dt in dts:
        ynew = get_ynew(dt, yold, r_exp, 1)
        ys_fd_exp.append(ynew)
        yold = ynew


    survival = y0[0]*s_erlang(t_eval, r, k)
    pdf = p_erlang(t_eval, r, k)

    fig, axs = plt.subplots(2, 1, figsize=(15, 15), sharex=True)

    pdf_idx = 0
    survival_idx = 1

    exp_label = 'Exp(x; r={:4.1e})'.format(r_exp)
    erlang_label = 'Erlang(x; r={:4.1e},k={})'.format(r, k)

    axs[pdf_idx].plot(t_eval, pdf, 'c-', label='PDF {}'.format(erlang_label))
    axs[pdf_idx].plot(t_eval, pdf_exp, 'c--', label='PDF {}'.format(exp_label))

    axs[survival_idx].plot(t_eval, survival_exp, 'm--', label='{}'.format(exp_label))
    axs[survival_idx].plot(t_eval, survival, 'm-', label='{}'.format(erlang_label) )

    axs[survival_idx].plot(t_eval, np.sum(sol_lct.y, axis=0), marker='o', linewidth=0, markeredgecolor='c',
                markerfacecolor='none', label=r'$I$')

    for istage in range(k):
        axs[survival_idx].plot(t_eval, sol_lct.y[istage], ':', label='$I_{}$'.format(istage + 1))



    # axs[1].plot(sol_exp.t, sol_exp.y[0], 'bx', label='Exp survival as IVP')
    # axs[1].plot(sol_exp.t, ys_fd_exp, color='b', marker='o', linewidth=0, 
                # markerfacecolor='none', label='Exp survival as FD')
    # axs[1].plot(t_eval, survival_exp, 'm--', label='Survival Exp')
    # axs[1].plot(t_eval, survival, 'm-', label='Survival Erlang(r={:4.1e},k={})'.format(r, k))

    # axs[1].plot(t_eval, sol_lct.y[-1], 'c-', label='LCT')

    # axs[1].plot(sol.t, sol.y[0], 'mx', label='Erlang survival as IVP')
    # ax.plot(t_eval, ys, 'ro')

    # Styling
    font_size = 16
    for ax in axs:
        ax.legend(fontsize=font_size) 

        ax.tick_params(axis='x', labelsize=font_size)
        ax.tick_params(axis='y', labelsize=font_size)


    axs[-1].set_xlabel('x', fontsize=font_size)
    axs[pdf_idx].set_ylabel('Density', fontsize=font_size)
    axs[survival_idx].set_ylabel('Individuals', fontsize=font_size)


    plt.show()


if __name__ == '__main__':
    example()