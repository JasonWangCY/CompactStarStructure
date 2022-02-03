# -------------------------------------------------------------------
# equations that describe hydrostatic equilibrium
# Inputs
#   - rho: density
#   - r: radius
#   - M: enclosed mass of the star/spherical body
# -------------------------------------------------------------------

from StellarStructure.constants import CGS_unit as cg
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Normal matter
def EOS_rel(rho, m_dm):
    """
    Equation of state for relativistic gas: P(rho)
    """
    return cg.K_rel * (rho * cg.Y_e)**(4/3)

def EOS_nr(rho, m_dm):
    """
    Equation of state for non-relativistic gas: P(rho)
    """
    return cg.K_nr * (rho * cg.Y_e)**(5/3)

def EOS(rho, m_dm):
    """
    Combined equation of state for WD
    """
    # P_nr = EOS_nr(rho, m_dm)
    # P_rel = EOS_rel(rho, m_dm)
    # sqrt_term = np.sqrt( 1/P_nr/P_nr + 1/P_rel/P_rel )
    # return 1/sqrt_term

    x = cg.x_const * rho**(1/3)
    return cg.K_exact * ( x*np.sqrt(1+x*x)*(2/3*x*x-1) + np.log(x+np.sqrt(1+x*x)) )

def EOS_deriv(rho, m_dm):
    """
    First derivative of the combined equation of state
    """
    # P_nr = EOS_nr(rho, m_dm)
    # P_rel = EOS_rel(rho, m_dm)
    # P = EOS(rho, m_dm)
    # return P/rho * ( 5/3*(P/P_nr)**2 + 4/3*(P/P_rel)**2 )

    x = cg.x_const * rho**(1/3)
    term1 = np.sqrt(1+x*x) * (2/3*x*x-1)
    term2 = x*x/np.sqrt(1+x*x) * (2/3*x*x-1)
    term3 = x*np.sqrt(1+x*x)*4/3*x
    term4 = ( 1+x/np.sqrt(1+x*x) ) / (x + np.sqrt(1+x*x))
    return cg.K_exact*(term1+term2+term3+term4) * cg.x_const * rho**(-2/3) / 3
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Dark matter
def EOS_rel_dm(rho, m_dm):
    """
    Equation of state for relativistic dark matter: P(rho)
    """
    return cg.K_rel_dm * (rho)**(4/3)

def EOS_nr_dm(rho, m_dm):
    """
    Equation of state for non-relativistic dark matter: P(rho)
    """
    return cg.K_nr_dm * (rho)**(5/3) / m_dm

def EOS_dm(rho, m_dm):
    """
    Combined equation of state for WD
    """
    P_nr = EOS_nr_dm(rho, m_dm)
    P_rel = EOS_rel_dm(rho, m_dm)
    sqrt_term = np.sqrt( 1/P_nr/P_nr + 1/P_rel/P_rel )
    return 1/sqrt_term

def EOS_deriv_dm(rho, m_dm):
    """
    First derivative of the combined equation of state
    """
    P_nr = EOS_nr_dm(rho, m_dm)
    P_rel = EOS_rel_dm(rho, m_dm)
    P = EOS_dm(rho, m_dm)
    return P/rho * ( 5/3*(P/P_nr)**2 + 4/3*(P/P_rel)**2 )
# ---------------------------------------------------------------------------


def dMdr(r, M, M_dm, rho, m_dm):
    """
    Mass of a spherical body
    """
    return 4 * cg.pi * r*r *rho

def drhodr_Newt(r, M, M_dm, rho, EOS, EOS_deriv, m_dm):
    """
    Hydrostatic equilibrium for a star modelled as a Newtonian body
    """
    M_tot = M + M_dm
    return -cg.G * M_tot * rho / r / r / EOS_deriv(rho, m_dm)

def drhodr_TOV(r, M, M_dm, rho, EOS, EOS_deriv, m_dm):
    """
    Tolman–Oppenheimer–Volkoff equation for a star modelled using general relativity
    """
    P = EOS(rho, m_dm)
    M_tot = M + M_dm
    return -cg.G*M_tot*rho/r/r * (1+P/rho/cg.c/cg.c) * \
        (1 + 4*cg.pi*r*r*r*P/M_tot/cg.c/cg.c) \
            / (1 - 2*cg.G*M_tot/r/cg.c/cg.c) / EOS_deriv(rho, m_dm)


def RK4(r0, M0, rho0, M_dm0, rho_dm0, h, f, g, m_dm):
    """
    Runge-Kutta method to 4th order to solve two simultaneous ODEs w.r.t. x
    Inputs:
        - r0, M0, rho0: initial conditions
        - h: step size
        - EOS: equation of state function
        - f(x, y, z) & g(x, y, z): functions for the two ODEs
            eg. f = dMdr(r, M, M_dm, rho) & g = dPdr(r, M, M_dm, rho, EOS_deriv)
    """
    x = r0
    y = M0
    z = rho0
    y2 = M_dm0
    z2 = rho_dm0
    x_list = [x / 6.9634e10]
    y_list = [y / 1.9885e33]
    z_list = [z * 1000]
    P_list = [EOS(z, m_dm) * 0.1]
    y2_list = []
    z2_list = []
    P2_list = []
    if m_dm != 0:
        y2_list = [y2 / 1.9885e33]
        z2_list = [z2 * 1000]
        P2_list = [EOS_dm(z2, m_dm) * 0.1]

    while z.real > 1e-5:
        k0 = f(x, y, y2, z, m_dm)
        k1 = f(x+h/2, y+h/2*k0, y2, z, m_dm)
        k2 = f(x+h/2, y+h/2*k1, y2, z, m_dm)
        k3 = f(x+h, y+h*k2, y2, z, m_dm)

        l0 = g(x, y, y2, z, EOS, EOS_deriv, m_dm)
        l1 = g(x+h/2, y, y2, z+h/2*l0, EOS, EOS_deriv, m_dm)
        l2 = g(x+h/2, y, y2, z+h/2*l1, EOS, EOS_deriv, m_dm)
        l3 = g(x+h, y, y2, z+h*l2, EOS, EOS_deriv, m_dm)

        if m_dm != 0:
            kk0 = f(x, y, y2, z2, m_dm)
            kk1 = f(x+h/2, y, y2+h/2*kk0, z2, m_dm)
            kk2 = f(x+h/2, y, y2+h/2*kk1, z2, m_dm)
            kk3 = f(x+h, y, y2+h*kk2, z2, m_dm)

            ll0 = g(x, y, y2, z2, EOS_dm, EOS_deriv_dm, m_dm)
            ll1 = g(x+h/2, y, y2, z2+h/2*ll0, EOS_dm, EOS_deriv_dm, m_dm)
            ll2 = g(x+h/2, y, y2, z2+h/2*ll1, EOS_dm, EOS_deriv_dm, m_dm)
            ll3 = g(x+h, y, y2, z2+h*ll2, EOS_dm, EOS_deriv_dm, m_dm)

        x = x+h
        y = y + 1/6 * h * (k0 + 2*k1 + 2*k2 + k3)
        z = z + 1/6 * h * (l0 + 2*l1 + 2*l2 + l3)
        P = EOS(z, m_dm)

        temp_x = x / 6.9634e10      # convert r to units of solar radius
        temp_y = y / 1.9885e33      # convert M to units of solar mass
        temp_z = z * 1000           # convert rho to units of kg/m^3
        temp_P = P * 0.1            # convert P to units of Pa
        x_list.append(temp_x)
        y_list.append(temp_y)
        z_list.append(temp_z)
        P_list.append(temp_P)

        if m_dm != 0:
            y2 = y2 + 1/6 * h * (kk0 + 2*kk1 + 2*kk2 + kk3)
            z2 = z2 + 1/6 * h * (ll0 + 2*ll1 + 2*ll2 + ll3)
            P2 = EOS_dm(z2, m_dm)

            temp_y = y2 / 1.9885e33      # convert M to units of solar mass
            temp_z = z2 * 1000           # convert rho to units of kg/m^3
            temp_P = P2 * 0.1            # convert P to units of Pa
            y2_list.append(temp_y)
            z2_list.append(temp_z)
            P2_list.append(temp_P)
    
    return x_list, y_list, z_list, P_list, y2_list, z2_list, P2_list


def plot_graph(r_list, M_list, rho_list, P_list):
    plt.rcParams['figure.figsize'] = [10, 10]

    fig = plt.figure()
    ax3 = fig.add_subplot(313)
    ax2 = fig.add_subplot(312)
    ax1 = fig.add_subplot(311)

    ax1.set_title(r'Pressure $P$ against radius $r$')
    ax1.plot(r_list, P_list)
    ax1.set_ylabel(r'$P$ (Pa)')

    ax2.set_title(r'Mass $M$ against radius $r$')
    ax2.plot(r_list, M_list)
    ax2.set_ylabel(r'$M$ ($M_\odot$)')

    ax3.set_title(r'Density $\rho$ against radius $r$')
    ax3.plot(r_list, rho_list)
    ax3.set_ylabel(r'$\rho$ (kg/$m^3$)')
    ax3.set_xlabel(r'$r$ ($R_\odot$)')

    plt.tight_layout()

    ax1.tick_params(axis='both', which='major', direction = 'in')
    ax2.tick_params(axis='both', which='major', direction = 'in')
    ax3.tick_params(axis='both', which='major', direction = 'in')

    # ax1.ticklabel_format(style = 'sci', axis='y', scilimits=(0,0))
    # ax2.ticklabel_format(style = 'sci', axis='y', scilimits=(0,0))
    # plt.setp(ax1.get_xticklabels(), visible=False)
    # plt.setp(ax2.get_xticklabels(), visible=False)

    plt.show()
    # plt.savefig('RK4_TOV.png', facecolor='white', transparent=False)


if __name__ == '__main__':
    # python -m StellarStructure.functions
    print('This is for functions to be imported only.')