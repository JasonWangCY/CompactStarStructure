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


def EOS_rel(rho):
    """
    Equation of state for relativistic gas: P(rho)
    """
    return cg.K_rel * (rho * cg.Y_e)**(4/3)

def EOS_nr(rho):
    """
    Equation of state for non-relativistic gas: P(rho)
    """
    return cg.K_nr * (rho * cg.Y_e)**(5/3)

def EOS(rho):
    """
    Combined equation of state for WD
    """
    P_nr = EOS_nr(rho)
    P_rel = EOS_rel(rho)
    sqrt_term = np.sqrt( 1/P_nr/P_nr + 1/P_rel/P_rel )
    return 1/sqrt_term

def EOS_deriv(rho):
    """
    First derivative of the combined equation of state
    """
    P_nr = EOS_nr(rho)
    P_rel = EOS_rel(rho)
    P = EOS(rho)
    return P/rho * ( 5/3*(P/P_nr)**2 + 4/3*(P/P_rel)**2 )

def dMdr(r, M, rho):
    """
    Mass of a spherical body
    """
    return 4 * cg.pi * r*r *rho

def drhodr_Newt(r, M, rho):
    """
    Hydrostatic equilibrium for a star modelled as a Newtonian body
    """
    return -cg.G * M * rho / r / r / EOS_deriv(rho)

def drhodr_TOV(r, M, rho):
    """
    Tolman–Oppenheimer–Volkoff equation for a star modelled using general relativity
    """
    P = EOS(rho)
    return -cg.G*M*rho/r/r * (1+P/rho/cg.c/cg.c) * \
        (1 + 4*cg.pi*r*r*r*P/M/cg.c/cg.c) \
            / (1 - 2*cg.G*M/r/cg.c/cg.c) / EOS_deriv(rho)


def RK4(r0, M0, rho0, h, EOS, f, g):
    """
    Runge-Kutta method to 4th order to solve two simultaneous ODEs w.r.t. x
    Inputs:
        - r0, M0, rho0: initial conditions
        - h: step size
        - EOS: equation of state function
        - f(x, y, z) & g(x, y, z): functions for the two ODEs
            eg. f = dMdr(r, M, rho) & g = dPdr(r, M, rho)
    """
    x = r0
    y = M0
    z = rho0
    x_list = [x / 6.9634e10]
    y_list = [y / 1.9885e33]
    z_list = [z * 1000]
    P_list = [EOS(z) * 0.1]

    while z.real > 1e-5:
        k0 = f(x, y, z)
        k1 = f(x+h/2, y+h/2*k0, z)
        k2 = f(x+h/2, y+h/2*k1, z)
        k3 = f(x+h, y+h*k2, z)

        l0 = g(x, y, z)
        l1 = g(x+h/2, y, z+h/2*l0)
        l2 = g(x+h/2, y, z+h/2*l1)
        l3 = g(x+h, y, z+h*l2)

        x = x+h
        y = y + 1/6 * h * (k0 + 2*k1 + 2*k2 + k3)
        z = z + 1/6 * h * (l0 + 2*l1 + 2*l2 + l3)
        P = EOS(z)

        temp_x = x / 6.9634e10      # convert r to units of solar radius
        temp_y = y / 1.9885e33      # convert M to units of solar mass
        temp_z = z * 1000           # convert rho to units of kg/m^3
        temp_P = P * 0.1            # convert P to units of Pa
        x_list.append(temp_x)
        y_list.append(temp_y)
        z_list.append(temp_z)
        P_list.append(temp_P)
    
    return x_list, y_list, z_list, P_list


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