"""
Solve for the structure of a white dwarf.
Constants are in cgs units.
"""

from StellarStructure.SolveEOS import SolveEOS
import StellarStructure.functions as fn
from StellarStructure.constants import CGS_unit as cg
import matplotlib.pyplot as plt
import numpy as np


def single_test(rho0=1e9):
    """
    Test for a specific value of rho to model structure of a White Dward (WD)
    """
    r0 = 0.0001
    h = 10000

    M0 = 4/3 * cg.pi * r0**3 * rho0
    eos = SolveEOS(r0, h)
    r_list, M_list, rho_list, P_list = eos.solve_TOS(M0, rho0, TOV=False, graph_flag=True)


def MR_relation(rho_low=3, rho_high=13):
    r0 = 0.0001
    rho0_list = np.logspace(rho_low, rho_high, 50)
    h = 10000
    totR_list = []
    totM_list = []
    totR_list_TOV = []
    totM_list_TOV = []

    eos = SolveEOS(r0, h)

    for rho0 in rho0_list:
        M0 = 4/3 * cg.pi * r0**3 * rho0

        r_list, M_list, rho_list, P_list = eos.solve_TOS(M0, rho0, TOV=False, graph_flag=False)       
        total_R = r_list[-1].real
        total_M = M_list[-1].real

        r_list_TOV, M_list_TOV, rho_list_TOV, P_list_TOV = eos.solve_TOS(M0, rho0, TOV=True, graph_flag=False) 
        total_R_TOV = r_list_TOV[-1].real
        total_M_TOV = M_list_TOV[-1].real

        totR_list.append(total_R)
        totM_list.append(total_M)
        totR_list_TOV.append(total_R_TOV)
        totM_list_TOV.append(total_M_TOV)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    ax.plot(totM_list, totR_list, c='b', label='Newtonian')
    ax.plot(totM_list_TOV, totR_list_TOV, c='r', label='TOV')
    ax.set_title('Mass-radius relationship of WD')
    ax.set_ylabel(r'$r$ ($R_\odot$)')
    ax.set_xlabel(r'$M$ ($M_\odot$)')
    ax.tick_params(axis='both', which='major', direction = 'in')
    ax.legend()
    plt.show()


def run():
    """
    Change this function to choose which to run
    """
    # single_test()
    MR_relation()

if __name__ == '__main__':
    run()