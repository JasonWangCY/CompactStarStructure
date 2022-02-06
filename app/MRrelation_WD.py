"""
Solve for the structure of a white dwarf.
Constants are in cgs units.
"""

from StellarStructure.SolveEOS import SolveEOS
from StellarStructure.constants import CGS_unit as cg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def MR_relation(rho_low=3, rho_high=13, rho_dm0=1e10, m_dm=380):
    """
    Generate the mass-radius relationship of a white dwarf
    """
    r0 = 0.0001
    M_dm0 = 4/3 * cg.pi * r0**3 * rho_dm0
    rho0_list = np.logspace(rho_low, rho_high, 50)
    h = 10000
    # totR_list = []
    # totM_list = []
    totR_list_TOV = []
    totM_list_TOV = []

    eos = SolveEOS(r0, h, m_dm)

    for rho0 in rho0_list:
        M0 = 4/3 * cg.pi * r0**3 * rho0

        # r_list, M_list, rho_list, P_list, r_dm_list, M_dm_list, rho_dm_list, P_dm_list\
        #      = eos.solve_TOS(M0, rho0, M_dm0, rho_dm0, TOV=False, graph_flag=False)       
        # total_R = r_list[-1].real
        # total_M = M_list[-1].real

        r_list_TOV, M_list_TOV, rho_list_TOV, P_list_TOV, r_dm_list_TOV, M_dm_list_TOV, rho_dm_list_TOV, P_dm_list_TOV\
             = eos.solve_TOS(M0, rho0, M_dm0, rho_dm0, TOV=True, graph_flag=False) 
        total_R_TOV = r_list_TOV[-1].real
        total_M_TOV = M_list_TOV[-1].real

        # totR_list.append(total_R)
        # totM_list.append(total_M)
        totR_list_TOV.append(total_R_TOV)
        totM_list_TOV.append(total_M_TOV)

    plt.rcParams['figure.figsize'] = [10, 10]
    fig, ax = plt.subplots()
    # ax.plot(totM_list, totR_list, c='b', label='Newtonian')
    ax.plot(totM_list_TOV, totR_list_TOV, c='r', label='TOV')
    ax.set_title('Mass-radius relationship of WD')
    ax.set_ylabel(r'$r$ ($R_\odot$)')
    ax.set_xlabel(r'$M$ ($M_\odot$)')
    ax.tick_params(axis='both', which='major', direction = 'in')
    ax.legend()
    plt.show()

    df = pd.DataFrame({"totM_list_TOV": totM_list_TOV, "totR_list_TOV": totR_list_TOV})
    df.to_csv("app/output/MRrelation_exact_WD_dm380_rho1e10.csv", index=False)


if __name__ == '__main__':
    MR_relation()