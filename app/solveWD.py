"""
Solve for the structure of a white dwarf.
Constants are in cgs units.
"""

from StellarStructure.SolveEOS import SolveEOS
from StellarStructure.constants import CGS_unit as cg


def single_test(rho0=1e9):
    """
    Test for a specific value of rho to model structure of a White Dward (WD)
    """
    r0 = 0.0001
    h = 10000

    M0 = 4/3 * cg.pi * r0**3 * rho0
    eos = SolveEOS(r0, h)
    r_list, M_list, rho_list, P_list = eos.solve_TOS(M0, rho0, TOV=False, graph_flag=True)


if __name__ == '__main__':
    single_test()