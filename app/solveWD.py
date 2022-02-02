"""
Solve for the structure of a white dwarf.
Constants are in cgs units.
"""

from StellarStructure.SolveEOS import SolveEOS
from StellarStructure.constants import CGS_unit as cg
import pandas as pd

def single_test(rho0=1e9):
    """
    Test for a specific value of rho to model structure of a White Dward (WD)
    """
    r0 = 0.0001
    h = 1000000

    M0 = 4/3 * cg.pi * r0**3 * rho0
    eos = SolveEOS(r0, h)
    r_list, M_list, rho_list, P_list = eos.solve_TOS(M0, rho0, TOV=True, graph_flag=True)

    df = pd.DataFrame({"r_list": r_list, "M_list": M_list, "rho_list": rho_list, "P_list": P_list})
    df.to_csv("app/output/single_test.csv", index=False)



if __name__ == '__main__':
    single_test()