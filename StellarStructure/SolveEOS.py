# -------------------------------------------------------------------
# solve using Runge-Kutta method to 4th order 
# -------------------------------------------------------------------

import StellarStructure.functions as fn

class SolveEOS:
    """
    Solve the hydrostatic equations with RK4 and the corresponding equation of state
    """
    def __init__(self, r0, h, m_dm):
        """
        Inputs:
            - r0, rho0, M0: initial conditions
        """
        self.r0 = r0
        self.m_dm = m_dm
        self.h = h

    def solve_TOS(self, M0, rho0, M_dm0, rho_dm0, TOV=True, graph_flag=True):
        """
        Input flags:
            - TOV: Newtonian method (default as False)
                --> if True, use simple Newtonian version of dP/dr
                --> if False, use TOV version of dP/dr
        """
        if TOV:
            drhodr = fn.drhodr_TOV
        else:
            drhodr = fn.drhodr_Newt
        dMdr = fn.dMdr

        r_list, M_list, rho_list, P_list, r_dm_list, M_dm_list, rho_dm_list, P_dm_list \
            = fn.RK4(self.r0, M0, rho0, \
            M_dm0, rho_dm0, self.h, dMdr, drhodr, self.m_dm)

        if graph_flag:
            fn.plot_graph(r_list, M_list, rho_list, P_list)
            if self.m_dm != 0:
                fn.plot_graph(r_dm_list, M_dm_list, rho_dm_list, P_dm_list)

        return r_list, M_list, rho_list, P_list, r_dm_list, M_dm_list, rho_dm_list, P_dm_list

if __name__ == '__main__':
    # python -m StellarStructure.SolveEOS
    solve = SolveEOS(0.0001, 10000)
    print(solve.__dict__)