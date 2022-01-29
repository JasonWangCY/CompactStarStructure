# -------------------------------------------------------------------
# solve using Runge-Kutta method to 4th order 
# -------------------------------------------------------------------

import StellarStructure.functions as fn

class SolveEOS:
    """
    Solve the hydrostatic equations with RK4 and the corresponding equation of state
    """
    def __init__(self, r0, h):
        """
        Inputs:
            - r0, rho0, M0: initial conditions
        """
        self.r0 = r0
        self.h = h

    def solve_TOS(self, M0, rho0, TOV=True, graph_flag=True):
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

        r_list, M_list, rho_list, P_list = fn.RK4(self.r0, M0, rho0, self.h, fn.EOS, dMdr, drhodr)

        if graph_flag:
            fn.plot_graph(r_list, M_list, rho_list, P_list)

        return r_list, M_list, rho_list, P_list

if __name__ == '__main__':
    # python -m StellarStructure.SolveEOS
    solve = SolveEOS(0.0001, 10000)
    print(solve.__dict__)