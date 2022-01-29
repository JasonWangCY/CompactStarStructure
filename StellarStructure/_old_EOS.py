"""
The following equations of state are no longer needed, but kept for the sake of future reference.
"""
from StellarStructure.constants import SI_unit as cs
import numpy as np


def EOS_rel(rho):
    """
    Equation of state for relativistic gas modelling the star as a polytrope (n=3)
    See https://www.astro.princeton.edu/~gk/A403/polytrop.pdf
    """
    K = 0.3639 * cs.G * (cs.M)**(2/3)         # polytropic constant
    Lambda = 4/3                        # polytropic index (lambda = 1+1/n)

    return K * rho**(Lambda)

def EOS_rel_solve(P):
    K = 0.3639 * cs.G * (cs.M)**(2/3)         # polytropic constant
    Lambda = 4/3                        # polytropic index (lambda = 1+1/n)

    return (P/K)**(1/Lambda)


def EOS_deg_rel(rho):
    """
    Equation of state for relativistic degenerate electron gas
    See https://www.ucolick.org/~woosley/ay112-14/lectures/lecture5.4x.pdf (WRONG?!?!)
    See https://scholarworks.wmich.edu/cgi/viewcontent.cgi?article=5320&context=masters_theses (p.39)
    """
    return cs.C4 * (0.5*rho)**(4/3)

def EOS_deg_rel_solve(P):
    return (P / cs.C4)**(3/4) * 2


def EOS_deg(rho):
    """
    Equation of state for non-relativistic degenerate electron gas
    See https://www.ucolick.org/~woosley/ay112-14/lectures/lecture5.4x.pdf (WRONG?!?!)
    See https://scholarworks.wmich.edu/cgi/viewcontent.cgi?article=5320&context=masters_theses (p.39)
    """
    return cs.C5 * (0.5*rho)**(5/3)

def EOS_deg_solve(P):
    return (P / cs.C5)**(3/5) * 2


def EOS_deg_comp(rho):
    """
    Comprehensive equation of state for ideal degenerate electron gas
    """
    a = cs.C2 * (cs.C3 * rho)**(1/3)
    return cs.C1 * ( a * np.sqrt(1+a*a) * (2/3*a*a-1) + np.log(a + np.log(1+a*a)) )

def EOS_deg_comp_solve(x, P):
    """
    Needs to be solved using fsolve (accuracy not ideal)
    """
    a = cs.C2 * (cs.C3 * x)**(1/3)
    return cs.C1 * ( a * np.sqrt(1+a*a) * (2/3*a*a-1) + np.log(a + np.log(1+a*a)) ) - P