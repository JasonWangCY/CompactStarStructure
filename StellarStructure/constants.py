# -------------------------------------------------------------------
# a module that contains useful constant values/parameters to be imported
# -------------------------------------------------------------------
import numpy as np

class SI_unit:
    pi = 3.1416
    AU = 1.4960e11              # in m
    solar_mass = 1.9885e30      # in kg
    solar_radius = 6.9634e8     # in m 
    G = 6.6741e-11              # grav constant in m^3/kg/s^2
    c = 2.9979e8                # speed of light in m/s

    m_e = 9.1094e-31            # mass of electron in kg
    h_bar = 1.0546e-34          # in m^2 kg/s
    h = 6.626e-34      
    y_e = 0.5                   # for carbon, oxygen, neon etc
    m_H = 1.67e-27              # mass of a hydrogen-1 atom (in kg)

    """Parameter for EOS_rel"""
    M = 1.9885e30               # mass of the Sun (in kg)

    """Constants for EOS_deg_comp"""
    C1 = 1.8005e22              # 2 * m_e**4 * c**5 / (16 * pi*pi * h_bar**3)
    C2 = 3.6618e21              # 1/(m_e*c)
    C3 = 1.0398e-74             # 6 * pi*pi * h_bar**3 * y_e / 2 / m_H

    """Constants for EOS_deg"""
    C4 = 2*pi*c/(3*h**3) * (3*h**3/(8*pi*m_H))**(4/3)       # relativistic
    # C4 = (h*c/8) * (3/pi)**(1/3) / m_H**(4/3)
    # print(C4)
    # C4 = 1.2444e6
    C5 = 8*pi*c/(15*h**3) * (3*h**3/(8*pi*m_H))**(5/3)      # non-relativistic

class CGS_unit:
    pi = np.pi
    c = 2.9979e8 * 1e2
    Y_e = 0.5
    K_nr = 1.0036e13
    K_rel = 1.2435e15
    G = 6.6741e-11 * 1e3

if __name__ == '__main__':
    cgs = CGS_unit
    v = cgs.__dict__
    v_list = list(cgs.__dict__.keys())
    for key in v_list[1:-3]:
        print(key, '=', v[key])