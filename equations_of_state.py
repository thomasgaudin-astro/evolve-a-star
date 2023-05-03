import pandas as pd
import numpy as np
from math import exp, log10, log
from scipy.interpolate import RectBivariateSpline
from mpmath import polylog

# constants
h = 6. 626176e-27 #erg/s #planck's const
me = 9.1094e-28 #electron mass (g)
kB = 1.3807e-16 #cm2 g s-2 K-1 boltzmann const

a = 7.56471e-15 #radiation constant in ergs.cm^-3.K^-4
mp = 1 #proton mass in amu
amu_to_g = 1.66054e-24 #1 amu in g
mp_in_g = mp*amu_to_g #convert mp from amu to g ****************************
ma = mp_in_g

#calculate new abundances
X2 = calc_new_abund(Ai, X, r_x1, r_x2, t1, t2, rho)
Y2 = calc_new_abund(Ai, Y, r_x1, r_x2, t1, t2, rho)

"""Number density needed to calculate psi"""
def number_density(X,rho):
    NA = 6.02e23 #avogadro's number = 6.02e23 mol^-1,  beacause 1mol = 6.02e23 particles
    mu_e = 2/(1+X)
    n = rho*NA/mu_e
    return n

"""ψ to be used for the general equation of state"""
def psi(T,X,rho):
    ne = number_density(X,rho)
    ps = np.log( ne*h**3 / (2*(2*np.pi*me*kB*T)**3/2))
    return ps

# use rho to make guess for X, Y and Z (carbon)
File = np.loadtxt("../Data/structure_00000.txt", dtype=float) #load the structure file
carbon_mf = np.mean(File[:,32])

# the equation of state...pressure
"""Need to provide Temp, H and He, abundance. Z is the average carbon mass function from the
structure file"""
def pressure(T,X,rho, Y, Z=carbon_mf):
    Psi = psi(T,X,rho)
    mu_I, mu_e = (X/1 + Y/4 + Z/14)**-1 , 2/(1+X) 
    mu = (1/mu_I + 1/mu_e)**-1 #"""mean molecular weight """ 
    integral = -3*kB*T*np.sqrt(np.pi/2) * polylog(5/2, -np.exp(1)**ps) / ((1/kB*me*T)**(3/2))
    pres = (a*T_init**4)/3 + rho*kB*T/(mu*ma) + 8*np.pi/(3*h**3) * integral
    return pres
