import pandas as pd
import numpy as np
from scipy.interpolate import RectBivariateSpline

#constants
a = 7.56471e-15 #erg cm-3 K-4
G = 6.67e-8 #units
ma = 1.66e-24 #g
NA = 6.022e23 #atoms mol^-1
h = 6.63e-27
kb = 1.38e-16 #cgs
kb_ev = 8.62e-5 #ev K^-1
m_He = 4.0026 * ma #g
c = 3.0e10 #cm s^-1

def delta_calc(P, T, X, Y, Z, old_rho):
    
    rho = rho_calc(P, T, X, Y, Z, old_rho)
   
    T2 = T + (0.01 * T)
    
    rho_2 = rho_calc(P, T2, X, Y, Z, old_rho)
    
    delta = - ((rho2 - rho)/(T2-T)) *(T2/rho_2)
    
    return delta

def alpha_calc(P, T, X, Y, Z, old_rho):
    
    rho = rho_calc(P, T, X, Y, Z, old_rho)
   
    P2 = P + (0.01 * P)
    
    rho_2 = rho_calc(P2, T, X, Y, Z, old_rho)
    
    alpha = ((rho2 - rho)/(P2-P)) *(P2/rho_2)
    
    return alpha

def calc_beta(P, T, X, Y, Z, old_rho):
    
    mu_I, mu_e = (X/1 + Y/4 + Z/14)**-1 , 2/(1+X) 
    mu = (1/mu_I + 1/mu_e)**-1
    
    rho = rho_calc(P, T, X, Y, Z, old_rho)
    
    P_ion = (rho / (mu * ma) ) * kb * T
    
    beta = P_ion / P
    
    return beta
    
def nabla_ad_calc(P, T, X, Y, Z, old_rho):
    
    beta = calc_beta(P, T, X, Y, Z, old_rho)
    
    Nabla_ad = (1 + ((1-beta)*(4+beta))/(beta**2)) / ((5/2) + (4*(1-beta)*(4+beta))/(beta**2))
    
    return Nabla_ad

def cp_calc(P, T, X, Y, Z, old_rho):
    
    delta = delta_calc(P, T, X, Y, Z, old_rho)
    rho = rho_calc(P, T, X, Y, Z, old_rho)
    Nabla_ad = nabla_ad_calc(P, T, X, Y, Z, old_rho)
    
    cp = (P * delta) / (T * rho * nabla_ad)
    
    return cp

def calc_cv(P, T, X, Y, Z, old_rho):
    
    delta = delta_calc(P, T, X, Y, Z, old_rho)
    alpha = alpha_calc(P, T, X, Y, Z, old_rho)
    
    rho = rho_calc(P, T, X, Y, Z, old_rho)
    cp = cp_calc(P, T, X, Y, Z, old_rho)
    
    cv = cp - (P *(delta**2)) / (rho*T*alpha)
    
    return cv

def calc_logR(T, rho):
    
    T6 = 1e-6 * T #kelvin
    
    R = rho / (T6**3) 
    
    return log10(R)

def load_opacities(tab):
    
    opacities = pd.read_csv(tab, header=0, index_col=0)
    
    opacity = opacities.fillna(value=0)
    opacity = 10**(opacity.loc[:7.5,:"-1.5"])
    
    log_T = 10**(np.array(opacity.index))
    
    log_R = (np.zeros(len(opacity.columns)))

    for val, ind in zip(opacity.columns, range(len(opacity.columns))):
        log_R[ind] = 10**float(val)
        
    return opacity, log_T, log_R

def kappa_calc(P, T, X, Y, Z, old_rho, opacity_tab):
    
    rho = rho_calc(P, T, X, Y, Z, old_rho)
    
    opacities, log_T, log_R = load_opacities(opacity_tab)
    
    kap_func = RectBivariateSpline(log_T, log_R, opacities)
    
    logT = log10(T)
    
    logR = calc_logR(T, rho)
    
    kappa = kap_func(logT, logR)[0][0]
    
    return kappa

def nabla_rad_calc(L, P, T, M, X, Y, Z, old_rho, opacity_tab):
    
    kappa = kappa_calc(P, T, X, Y, Z, old_rho, opacity_tab)
    
    Nabla_rad = (3 * kappa * L * P) / (16 * pi * a * c * G * M * (T**4))
    
    return Nabla_rad

def Nabla_calc(L, P, T, M, X, Y, Z, opacity_tab):
    
    Nabla_ad = nabla_ad_calc(P, T, X, Y, Z)
    Nabla_rad = nabla_rad_calc(L, P, T, M, X, Y, Z, opacity_tab)
    
    Nabla = min(Nabla_ad, Nabla_rad)
    
    return Nabla
