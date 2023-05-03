import pandas as pd
import numpy as np
from math import exp, log10
from scipy.interpolate import RectBivariateSpline

def calc_delta(rho2, rho1, T2, T1):
    
    delta = - ((rho2 - rho1)/(T2-T1)) *(T2/rho2)
    
    return delta

def calc_alpha(rho2, rho1, P2, P1):
    
    alpha = ((rho2 - rho1)/(P2-P1)) *(P2/rho2)
    
    return alpha

def calc_cp(P, T, rho, nabla_ad, delta):
    
    cp = (P * delta) / (T * rho * nabla_ad)
    
    return cp

def calc_cv(P, T, rho, delta, alpha, cp):
    
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

def calc_kappa(opacity_tab, T, P):
    
    rho = calc_rho(T, P)
    
    opacities, log_T, log_R = load_opacities(opacity_tab)
    
    kap_func = RectBivariateSpline(log_T, log_R, opacities)
    
    logT = log10(T)
    
    logR = calc_logR(T, rho)
    
    kappa = kap_func(logT, logR)[0][0]
    
    return kappa

def calc_Nabla_ad(P1, P2, T1, T2):
    
    Nabla_ad = ((T2 - T1) / (P2 - P1)) * (P2 / T2)
    
    return Nabla_ad

def calc_Nabla_rad(kappa, L, P, M, T):
    
    a = 7.56471e-15 #erg cm-3 K-4
    
    G = 6.67e-8 #units
    
    Nabla_rad = (3 * kappa * L * P) / (16 * pi * a * c * G * M * (T**4))
    
    return Nabla_rad

def calc_Nabla(Nabla_ad, Nabla_rad):

    Nabla = min(Nabla_ad, Nabla_rad)
    
    return Nabla
