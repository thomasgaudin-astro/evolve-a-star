import pandas as pd
import numpy as np

def load_opacities(tab):
    
    opacities = pd.read_csv(tab, header=0, index_col=0)
    
    opacity = opacities.fillna(value=0)
    opacity = 10**(opacity.loc[:7.5,:"-1.5"])
    
    log_T = 10**(np.array(opacity.index))
    
    log_R = (np.zeros(len(opacity.columns)))

    for val, ind in zip(opacity.columns, range(len(opacity.columns))):
        log_R[ind] = 10**float(val)
        
    return opacity, log_T, log_R

def calc_beta(T, rho):
    
    beta = ((rho * kb * T) / (mu * mH)) * (1 / (((rho * kb * T) / (mu * mH)) + (a / 3)*(T**4)))
    
    return beta

def calc_Nabla_ad(P1, P2, T1, T2):
    
    Nabla_ad = ((T2 - T1) / (P2 - P1)) * (P2 / T2)
    
    return Nabla_ad

def calc_Nabla_rad(kappa, L, P, M, T):
    
    Nabla_rad = (3 * kappa * L * P) / (16 * pi * a * c * G * M * (T**4))
    
    return Nabla_rad

def calc_Nabla(Nabla_ad, Nabla_rad):

    Nabla = min(Nabla_ad, Nabla_rad)
    
    return Nabla
