import numpy as np
import pandas as pd

def eps_ppc(X, T, rho):
    
    T9 = 1e-9 * T 
    
    eps_p = ((2.4e4 * rho * (X**2)) / (T9**(2/3))) * exp(-3.88 / (T9**(1/3)))
             
    return eps_p
  
def eps_cno(X, Z, T, rho):
    
    T9 = 1e-9 * T 
    
    eps_c = ((4.4e25 * rho * X * Z) / (T9**(2/3))) * exp(-15.228 / (T9**(1/3)))
    
    return eps_c
