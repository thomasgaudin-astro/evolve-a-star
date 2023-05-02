import numpy as np
import pandas as pd

def calc_zeta(Zi, xi, Ai):
    #Zi - proton number
    #xi - abundance
    #Ai - Reduced Mass
    
    zeta_i = []
    
    for Z, x, A in zip(Zi, xi, Ai):
        
        zet = ( ((Zi) * (Zi + 1)) / Ai ) * xi
    
        zeta_i.append(zet)
        
    zeta = sum(zeta)
    
    return zeta

def calc_f(Z1, Z2, zeta, rho, T): 
    
    T6 = 1e-6 * T 
    
    zet_rho = (zeta * rho) / (T6**3)
    
    f = exp(0.188 * Z1 * Z2 * (zet_rho**(1/2)))
    
    return f

def eps_ppc(X, T, rho):
    
    T9 = 1e-9 * T 
    
    eps_p = ((2.4e4 * rho * (X**2)) / (T9**(2/3))) * exp(-3.88 / (T9**(1/3)))
             
    return eps_p
  
def eps_cno(X, Z, T, rho):
    
    T9 = 1e-9 * T 
    
    eps_c = ((4.4e25 * rho * X * Z) / (T9**(2/3))) * exp(-15.228 / (T9**(1/3)))
    
    return eps_c

def eps_3alph(T, Y, rho, Z_alph, Z_Be, f1, f2):
    
    T9 = 1e-9 * T 
    
    rho_y_T = ( ( (rho**2) * (Y**3) ) / (T9**3) )
    
    eps_3a = (5.08e8)*f1*f2*rho_y_T*exp(-4.4027/T9)
    
    return eps_3a

    
