import numpy as np
import pandas as pd
from math import exp

def reduced_mass(num):
    
    A = num/(1+num)
    
    return A

Zi_1 = [1, 2, 6]
xi_1 = [X, Y, Z]
Ai_1 = [reduced_mass(1), reduced_mass(2), reduced_mass(6)]

Zi_2 = [1, 2, 4, 6]
xi_2 = [X, Y, Be, Z]
Ai_2 = [reduced_mass(1), reduced_mass(2), reduced_mass(4), reduced_mass(6)]

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

def calc_Be_abund(f_alpha_alph, T, rho, A_alph, A_Be, Y, chi_alph=-91.78e-3):
    
    ma = 1.66e-24 #g
    NA = 6.022e23 #atoms mol^-1
    h = 6.63e-27
    kb = 1.38e-16 #cgs
    kb_ev = 8.62e-5 #ev K^-1
    m_He = 4.0026 * ma #g
    
    x = (2*pi*(m_He/2)*kb*T)/(h**2)
    
    x_Be = ((A_alph**2)/A_Be) * (1/((Y**2)*rho*NA)) * (1/f_alpha_alph) * (x**(-3/2)) * exp(chi_alph/(kb_ev*T))
    
    return x_Be

def eps_ppc(X, T, rho):
    
    T9 = 1e-9 * T 
    
    eps_p = ((2.4e4 * rho * (X**2)) / (T9**(2/3))) * exp(-3.88 / (T9**(1/3)))
             
    return eps_p
  
def eps_cno(X, Z, T, rho):
    
    T9 = 1e-9 * T 
    
    eps_c = ((4.4e25 * rho * X * Z) / (T9**(2/3))) * exp(-15.228 / (T9**(1/3)))
    
    return eps_c

def eps_3alph(T, Y, rho, f1, f2):
    
    T9 = 1e-9 * T 
    
    rho_y_T = ( ( (rho**2) * (Y**3) ) / (T9**3) )
    
    eps_3a = (5.08e8)*f1*f2*rho_y_T*exp(-4.4027/T9)
    
    return eps_3a

def eps_grav(cp, T2, P2, P1, T1, t2, t1, nabla_ad):
    
    eps_g = - (cp*T) ((1/T2)*((T2-T1)/(t2-t1)) - (nabla_ad/P2)((P2-P1)/(t2-t1)))
    
    return eps_g
    
def eps_nuc(eps_p, eps_c, eps_3a):
    
    eps_n = sum(eps_p, eps_c, eps_3a)
    
    return eps_n

def calc_reac_rate(eps_x, Qx, rho):
    
    r_x = (eps_x * rho) / Qx
    
    return r_x

def calc_new_abund(Ai, X1, r_x1=[], r_x2=[], t1, t2, rho):
    #r_x1 - list of creation reaction rates
    #r_x2 - list of destruction reaction rates
    #Ai - species reduced atomic mass
    #X1 - Original abundance
    #t2, t1 - time (seconds)
    
    NA = 6.022e23 #atoms mol^-1
    
    dXi_dt = (Ai / (rho*NA)) * (sum(r_x1)-sum(r_x2))
    
    X2 = X1 + (dXi_dt * (t2-t1))
    
    return X2
