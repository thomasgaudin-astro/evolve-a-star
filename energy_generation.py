import numpy as np
import pandas as pd
from math import exp

ma = 1.66e-24 #g
NA = 6.022e23 #atoms mol^-1
h = 6.63e-27
kb = 1.38e-16 #cgs
kb_ev = 8.62e-5 #ev K^-1
m_He = 4.0026 * ma #g

def reduced_mass(num):
    
    A = num/(1+num)
    
    return A

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

def calc_f(P, T, X, Y, Z, old_rho, Z1, Z2, zeta): 
    
    rho = rho_calc(P, T, X, Y, Z, old_rho)
    
    T6 = 1e-6 * T 
    
    zet_rho = (zeta * rho) / (T6**3)
    
    f = math.exp(0.188 * Z1 * Z2 * (zet_rho**(1/2)))
    
    return f

def calc_Be_abund(P, T, X, Y, Z, old_rho, A_alph, A_Be, chi_alph=-91.78e-3):
    
    rho = rho_calc(P, T, X, Y, Z, old_rho)
    
    #without Be
    Zi_1 = [1, 2, 6]
    xi_1 = [X, Y, Z]
    Ai_1 = [reduced_mass(1), reduced_mass(2), reduced_mass(6)]
    
    zeta_no_Be = calc_zeta(Zi_1, xi_1, Ai_1)
    
    f_alph_alph = calc_f(P, T, X, Y, Z, old_rho, 2, 2, zeta_no_Be)
    
    x = (2*pi*(m_He/2)*kb*T)/(h**2)
    
    x_Be = ((A_alph**2)/A_Be) * (1/((Y**2)*rho*NA)) * (1/f_alpha_alph) * (x**(-3/2)) * exp(chi_alph/(kb_ev*T))
    
    return x_Be

def eps_ppc(P, T, X, Y, Z, old_rho):
    
    T9 = 1e-9 * T 
    
    rho = rho_calc(P, T, X, Y, Z, old_rho)
    
    eps_p = ((2.4e4 * rho * (X**2)) / (T9**(2/3))) * exp(-3.88 / (T9**(1/3)))
             
    return eps_p
  
def eps_cno(P, T, X, Y, Z, old_rho):
    
    T9 = 1e-9 * T 
    
    rho = rho_calc(P, T, X, Y, Z, old_rho)
    
    eps_c = ((4.4e25 * rho * X * Z) / (T9**(2/3))) *exp(-15.228 / (T9**(1/3)))
    
    return eps_c

def eps_3alph(P, T, X, Y, Z, old_rho):
    
    T9 = 1e-9 * T 
    
    rho = rho_calc(P, T, X, Y, Z, old_rho)
    
    #without Be
    Zi_1 = [1, 2, 6]
    xi_1 = [X, Y, Z]
    Ai_1 = [reduced_mass(1), reduced_mass(2), reduced_mass(6)]
    
    zeta_no_Be = calc_zeta(Zi_1, xi_1, Ai_1)
    
    f_alph_alph = calc_f(P, T, X, Y, Z, old_rho, 2, 2, zeta_no_Be)
    
    Be = calc_Be_abund(P, T, X, Y, Z, old_rho, reduced_mass(2), reduced_mass(4), chi_alph=-91.78e-3)
    
    #with Be
    Zi_2 = [1, 2, 4, 6]
    xi_2 = [X, Y, Be, Z]
    Ai_2 = [reduced_mass(1), reduced_mass(2), reduced_mass(4), reduced_mass(6)]
    
    zeta_Be = calc_zeta(Zi_1, xi_1, Ai_1)
    
    f_alph_Be = calc_f(P, T, X, Y, Z, old_rho 2, 4, zeta_Be)
    
    rho_y_T = ( ( (rho**2) * (Y**3) ) / (T9**3) )
    
    eps_3a = (5.08e8)*f1*f2*rho_y_T*exp(-4.4027/T9)
    
    return eps_3a

def eps_grav(P, T, X, Y, Z, old_rho, P_old, T_old, time_step):
    
    cp = cp_calc(P, T, X, Y, Z, old_rho)
    nabla_ad = nabla_ad_calc(P, T, X, Y, Z, old_rho)
    
    eps_g = - (cp*T) ((1/)*((T-T_old)/(time_step)) - (nabla_ad/P)((P-P_old)/(time_step)))
    
    return eps_g
    
def eps_nuc(P, T, X, Y, Z, old_rho):
    
    eps_p = eps_ppc(P, T, X, Y, Z, old_rho)
    eps_c = eps_cno(P, T, X, Y, Z, old_rho)
    eps_3a = eps_3alph(P, T, X, Y, Z, old_rho)
    
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
