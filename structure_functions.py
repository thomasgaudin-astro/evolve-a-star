import pandas as pd
import numpy as np
from scipy.interpolate import RectBivariateSpline
from mpmath import polylog
from scipy import interpolate

# constants
h = 6.626176e-27 #erg/s #planck's const
me = 9.1094e-28 #electron mass (g)
kB = 1.3807e-16 #cm2 g s-2 K-1 boltzmann const

a = 7.56471e-15 #radiation constant in ergs.cm^-3.K^-4
mp = 1 #proton mass in amu
amu_to_g = 1.66054e-24 #1 amu in g
mp_in_g = mp*amu_to_g #convert mp from amu to g ****************************
ma = mp_in_g

#constants
G = 6.67e-8 #units
NA = 6.022e23 #atoms mol^-1
kb_ev = 8.62e-5 #ev K^-1
m_He = 4.0026 * ma #g
c = 3.0e10 #cm s^-1

# file to interpolate values for integral in equation of state
fermi_dirac = np.loadtxt("../Data/fermi_dirac.txt",skiprows=1)
x = fermi_dirac[0][1:]
y = fermi_dirac[:,0][1:]  #these are Psi values  
z = fermi_dirac[1:][:,1:]
fermi_interp = interpolate.interp2d(x, y, z, kind='cubic')

"""Number density needed to calculate psi"""
def number_density(X,rho):
    NA = 6.02e23 #avogadro's number = 6.02e23 mol^-1,  beacause 1mol = 6.02e23 particles
    mu_e = 2/(1+X)
    n = rho*NA/mu_e
#     print('n =',n, ', rho =', rho)
    return n

"""ψ to be used for the general equation of state"""
def psi(T,X,rho):
    ne = number_density(X,rho)
#     print('ne =',ne)
    ps = np.log( ne*h**3 / ( 2* (2*np.pi*me*kB*T)**(3/2) ) )
    return ps

# use rho to make guess for X, Y and Z (carbon)
#File = np.loadtxt("../Data/structure_00000.txt", dtype=float) #load the structure file
#carbon_mf = np.mean(File[:,32])

# the equation of state...Pressure
#Need to provide Temp, H and He, abundance. Z is the average carbon mass function from the structure file
# make a guess for rho & calculate pressure
def pressure(T,X,rho, Y, Z):
    Psi = psi(T,X,rho)
    mu_I, mu_e = (X/1 + Y/4 + Z/14)**-1 , 2/(1+X) 
    mu = (1/mu_I + 1/mu_e)**-1 #"""mean molecular weight """ 

    nu = 3/2
    integral = ( (2*me*kB*T)**(5/2) )/(2*me) * fermi_interp(nu,Psi)
    pres = (a*T**4)/3 + rho*kB*T/(mu*ma) + 8*np.pi/(3*h**3) * integral
    
    if Psi < -2:
        # ignore the degeneracy part if Psi is very negative
        pres = (a*T**4)/3 + rho*kB*T/(mu*ma) 
    
    return pres


# rho_calc function
#------------------------------------------------------------------------------------------------------------------------------------------------
# function that uses a rho guess to calculate pressure, compares that pressure to henyey pressure by taking difference.
# it perturbs that rho, repeats the above process, then finds the delta_rho (rho guess should be changed by).
def rho_calc(P,T, X, Y, Z, rho):
    mu_I, mu_e = (X/1 + Y/4 + Z/14)**-1 , 2/(1+X) 
    mu = (1/mu_I + 1/mu_e)**-1
    
    Psi = psi(T,X,rho)
    P_hen = P
    
    # perturb rho iteratively until the pressure difference is approximately 0 
    diff = 1000
    while not np.isclose(diff, 0, atol = 0.0):
        P_rho = float(pressure(T,X,rho, Y, Z).real)

        diff = P_hen - P_rho
        
        rho_pert = rho*1.05 # perturb the old rho and use it to calculate a new pressure

        #calculate pressure again using new rho. Redefine P_rho
        P_rho_pert = float(pressure(T,X,rho_pert, Y, Z).real)

        del_rho = (rho_pert - rho)/(P_rho - P_rho_pert) * (P_rho - P_hen)
                
        rho += del_rho
                    
    return rho


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
    
    P_ion = (rho / (mu * ma) ) * kB * T
    
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
