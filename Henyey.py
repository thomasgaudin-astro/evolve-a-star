import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd


# Constants
G = 1
A_RAD = 1
C = 1
M_SUN = 1
M_T = 1
SIGMA_SB = 1
DEL_T = 1


# Dummy functions for various parameters
def rho_calc(pre, tem, x, y, z, rhopre):
    return()

def eps_nuc_calc(pre, tem, x, y, z, rhopre):
    return()

def eps_gra_calc(pre, tem, x, y, z, rhopre, prepre, tempre):
    return()

def cp_calc(pre, tem, x, y, z, rhopre):
    return()

def delta_calc(pre, tem, x, y, z, rhopre):
    return()

def nabla_calc(pre, lum, tem, mas, x, y, z, rhopre, filename):
    return()

def kappa_calc(pre, tem, x, y, z, rhopre, filename):
    return()

def calc_new_abund(x, y, z, rhopre, eps_p, eps_c, eps_3a, time_step):
    return()


# Function to calculate derivative
def derivative_calc(param, pre_param):
    deriv = (param-pre_param)/DEL_T
    return(deriv)


# Functions to calculate A values
# Pressure at the surface of the star
def BP(re, Pe, Te, x, y, z, rhopre):
    kappa_s = kappa_calc(Pe, Te, x, y, z, rhopre, './opacity_table.csv')
    bP = (2/3) * (G*M_T*M_SUN)/(kappa_s*re**2)
    return(bP)

# Temperature at the surface of the star
def BT(re, Le):
    bT = (Le / (4*np.pi*SIGMA_SB*re**2))**(1/4)
    return(bT)

# Radius disparity measure at shell j > 0
def Aj1(rj, Pj, Tj, Mj, rj1, Pj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1):
    rhoj = rho_calc(Pj, Tj, x, y, z, x1, y1, z1, rhopre)
    rhoj1 = rho_calc(Pj1, Tj1, x, y, z, x1, y1, z1, rhopre1)
    aj1 = (rj1-rj)/(Mj1-Mj) - 1/(8*np.pi)*(1/(rj1**2*rhoj1) + 1/(rj**2*rhoj))
    return(aj1)

# Pressure disparity measure at shell j > 0
def Aj2(rj, Pj, Mj, rj1, Pj1, Mj1):
    aj2 = (Pj1-Pj)/(Mj1-Mj) + G/(8*np.pi)*((Mj1)/(rj1**4) + (Mj)/(rj**4))
    return(aj2)

# Luminosity disparity measure at shell j > 0
def Aj3(Pj, Lj, Tj, Mj, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1, Pprej, Tprej, Pprej1, Tprej1):
    eps_nucj = eps_nuc_calc(Pj, Tj, x, y, z, rhopre)
    eps_nucj1 = eps_nuc_calc(Pj1, Tj1, x1, y1, z1, rhopre1)
    cpj = cp_calc(Pj, Tj, x, y, z, rhopre)
    cpj1 = cp_calc(Pj1, Tj1, x1, y1, z1, rhopre1)
    del_Tj = derivative_calc(Tj, Tprej)
    del_Tj1 = derivative_calc(Tj1, Tprej1)
    deltaj = delta_calc(Pj, Tj, x, y, z, rhopre)
    deltaj1 = delta_calc(Pj1, Tj1, x1, y1, z1, rhopre1)
    rhoj = rho_calc(Pj, Tj, x, y, z, rhopre)
    rhoj1 = rho_calc(Pj1, Tj1, x1, y1, z1, rhopre1)
    del_Pj = derivative_calc(Pj, Pprej)
    del_Pj1 = derivative_calc(Pj1, Pprej1)
    aj3 = (Lj1-Lj)/(Mj1-Mj) - 1/2*(eps_nucj1 + eps_nucj) + 1/2*(cpj1*del_Tj1 + cpj*del_Tj) - 1/2*(deltaj1/rhoj1*del_Pj1 + deltaj/rhoj*del_Pj)
    return(aj3)

# Temperature disparity measure at shell j > 0
def Aj4(rj, Pj, Lj, Tj, Mj, rj1, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1):
    nablaj = nabla_calc(Pj, Lj, Tj, Mj, x, y, z, rhopre, './opacity_table.csv')
    nablaj1 = nabla_calc(Pj1, Lj1, Tj1, Mj1, x1, y1, z1, rhopre1, './opacity_table.csv')
    aj4 = (Tj1-Tj)/(Mj1-Mj) + G/(8*np.pi)*((Tj1*Mj1*nablaj1)/(Pj1*rj1**4) + (Tj*Mj*nablaj)/(Pj*rj**4))
    return(aj4)

# Radius disparity measure at core
def A01(pre0, tem0, rad1, mas1, x, y, z, rhopre):
    rho0 = rho_calc(pre0, tem0, x, y, z, rhopre)
    a01 = rad1 - ((3*mas1)/(4*np.pi*rho0))**(1/3)
    return(a01)

# Pressure disparity measure at core
def A02(pre0, tem0, pre1, mas1, x, y, z, rhopre):
    rho0 = rho_calc(pre0, tem0, x, y, z, rhopre)
    a02 = pre1 - pre0 + (3*G)/(8*np.pi)*((4*np.pi*rho0)/3)**(4/3)*mas1**(2/3)
    return(a02)

# Luminosity disparity measure at core
def A03(pre0, tem0, lum1, mas1, x, y, z, rhopre, pre0pre, tem0pre):
    eps_nuc0 = eps_nuc_calc(pre0, tem0, x, y, z, rhopre)
    eps_gra0 = eps_gra_calc(pre0, tem0, x, y, z, rhopre, pre0pre, tem0pre)
    a04 = lum1 - (eps_nuc0+eps_gra0)*mas1
    return(a04)

# Temperature disparity measure at core
def A04(pre0, tem0, tem1, mas1, x, y, z, rhopre, pre0pre, tem0pre):
    rho0 = rho_calc(pre0, tem0, x, y, z, rhopre)
    kappa0 = kappa_calc(pre0, tem0, x, y, z, rhopre, './opacity_table.csv')
    eps_nuc0 = eps_nuc_calc(pre0, tem0, x, y, z, rhopre)
    eps_gra0 = eps_gra_calc(pre0, tem0, x, y, z, rhopre, pre0pre, tem0pre)
    a03 = tem1**4 - tem0**4 + (kappa0*(eps_nuc0+eps_gra0))/(2*A_RAD*C)*(3/(4*np.pi))**(2/3)*rho0**(4/3)*mas1**(2/3)
    return(a03)


# Partial derivative finders
# Partial derivatives of BP
def partial_BP(re, Pe, Te, delre, delPe, delTe, x, y, z, rhopre):
    bP = BP(re, Pe, Te, x, y, z, rhopre)
    part_re = (bP(re+delre, Pe, Te, x, y, z, rhopre)-bP)/delre
    part_Pe = (bP(re, Pe+delPe, Te, x, y, z, rhopre)-bP)/delPe
    part_Le = 0
    part_Te = (bP(re, Pe, Te+delTe, x, y, z, rhopre)-bP)/delTe
    return(np.array([part_re, part_Pe, part_Le, part_Te]))

# Partial derivatives of BT
def partial_BT(re, Le, delre, delLe):
    bT = BT(re, Le)
    part_re = (BT(re+delre, Le)-bT)/delre
    part_Pe = 0
    part_Le = (BT(re, Le+delLe)-bT)/delLe
    part_Te = 0
    return(np.array([part_re, part_Pe, part_Le, part_Te]))

# Partial derivatives of Aj1
def partial_Aj1(rj, Pj, Tj, Mj, rj1, Pj1, Tj1, Mj1, delrj, delPj, delTj, delrj1, delPj1, delTj1, x, y, z, x1, y1, z1, rhopre, rhopre1):
    aj1 = Aj1(rj, Pj, Tj, Mj, rj1, Pj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)
    part_rj = (Aj1(rj+delrj, Pj, Tj, Mj, rj1, Pj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delrj
    part_Pj = (Aj1(rj, Pj+delPj, Tj, Mj, rj1, Pj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delPj
    part_Lj = 0
    part_Tj = (Aj1(rj, Pj, Tj+delTj, Mj, rj1, Pj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delTj
    part_rj1 = (Aj1(rj, Pj, Tj, Mj, rj1+delrj1, Pj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delrj1
    part_Pj1 = (Aj1(rj, Pj, Tj, Mj, rj1, Pj1+delPj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delPj1
    part_Lj1 = 0
    part_Tj1 = (Aj1(rj, Pj, Tj, Mj, rj1, Pj1, Tj1+delTj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delTj1
    return(np.array([part_rj, part_Pj, part_Lj, part_Tj, part_rj1, part_Pj1, part_Lj1, part_Tj1]))

# Partial derivatives of Aj2
def partial_Aj2(rj, Pj, Mj, rj1, Pj1, Mj1, delrj, delPj, delrj1, delPj1):
    aj1 = Aj2(rj, Pj, Mj, rj1, Pj1, Mj1)
    part_rj = (Aj2(rj+delrj, Pj, Mj, rj1, Pj1, Mj1)-aj1)/delrj
    part_Pj = (Aj2(rj, Pj+delPj, Mj, rj1, Pj1, Mj1)-aj1)/delPj
    part_Lj = 0
    part_Tj = 0
    part_rj1 = (Aj2(rj, Pj, Mj, rj1+delrj1, Pj1, Mj1)-aj1)/delrj1
    part_Pj1 = (Aj2(rj, Pj, Mj, rj1, Pj1+delPj1, Mj1)-aj1)/delPj1
    part_Lj1 = 0
    part_Tj1 = 0
    return(np.array([part_rj, part_Pj, part_Lj, part_Tj, part_rj1, part_Pj1, part_Lj1, part_Tj1]))

# Partial derivatives of Aj3
def partial_Aj3(Pj, Lj, Tj, Mj, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1, Pprej, Tprej, Pprej1, Tprej1, delPj, delLj, delTj, delPj1, delLj1, delTj1):
    aj1 = Aj3(Pj, Lj, Tj, Mj, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1, Pprej, Tprej, Pprej1, Tprej1)
    part_rj = 0
    part_Pj = (Aj3(Pj+delPj, Lj, Tj, Mj, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1, Pprej, Tprej, Pprej1, Tprej1)-aj1)/delPj
    part_Lj = (Aj3(Pj, Lj+delLj, Tj, Mj, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1, Pprej, Tprej, Pprej1, Tprej1)-aj1)/delLj
    part_Tj = (Aj3(Pj, Lj, Tj+delTj, Mj, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1, Pprej, Tprej, Pprej1, Tprej1)-aj1)/delTj
    part_rj1 = 0
    part_Pj1 = (Aj3(Pj, Lj, Tj, Mj, Pj1+delPj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1, Pprej, Tprej, Pprej1, Tprej1)-aj1)/delPj1
    part_Lj1 = (Aj3(Pj, Lj, Tj, Mj, Pj1, Lj1+delLj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1, Pprej, Tprej, Pprej1, Tprej1)-aj1)/delLj1
    part_Tj1 = (Aj3(Pj, Lj, Tj, Mj, Pj1, Lj1, Tj1+delTj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1, Pprej, Tprej, Pprej1, Tprej1)-aj1)/delTj1
    return(np.array([part_rj, part_Pj, part_Lj, part_Tj, part_rj1, part_Pj1, part_Lj1, part_Tj1]))

# Partial derivatives of Aj4
def partial_Aj4(rj, Pj, Lj, Tj, Mj, rj1, Pj1, Lj1, Tj1, Mj1, delrj, delPj, delLj, delTj, delrj1, delPj1, delLj1, delTj1, x, y, z, x1, y1, z1, rhopre, rhopre1):
    aj1 = Aj4(rj, Pj, Lj, Tj, Mj, rj1, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)
    part_rj = (Aj4(rj+delrj, Pj, Lj, Tj, Mj, rj1, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delrj
    part_Pj = (Aj4(rj, Pj+delPj, Lj, Tj, Mj, rj1, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delPj
    part_Lj = (Aj4(rj, Pj, Lj+delLj, Tj, Mj, rj1, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delLj
    part_Tj = (Aj4(rj, Pj, Lj, Tj+delTj, Mj, rj1, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delTj
    part_rj1 = (Aj4(rj, Pj, Lj, Tj, Mj, rj1+delrj1, Pj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delrj1
    part_Pj1 = (Aj4(rj, Pj, Lj, Tj, Mj, rj1, Pj1+delPj1, Lj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delPj1
    part_Lj1 = (Aj4(rj, Pj, Lj, Tj, Mj, rj1, Pj1, Lj1+delLj1, Tj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delLj1
    part_Tj1 = (Aj4(rj, Pj, Lj, Tj, Mj, rj1, Pj1, Lj1, Tj1+delTj1, Mj1, x, y, z, x1, y1, z1, rhopre, rhopre1)-aj1)/delTj1
    return(np.array([part_rj, part_Pj, part_Lj, part_Tj, part_rj1, part_Pj1, part_Lj1, part_Tj1]))

# Partial derivatives of A01
def partial_A01(pre0, tem0, rad1, mas1, delpre0, deltem0, delrad1, x, y, z, rhopre):
    a01 = A01(pre0, tem0, rad1, mas1, x, y, z, rhopre)
    partial_pre0 = (A01(pre0+delpre0, tem0, rad1, mas1, x, y, z, rhopre)-a01)/delpre0
    partial_tem0 = (A01(pre0, tem0+deltem0, rad1, mas1, x, y, z, rhopre)-a01)/deltem0
    partial_rad1 = (A01(pre0, tem0, rad1+delrad1, mas1, x, y, z, rhopre)-a01)/delrad1
    partial_pre1 = 0
    partial_lum1 = 0
    partial_tem1 = 0
    return(np.array([partial_pre0, partial_tem0, partial_rad1, partial_pre1, partial_lum1, partial_tem1]))

# Partial derivatives of A02
def partial_A02(pre0, tem0, pre1, mas1, delpre0, deltem0, delpre1, x, y, z, rhopre):
    a02 = A02(pre0, tem0, pre1, mas1, x, y, z, rhopre)
    partial_pre0 = (A02(pre0+delpre0, tem0, pre1, mas1, x, y, z, rhopre)-a02)/delpre0
    partial_tem0 = (A02(pre0, tem0+deltem0, pre1, mas1, x, y, z, rhopre)-a02)/deltem0
    partial_rad1 = 0
    partial_pre1 = (A02(pre0, tem0, pre1+delpre1, mas1, x, y, z, rhopre)-a02)/delpre1
    partial_lum1 = 0
    partial_tem1 = 0
    return(np.array([partial_pre0, partial_tem0, partial_rad1, partial_pre1, partial_lum1, partial_tem1]))

# Partial derivatives of A03
def partial_A03(pre0, tem0, lum1, mas1, delpre0, deltem0, dellum1, x, y, z, rhopre, pre0pre, tem0pre):
    a03 = A03(pre0, tem0, lum1, mas1, x, y, z, rhopre, pre0pre, tem0pre)
    partial_pre0 = (A03(pre0+delpre0, tem0, lum1, mas1, x, y, z, rhopre, pre0pre, tem0pre)-a03)/delpre0
    partial_tem0 = (A03(pre0, tem0+deltem0, lum1, mas1, x, y, z, rhopre, pre0pre, tem0pre)-a03)/deltem0
    partial_rad1 = 0
    partial_pre1 = 0
    partial_lum1 = (A03(pre0, tem0, lum1+dellum1, mas1, x, y, z, rhopre, pre0pre, tem0pre)-a03)/dellum1
    partial_tem1 = 0
    return(np.array([partial_pre0, partial_tem0, partial_rad1, partial_pre1, partial_lum1, partial_tem1]))

# Partial derivatives of A04
def partial_A04(pre0, tem0, tem1, mas1, delpre0, deltem0, deltem1, x, y, z, rhopre, pre0pre, tem0pre):
    a04 = A04(pre0, tem0, tem1, mas1, x, y, z, rhopre, pre0pre, tem0pre)
    partial_pre0 = (A04(pre0+delpre0, tem0, tem1, mas1, x, y, z, rhopre, pre0pre, tem0pre)-a04)/delpre0
    partial_tem0 = (A04(pre0, tem0+deltem0, tem1, mas1, x, y, z, rhopre, pre0pre, tem0pre)-a04)/deltem0
    partial_rad1 = 0
    partial_pre1 = 0
    partial_lum1 = 0
    partial_tem1 = (A04(pre0, tem0, tem1+deltem1, mas1, x, y, z, rhopre, pre0pre, tem0pre)-a04)/deltem1
    return(np.array([partial_pre0, partial_tem0, partial_rad1, partial_pre1, partial_lum1, partial_tem1]))


# Main function to gather information together
def henyey():
    # Choosing what structure file to read-in and reading it in (reversing order so core is first)
    pre_num = 0
    structure_df = pd.read_csv(f'structure_{pre_num:05}.txt',sep='\s+')[::-1]

    # Converting relevant structure columns to arrays
    pre_mass_array = structure_df['M_r'].to_numpy()
    pre_radius_array = structure_df['r'].to_numpy()
    pre_pressure_array = structure_df['P'].to_numpy()
    pre_luminosity_array = structure_df['L_r'].to_numpy()
    pre_temperature_array = structure_df['T'].to_numpy()
    pre_rho_array = structure_df['rho'].to_numpy()
    pre_X_array = structure_df['X'].to_numpy()
    pre_Y_array = structure_df['Y'].to_numpy()
    pre_Z_array = 1 - pre_X_array - pre_Y_array
    pre_eps_pp_array = structure_df['epsilon_pp'].to_numpy()
    pre_eps_CNO_array = structure_df['epsilon_CNO'].to_numpy()
    pre_eps_3alpha_array = structure_df['epsilon_3alpha'].to_numpy()

    # Using values from the previous stage to calculate new abundances
    X_array, Y_array, Z_array = calc_new_abund(pre_X_array, pre_Y_array, pre_Z_array, pre_rho_array, pre_eps_pp_array, pre_eps_CNO_array, pre_eps_3alpha_array, DEL_T)

    # Initializing Henyey matrix H and vector A to appropriate sizes
    henyey_matrix = np.zeros((4*len(pre_mass_array)-2, 4*len(pre_mass_array)-2))
    henyey_vector = np.zeros(4*len(pre_mass_array)-2)
    
    # Initializing 'guesses' for the new values of radius, pressure, luminosity, and temperature
    radius_array = pre_radius_array*1.01
    pressure_array = pre_pressure_array*1.01
    luminosity_array = pre_luminosity_array*1.01
    temperature_array = pre_temperature_array*1.01

    # Running Henyey calculation
    while not np.isclose(np.sum(np.abs(henyey_vector)), 0):
        henyey_matrix = np.zeros((4*len(pre_mass_array)-2, 4*len(pre_mass_array)-2))
        henyey_vector = np.zeros(4*len(pre_mass_array)-2)

        for j, _ in enumerate(pre_mass_array):
            radj = radius_array[j]
            prej = pressure_array[j]
            lumj = luminosity_array[j]
            temj = temperature_array[j]
            masj = pre_mass_array[j]

            delradj = 0.05*radj
            delprej = 0.05*prej
            dellumj = 0.05*lumj
            deltemj = 0.05*temj

            xj = X_array[j]
            yj = Y_array[j]
            zj = Z_array[j]

            rhojpre = pre_rho_array[j]
            prejpre = pre_pressure_array[j]
            temjpre = pre_temperature_array[j]
            if j == (len(pre_mass_array)-1):
                henyey_vector[-2] = prej - BP(radj, prej, temj, xj, yj, zj, rhojpre)
                henyey_vector[-1] = temj - BT(radj, lumj)

                henyey_matrix[-2:-1,-4:] = partial_BP(radj, prej, temj, delradj, delprej, deltemj, xj, yj, zj, rhojpre)
                henyey_matrix[-1:,-4:] = partial_BT(radj, lumj, delradj, dellumj)
            else:
                radj1 = radius_array[j+1]
                prej1 = pressure_array[j+1]
                lumj1 = luminosity_array[j+1]
                temj1 = temperature_array[j+1]
                masj1 = pre_mass_array[j+1]

                delradj1 = 0.05*radj1
                delprej1 = 0.05*prej1
                dellumj1 = 0.05*lumj1
                deltemj1 = 0.05*temj1

                xj1 = X_array[j+1]
                yj1 = Y_array[j+1]
                zj1 = Z_array[j+1]

                rhoj1pre = pre_rho_array[j+1]
                prej1pre = pre_pressure_array[j+1]
                temj1pre = pre_temperature_array[j+1]
                if j ==0:
                    henyey_vector[4*j] = -A01(prej, temj, radj1, masj1, xj, yj, zj, rhojpre)
                    henyey_vector[4*j+1] = -A02(prej, temj, prej1, masj1, xj, yj, zj, rhojpre)
                    henyey_vector[4*j+2] = -A03(prej, temj, lumj1, masj1, xj, yj, zj, rhojpre, prejpre, temjpre)
                    henyey_vector[4*j+3] = -A04(prej, temj, temj1, masj1, xj, yj, zj, rhojpre, prejpre, temjpre)

                    henyey_matrix[4*j:4*j+1,:6] = partial_A01(prej, temj, radj1, masj1, delprej, deltemj, delradj1, xj, yj, zj, rhojpre)
                    henyey_matrix[4*j+1:4*j+2,:6] = partial_A02(prej, temj, prej1, masj1, delprej, deltemj, delprej1, xj, yj, zj, rhojpre)
                    henyey_matrix[4*j+2:4*j+3,:6] = partial_A03(prej, temj, lumj1, masj1, delprej, deltemj, dellumj1, xj, yj, zj, rhojpre, prejpre, temjpre)
                    henyey_matrix[4*j+3:4*j+4,:6] = partial_A03(prej, temj, lumj1, masj1, delprej, deltemj, dellumj, xj, yj, zj, rhojpre, prejpre, temjpre) 
                else:
                    henyey_vector[4*j] = -Aj1(radj, prej, temj, masj, radj1, prej1, temj1, masj1, xj, yj, zj, xj1, yj1, zj1, rhojpre, rhoj1pre)
                    henyey_vector[4*j+1] = -Aj2(radj, prej, masj, radj1, prej1, masj1)
                    henyey_vector[4*j+2] = -Aj3(prej, lumj, temj, masj, prej1, lumj1, temj1, masj1, xj, yj, zj, xj1, yj1, zj1, rhojpre, rhoj1pre, prejpre, temjpre, prej1pre, temj1pre)
                    henyey_vector[4*j+3] = -Aj4(radj, prej, lumj, temj, masj, radj1, prej1, lumj1, temj1, masj1, xj, yj, zj, xj1, yj1, zj1, rhojpre, rhoj1pre)

                    henyey_matrix[4*j:4*j+1,2+4*(j-1):10+4*(j-1)] = partial_Aj1(radj, prej, temj, masj, radj1, prej1, temj1, masj1, delradj, delprej, deltemj, delradj1, delprej1, deltemj1, xj, yj, zj, xj1, yj1, zj1, rhojpre, rhoj1pre)
                    henyey_matrix[4*j+1:4*j+2,2+4*(j-1):10+4*(j-1)] = partial_Aj2(radj, prej, masj, radj1, prej1, masj1, delradj, delprej, delradj1, delprej1)
                    henyey_matrix[4*j+2:4*j+3,2+4*(j-1):10+4*(j-1)] = partial_Aj3(prej, lumj, temj, masj, prej1, lumj1, temj1, masj1, xj, yj, zj, xj1, yj1, zj1, rhojpre, rhoj1pre, prejpre, temjpre, prej1pre, temj1pre, delprej, dellumj, deltemj, delprej1, dellumj1, deltemj1)
                    henyey_matrix[4*j+3:4*j+4,2+4*(j-1):10+4*(j-1)] = partial_Aj4(radj, prej, lumj, temj, masj, radj1, prej1, lumj1, temj1, masj1, delradj, delprej, dellumj, deltemj, delradj1, delprej1, dellumj1, deltemj1, xj, yj, zj, xj1, yj1, zj1, rhojpre, rhoj1pre)
        
        correction_vector = np.dot(np.linalg.inv(henyey_matrix), henyey_vector)*0.5
        correction_vector = np.insert(correction_vector, 0, 0)
        correction_vector = np.insert(correction_vector, 2, 0)

        radius_array = radius_array + correction_vector[::4]
        pressure_array = pressure_array + correction_vector[1::4]
        luminosity_array = luminosity_array + correction_vector[2::4]
        temperature_array = temperature_array + correction_vector[3::4]

    return()

# Main for argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'Final Project',
        description = 'Henyey method for evolution of a star'
        )
    args = parser.parse_args()

    henyey()
