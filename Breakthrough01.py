# %%
# Import packages
from pyapep.simsep import column
import numpy as np
# %% 
# Define the colum
# %%
A_cr = (1.57/2*1E-2)**2*np.pi
N_node = 11
c1 = column(0.3,A_cr, 3, N_node)

print(c1)
# %%
# isotherm function
# %%
par_Psat = [6.75580437e-03, 1.86780002e+01,
 2.57140000e+02, 2.34500000e+02]

par_isoAl = [5.62080479, 0.68058477, 16.61753254]

def Psat(T_K):
    [par0, par1, par2, par3] = par_Psat
    TT = T_K - 273.15
    term1 = par1- TT/par3
    term2 = TT/(par2+TT)
    P_H2O_pred = par0*np.exp(term1*term2)
    return P_H2O_pred

def iso_Al(P, T):
    P_sat = Psat(T)
    P_H2O = P[2]
    P_norm = P_H2O/P_sat
    Xmm, Cc, Kk = par_isoAl
    numo = Xmm*Cc*Kk*P_norm
    deno1 = 1 - Cc*P_norm
    deno2 = 1 - Cc*P_norm + Cc*Kk*P_norm
    q_H2O = numo/deno1/deno2
    #q_H2O = np.zeros(N_node)

    q_CH4 = np.zeros_like(P[0])
    q_CO2 = np.zeros_like(P[1])
    q_res = [q_CH4, q_CO2, q_H2O]
    return q_res

def iso_Al(P, T):
    q_CH4 = np.zeros_like(P[0])
    q_CO2 = np.zeros_like(P[1])
    q_H2O = np.zeros(N_node)
    q_res = [q_CH4, q_CO2, q_H2O]
    return q_res

void_frac = 0.4

c1.adsorbent_info(iso_Al,
                epsi = void_frac, 
                D_particle = 0.01,
                rho_s = 980)
print(c1)
# %%
# Gas property
# %%
c1.gas_prop_info(Mass_molar = [0.014, 0.044, 0.018],
                mu_viscosity= [0.005, 0.005, 0.005])
print(c1)
# %%
# mass_trans_info
c1.mass_trans_info(k_mass_transfer= [0, 0, 0.01],
                    a_specific_surf= 1,
                    D_dispersion = 1E-6)
print(c1)
# %%
# Thermal_info
# %%
c1.thermal_info(dH_adsorption=[0, 0, 5000],
                Cp_solid= 1000,
                Cp_gas = [300, 300, 300],
                h_heat_transfer= 0.2, k_conduct=0.0001,
                h_heat_ambient = 0.0,
                T_ambient = 298.15)
print(c1)
# %%
# Boundary conditions

# %%
# %%
# BT CASE I
# 0.967993 : 0.029938 : 0.002069
# BT CASE II
# 0.968067 : 0.029940 : 0.001996
# BT CASE III
# 0.968151 : 0.029943 : 0.001906

y_in = [0.967993, 0.029938, 0.002069]
u_sup_in = 0.05 # m/s
c1.boundaryC_info(P_outlet = 9.9, P_inlet= 10.1,
                T_inlet = 293.15, y_inlet = y_in,
                Cv_in = 0.1, Cv_out = 0.001, 
                Q_inlet = A_cr*u_sup_in,
                assigned_v_option=True)
print(c1)
# %%
# Initial conditions
# %%
P_init = 10.0*np.ones(N_node)
Tg_init = 293.15*np.ones(N_node)
y_init = [1*np.ones(N_node),
        np.zeros(N_node),
        np.zeros(N_node)]
q_init = [np.zeros(N_node),]*3

c1.initialC_info(P_initial = P_init,
                Tg_initial = Tg_init,
                Ts_initial = Tg_init,
                y_initial = y_init,
                q_initial = q_init)
print(c1)

# %%
# Run the simulation
# %%
c1.run_mamoen(4000,5, True)
# %%
c1.Graph(500, 0)
# %%
c1.Graph_P(100)
# %%
c1.run_mamoen(2,10, True)
# %%
c1.Graph(1, 0)
c1.Graph_P(1)
# %%
