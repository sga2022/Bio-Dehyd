# %% 
# Import pyAPEP
# %%
from pyapep import isofit
import numpy as np
import matplotlib.pyplot as plt

# %%
# Generate a data dummy data
# %%
qm_dum = 4.5 # mol/kg
K1_dum = 3.0 # 1/bar
K2_dum = 0.5 # 1/(bar^2)

P_dum = np.linspace(0, 2, 9)
T_dum = [280, 300, 305, 310, 320] # dummy T
dH_dum = 25000 # J/mol : heat of adsorption
R_gas = 8.3145
T_ref = T_dum[1] # Reference temperature = 300 K

q_dum = []

for T in T_dum:
    P_norm = np.exp(dH_dum/R_gas*(1/T - 1/T_ref))*P_dum
    err_dum = 0.05*np.random.rand(len(P_norm))
    q_tmp = qm_dum*(K1_dum*P_norm + 2*K2_dum*P_norm**2)/(
            1+ K1_dum*P_norm + K2_dum*P_norm**2)
    q_w_err = q_tmp + err_dum
    q_dum.append(q_w_err)
# %%
# Use isofit
# %%
P_list = [P_dum,]*len(q_dum)
res_diffT = isofit.fit_diffT(P_list,q_dum, T_dum, 1, tol = 2E-4)
iso_function, iso_param, zerr, fit_error, dH_found, Tref_found, theta_list = res_diffT
print("This data is fitted with {}.".format(zerr))

# %%
# Prediction for different T
# %%
P_pred = np.linspace(0,2.5, 51)
q_pred = []
for T in T_dum:
    q_tmp =iso_function(P_pred, T)
    q_pred.append(q_tmp)
# %% 
# Draw a graph
# %%

fig, ax = plt.subplots()

lines = ('-','--','-.',':', (0,(3,1,0.5,1,0.5)))
for qq_pr, li in zip(q_pred, lines): 
    ax.plot(P_pred, qq_pr,
            linestyle = li,
            color = 'k',
            linewidth = 1.7)

markers = ('o','^','D','v','p','sq')
for ma, qq,  T  in zip(markers, q_dum, T_dum):
    ax.plot(P_dum, qq, ma, markersize = 7.5,
            label = 'T = {0:.1f} K'.format(T),
            mfc = 'w', mec = 'k')
    
ax.legend(fontsize = 13)


ax.set_xlabel('pressure (bar)', fontsize = 13)
ax.set_ylabel('uptake (mol/kg)', fontsize = 13)
ax.grid(which='major', linestyle = ':')
fig.savefig('fig_test_T.png')


