# %%
# Data importing
# %%
from Dataimport import dewT, H2Oppm
import numpy as np
import matplotlib.pyplot as plt
# %%
# dewT in degree C
# H2Oppm in ppm (weight-based)

dewT_K = dewT + 273.15 # T in K (from degree C of dewT)
H2O_mol = H2Oppm*1E-1/8.31446/293.15 # mol/m^3 (at 20 dgr C)
molfr_H2O = H2Oppm/1E6 # mole fraction (mol/mol)
P_H2O = molfr_H2O # when 1 bar is applied
print(H2O_mol)
print(molfr_H2O)

plt.plot(dewT, molfr_H2O)
plt.plot(0, 0.006, 'o', markersize = 7)
#plt.ylabel('mole fraction (mol/mol)')
plt.ylabel('saturated vapor pressure (bar)')
plt.show()
# %%
# Model
# %%
def T2P_dew(TT, par):
    [par0, par1, par2, par3] = par
    term1 = par1- TT/par3
    term2 = TT/(par2+TT)
    P_H2O_pred = par0*np.exp(term1*term2)
    return P_H2O_pred
par_test = [6.1121E-3, 18.678, 257.14, 234.5]
Td_test = np.linspace(-80, 15, 48)
Pd_test = T2P_dew(Td_test, par_test)
plt.figure()
plt.plot(Td_test, Pd_test, 'o')
plt.plot(dewT, molfr_H2O)
### From this results, we can see this model fits the results very well
### Find the right model from this data
# %%
# Parameter estimation from the data
# %%
from scipy.optimize import minimize
def err(par):
    P_pred = T2P_dew(dewT, par)# predicted H2O satured vapor pressure
    diff = P_pred - molfr_H2O*1
    mse = np.mean(diff**2)
    return mse
opt_res = minimize(err, par_test)
print(opt_res)
# %%
print(opt_res.x)
par_opt = opt_res.x
T_dom = np.linspace(dewT[0], dewT[-1])
P_pred_opt = T2P_dew(T_dom, par_opt)
plt.figure()
plt.plot(dewT[::2], molfr_H2O[::2], 'ko')
plt.plot(T_dom, P_pred_opt)
# %%
# Pickle the model into the model
