# %%
# Importing isofit
# %%
from pyapep.isofit import best_isomodel
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# %%
# Activated Alumina data "iso_Al_01.csv"
df_Al = pd.read_csv('iso_Al_01.csv')
print(df_Al)
P_Al = np.array(df_Al.iloc[:,0].to_list())
q_Al = np.array(df_Al.iloc[:,1].to_list())

# %%
iso_res = best_isomodel(P_Al,q_Al)

# %%
print(iso_res)
iso_fun = iso_res[0]
P_ran = np.linspace(0,1,101)
q_test = iso_fun(P_ran)

plt.plot(P_Al, q_Al, 'o')
plt.plot(P_ran, q_test)
# %% 
# Def the 
def iso_Al(P, par):
    P_norm = P
    Xmm, Cc, Kk = par
    numo = Xmm*Cc*Kk*P_norm
    deno1 = 1 - Cc*P_norm
    deno2 = 1 - Cc*P_norm + Cc*Kk*P_norm
    q_res = numo/deno1/deno2
    return q_res

def isoerr_Al(par):
    q_pred = iso_Al(P_Al, par)
    diff = q_pred - q_Al
    mse = np.mean(diff**2)
    return mse

par_test = [5.884469, 0.659153666, 1]
#par_test = [15.58667, 5.884469, 0.659153666]
q_ran = iso_Al(P_ran, par_test)
plt.figure()
plt.plot(P_Al, q_Al, 'ko', markersize = 7)
plt.plot(P_ran, q_ran)
# %%
from scipy.optimize import minimize
opt_res = minimize(isoerr_Al, par_test, method = 'Nelder-mead')
#opt_res = minimize(isoerr_Al, par_test, method = 'Powell')

# %%
par_opt = opt_res.x
q_pred_res = iso_Al(P_ran, par_opt)
print(par_opt)
plt.figure()
plt.plot(P_Al, q_Al, 'ko', markersize = 7)
plt.plot(P_ran, q_pred_res)

# %%
