# %%
# Import packages
from pyapep.simsep import column
import numpy as np
# %% 
# Define the colum
# %%
A_cr = (1.57/2*1E-2)**2*np.pi
c1 = column(0.3,A_cr, 2, 101)
print(c1)
# %%
# isotherm function

par_Psat = 

def Psat():



def iso_Al(P, par):
    P_norm = P
    Xmm, Cc, Kk = par
    numo = Xmm*Cc*Kk*P_norm
    deno1 = 1 - Cc*P_norm
    deno2 = 1 - Cc*P_norm + Cc*Kk*P_norm
    q_res = numo/deno1/deno2
    return q_res
#c1.adsorbent_info()