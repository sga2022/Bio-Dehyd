#%% 
# Packages
#%%
import pandas as pd
import numpy as np
#%%
# Calibration data (dewpoints and concentration)
df = pd.read_csv('Dew2conc.csv')
#print(df)
# %%
dewT = np.array(df.iloc[:,0].to_list())
H2Oppm = np.array(df.iloc[:,1].to_list())
#print(dewT)
#print(H2Oppm)