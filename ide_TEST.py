import pyapep.isofit as isofit
import pyapep.simide as simide

# Find hydrogen isotherm (Opt. 1)
# Data import
P = [2, 3, 4, 5]
q = [1, 2, 3, 4]

# Find best isotherm function
H2_isotherm, par_H2, fn_type_H2, val_err_H2 = isofit.best_isomodel(P, q)

# Define nitrogen isotherm (Opt. 2)
# Data import
par_N2 = [2, 0.2, 0.0002]

par1 = [3, 0.1]
par2 = [2, 0.2, 0.0002]

def Lang(par, P, T):
    nume = par[0]*par[1]*P
    deno = 1 + par[1]*P
    q = nume/deno
    return q

def Quad(par, P, T):
    nume = par[0]*(par[1]*P + 2*par[2]*P**2)
    deno = 1 + par[1]*P + par[2]*P**2
    q = nume/deno
    return q

N2_isotherm = lambda P,T: Quad(par_N2, P, T)
H2_isotherm = lambda P,T: Lang(par1, P, T)

iso_list = [H2_isotherm, N2_isotherm]
iso_mix = lambda P,T : isofit.IAST(iso_list, P, T)

iso_mix_test = iso_mix([1,2], 300)
print(iso_mix_test)

CI1 = simide.IdealColumn(2, iso_mix, )

# Feed condition setting
P_feed = 8      # Feed presure (bar)
T_feed = 313.15    # Feed temperature (K)
y_feed = [3/4, 1/4] # Feed mole fraction (mol/mol)
CI1.feedcond(P_feed, T_feed, y_feed)

# Operating condition setting
P_high = 8 # High pressure (bar)
P_low  = 1 # Low pressure (bar)
CI1.opercond(P_high, P_low)
print(CI1)
# Simulation run
x_tail = CI1.runideal()

print(x_tail)       # Output: [x_H2, x_N2]
