# Contribution of ONE extra bosonic degree of freedom (g=1) to N_eff, standard radiation-era bookkeeping:
# rho_rad = (pi^2/30) [ g_gamma T_g^4 + (7/8) * 2 * N_eff * T_nu^4 ],  T_nu/T_g = (4/11)^(1/3) after e+e- annihilation.
# One boson with g=1 at temperature T_x adds Delta N_eff = (g/2) * (8/7) * (T_x/T_nu)^4 ... written out:
from fractions import Fraction
Tnu_over_Tg = (4/11)**(1/3)
def dneff(g, Tx_over_Tg):
    # rho_x = g*(pi^2/30)*Tx^4 ; one nu species (nu+nubar, g=2 fermionic) = (7/8)*2*(pi^2/30)*Tnu^4
    return g * Tx_over_Tg**4 / ((7/8)*2*Tnu_over_Tg**4)
print("g=1 extra boson sharing photon temperature after e+e-:  dN_eff =", dneff(1, 1.0))
print("g=1 extra boson decoupled with neutrinos (T_x = T_nu):   dN_eff =", dneff(1, Tnu_over_Tg))
print("Planck 2018 + BAO: N_eff = 2.99 +/- 0.17 (68%); SM 3.044 -> dN_eff/sigma:",
      dneff(1,1.0)/0.17, dneff(1,Tnu_over_Tg)/0.17)
