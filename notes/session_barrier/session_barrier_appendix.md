# Appendix — Barrier / Stability / Materials Thread

This appendix belongs to the same session write-up as Sections 0–8. Its purpose is practical rather than rhetorical: it collects the symbol dictionary, reduced-unit bookkeeping, artifact inventory, source-document roles, and the explicit non-claims / open theorem gaps that keep the thread readable after the fact.

The appendix is scoped to the **barrier-crossing / stability / damping / materials** chain developed in this session. It is not a complete glossary for the entire toy-model program.

---

## Appendix A. Symbol glossary

### A.1 Parent-theory and brane-reduction symbols

| Symbol | Meaning in this thread |
|---|---|
| \(x^M=(t,x,y,z,w)\) | Bulk spacetime coordinates of the parent \(4+1\) theory. |
| \(\mathbf X=(x,y,z,w)\) | Bulk spatial coordinate. |
| \(\mathbf x=(x,y,z)\) | Brane spatial coordinate. |
| \(w\) | Extra transverse / bulk coordinate. |
| \(\psi(\mathbf X,t)\) | Gauged GNLS matter/order parameter field. |
| \(\rho=|\psi|^2\) | Bulk matter density. |
| \(A_M\) | Localized \(4+1\) gauge potential. |
| \(F_{MN}\) | Gauge-field strength tensor. |
| \(Z(w)\) | EM localization profile appearing in the action. |
| \(Z_{\rm int}=\int Z(w)\,dw\) | Integrated localization weight controlling brane-effective EM coupling. |
| \(W(w)\) | Projection kernel defining what the brane observer measures. |
| \(\mu_0^{\rm eff}=\mu_0/Z_{\rm int}\) | Controlled brane-effective Maxwell coupling. |
| \(\eta_Q\) | Fixed topological puncture orientation carrying electric-charge sign. |
| \(e_*>0\) | Positive microscopic charge scale. |
| \(q_*=\eta_Q e_*\) | Signed microscopic branch charge. |
| \(q_{\rm eff}=q_*/\sqrt{Z_{\rm int}}\) | Observable brane charge after canonical zero-mode normalization. |
| \(j^A\) | Bulk number-current components, including \(j^w\). |
| \(S_{\rm leak}\) | Exact projected leakage source in brane continuity. |
| \(\varphi\) | Longitudinal brane velocity potential appearing in the exact longitudinal identity. |

### A.2 Mixed-sector / beyond-MHD symbols

| Symbol | Meaning in this thread |
|---|---|
| \(A_w\) | Extra gauge component that is suppressed only in the strict far-field brane limit. |
| \(F_{\mu w}\) | Mixed field-strength components. |
| \(E_w=-\partial_t A_w-\partial_w A_0\) | Transverse electric field. |
| \(C_a=F_{aw}\) | Mixed spatial field component. |
| \(J^w\) | Transverse electric current. |
| \(J^wE_w\) | Scalar-photon work channel that tracks energy export into the hidden sector. |
| \(v^w\) | Transverse velocity component driven by \(E_w\) when \(J^w\neq 0\). |
| \(h_{\rm sub}\) | Subscale / unresolved helicity, defined from projected minus resolved helicity. |
| \(\partial_t h_{\rm sub}\) | Helicity-export rate used as the spin-orientation diagnostic in Session II. |

### A.3 Moving-throat and bundle symbols

| Symbol | Meaning in this thread |
|---|---|
| \(R(\Omega,w,t)\) | Moving-throat shape field in the geometry lift. |
| \(\Sigma=r-R(\Omega,w,t)\) | Level-set representation of the throat surface. |
| \(a(t),L(t)\) | Old collective geometry variables, reinterpreted as low moments of the distributed throat. |
| \(q_{2m}\) | Grouped real \(P_2\) wall/support amplitudes in the moving-throat program. |
| \(K_*=K-C^2/\varpi^2\) | Effective wall stiffness after integrating out a stable support mode. |
| \(\Delta=\Omega_U^2\Omega_W^2-R^2\) | Internal one-port Maxwell/mixed determinant. |
| \(Q=G_U^2\Omega_W^2+2G_UG_WR+G_W^2\Omega_U^2\) | Static one-port numerator entering the conservative wall operator. |
| \(P=\Omega_U^2G_W+RG_U\) | Static outgoing-load factor in the one-port mixed bundle. |
| \(D_0=K_*-Q/\Delta\) | Static conservative wall operator in the one-port bundle. |
| \(N_0=P^2/\Delta^2\) | Static outgoing-transfer coefficient. |
| \(P_0=N_0/D_0\) | Static normalized outgoing prefactor that also enters the quadrupole-normalization chain. |
| \(m_{\hat 0}\) | Source-map factor in the outgoing-normalization product \(m_{\hat 0}^2P_0\). |

### A.4 Session-I stationary barrier symbols

| Symbol | Meaning in this thread |
|---|---|
| \(r\) | Reduced same-charge separation coordinate. |
| \(r_{\rm reg}\) | Short-distance regularization length used in the reduced barrier ansatz. |
| \(\lambda\) | Reduced localization / width parameter in the stress-test closures. |
| \(U\) | Logarithmic transfer-shape coordinate introduced when rigid-mouth factorization is relaxed. |
| \(V\) | Logarithmic dressing-leg coordinate introduced together with \(U\). |
| \(\epsilon_\eta\) | Dressing variable whose invariance is broken once \(U/V\) coupling is allowed. |
| \(\Xi_1\) | First-order dynamic barrier scalar; in Session I it matched the reduced \(U(r)\) branch. |
| \(\widetilde\Xi\) | Finite nonlinear rigid-mouth proxy carried in the stationary script. |
| \(a_U,a_V\) | Quadratic free-energy coefficients for the reduced \(U,V\) closure. |
| \(g_{UV}\) | Cross-coupling strength between \(U\) and \(V\). |
| \(s_0\) | Amplitude scale in the reduced transfer-source profile. |
| \(\sigma(z)\) | Axial mouth/core source profile, relaxed to a sign-changing multimode branch. |
| \(\mathfrak g[\sigma]\) | Mouth-bias functional of the axial source profile. |
| \(\mathcal R[\sigma]\) | Mixed-to-shell loading ratio used in the Family-1 source stress test. |
| \(\beta_Q,\beta_U,\beta_W\) | Primitive short-range source amplitudes entering the reduced barrier kernel. |
| \(\kappa\) | Yukawa inverse-range parameter appearing in \(e^{-2\kappa x}\)-type families. |
| \(V_{\rm eff}(r)\) | Reduced same-charge effective potential after the three constraint relaxations are turned on. |
| \(V_{\rm Coul}(r)=1/r\) | Pure Coulomb comparison branch. |

### A.5 Session-II dynamic and helicity symbols

| Symbol | Meaning in this thread |
|---|---|
| \(m_s\) | Reduced particle mass in the dynamic scattering problem. |
| \(\hbar_{\rm eff}\) | Reduced Planck-like scale entering the WKB formulas. |
| \(r_{\rm turn}\) | Classical turning point on the reduced potential branch. |
| \(V_{\rm peak}\) | Height of the reduced barrier peak. |
| \(r_{\rm peak}\) | Location of the reduced barrier peak. |
| \(E_{\rm sub}\) | Representative subbarrier energy used for the WKB comparison. |
| \(I_{\rm WKB}\) | WKB action integral on the reduced branch. |
| \(T_{\rm new},T_{\rm Coul}\) | Reduced tunneling weights for the lowered branch and pure Coulomb comparison. |
| \(\Omega_{ij}=-(q_s/m_s)F_{ij}\) | Canonical-vorticity identity used to define the magnetic/vortical sector. |
| aligned / anti-aligned | The two reduced spin closures compared in the helicity-export diagnostic. |
| \(\chi_\lambda=\lambda |\partial_r\ln V_{\rm eff}|\) | Gradient trigger for beyond-MHD / transverse-bypass behavior. |
| \(\lambda_{\rm th}(r_{\rm turn})\) | Threshold width extracted from the turning-point gradient. |

### A.6 Session-III stability symbols

| Symbol | Meaning in this thread |
|---|---|
| \(t_{\rm cross}\) | Characteristic crossing time through the lowered barrier region. |
| \(t_{\rm collapse}\) | Characteristic dressing-leg collapse time. |
| \(\mathcal S(E)=t_{\rm cross}/t_{\rm collapse}\) | Stability ratio used to define the Goldilocks window. |
| \(\mu_\eta\) | Wall-inertia scale for the dressing / geometry leg. |
| \(\chi_\lambda^{\rm peak}\) | Steepest logarithmic gradient of the reduced barrier branch. |
| \(E_{\rm safe,min}\) | Lower edge of the stable over-barrier window in the proton-proxy sweep. |
| \(v_{\rm safe,min}\) | Equivalent lower-edge speed. |
| \(\alpha\) | Ratio \(\mu_\eta/m_s\) in the heavy-throat scaling discussion. |

### A.7 Session-IV damping / heat-sink symbols

| Symbol | Meaning in this thread |
|---|---|
| \(\gamma_{\rm tot}\) | Total V-leg shedding rate. |
| \(\gamma_{\rm vac}\) | Vacuum-side shedding / radiation component. |
| \(\gamma_{\rm lattice}\) | Lattice-side shedding / heat-sink component. |
| \(\gamma_{\rm safe}(v_0)\) | Minimum total shedding that stabilizes a specific cold event. |
| \(\gamma_{\rm crit}\) | Envelope-level unconditional-stability threshold where the characteristic collapse time diverges. |
| \(E_{\rm diss}\) | Total dissipated dressing-leg energy during a crossing. |
| \(E_{\rm vac}\) | Vacuum-side share of dissipated energy. |
| \(E_{\rm lattice}\) | Lattice-heat share of dissipated energy. |
| \(v_0\) | Initial inward speed for the cold-event sweep; default value \(v_0=2.6\). |

### A.8 Session-V condensed-matter mapping symbols

| Symbol | Meaning in this thread |
|---|---|
| \(\lambda_{\rm ep}\) | Dimensionless electron-phonon coupling constant. |
| \(\omega_D\) | Debye frequency of the host material. |
| \(\zeta_{\rm ep}\) | Phenomenological proportionality constant in the damping map \(\gamma_{\rm lattice}^{\rm phys}=\zeta_{\rm ep}\lambda_{\rm ep}\omega_D\). |
| \(t_*\) | Physical time-unit conversion between reduced and laboratory time. |
| \(\lambda_{\rm phys}\) | Physical localization width entering the lattice mapping. |
| \(a_{\rm int}\) | Interstitial spacing used when identifying \(\lambda_{\rm phys}=a_{\rm int}/2\). |
| \(k_{\rm eff}\) | Effective harmonic-trap stiffness in the lattice mapping. |
| \(E_*\) | Physical energy scale used in the force-matching stiffness map. |
| \(T_1\) | Spin-lattice relaxation time. |
| \(\mathcal K_{\rm corr}\) | Korringa constant. |
| \(T_{\max}\) | Maximum temperature below which spin alignment survives the crossing time. |

---

## Appendix B. Reduced-unit and parameter ledger

### B.1 Carry-forward structural constants from the source stack

These constants and structural outputs were already frozen before the barrier session and were treated as carry-forward background rather than as newly refitted session parameters.

| Quantity | Carry-forward value / meaning | Source role |
|---|---|---|
| \(n\) | \(5\) | Frozen stiff-polytropic EOS exponent. |
| \(\kappa_\rho\) | \(1\) | Newtonian mass-dressing coefficient in the gravity-side ledger. |
| \(\kappa_{\rm add}\) | \(1/2\) | Added-mass topology coefficient for the \(w\)-uniform throat. |
| \(\kappa_{\rm PV}\) | \(3/2\) | Adiabatic one-DOF pressure-volume response coefficient. |
| \(\beta_{\rm 1PN}\) | \(3\) | Conservative 1PN precession ledger. |
| \(\alpha^2\) | \(3/4\) | Wake-mixing weight in the vector cross sector. |
| \(a_H\) | \(0\) | Helical admixture set to zero on the carried parity-even branch. |
| \(K_{\rm vec}\) | \(2/\pi^2\) | Vector-sector normalization fixed by the 1PN chain. |
| \(L/a\) | \(\approx 1.85\) | Preferred throat aspect ratio from the EM cavity/throat geometry branch. |
| \(\mu_0^{\rm eff}\) | \(\mu_0/Z_{\rm int}\) | Controlled brane-effective Maxwell coupling. |
| \(q_{\rm eff}\) | \(q_*/\sqrt{Z_{\rm int}}\) | Controlled brane-effective charge. |
| One-port admissibility | \(\Delta>0,\;D_0>0\) | Same bundle conditions used by the barrier audit and normalization chain. |
| Static mixed-kernel families | \(x^{-6},\;e^{-2\kappa x}/x^4,\;e^{-4\kappa x}/x^2\) | Barrier-audit restriction on the static same-charge corridor. |
| First linear dynamic mixed correction | phase lag / pumping, not conservative barrier lowering | Barrier-audit restriction on the linear outgoing channel. |

### B.2 Session-I stationary relaxed-constraint run

#### B.2.1 Parameter set

| Parameter | Value |
|---|---:|
| \(\lambda\) (`lam`) | 1.0 |
| \(r_{\rm reg}\) | 0.25 |
| \(E_{0,\rm amp}\) | 0.18 |
| \(\rho_0\) | 1.0 |
| \(\mu_w\) | 0.8 |
| \(q\) | 1.0 |
| \(a_U\) | 2.5 |
| \(a_V\) | 3.0 |
| \(g_{UV}\) | 0.95 |
| \(s_0\) | 0.9 |
| \(\epsilon_{\rm ref}\) | 0.3 |
| \(K_*\) | 4.0 |
| \(\Omega_U^2\) | 9.0 |
| \(\Omega_W^2\) | 16.0 |
| \(G_U\) | 1.0 |
| \(G_W\) | 1.25 |
| \(R_{\rm mix}\) | 1.35 |
| \(\beta_Q\) | 0.03 |
| \(\beta_{U0}\) | 0.15 |
| \(\beta_{W0}\) | 0.20 |
| \(\kappa\) | 1.0 |
| \(a_0\) | 2.2 |
| \(b_0\) | -0.6 |
| \(r_\sigma\) | 0.8 |
| \(\xi_R\) | 0.9 |
| \(\eta_{\rm leak}\) | 0.03 |
| \(\eta_{UV}\) | 0.22 |
| \(kk_{\rm amp}\) | 0.0 |
| \(r_{\rm F1}\) | 1.77799353547498 |

#### B.2.2 Headline outputs

| Observable | Value |
|---|---:|
| \(\Delta\) | 142.17750000 |
| \(D_0\) | 3.76481862 |
| \(r_{\rm eval}\) | 1.00217028 |
| \(\Xi_1(r_{\rm eval})\) | 0.14313458 |
| \(\widetilde\Xi(r_{\rm eval})\) | 0.14352690 |
| \(U(r_{\rm eval})\) | 0.14313458 |
| \(V(r_{\rm eval})\) | -0.03619791 |
| \(\epsilon_\eta(r_{\rm eval})\) | 0.28933482 |
| \(S_{\rm leak}(r_{\rm eval})\) | 0.03663918 |
| \(J^wE_w(r_{\rm eval})\) | 0.02108684 |
| \(\mathfrak g[\sigma](r_{\rm eval})\) | 0.82823667 |
| \(\mathcal R[\sigma](r_{\rm eval})\) | 0.21677037 |
| \(\sigma_{\min}\) | -0.08979545 |
| sign-changing branch? | True |
| strongest softening point \(r_{\rm soft}\) | 0.18000000 |
| \(V_{\rm eff}(r_{\rm soft})\) | 1.74701126 |
| \(V_{\rm Coul}(r_{\rm soft})\) | 5.55555556 |
| \(V_{\rm eff}/V_{\rm Coul}\) | 0.31446203 |
| UV energy drop at \(r_{\rm soft}\) | 0.21064278 |
| bulk work at \(r_{\rm soft}\) | 1.51632107 |

### B.3 Session-II dynamic scattering / helicity run

#### B.3.1 Additional dynamic parameters

| Parameter | Value |
|---|---:|
| \(m_s\) | 1.0 |
| \(\hbar_{\rm eff}\) | 1.0 |
| \(r_0\) | 5.0 |
| \(r_{\rm contact}\) | 0.18 |
| \(E_{\rm sub}\) | 2.5 |
| `cross_factor` | 1.02 |
| \(\xi_{\rm spin}\) | 0.4 |
| \(C_{0,\rm spin}\) | 0.05 |
| \(r_{\rm spin}\) | 0.8 |
| \(r_{\rm core}\) | 0.15 |
| \(\eta_h\) | 0.25 |
| \(\Omega_{0,\rm sub}\) | 0.7 |
| \(dt\) | 0.0001 |
| \(t_{\max}\) | 10.0 |

#### B.3.2 Headline outputs

| Observable | Value |
|---|---:|
| \(r_{\rm peak}\) | 0.23944389 |
| \(V_{\rm peak}\) | 3.42933112 |
| \(V_{\rm eff}(r_0=5)\) | 0.19999794 |
| Coulomb \(V(r_0)\) | 0.20000000 |
| \(v_{\rm crit,new}\) | 2.54139063 |
| Coulomb contact speed to \(r_{\rm contact}\) | 3.27278339 |
| \(v_{0,\rm sub}\) | 2.14476202 |
| \(r_{\rm turn,new}\) | 0.39096144 |
| \(r_{\rm turn,Coul}\) | 0.40000141 |
| inner turning point | 0.19039548 |
| \(I_{\rm new}\) | 0.19744614 |
| \(I_{\rm Coul}\) | 0.30222297 |
| \(T_{\rm new}\) | \(6.73752615\times 10^{-1}\) |
| \(T_{\rm Coul}\) | \(5.46377065\times 10^{-1}\) |
| \(T_{\rm new}/T_{\rm Coul}\) | 1.23312756 |
| fusion-probability increase | 23.3128% |
| \(\Xi_1(r_{\rm turn})\) | 0.34437471 |
| \(\lambda_{\rm th}(r_{\rm turn})\) | 0.42826825 |
| above-threshold demo \(v_{0,\rm cross}\) | 2.59221845 |
| Coulomb turning point at same \(v_0\) | 0.28091705 |
| aligned \(\max(\partial_t h_{\rm sub})\) | 281.79830789 |
| anti-aligned \(\max(\partial_t h_{\rm sub})\) | 56.96878122 |
| peak helicity-export ratio | 4.94653917 |
| aligned \(h_{\rm sub,final}\) | 20.58070146 |
| anti-aligned \(h_{\rm sub,final}\) | 5.00843357 |
| integrated helicity-export ratio | 4.10920923 |
| scanned \(\lambda_{\rm th}\) range | [0.40673709, 1.06949146] |
| scanned \(\Xi_1\) range | [0.25095422, 0.53817934] |
| scanned WKB-multiplier range | [1.18016972, 1.31627906] |

### B.4 Session-III proton-proxy stability run

#### B.4.1 Stability-sweep additions

| Parameter | Value |
|---|---:|
| \(E_{\rm sub,reference}\) | 2.5 |
| proton mass ratio | 1836.15267343 |
| \(v\)-multiplier min | 1.001 |
| \(v\)-multiplier max | 5.0 |
| number of energy points | 350 |
| number of dynamic samples | 24 |
| \(dt_{\rm dynamic}\) | 0.002 |
| \(t_{\max,\rm dynamic}\) | 500.0 |
| `lambda_choice` | trigger |

#### B.4.2 Headline outputs

| Observable | Value |
|---|---:|
| trigger width \(\lambda_{\rm th}(r_{\rm turn})\) | 0.42826825 |
| chosen \(\lambda_{\rm eff}\) | 0.42826825 |
| proton-proxy \(m_s\) | 1836.15267343 |
| proton-proxy \(\mu_\eta\) | 1836.15267343 |
| proton-proxy threshold speed \(v_{\rm crit,p}\) | 0.05930851 |
| steepest-gradient location | 0.18000000 |
| \(\chi_\lambda^{\rm peak}\) | 21.73204372 |
| \(t_{\rm collapse}\) | 9.43066476 |
| analytic \(E_{\rm safe,min}\) | 5.32265943 |
| analytic \(v_{\rm safe,min}\) | 0.07469791 |
| lower edge in threshold units | \(1.25948037\,v_{\rm crit,p}\) |
| numeric lower edge on scan | 5.36393605 |
| numeric upper edge on scan | 80.93332737 |
| aligned min transit | 0.20400000 |
| aligned max transit | 4.05400000 |
| raw-width \(\chi_\lambda^{\rm peak}\) | 50.74399964 |
| raw-width \(t_{\rm collapse}\) | 6.17163516 |
| raw-width \(E_{\rm safe,min}\) | 27.53273095 |
| reported main safe window | [5.32265943, 80.93332737] |
| equivalent speed window | [0.07469791, 0.29654256] |

### B.5 Session-IV damped V-leg run

#### B.5.1 Damping-specific additions

| Parameter | Value |
|---|---:|
| \(v_{0,\rm cold}\) | 2.6 |
| \(\mu_\eta\) | 1.0 |
| vacuum fraction | 0.25 |
| `gamma_scan_factor_max` | 1.4 |
| number of \(\gamma\) points | 321 |
| `lambda_choice` | model |

#### B.5.2 Headline outputs

| Observable | Value |
|---|---:|
| \(r_{\rm peak}\) | 0.23944389 |
| \(V_{\rm peak}\) | 3.42933112 |
| \(V_{\rm eff}(r_0)\) | 0.19999794 |
| \(v_{\rm crit}\) | 2.54139063 |
| chosen cold speed \(v_0\) | 2.60000000 |
| cold-event energy \(E_{\rm cold}\) | 3.57999794 |
| lowered-branch status | contact |
| pure Coulomb status | turn |
| pure Coulomb turning radius at \(v_0=2.6\) | 0.27933174 |
| actual lowered-branch time to contact | 2.11290000 |
| width used in \(t_{\rm cross}\) | 1.00000000 |
| trigger-width cross-check | 0.42826825 |
| \(\chi_\lambda^{\rm peak}\) | 50.74399964 |
| undamped \(t_{\rm collapse,0}\) | 0.14402764 |
| characteristic \(t_{\rm cross}(v_0=2.6)\) | 1.82169718 |
| undamped stability ratio \(\mathcal S(0)\) | 12.64824697 |
| \(\gamma_{\rm crit}\) | 6.94311167 |
| \(\gamma_{\rm safe}(v_0=2.6)\) | 6.39417302 |
| \(\gamma_{\rm vac}/\gamma_{\rm tot}\) | 0.25 |
| \(\gamma_{\rm lattice}/\gamma_{\rm tot}\) | 0.75 |
| \(\gamma_{\rm vac}\) at \(\gamma_{\rm safe}\) | 1.59854325 |
| \(\gamma_{\rm lattice}\) at \(\gamma_{\rm safe}\) | 4.79562976 |
| \(E_{\rm diss}\) at \(\gamma_{\rm safe}\) | 0.01033460 |
| \(E_{\rm vac}\) at \(\gamma_{\rm safe}\) | 0.00258365 |
| \(E_{\rm lattice}\) at \(\gamma_{\rm safe}\) | 0.00775095 |
| \(E_{\rm store,final}\) at \(\gamma_{\rm safe}\) | 0.00960157 |
| \(\gamma_{\rm vac}\) at \(\gamma_{\rm crit}\) | 1.73577792 |
| \(\gamma_{\rm lattice}\) at \(\gamma_{\rm crit}\) | 5.20733375 |
| \(E_{\rm lattice}\) at \(\gamma_{\rm crit}\) | 0.00770830 |

### B.6 Session-V material-mapping inputs and outputs

#### B.6.1 Inputs carried from earlier reduced runs

| Quantity | Value |
|---|---:|
| required \(\gamma_{\rm lattice}\) | 4.79562976 |
| reduced \(t_{\rm cross}\) | 1.82169718 |
| reduced \(\lambda_{\rm th}\) | 0.42826825 |
| reduced \(r_{\rm turn}\) | 0.39096144 |
| reduced \(V_{\rm eff}(r_{\rm turn})\) | 2.50000000 |

#### B.6.2 Headline mapping formulas and benchmark numbers

| Mapping output | Result |
|---|---|
| \((\lambda_{\rm ep}\omega_D)_{\min}\) | \(4.79562976/(t_*\zeta_{\rm ep})\) |
| equivalent threshold | \(8.73618521011608/(t_{\rm cross}^{\rm phys}\zeta_{\rm ep})\) |
| turnover product | \(\lambda_{\rm ep}\omega_D\,t_{\rm cross}^{\rm phys}\ge 8.73618521011608/\zeta_{\rm ep}\) |
| harmonic-trap geometry ratio | \(\chi_{\lambda,\rm lattice}=2\lambda_{\rm phys}/r_{\rm phys}\) |
| threshold geometry condition | \(r_{\rm phys}\le 2\lambda_{\rm phys}\) |
| \(r_{\rm turn}^{\rm phys}/\lambda_{\rm phys}\) | 0.912889153001653 |
| \(\chi_{\lambda,\rm lattice}(r_{\rm turn})\) | 2.19084649 |
| force-matched stiffness | \(k_{\rm eff,req}=2.7385581171381\,E_*[\mathrm{eV}]/\lambda_{\rm phys}^2[\mathrm{\AA}^2]\) |
| interstitial-spacing form | \(k_{\rm eff,req}=10.9542324685524\,E_*[\mathrm{eV}]/a_{\rm int}^2[\mathrm{\AA}^2]\) |
| \(T_{\max}\) | \(\mathcal K_{\rm corr}/t_{\rm cross}^{\rm phys}\) |
| reduced-unit form of \(T_{\max}\) | \(0.548938655106224\,\mathcal K_{\rm corr}/t_*\) |

### B.7 Reduced-unit conventions that stayed unresolved

The session intentionally left several physical conversion factors unresolved. They are recorded here because many later interpretation mistakes come from forgetting them.

| Quantity | Status in this session |
|---|---|
| \(t_*\) | Unfixed reduced-to-physical time conversion. |
| \(\lambda_{\rm phys}\) | Unfixed physical interpretation of the reduced width \(\lambda\). |
| \(E_*\) | Unfixed physical energy scale used in the force-matching stiffness map. |
| \(\zeta_{\rm ep}\) | Unfixed proportionality between phenomenological lattice damping and \(\lambda_{\rm ep}\omega_D\). |
| \(\mathcal K_{\rm corr}\) | Material-specific Korringa constant, not inserted numerically in the session. |
| candidate-host elastic / electronic data | Not yet screened against the derived thresholds. |

---

## Appendix C. Artifact index

### C.1 Session write-up files

| File | Role |
|---|---|
| `session_barrier_thread_toc.md` | Session table of contents used to organize the thread. |
| `session_barrier_sections_0_2.md` | Draft write-up of Sections 0–2. |
| `session_barrier_sections_3_5.md` | Draft write-up of Sections 3–5. |
| `session_barrier_sections_6_8.md` | Draft write-up of Sections 6–8. |
| `session_barrier_writeup_0_8_clean.md` | Clean merged main write-up without appendix. |

### C.2 Core source-stack files used in the write-up

| File | Role in the barrier thread |
|---|---|
| `4d_summary.md` | Parent \(4+1\) action, projection vs reduction, leakage identity, charge ontology, brane Maxwell hook. |
| `4d_em_fields_summary.md` | Localized Maxwell sector, \(q_{\rm eff}\), \(\mu_0^{\rm eff}\), KK/Yukawa corrections, current-conservation structure. |
| `4d_plasma_summary.md` | Beyond-MHD mixed channels, \(E_w\), \(C_a\), \(J^w\), projected leakage, open-system interpretation. |
| `moving_throat_pde_program_compact.md` | Claim-status firewall, grouped real \(P_2\) program status, open branch-realization bottleneck. |
| `barrier_audit_full.md` | Same-charge static/dynamic mixed-kernel restrictions and corridor narrowing. |

### C.3 Session-I computational artifacts

| File | Role |
|---|---|
| `stress_test_relaxed_constraints_sympy.py` | SymPy stress test for relaxed stationary constraints. |
| `stress_test_relaxed_constraints_report.txt` | Stationary run report with formulas, parameters, and outputs. |
| `relaxed_constraints_veff.png` | Plot of the reduced \(V_{\rm eff}(r)\) against Coulomb. |
| `relaxed_constraints_diagnostics.png` | Leakage / \(\Xi_1\) / mouth-response diagnostics for Session I. |

### C.4 Session-II computational artifacts

| File | Role |
|---|---|
| `dynamic_scattering_helicity_lambda_sympy.py` | Dynamic scattering and helicity continuation. |
| `dynamic_scattering_helicity_lambda_report.txt` | Turning-point, WKB, helicity-export, and \(\lambda_{\rm th}\) report. |
| `dynamic_scattering_potential.png` | Potential plot for the dynamic continuation. |
| `dynamic_scattering_trajectories.png` | Trajectory comparison for lowered branch vs Coulomb. |
| `dynamic_helicity_export.png` | Helicity-export comparison for aligned vs anti-aligned closures. |
| `lambda_threshold_curve.png` | \(\lambda_{\rm th}\) / turning-point trigger curve. |

### C.5 Session-III computational artifacts

| File | Role |
|---|---|
| `proton_proxy_stability_sweep_sympy.py` | Proton-proxy Goldilocks-window sweep. |
| `proton_proxy_stability_report.txt` | Stability-ratio report and safe-window summary. |
| `proton_proxy_stability_timescales.png` | Crossing-vs-collapse timescale plot. |

### C.6 Session-IV computational artifacts

| File | Role |
|---|---|
| `damped_shedding_cold_sweep_sympy.py` | Damped V-leg / cold-sweep script. |
| `damped_shedding_report.txt` | Damped-shedding report with \(\gamma_{\rm safe}\), \(\gamma_{\rm crit}\), and heat partition. |
| `damped_shedding_timescales.png` | Timescale comparison under damping. |
| `damped_shedding_heat_partition.png` | Vacuum-vs-lattice heat partition plot. |

### C.7 Session-V computational artifacts

| File | Role |
|---|---|
| `material_mapping_condensed_matter_sympy.py` | Condensed-matter parameter map script. |
| `material_mapping_condensed_matter_report.txt` | Threshold equations for \(\lambda_{\rm ep}\omega_D\), \(k_{\rm eff}\), and \(T_{\max}\). |

### C.8 Workspace note on expected-but-missing plotting artifacts

The current workspace snapshot did **not** contain the three material-mapping plot files that were provisionally named in the early table of contents (`material_mapping_lambdaep_vs_Debye.png`, `material_mapping_keff_vs_lambda.png`, `material_mapping_Tmax_vs_Korringa.png`). The session thread therefore has script-plus-report coverage for Session V, but not companion plot files under those exact names in the present directory snapshot.

---

## Appendix D. Source-document map

### D.1 Parent-theory and reduction backbone

| Source file | What it contributed to this session |
|---|---|
| `4d_summary.md` | Provided the exact parent GNLS + localized Maxwell + geometry action, the projection/reduction distinction, the exact projected leakage source, the exact longitudinal identity, and the corrected charge ontology used throughout the session. |
| `4d_em_fields_summary.md` | Supplied the localized \(4+1\) Maxwell equations, the controlled brane reduction to \(\mu_0^{\rm eff}=\mu_0/Z_{\rm int}\), the canonical \(q_{\rm eff}\) formula, and the KK/Yukawa packaging that underlies the Coulomb-plus-short-range view of the EM sector. |
| `4d_plasma_summary.md` | Made the beyond-MHD mixed sector explicit: \(A_w\), \(F_{\mu w}\), \(E_w\), \(C_a\), \(J^w\), projected leakage, and the interpretation of brane non-ideality as conservative export into hidden or higher-mode channels. |

### D.2 Moving-throat and barrier-audit backbone

| Source file | What it contributed to this session |
|---|---|
| `moving_throat_pde_program_compact.md` | Supplied the session’s status firewall: what counts as exact, exact within closure, reduced/controlled, effective closure, or open. Also stated explicitly that the main remaining theorem gap is **branch realization** on the completed moving-throat PDE. |
| `barrier_audit_full.md` | Narrowed the same-charge corridor before the session began: the static one-port mixed bundle only renormalizes short-range families, and the first linear outgoing correction is phase lag / pumping rather than conservative long-range barrier lowering. |

### D.3 Session-generated reports

| Source file | What it contributed to this session |
|---|---|
| `stress_test_relaxed_constraints_report.txt` | Recorded the stationary relaxed-constraint closure, the explicit \(S_{\rm leak}\) and \(J^wE_w\) channels, the \(U/V\) drain, the sign-changing source branch, and the \(V_{\rm eff}/V_{\rm Coul}\approx 0.314\) near-contact softening. |
| `dynamic_scattering_helicity_lambda_report.txt` | Recorded the dynamic continuation: barrier peak, threshold speed, turning points, WKB enhancement, aligned-vs-anti-aligned helicity export, and the \(\lambda_{\rm th}\) trigger curve. |
| `proton_proxy_stability_report.txt` | Recorded the Goldilocks timescale sweep under proton-proxy mass/inertia scaling, including \(t_{\rm cross}\), \(t_{\rm collapse}\), and the one-sided safe window. |
| `damped_shedding_report.txt` | Recorded the damped V-leg closure, \(\gamma_{\rm safe}\), \(\gamma_{\rm crit}\), and the split of dissipated energy into vacuum and lattice channels. |
| `material_mapping_condensed_matter_report.txt` | Recorded the condensed-matter thresholds for lattice turnover, harmonic-trap stiffness, and Korringa-limited temperature ceiling. |

### D.4 How the source map is supposed to be read

The session used the core source stack to constrain **what kind of branch was even allowed** and then used the session reports to explore one reduced barrier/stability/materials branch numerically. The source stack did not “prove” the final session outputs directly. It defined the allowed ontology, the admissible closures, and the already-existing restrictions that the session outputs were expected to respect.

---

## Appendix E. Non-claims and open theorem gaps

This appendix is deliberately blunt. It records what the session **did not** show, because many later misunderstandings come from reading a reduced closure as if it were a theorem of the full moving-throat PDE.

### E.1 The session did not solve the full moving-throat PDE

The compact moving-throat program still treats the actual branch realization problem as open. The session’s barrier/stability/materials chain therefore lives inside a hierarchy of reduced closures and phenomenological completions. It is strongest as a **reduced / controlled reduction** plus **effective closure** story, not as a completed theorem of the fully solved bulk/interface PDE.

### E.2 The session did not discover a new long-range same-charge attraction

The barrier audit already constrained the static same-charge mixed bundle to short-range families. Nothing in Sessions I–V overturned that. The stationary reduction lowered the near-contact barrier, but it did so through short-range families, leakage/export channels, dressing-leg drain, and source compensation. The dynamic audit likewise did not turn the first linear outgoing correction into a new conservative long-range force. The live corridor remained a **short-range/open-system bypass**, not a replacement of Coulomb by a hidden attractive \(1/r\)-law.

### E.3 The stationary relaxed-constraint branch is still closure-dependent

The Session-I result depends on a specific phenomenological choice for the \(U/V\) free energy, a specific scaling of the leakage channel, and a specific sign-changing source completion. The fact that the reduced \(V_{\rm eff}(r)\) dropped to about \(0.314\) of Coulomb in the plotted near-contact window is therefore a result **inside that declared closure**, not a branch-independent theorem.

### E.4 The helicity-export preference is not yet a spin theorem

The aligned-vs-anti-aligned comparison in Session II established a strong preference for aligned spins on the chosen reduced mixed-sector closure, because aligned spins exported far more subscale helicity. That is a useful reduced diagnostic, but it is not yet a full theorem of the microscopic spin-support structure of the parent theory.

### E.5 The collapse and damping times are characteristic closures, not solved geometry times

The \(t_{\rm cross}\), \(t_{\rm collapse}\), and \(t_{\rm collapse}^{\rm damped}\) formulas were characteristic estimates chosen to stress-test the branch. In particular, the finite \(\gamma_{\rm crit}\) of Session IV comes from the adopted **envelope** closure. A raw damped inverted-oscillator equation would not, by itself, generate the same finite unconditional-stability threshold. So \(\gamma_{\rm crit}\) should be read as a property of the adopted damped-envelope reduction, not as a closure-free theorem.

### E.6 The proton-proxy window did not prove generic heavy-particle stability

Session III showed that a one-sided stable window opens when the heavy-throat scaling is combined with the trigger-width choice inherited from Session II. It did **not** show that proton-like defects are generically stable simply because they are heavy. The safe lower edge remained highly sensitive to the confinement-width / \(\chi_\lambda\) choice.

### E.7 The condensed-matter map did not yet prove any specific candidate material

Session V stopped at threshold equations. It did not insert real Palladium, Titanium, hydride, deuteride, alloy, or stressed-lattice data and did not prove that any particular host already satisfies all three thresholds simultaneously. The material screen is still a next step, not a completed result.

### E.8 The session did not prove a room-temperature reactor theorem

A few higher-level interpretations suggested that very fast crossing could make spin scrambling irrelevant in ordinary laboratory conditions. The session itself did **not** prove that. It produced the symbolic ceiling
\[
T_{\max}=\frac{\mathcal K_{\rm corr}}{t_{\rm cross}^{\rm phys}},
\]
but it left both \(t_*\) and the material-specific Korringa constant unfixed. So any room-temperature or hot-operation claim remains interpretive until the physical unit map and material constants are inserted.

### E.9 The stiffness map contains an extra force-matching assumption

The \(\chi_\lambda\) criterion alone only fixes the geometry ratio \(r_{\rm phys}/\lambda_{\rm phys}\) for a harmonic trap. The explicit \(k_{\rm eff}\) formula required one more assumption: matching the absolute harmonic-trap force to the reduced-model barrier force at the turning point. That formula is therefore useful, but it is not a pure consequence of \(\chi_\lambda\) by itself.

### E.10 The most important remaining theorem gaps

The remaining gaps are now fairly narrow.

1. **Physical-unit calibration:** determine \(t_*\), \(\lambda_{\rm phys}\), \(E_*\), and \(\zeta_{\rm ep}\).
2. **Candidate-material screening:** test real materials against the derived lattice-turnover, stiffness, and Korringa thresholds.
3. **True branch realization:** determine whether the completed moving-throat PDE actually returns the isotropic / grouped-\(P_2\) / mixed-sector branch that the reduced barrier scripts assumed.
4. **Placement/orbit-lock kill test:** determine whether the actual branch survives the static placement / orbit-lock constraints highlighted by the barrier audit.
5. **Higher-fidelity damping/export closure:** replace the phenomenological V-leg envelope law with a more microscopic wall/support/mixed-sector export calculation.

### E.11 The cleanest falsifier after this appendix

The narrowest next falsifier is probably not “build a full reactor model.” It is:

1. calibrate reduced units tightly enough to turn the Session-V threshold equations into real laboratory inequalities,
2. test whether any plausible host material satisfies them,
3. and, only if that screen survives, return to the moving-throat PDE to see whether the actual branch realization preserves the same corridor once the hidden channels are no longer packaged phenomenologically.

That order keeps the project faithful to its established methodology: derive the tightest kill test first, then go deeper only if the target remains alive.
