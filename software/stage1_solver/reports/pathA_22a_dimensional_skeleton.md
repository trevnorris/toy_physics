# PathA 22a dimensional skeleton

## Result

- Homogeneity: `HOMOGENEITY_PASS`.
- P0 finding: `P0_DIMENSIONLESS_AFTER_EXPLICIT_FREQUENCY_NORMALIZATION`. The actual patha_extraction.py mixed-port formulas derive N0=P^2/Delta^2 with dimension M/L, so faithful N0/D0 has dimension T^2. Multiplying by the explicit (c_s/a)^2 frequency-normalization factor makes P0 dimensionless. With that normalized P0, R_norm is homogeneous; without it, the gate would fail.
- Headline: `TUNABILITY_CHANNEL_PRESENT`.
- Tunability channels: >= 2 (`chi_Q / S_port`, `lambda_gamma=c/c_s`).
- Strict class-(d) free knobs: `chi_Q / S_port`.
- Strict-ledger disclosure: The strict class labels under-count tunability: lambda_gamma and the source-map residual cluster are underived class-(c) quantities in a proof ledger, but become calibration knobs under a calibrate-predict methodology.
- Scoped next step: Before deriving broad G/mass/source machinery, derive the minimal combination xi=mhat0^2*S_port/[G*c_s^5/(a^5*c^5)], with chi_Q/S_port fixed only by actual outgoing compact DtN matching and lambda_gamma fixed only by brane zero-mode reduction; otherwise both remain calibration channels.

## Dimension table

| ingredient | dimension | SymPy monomial | derivation | provenance |
| --- | --- | --- | --- | --- |
| `G_3D` | `L^3 T^-2 M^-1` | `L**3/(M*T**2)` | Newton force dimension in observed 3-space: [G]=L^3 T^-2 M^-1. | software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:31; software/stage1_solver/reports/pathA_19_dimensional_foundation.md |
| `c_s` | `L T^-1` | `L/T` | EOS sound-speed law c_s^2=5 K rho0^4/m_GNLS. | software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:59; dimensional_check.py pathA_19_F3 |
| `c` | `L T^-1` | `L/T` | Observed brane light speed. | pathA_20/20b carried as lambda_gamma=c/c_s. |
| `a` | `L` | `L` | Throat/core length scale. | research/pde/paper/pde.tex:2061 |
| `mhat0` | `L^-1 T^-1 M^-1/2` | `1/(L*sqrt(M)*T)` | Source-map factor with [mhat0]=L^-1 T^-1 M^-1/2, so mhat0^2 matches the full target dimension. | software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:35 |
| `K,B0,Z0,D0` | `L^-1 T^-2 M` | `M/(L*T**2)` | Reduced static denominator compiler: D0=K-B0-Z0; Z0 also derives to reduced stiffness. | software/stage1_solver/src/stage1_solver/patha_extraction.py:397-445,473-482; research/pde/paper/pde.tex:1849-1872 |
| `N0,N2,N4 from code formulas` | `L^-1 M` | `M/L` | N0=P^2/Delta^2 derives to K*T^2=M/L; N2 and N4 derive to M/L*T^2 and M/L*T^4. | software/stage1_solver/src/stage1_solver/patha_extraction.py:397-445 |
| `P0_faithful=N0/D0` | `T^2` | `T**2` | Faithful code-formula ratio before normalization; dimension T^2. | research/pde/paper/pde.tex:2018-2026; software/stage1_solver/src/stage1_solver/patha_extraction.py:482 |
| `P0=(c_s/a)^2*N0/D0` | `1` | `1` | Dimensionless static outgoing prefactor after explicit frequency normalization. | software/stage1_solver/src/stage1_solver/patha_extraction.py:397-445,473-482; research/pde/paper/pde.tex:2018-2026 |
| `S_port == chi_Q` | `1` | `1` | Dimensionless outgoing-normalization scalar; S_port is the code slot in observable_residuals/extract_branch. | software/stage1_solver/src/stage1_solver/patha_extraction.py:526-544,609-624; research/pde/paper/pde.tex:1980-1998 |
| `Gamma5` | `T^5` | `T**5` | Gamma5=chi_Q P0 a^5/(27 c_s^5); the radiation time^5 factor is explicit, not hidden in P0. | research/pde/paper/pde.tex:2053-2061 |
| `54*G*c_s^5/(5*a^5*c^5)` | `L^-2 T^-2 M^-1` | `1/(L**2*M*T**2)` | Target dimension is M^-1 L^-2 T^-2. | research/pde/paper/pde.tex:2082 |

## Homogeneity checks

- `D0=K-B0-Z0 reduced static denominator`: **CONSISTENT** (expected `L^-1 T^-2 M`, actual `L^-1 T^-2 M`). Uses the reduced coefficient compiler contract in patha_extraction.py and pde.tex. terms=K:L^-1 T^-2 M, B0:L^-1 T^-2 M, Z0:L^-1 T^-2 M
- `mixed-port Delta=Omega_U^2*Omega_W^2-R^2`: **CONSISTENT** (expected `T^-4`, actual `T^-4`). Uses pde.tex building-block dimensions Omega_U,Omega_W:T^-1 and R:T^-2. terms=Omega_U^2*Omega_W^2:T^-4, R^2:T^-4
- `mixed-port P=Omega_U^2*g_W+R*g_U`: **CONSISTENT** (expected `L^-1/2 T^-4 M^1/2`, actual `L^-1/2 T^-4 M^1/2`). Uses [g_U]^2=[g_W]^2=[K]*T^-2 from the same building-block dictionary. terms=Omega_U^2*g_W:L^-1/2 T^-4 M^1/2, R*g_U:L^-1/2 T^-4 M^1/2
- `Z0=Q/Delta derived from mixed-port formula`: **CONSISTENT** (expected `L^-1 T^-2 M`, actual `L^-1 T^-2 M`). Cross-check: Z0 remains a reduced stiffness coefficient.
- `Z2 derived from mixed-port formula`: **CONSISTENT** (expected `L^-1 M`, actual `L^-1 M`). Cross-check: Z2=K*T^2 matches the existing reduced_M_B2_Z2_N2 dictionary entry.
- `Z4 derived from mixed-port formula`: **CONSISTENT** (expected `L^-1 T^2 M`, actual `L^-1 T^2 M`). Cross-check: Z4=K*T^4 matches the existing reduced_B4_Z4_N4 dictionary entry.
- `N0=P^2/Delta^2 derived from code formula`: **CONSISTENT** (expected `L^-1 M`, actual `L^-1 M`). This is derived from patha_extraction.py, not asserted from the old reduced_K_B0_Z0_N0 dictionary bucket.
- `N2 derived from code formula`: **CONSISTENT** (expected `L^-1 T^2 M`, actual `L^-1 T^2 M`). The N tower is shifted by T^2 relative to the old asserted N2 bucket.
- `N4 derived from code formula`: **CONSISTENT** (expected `L^-1 T^4 M`, actual `L^-1 T^4 M`). The N tower is shifted by T^2 relative to the old asserted N4 bucket.
- `faithful P0=N0/D0 before frequency normalization`: **CONSISTENT** (expected `T^2`, actual `T^2`). The actual code formula gives P0 with T^2 before applying the explicit frequency-normalization factor.
- `P0 frequency-normalization factor (c_s/a)^2`: **CONSISTENT** (expected `T^-2`, actual `T^-2`). The paper's dimensionless P0 requires this explicit T^-2 factor; no a^5/c_s^5 radiation factor is hidden here.
- `normalized P0=(c_s/a)^2*N0/D0 dimension`: **CONSISTENT** (expected `1`, actual `1`). This is the P0 used by the R_norm homogeneity gate and by Gamma5.
- `mhat0 source-map dimension`: **CONSISTENT** (expected `L^-1 T^-1 M^-1/2`, actual `L^-1 T^-1 M^-1/2`). Source-map dimensional table carried by pathA_21b; this group rechecks the monomial.
- `3D GR target 54*G*c_s^5/(5*a^5*c^5)`: **CONSISTENT** (expected `L^-2 T^-2 M^-1`, actual `L^-2 T^-2 M^-1`). Uses 3D Newton G=L^3 T^-2 M^-1; numerical factors are dimensionless.
- `S_port / chi_Q outgoing-normalization scalar`: **CONSISTENT** (expected `1`, actual `1`). S_port is the code slot occupied by the paper's chi_Q.
- `R_norm=mhat0^2*S_port*P0 - 54*G*c_s^5/(5*a^5*c^5)`: **CONSISTENT** (expected `L^-2 T^-2 M^-1`, actual `L^-2 T^-2 M^-1`). This is the actual homogeneity gate; it is not an x==x self-check. terms=mhat0^2*S_port*P0:L^-2 T^-2 M^-1, GR_target:L^-2 T^-2 M^-1

## P0 normalization reading

- Faithful `N0/D0` dimension before normalization: `T^2`.
- Required normalization: `(c_s/a)^2` with dimension `T^-2`.
- Normalized `P0` dimension: `1`.
- The outgoing radiation carrier is explicit in `Gamma5=chi_Q*P0*a^5/(27*c_s^5)`; it is not hidden inside `P0`.

## Dimensionless audit

- `G = (a*c_s^2/m_GNLS) * g_G = (5*a*K*rho0^4/m_GNLS^2) * g_G`
- `mhat0 = (c_s/(a^2*sqrt(m_GNLS))) * g_mhat`
- `c = lambda_gamma*c_s`
- `xi*S_port*P0 = P0*chi_Q*g_mhat**2*lambda_gamma**5/g_G`.
- All g_* factors are dimensionless and underived unless separately fixed.

## Knob ledger

| knob/residual | class | role | provenance | verdict effect |
| --- | --- | --- | --- | --- |
| `sigma_Q^can=4*a^5/(27*c_s^5)` | (a) fixed-by-prior-derivation | Canonical compact outgoing fingerprint factor; dimensional T^5 carrier for Gamma5. | research/pde/paper/pde.tex:1980-1993,2061 | Not tunable once the canonical compact outgoing branch is selected. |
| `grouped signature (1,1/2,-1)` | (a) fixed-by-prior-derivation | Weak-axisymmetric grouped transport signature. | research/pde/paper/pde.tex:92,2371-2384 | Does not affect the isotropic homogeneity gate. |
| `lambda_gamma=c/c_s` | (c)/(d) underived residual; TRUE FREE calibration knob under calibrate-predict | Observed brane cone versus sound cone ratio; enters xi as lambda_gamma^5 and is equally unpinned by the current reduction. | software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md; research/pde/paper/pde.tex:2551 | A future brane zero-mode reduction may pin it; under calibrate-predict it is a tunability channel, not a prediction. |
| `chi_Q / S_port` | (d) TRUE FREE calibration knob | Outgoing-normalization scalar. Code S_port multiplies mhat0^2*P0 in exactly the paper's chi_Q slot. | software/stage1_solver/src/stage1_solver/patha_extraction.py:526-544,609-624; research/pde/paper/pde.tex:1980-1998,2071-2082 | TUNABILITY_CHANNEL_PRESENT unless the actual outgoing compact DtN derivation independently fixes chi_Q; Stage 104/105 canonical sigma_Q choice does not prove the Path-A branch is canonical. |
| `P0 branch value=N0/D0` | (b) branch-determined (target-blind) | Static outgoing prefactor from reduced overlap data after the branch solve. | software/stage1_solver/src/stage1_solver/patha_extraction.py:397-445,473-482; research/pde/paper/pde.tex:1849-1872,2018-2026 | Not a calibration knob if the profile/overlaps are solved target-blind. |
| `R0/a` | (b) branch-determined (target-blind) | Finite throat radius profile scale in the branch data. | software/stage1_solver/src/stage1_solver/patha_extraction.py:25-37,64-75; research/pde/paper/pde.tex:2551 | Conditional on solving the stationary branch; not automatically tunable. |
| `L/a` | (b) branch-determined (target-blind) | Wall/worldtube length ratio used by finite-throat profiles. | software/stage1_solver/src/stage1_solver/patha_extraction.py:40-53,691-700 | Geometry/domain data; only tunable if later allowed as calibration. |
| `W_eff/a` | (c) underived residual | Effective kernel/support width entering the profile and force normalization. | software/stage1_solver/decisions/13_emergent_constants_derivation.md:433-445 | Pure-number closure needs the branch/kernel form. |
| `Theta_Q` | (c) underived residual | Quadrupole source/shape residual not closed by dimensional analysis. | software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:72-76 | Carried to the minimal xi/source-map derivation. |
| `J-selector / flux law` | (c) underived residual | Mass/source bridge rate selector inherited from the pathA_20/pathA_21 chain. | software/stage1_solver/decisions/13_emergent_constants_derivation.md:407,433-445 | Affects mhat/G forms, not dimensional homogeneity. |
| `alpha_J` | (c) underived residual | Dimensionless mass-bridge coefficient, including the h versus hbar/2pi convention. | software/stage1_solver/decisions/13_emergent_constants_derivation.md:407; software/stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md | Pending source/boundary/Hamiltonian derivation; not set to one. |
| `branch-kernel choices` | (c) underived residual | Profile/kernel choices that determine overlap integrals and the branch packet. | software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:72-76 | Indeterminate until forms are derived or branch solve fixes them target-blind. |
| `g_G` | (c) underived residual | Arbitrary dimensionless multiplier in G=(a*c_s^2/m_GNLS)*g_G. | software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:63-71 | Must remain symbolic in xi; dimensional analysis cannot set it to one. |
| `g_mhat` | (c) underived residual | Arbitrary dimensionless multiplier in mhat0=(c_s/(a^2*sqrt(m_GNLS)))*g_mhat. | software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:63-71 | Must remain symbolic in xi; source-map form is the scoped next step. |

## Calibrate-Predict Disclosure

- Tunability lower bound: `>= 2` channels: `chi_Q / S_port`, `lambda_gamma=c/c_s`.
- `chi_Q/S_port` is expected to be pinned only by an actual outgoing compact DtN derivation.
- `lambda_gamma=c_gamma/c_s` is expected to be pinned only by the brane zero-mode reduction and enters the reduction as `lambda_gamma^5`.
- Both are underived branch/reduction quantities; both become free knobs under calibrate-predict.
- Strict class-(c) source-map residuals that also become tunable under calibrate-predict: `Theta_Q`, `alpha_J`, `W_eff/a`, `branch-kernel choices`.

## Known Inconsistency

- `KNOWN_PDE_TEX_INCONSISTENCY`: eq:outgoing-BT-target forces dimensionful mhat with [mhat]=L^-1 T^-1 M^-1/2, while eq:outgoing-natural-source-map writes mhat=1+O(a^2/r^2), which is dimensionless.
- Code/harness reading: `[mhat]=L^-1 T^-1 M^-1/2`.
- Required resolution: The eventual source-map derivation must reconcile this; pathA_22a uses the dimensionful reading.

## Negative controls

- Planted dimensional mismatch: `CAUGHT_DIMENSIONAL_MISMATCH`; expected `HOMOGENEITY_FAILURE`, actual `HOMOGENEITY_FAILURE`.
- Planted unresolved dimensionless factor: `PRESERVED_UNRESOLVED_FACTOR`; expected `TUNABILITY_CHANNEL_PRESENT or INDETERMINATE_NEEDS_FORMS`, actual `TUNABILITY_CHANNEL_PRESENT`.
- Reachability: homogeneity failure -> `HOMOGENEITY_FAILURE`; tunability -> `TUNABILITY_CHANNEL_PRESENT`; no class-(d) but residual forms -> `INDETERMINATE_NEEDS_FORMS`.
- Wrong `N0` stiffness assertion: `CAUGHT_WRONG_N0_ASSERTION`.

## Residuals

- chi_Q/S_port is not canonically fixed by the current code path; S_port defaults to 1.0 as a convention.
- lambda_gamma enters as lambda_gamma^5 and remains unpinned after pathA_20b; it is a second tunability channel under calibrate-predict.
- Theta_Q, alpha_J, W_eff/a, and branch-kernel source-map residuals are strict class-(c) gaps but become tunable under calibrate-predict.
- g_G, g_mhat, alpha_J, J-selector, W_eff/a, Theta_Q, and branch-kernel forms remain underived.
- pde.tex is internally inconsistent on mhat; the code uses the dimensionful eq:outgoing-BT-target reading.
- Faithful N0/D0 has dimension T^2; P0 is dimensionless only after the explicit (c_s/a)^2 frequency normalization.

