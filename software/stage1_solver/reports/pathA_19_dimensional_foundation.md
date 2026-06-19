# Path-A 19 dimensional foundation reference

## Verdict

- **Base set:** `RETAIN_L_T_M`.
- **Reason:** m_GNLS is an explicit action parameter and no action-level derivation ties m_defect to inflow. The LT projection is a natural-unit representation only.
- **Mass split:** `m_GNLS` is the constituent mass in the GNLS action; `m_defect` is a throat branch mass and is not derived from inflow in this step.
- **No normalization change:** this report does not derive or modify any frozen normalization factor.

## F1 mass fork

| Quantity | Status | Dimension/result | Source |
|---|---|---|---|
| `m_GNLS` | explicit action parameter | `M` retained | `part01_parent_geometry.tex:174-219`; `pde.tex:326-406` |
| `m_defect` | not action-derived here | blocked by `INFLOW_MASS_SOURCE_MISSING` | `brane_bulk_ontology.tex:1267-1297`, `1998-2039` |
| `hbar*J/c_gamma^2` | dimensional conversion only | `M` | harness check `conditional defect rest-frequency conversion` |

The current action contains `m_GNLS` in the kinetic operator, current, Madelung velocity, Euler equation, and sound-speed law. The ontology supplies drainage/volume-deficit scaling for defect mass, but not a boundary source, Noether charge, or Hamiltonian energy theorem tying `m_defect` to the inflow rate. Therefore the mass-emergence hypothesis is rejected for this foundation gate and `{L,T,M}` is retained.

## F2 frame-tagged flux and a pin

| Flux | Frame | Dimension | Interpretation |
|---|---|---|---|
| `J_bulk = int rho v.dSigma_3` | 4D bulk closed 3-surface | `T^-1` | number flux/rate |
| `J_brane = int rho_3 v.dS_2` | 3D brane reduction | `T^-1` | number flux/rate after transverse reduction |
| `Q_vol,bulk = rho^-1 J` | 4D bulk | `L^4 T^-1` | bulk four-volume flux |
| `Q_vol,brane = rho_3^-1 J` | 3D brane | `L^3 T^-1` | brane volume flux |
| `m_GNLS J` | action constituent mass flux | `M T^-1` | mass-per-particle times number rate |

Gauss shape-independence holds only in a region with no enclosed source or leakage. Projection creates `S_leak` (`part01_parent_geometry.tex:298-330`; `pde.tex:512-539`), and the throat bottom is explicitly open/closed/connected pending microscopic input (`brane_bulk_ontology.tex:1998-2039`). No-net-accretion is therefore carried as `NO_NET_ACCRETION_BC_UNDERIVED`.

`a` is a mouth-radius collective moment, not a fundamental coordinate: `a0=R0(0)` and `a(t)` is the mouth average (`part01_parent_geometry.tex:447-510`; `pde.tex:633-648`). The conserved rate `J` is the better invariant label; `a` remains branch geometry consumed by downstream scale maps.

## F3 pins, healing length, and dictionary

- `4` pins on `3` base dimensions leave one relation: `a*c_s0*m_GNLS/hbar = 1`, i.e. `a = hbar/(m_GNLS*c_s0)`.
- GNLS core balance gives `xi_h = sqrt(hbar^2/(2*m_GNLS*h0)); h0=(5K/4)*rho0^4=(m_GNLS*c_s0^2)/4`.
- **Derived healing scale:** `xi_h = sqrt(2)*hbar/(m_GNLS*c_s0)`.
- Consequence: If a is identified with the GNLS healing core, then a/xi_h is a convention/branch factor; the raw four pins correspond to a=xi_h/sqrt(2).

| Quantity | Current status |
|---|---|
| `hbar` | independent action constant, dimension `M L^2 T^-1` |
| `m_GNLS` | independent action constant, dimension `M` |
| `K` | independent EOS constant after choosing the EOS, dimension `M L^18 T^-2` |
| chosen state `rho0` | independent branch/state datum, dimension `L^-4` in 4D bulk |
| `c_s0` | derived from `c_s0^2=5K rho0^4/m_GNLS` |
| `xi_h` | derived from GNLS kinetic/enthalpy balance |
| `a` | derived/branch collective geometry if identified with a core scale; otherwise an input branch moment |
| `m_defect` | not emergent in this step; blocked by `INFLOW_MASS_SOURCE_MISSING` |

Dictionary confirmation: the harness 4D action dictionary is homogeneous for the GNLS, gauge coupling, Maxwell sector with explicit `c` factors, wall action, current, and flux claims. The observed 3D GR target remains flagged as a downstream conversion problem, not a pathA_19 base-system change.

## F4 paper-prose reconciliation

| Source | Statement | Classification | Note |
|---|---|---|---|
| `part01_parent_geometry.tex:140-144`, `174-203` | `rho=|psi|^2`, GNLS action/EOS/sound speed | AGREES | Implies 4D `rho=L^-4`, `psi=L^-2`, `K=M L^18 T^-2`. |
| `part01_parent_geometry.tex:213-219`, `268-291` | current, velocity, Euler/vorticity identities contain `m` | AGREES | Supports `m_GNLS` as action content. |
| `part01_parent_geometry.tex:298-330` | normalized projection and leakage source | AGREES/AMBIGUOUS | Agrees on open-system projection; normalized kernel keeps dimensions distinct from integrated 3D reduction. |
| `part01_parent_geometry.tex:447-510` | `a0=R0(0)`, `a(t)` mouth average | AGREES | `a` is a length and a collective moment. |
| `pde.tex:326-352`, `396-406` | parent action/EOS/current | AGREES | Same 4D parent dictionary as the harness. |
| `pde.tex:512-539` | projected continuity has `S_leak` | AGREES | Blocks unqualified no-net-accretion. |
| `pde.tex:633-648`, `1074-1118` | `a,L` are collective moments and reduced wall coordinates | AGREES | Supports reassessing the `a=1` pin. |
| `brane_bulk_ontology.tex:1267-1297`, `1967-1975` | mass as drainage/volume deficit; charge as vorticity flux | AMBIGUOUS | Physical scaling prose, not an action-level mass theorem. |
| `brane_bulk_ontology.tex:1998-2039` | bottom open/closed/connected | AGREES | Carries no-net-accretion as a gap. |
| `em_fields.tex:1717-1721` | `rho0` as `kg m^-3` mass density | WRONG-3D-CONVENTION | Legacy 3D/SI prose; not the 4D number-density harness dictionary. |
| `em_fields.tex:1723-1726`, `1738-1752` | `c_s`, `a`, `L`, circulation dimensions | AGREES | Correct in the stated SI/3D frame. |
| `em_fields.tex:1728-1736`, `1782-1786` | pressure/enthalpy per mass in SI units | WRONG-3D-CONVENTION | Useful prose but not the 4D action density/enthalpy dictionary. |
| `em_fields.tex:1757-1779` | `V=pi a^2L`, `m_G=kappa_m rho0 V`, `q` units | WRONG-3D-CONVENTION | Legacy 3D throat-volume bookkeeping; not a derivation of `m_defect` from `J`. |

## F5 gaps carried forward

- `INFLOW_MASS_SOURCE_MISSING`: BLOCKS_MASS_EMERGENCE. Retain base dimensions {L,T,M}; record J as a conserved rate label only. A positive m_defect result must be derived in pathA_21 from a defect source, boundary energy, or Hamiltonian/Noether charge. Source: part01_parent_geometry.tex:174-219 and pde.tex:326-406 keep m_GNLS in the exact action/current; brane_bulk_ontology.tex:1267-1297 gives drainage scaling, not a boundary/source/Noether/Hamiltonian derivation.
- `NO_NET_ACCRETION_BC_UNDERIVED`: CARRIED_FORWARD. Gauss flux is surface-independent only in no-source/no-leakage regions. No-net-accretion must be supplied as a boundary condition before using recirculation as a physical closure. Source: part01_parent_geometry.tex:298-330 and pde.tex:512-539 give projected leakage; brane_bulk_ontology.tex:1998-2039 leaves the throat bottom open/closed/connected.
- `A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT`: CARRIED_FORWARD. Use J and explicit branch data as invariant scale-map inputs; treat a as derived geometry. No normalization value is changed in pathA_19. Source: part01_parent_geometry.tex:447-510 and pde.tex:633-648,1074-1118 define a(t) as a mouth-average collective moment.
- `EOS_FROM_GNLS_FACTOR`: pathA_20 must keep the factor `h0=(m_GNLS*c_s0^2)/4` and the derived healing scale `sqrt(2)*hbar/(m_GNLS*c_s0)`.
- `M_TO_G_UNIFICATION`: pathA_21 must derive any defect-mass/back-reaction relation; pathA_19 does not prove it.
- `SCALE_MAP_INPUTS`: pathA_22 must consume `J`, branch geometry `a(J,branch)` if derived, `rho0`, `K`, `m_GNLS`, `hbar`, and the 3D reduction/conversion factors.

## Algebraic harness summary

- Algebraic checks: 17 consistent, 0 inconsistent, 17 total.
- Acceptance status: `PASS_WITH_NAMED_RESIDUALS`.

Flagged algebraic residuals:
- `formal_4D_R_norm_target_not_dimensionless_without_conversion`: FLAGGED_EXISTING_pathA_18_GAP; actual `L^-1 T^-2 M^-1`, factor needed `L T^2 M`. Preserve pathA_18 behavior and do not repair R_norm in pathA_19; the conversion belongs to pathA_22.
- `observed_3D_GR_target_not_dimensionless_without_conversion`: FLAGGED_FOR_pathA_21_pathA_22; actual `L^-2 T^-2 M^-1`, factor needed `L^2 T^2 M`. Do not interpret dimensional matching as a derived B2c normalization; the conversion belongs to pathA_21/pathA_22.
- `LT_R_norm_gate_fails_without_new_conversion_factor`: REJECTS_TRUE_LT_BASE; actual `L^-2 T^-2`, factor needed `L^2 T^2`. {L,T} is only a natural-unit representation in this step, not the base dimensional system.

