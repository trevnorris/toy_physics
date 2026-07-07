# ledger_stage004_gnls_action_dimensional_foundation

## Status

EARNED (dimensional foundation) -- `RETAIN_L_T_M`, `PASS_WITH_NAMED_RESIDUALS`.

The exact-symbolic claims of this stage (the 4D action dimensional dictionary,
the derived sound-speed/healing-length anchors, the pin null-relation, the flux
dimensions, and the base-system verdict) are closed. The honest landing carries
**named residuals** -- a labeled non-derivation (`m_defect` is NOT emergent here)
and three deferred `{L,T}`-conversion gaps -- which are first-class and are NOT
softened or repaired in this stage (they are the obligations of pathA_21/22).

This is the first stage of **Part I (The medium)**: it fixes the dimensional
substrate every later sector reduces from.

## Purpose

Establish the dimensional foundation of the medium: the field content and
`{L, T, M}` dimensional dictionary of the GNLS parent action, the derived
sound-speed and healing-length anchors, the pin null-relation
`a = hbar/(m_GNLS c_s0)`, and the base-system verdict `RETAIN_L_T_M` -- the base
dimensional system is `{L, T, M}`; the `{L, T}` projection is a natural-unit
*representation only*, not a base-system change. The honest carried gap: defect
(gravitational) mass `m_defect` is not shown emergent at this gate.

## Provenance

Earned source (provenance, not content): `software/stage1_solver/reports/
pathA_19_dimensional_foundation.md`, verdict `RETAIN_L_T_M` /
`PASS_WITH_NAMED_RESIDUALS`. The SymPy logic reshaped here was extracted from the
shared harness `software/stage1_solver/src/stage1_solver/dimensional_check.py`
(`run_patha19_foundation`, `--patha19-foundation`; 14 foundation + 3
`{L,T}`-representation checks = 17, already passing). The derivation below is
inlined so a reader never needs to open that report or the harness.

Script-backed status: the exact residuals are checked by
`scripts/ledger_stage004_gnls_action_dimensional_foundation_sympy_audit.py` and
independently by
`mathematica/ledger_stage004_gnls_action_dimensional_foundation_mathematica_audit.wl`.
The `.wl` is a genuinely independent route (native `UnitDimensions`/`NullSpace`
construction), not a transliteration of the `.py`. Both scripts include
able-to-fail mutation probes (dimensional-firewall ablations + a pin-corruption
break) and derive -- do not hardcode -- the dictionary, the pin relation, the
healing scale, and the verdict.

## 0. Why needed

Every downstream sector (gravity, light, charge, magnetism) consumes the
medium's dimensional dictionary, its sound speed `c_s0`, and its healing scale.
Before any dynamics, the dimensional substrate must be pinned and its base
system settled: is mass an *emergent* quantity (so the base system could reduce
to `{L, T}`), or an *explicit action constant* (so `{L, T, M}` is retained)?
This stage answers that question and hands the derived anchors forward.

## 1. The GNLS parent action and field content

The parent medium is a Gross-Pitaevskii / nonlinear-Schrodinger (GNLS) condensate
in the 4D bulk. The order parameter is `psi`, with number density
`rho = |psi|^2`. The action density carries the Schrodinger kinetic term, the
gradient (quantum-pressure) term, and the equation-of-state term:

```text
L ~ i hbar psi* d_t psi  -  (hbar^2 / 2 m_GNLS) |grad psi|^2  -  U(rho),
U(rho) = (K/4) rho^5,   so   P = K rho^5   (the imposed EOS closure).
```

`m_GNLS` is the constituent (action) mass; `hbar` and `K` are the other action
constants. These are POSTULATED action content, labeled as such.

## 2. The 4D-bulk dimensional dictionary (two tiers)

Dimensions are carried as exact `{L, T, M}` exponent triples. The dictionary has
two tiers, kept distinct (matching the source F3 independent-vs-derived table).

**(i) Primitive input dimensions** -- posted, not derived from nothing; their
only obligation is a consistency check against their action usage:

```text
[hbar]    = M L^2 T^-1     (independent action constant)
[m_GNLS]  = M              (independent action constant)
[rho0]    = L^-4           (chosen state datum, 4D-bulk number density)
```

**(ii) Derived-by-composition dimensions** -- forced by dimensional homogeneity
of the action terms:

```text
[psi]  = L^-2         from the kinetic vs gradient term equality
[rho]  = L^-4         = [psi^2]
[K]    = M L^18 T^-2  from (K/4) rho^5 sharing the kinetic-density dimension
```

Homogeneity is asserted across the full 4D action dictionary -- the 14 harness
foundation checks (GNLS kinetic + gradient + EOS/enthalpy density, the gauge
coupling `q A_i/hbar`, the Maxwell sector with explicit `c` factors, the wall
action, the current `j = rho v`, the flux, and the **bulk continuity equation**
`d_t rho + div(rho v) = 0`) plus the 3 `{L,T}`-representation checks
(`local GNLS`, `local Maxwell`, `local wall`, evaluated in the M-projected
`{L, T}` representation; these feed the verdict discriminator in Section 6).

## 3. Derived velocity and healing anchors

The sound speed follows from the EOS (its full derivation is owned by stage005 /
pathA_20; here it appears as a consequence of the dictionary and is dimensionally
confirmed):

```text
c_s0^2 = 5 K rho0^4 / m_GNLS,     [c_s0] = L T^-1.
```

The GNLS core balance fixes the enthalpy scale and the healing length:

```text
h0 = (5K/4) rho0^4 = (m_GNLS c_s0^2)/4,
xi_h = sqrt( hbar^2 / (2 m_GNLS h0) ) = sqrt(2) hbar / (m_GNLS c_s0),   [xi_h] = L.
```

`xi_h` is DERIVED from `h0` (not hardcoded). Carried forward to pathA_20 as
`EOS_FROM_GNLS_FACTOR`.

## 4. The pin null-relation

Impose the four natural pins `{a = 1, c_s0 = 1, hbar = 1, m_GNLS = 1}` on the
three base dimensions `{L, T, M}`. The pin exponent matrix has **rank 3, nullity
1**: exactly one dimensionless monomial survives. Its null vector gives

```text
a c_s0 m_GNLS / hbar = 1,     i.e.   a = hbar / (m_GNLS c_s0).
```

Consequence: if `a` is identified with the GNLS healing core, then
`a / xi_h = 1/sqrt(2)` (a convention/branch factor); the raw four pins correspond
to `a = xi_h / sqrt(2)`. Because `a` is fixed by the pin choice (a mouth-radius
collective moment), not by a base dimension, it is carried as
`A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT`; the conserved rate `J` is the invariant
label.

## 5. Flux dimensions

```text
[J_bulk]  = [J_brane]   = T^-1        (number rate)
[Q_vol,bulk] = rho^-1 J = L^4 T^-1
[Q_vol,brane]           = L^3 T^-1
[m_GNLS J]              = M T^-1       (mass-per-particle times number rate)
```

Gauss shape-independence holds only in a no-source/no-leakage region; the
projected continuity carries a leakage source `S_leak`, so no-net-accretion is a
boundary condition, not derived here (`NO_NET_ACCRETION_BC_UNDERIVED`).

## 6. The verdict `RETAIN_L_T_M` (two-conjunct)

The base-system verdict is COMPUTED from the conjunction of two conditions:

- **(a) the `{L,T}`-rejection gate.** Under the full `{L, T, M}` dictionary every
  Section-2 check is a zero residual (consistent). Under the M-dropped `{L, T}`
  representation, at least one gate is non-dimensionless -- the three flagged
  targets need `M`-bearing conversion factors:

  ```text
  formal_4D_R_norm    : actual L^-1 T^-2 M^-1,  factor needed  L T^2 M
  observed_3D_GR      : actual L^-2 T^-2 M^-1,  factor needed  L^2 T^2 M
  LT_R_norm gate      : actual L^-2 T^-2,       factor needed  L^2 T^2   (REJECTS_TRUE_LT_BASE)
  ```

  So `{L, T}` cannot be the base system; it is a natural-unit representation only.

- **(b) the mass fork.** `m_GNLS` remains an explicit action constant of
  dimension `M`, while `m_defect` is NOT derived here
  (`INFLOW_MASS_SOURCE_MISSING`). `M` is therefore genuinely load-bearing and
  irreducible: `m_defect_derived_here = False`.

Both conjuncts hold, so the verdict is `RETAIN_L_T_M`. The verdict is able-to-fail
on both: restoring `M` (or forcing the LT gates dimensionless) changes conjunct
(a); flipping `m_defect_derived_here = True` (counterfactual) removes conjunct
(b). Neither is hardcoded.

## 7. Carried gaps (first-class, not repaired)

- `INFLOW_MASS_SOURCE_MISSING` (BLOCKS_MASS_EMERGENCE): `m_defect` must be derived
  in pathA_21 from a defect source, boundary energy, or Noether/Hamiltonian
  charge. The candidate bridge `hbar J / c_gamma^2 = M` is a **dimensional
  conversion only**, not a mass theorem.
- `LT_R_norm_gate_fails_without_new_conversion_factor` (REJECTS_TRUE_LT_BASE),
  `formal_4D_R_norm_target_not_dimensionless_without_conversion`,
  `observed_3D_GR_target_not_dimensionless_without_conversion`: the conversions
  belong to pathA_21/22; not repaired here (pathA_18 behavior preserved).
- `A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT`, `EOS_FROM_GNLS_FACTOR`,
  `NO_NET_ACCRETION_BC_UNDERIVED`, `M_TO_G_UNIFICATION`, `SCALE_MAP_INPUTS`:
  carried forward verbatim as provenance.

## What this achieves physically

The medium's dimensional substrate is fixed: a `{L, T, M}` base system with
`hbar`, `m_GNLS`, `K` as independent action constants and `rho0` a chosen state
datum, from which `c_s0`, `xi_h`, and the pin geometry `a` follow. Mass is an
explicit action parameter, not an emergent quantity -- the first, load-bearing
scoping decision of the whole program. Every later Part inherits this dictionary.

## What is still missing

Defect-mass emergence (`m_defect` from inflow), the `{L,T}` conversion factors,
and the no-net-accretion boundary condition are all deferred (pathA_21/22),
carried as named residuals. This stage settles the base system and the anchors;
it does not close those gaps.

## Next step

Stage005 (I-2) derives the sound speed `c_s^2 = 5 K rho^4 / m` proper from the
EOS and establishes the light/sound ratio `lambda_gamma = c_gamma/c_s` as a free
calibration input (`C_GAMMA_RATIO_UNDERDETERMINED`), consuming this stage's
dictionary and `EOS_FROM_GNLS_FACTOR` handoff.
