⚠ **WARNING — this file holds S9 result values and expected verdicts.** If you are building or
reviewing the comparator, do not read it. That is a warning, not a control.

# S9 pilot adjudication

## 1 · Status header

**Status: ORCHESTRATOR-ONLY. ⛔ Never readable by a comparator builder. Applied only after the generic
comparator is frozen. Folds into the step record after the pilot closes.**

## 2 · Expected per-row verdicts

Today's comparator over S9's committed cross-engine pairs gives the A1 differential target below.

| quantity | Mathematica tag | SymPy tag | verdict |
|---|---|---|---|
| `factored_determinant` | `WL_S9_DETERMINANT` | `PY_S9_MAIN_DET_M_FACTORED` | `AGREE` |
| `full_root_set` | `WL_S9_ROOT_MULTISET` | `PY_S9_MAIN_ROOT_MULTISET` | `AGREE` |
| `transverse_multiplicity` | `WL_S9_ROOT2_E2` | `PY_S9_MAIN_ROOT_E2` | `AGREE` |
| `transverse_speed_squared` | `WL_S9_CANDIDATE_SPEED_SQUARED1` | `PY_S9_MAIN_SPEED_SQUARED_CANDIDATES` | `AGREE` |
| `dispersion_scaling_residual_flexural` | `WL_S9_X7_ROOT1_SCALING_RESIDUAL` | `PY_S9_X7_ROOT_SCALING_RESIDUAL` | `AGREE` |
| `inertia_coefficient_dimension` | `WL_S9_RHO_DIMENSION` | `PY_S9_MAIN_DIM_PRIMARY_INERTIA` | `AGREE` |
| `stiffness_coefficient_dimension` | `WL_S9_MU_DIMENSION` | `PY_S9_MAIN_DIM_PRIMARY_STIFFNESS` | `AGREE` |
| `coefficient_dimension_difference` | `WL_S9_MU_RHO_DIMENSION_DIFFERENCE` | `PY_S9_MAIN_DIM_STIFFNESS_MINUS_INERTIA` | `AGREE` |
| `implied_speed_dimension` | `WL_S9_SPEED_SQUARED_IMPLIED_DIMENSION` | `PY_S9_MAIN_DIM_SPEED_FROM_EXPRESSION` | `AGREE` |
| `speed_dimension_difference` | `WL_S9_SPEED_SQUARED_DIMENSION_DIFFERENCE` | `PY_S9_MAIN_DIM_SPEED_DIFFERENCE` | `AGREE` |
| `dynamical_matrix_route_residual` | `WL_S9_DYNAMICAL_MATRIX_RESIDUAL` | `PY_S9_MAIN_ROUTE_OPERANDS_AND_RESIDUAL` | `AGREE` |
| `bare_field_coefficient_dimension` | `WL_S9_X8_MU_G_DIMENSION` | `PY_S9_X8_DIM_COEFFICIENTS` | `AGREE` |

## 3 · Load-bearing declarations

Literal stdout from `reduction/measurements/declaration_load_ablation.py S9`:

```
STEP=S9 CONFIG=/var/projects/toy_physics/research/pde_ledger_v3/reduction/checks_S9.yaml
BASELINE CROSS_ENGINE: agree=12

DECLARED_NAMING_EXCEPTIONS=6 DECLARED_SYMBOL_IDENTITIES=1

DROP naming[rho_z]  CROSS_ENGINE: agree=12
  changed_rows=0
DROP naming[mu_F]  CROSS_ENGINE: agree=11 disagree=1
  changed_rows=1
  transitions={'AGREE->DISAGREE': 1}
  rows=dispersion_scaling_residual_flexural
DROP naming[mu_G]  CROSS_ENGINE: agree=12
  changed_rows=0
DROP naming[lambda_rho]  CROSS_ENGINE: agree=12
  changed_rows=0
DROP naming[lambda_mu]  CROSS_ENGINE: agree=12
  changed_rows=0
DROP naming[lambda_scale]  CROSS_ENGINE: agree=11 naming_mismatch=1
  changed_rows=1
  transitions={'AGREE->NAMING_MISMATCH': 1}
  rows=dispersion_scaling_residual_flexural
DROP identity[py:omega2]  CROSS_ENGINE: agree=11 disagree=1
  changed_rows=1
  transitions={'AGREE->DISAGREE': 1}
  rows=factored_determinant
```

## 4 · S9 result values the plan must not carry

| removed item | repository-recorded value | source |
|---|---|---|
| computed-value outcome | No S9 value moves. | `steps/S9_light_requires_shear.md` |
| S9 baseline aggregate | `12/12 agree, 0 disagree` | current comparator over the committed engine outputs |
| S10 verdict cells | `545/690`; `145 no verdict`; `11 false DISAGREE` | `steps/S10_two_transverse_photons.md` and the prior plan measurement |
| squared-speed/result form and old L3 exposure | `c² = μ_R/ρ_br`; the removed L3 cases applied coefficient `17`, a sign inversion, and fabricated dimensionless `E2=7`, and each reported `L3 VERDICT: PASS` | committed S9 engine outputs and `steps/S9_light_requires_shear.md` |
| propagating transverse mode count / expected witness multiplicity | `E2 = 2` | committed S9 engine outputs and `steps/S9_light_requires_shear.md` |
| concrete symbol identity | `omega2 = omega**2` | `reduction/checks_S9.yaml` |
| root-set cardinality | `3` | `reduction/checks_S9.yaml` and committed S9 engine outputs |

## 5 · Held-out mutants

⛔ **EMPTY.** Held-out mutants must be authored only after the comparator interface is frozen, by a party
that is not the builder.

## 6 · Probe points / seeds

⛔ **EMPTY.** Draw after implementation and re-draw for every adjudication run.

## 7 · Control structural expectations

| control | what must move |
|---|---|
| X1 · inertial-coefficient control | the determinant and propagating-root location must change |
| X2 · stiffness-coefficient control | the determinant and propagating-root location must change |
| X3 · divergence-only action | the longitudinal/transverse response must swap: a longitudinal root appears and the transverse family changes |
| X4 · gradient-elastic action | a longitudinal root must appear in addition to the transverse family |
| X5 · stiffness-sign control | the sign of an emitted squared-frequency root must change |
| X6 · anisotropic-inertia control | the transverse multiplicity must reduce and a coincident root must split |
| X7 · flexural control | the root-scaling response must change |
| X8 · bare-field control | a root location and the root-scaling response must change |
