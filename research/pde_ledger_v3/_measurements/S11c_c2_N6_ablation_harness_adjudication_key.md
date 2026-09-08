# S11c-c2 N6 ablation-harness ADJUDICATION KEY — ⛔ ORCHESTRATOR-ONLY (never handed to the blind builder)

Expected cones per cleared knife (directive `S11c_c2_N6_ablation_harness_directive.md`, both round-2 legs SOUND),
for adjudicating the harness build + transcripts. ⛔ Do NOT copy any expected value into the directive or a harness.
Baselines being made durable: `R_N6`=18/288 nonzero (non-RHO4 cases; MATERIAL_ADVECTED.RHO4 = 0/288);
`SPLIT_CHECK`=0 (exact affine identity); `CARRIER_BRIDGE_RESIDUAL`=0; `R_cov`=no-nonzero (δ≈2.6e-22, SymPy covariance).
Primes: SymPy `{1000000009, 998244353, 1000000033}`; WL `{1000000009, 998244353, 1004535809}` (blind, 2 shared/1 disjoint).
⭐ Global: every NO-OP (identity patch) diff MUST be identically 0 (no fabricated bite); every FORM knife must match
its moves/leaves cone; every ×2 companion must be a coefficient multiple of an existing term (not a new/vanished
structural coupling). Pinned case = `LAB_HELD × RHOBR_CONSTANT` for the SymPy harnesses; the WL harness runs all 4.

## Harness 1 — WL blind engine  (certified: WL_S11CC2_N6RC_R_N6, N6RC_SPLIT_CHECK, N6COV_R_COV, guards)
- **K_carrier** (`materialNormalKnife 0→1`, FORM): MOVES `CARRIER_MATERIAL`/`CARRIER_BRIDGE_RESIDUAL`/`CARRIER_CHANNEL`/
  `MATERIAL_OPERAND`/`R_N6`/`SPLIT_CHECK`; LEAVES `CARRIER_EULERIAN`, `SOURCE_EULERIAN/MATERIAL`, `SOURCE_BRIDGE_RESIDUAL`,
  `SOURCE_ACTUAL/PREDICTED`, `R_COV`, `FROZEN_PHI`. ⛔ `R_COV_INCREMENT` (reads `cm` 896-900) + `ACTUAL_CONTROL_PARAMETERS`
  (embeds the knob 977-980) are NOT in DEAD (they move). ×2 = `materialNormalKnife 2`.
- **K_junk** (`actualJunkCoefficient 0→1`, FORM; gate `actualJunkCase=MATERIAL_ADVECTED×RHO4`): MOVES `R_COV`/
  `SOURCE_ACTUAL`/`R_COV_INCREMENT`/`R_COV_CONTROL_DELTA`/`R_N6` in the junk case; LEAVES `CARRIER_*`,
  `SOURCE_PREDICTED`, `FROZEN_PHI`, `PHI_DOMAIN_CENSUS`. ×2 = `actualJunkCoefficient 2`.
- **K_split_route** (`sourceChannel` final `False→True`, FORM): admits the affine sig-0 `−C_M·p` family into the source
  channel → MOVES `SPLIT_CHECK` off the identity; LEAVES `CARRIER_*`, `SOURCE_EULERIAN/MATERIAL`, `R_COV`. ×2 = keep
  `False`, replace `es-ms`→`2(es-ms)`.

## Harness 2 — covariance  (certified: N6COV_R_COV, SOURCE_ACTUAL/PREDICTED, R_COV_INCREMENT, R_COV_CONTROL_DELTA, FROZEN_PHI, PHI_DOMAIN_CENSUS)
- **K_junk** (`ACTUAL_JUNK 0→1`, FORM): MOVES `R_COV`(4 cols)/`SOURCE_CONTROL_DELTA`; LEAVES `SOURCE_PREDICTED`,
  `FROZEN_PHI`, `PHI_DOMAIN_CENSUS`. ×2 = `ACTUAL_JUNK 2`.
- **K_circular** (run block, `mu_pred=mu_actual`, FORM): COLLAPSES `SOURCE_PREDICTED`→`SOURCE_ACTUAL` (V_E≡V_M) ⇒
  `R_cov`≡0 by construction (certifies the shipped non-circular prediction is load-bearing); LEAVES `SOURCE_ACTUAL`,
  `FROZEN_PHI`, `PHI_DOMAIN_CENSUS`. ×2 = `mu_pred=[2*t for t in predicted_amplitude(...)]`.
- **K_rank** (`image_of` rank≥2→parent image, FORM): MOVES `R_COV`(→18 in R_COV_INCREMENT); LEAVES `SOURCE_ACTUAL`
  (⛔ `FROZEN_PHI`/`PHI_DOMAIN_CENSUS` move — NOT in DEAD). ×2 = `(2 if rank≥2 else 1)*total_derivative(...)`.
- Coefficient companion `ACTUAL_A_RHO 1→2`: MOVES `R_COV`/`R_COV_CONTROL_DELTA`(84 cols); `SOURCE_PREDICTED` fixed.

## Harness 3 — diagnostic  (certified: REP_INVARIANCE_{EULERIAN,MATERIAL}_OPERAND + RESIDUAL [=R_N6], SLOT/CLOSURE guards, + TILT & N4 control triples)
- **K_EW_rowdrop** (`rows.pop('E_W')` after flatten, FORM): MOVES `m_rows`→`m_coeff`→`M`→`REP_INVARIANCE_RESIDUAL`
  (+ TILT/factory routes share face_factory); LEAVES `REP_INVARIANCE_EULERIAN_OPERAND` (imported `e_coeff`),
  `MU_RECONSTRUCTION_{IMPORTED,NATIVE,RESIDUAL}` (from `constitutive`). Companion `K_ewsign` (`'E_W':-folded['E_W']`)
  is COEFFICIENT. ×2 = `'E_W':2*folded['E_W']`.
- **K_slotdrop** (drop `delta_p_plus` from the `slots` tuple, FORM): MOVES the coeff extractions (797/809) → `E`/`M`/
  `R_N6`/guards; LEAVES `MU_RECONSTRUCTION_*`, `MU_AMPLITUDE`, `FACE_VELOCITY` (both routes). ×2 = local
  `pressure_coefficients` shadow scaling only the `delta_p_plus` column, all four slots kept.
- ⛔ No advection knife (dropped). The engine's own `CONTROL_INDEPENDENCE {BASE,CORRUPTED,RESIDUAL}` for `TILT` +
  `N4_ADVECTION` are emitted (nonzero at the pinned RHOBR case) as corroboration, NOT the harness's own cert.

## Harness 4 — reconcile  (certified: N6RC_R_N6, SPLIT_CHECK, CARRIER/SOURCE bridges + channels + operands)
- **K_normal** (scoped `build_material_face_source` wrapper: `components[0]+=a.grad_W[1]`, renormalize; restored in
  `finally`, FORM): MOVES `CARRIER_BRIDGE_RESIDUAL`(0→4)/`CARRIER_CHANNEL`/`C_M`/`R_N6`; LEAVES `SOURCE_EULERIAN/MATERIAL`,
  `SOURCE_BRIDGE_RESIDUAL`, `EULERIAN_OPERAND` (carrier-only: `m_v`=`build_material_velocity` runs BEFORE the wrapper).
  ⛔ NOT global, NOT Eulerian a:850, NOT `material_inverse_transpose` a:694. ×2 = `components[0]=2*components[0]`.
- **K_source_route** (`source,_=n.build_increment(...,slots)`, FORM): injects sig-0 `−C·p` → MOVES `SPLIT_CHECK`(0→12);
  LEAVES `CARRIER_*`. ×2 = keep `closed_response`, scale its `ds` arg ×2.
- **K_operand_swap** (`ms[s].get(...)`→`es[s].get(...)` at 211-212, `es−es` collapse, FORM): MOVES `SOURCE_CHANNEL`/
  `CROSS_CHANNEL`/`SPLIT_CHECK`(0→18-ish); LEAVES `CARRIER_*`. ⛔ the whole-operand ms↔es swap is `−ds` = ×−1 coefficient
  (excluded). ×2 = `2*ms[s].get(...)`.
