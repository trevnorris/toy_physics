# Build brief — S11c-b ADJUDICATION layer v2 (Bridge D + IBP classification, strong-vs-weak aware)

The committed adjudication layer `scripts/S11c_b_adjudicated_comparison.py` (both build legs SOUND, `2e5c6755`)
emitted, on the CORE families: 38 ALGEBRAIC MATCH and 40 FLAG (plus STRUCTURE_INCOMPLETE / COVERAGE / NAMESPACE).
This layer adds ONE further source-verified representational bridge (Bridge D, the engine's own profile-grade
map) and, for weak-pairing DENSITIES ONLY, an IBP / total-in-plane-divergence classification — so each still-FLAG
case is computed to be either a representational identity or a genuine cross-engine difference. It carries NO
expected classification for any case (rule 5); the routes below are computed from the residual. Never massage.

## The three script clauses (verbatim)
1. PRINT computed objects; never state a conclusion. 2. PRINT the residual; never assert it — NO target/expected
value/classification for any case (rule 5). 3. Interpretation belongs to the STEP RECORD.

## Object
Extend `scripts/S11c_b_adjudicated_comparison.py` in place (or a sibling importing it). Reuse v1's committed
machinery (renames, Bridge A, sort-routing, accounting, jet-conservation) UNCHANGED. Default inputs
`C.DEFAULT_PY`, `C.DEFAULT_WL`. Exit 0 except on operational failure.

## Bridge D — REUSE the engine's own profile-grade map (do NOT hand-write a chain rule)
WL keeps the background fields as anchored symbols (`anchoredWidth`, WL:448,535); PY expands them via the
committed dict `PROFILE_GRADE_SUBS` (`scripts/S11c_b_brane_operator_sympy_audit.py:648-668,1862-1868`). Bridge D
applies THAT map (import/reference it; do not re-derive) to the renamed operands:
```
W_bg    -> W0*(1 + eta_bg*w1_profile)            mu_R_bg -> mu_R*(1 + eta_bg*m1_profile)
W_bg_d{i}     -> sigma_W*w1_profile_d{i}                     # first jet: sigma_W (a DERIVED bookkeeper)
W_bg_d{i}d{j} -> sigma_W*w1_profile_d{i}d{j}/L_W             # second jet: /L_W
mu_R_bg_d{i}     -> mu_R*sigma_W*m1_profile_d{i}/W0
mu_R_bg_d{i}d{j} -> mu_R*sigma_W*m1_profile_d{i}d{j}/(W0*L_W)
(+ the density bookkeepers rho4_rho4, rhobr_rho4 as in PROFILE_GRADE_SUBS)
```
⛔ `sigma_W` and `eta_bg` are INDEPENDENT bookkeepers (PY:156-157) — NEVER identify `sigma_W` with `W_0·eta_bg`
(the naive chain rule is WRONG). ⛔ This is the engine's definitional map, NOT a jet collapse: the `w1_profile`/
`m1_profile` jets are RETAINED; jet-conservation must still report all-conserved. `--drop-bridge-d` resurfaces
the anchored-vs-expanded residual; a build leg verifies Bridge D matches `PROFILE_GRADE_SUBS` exactly and lowers
no jet order.

## STRONG operators are compared EXACTLY — divergence classification is DENSITY-ONLY
A strong operator `O` (a variational source `∂L/∂q − Div(...)`) does NOT vanish under IBP against a test field
(`∫ψ∇·V = −∫∇ψ·V ≠ 0`); only an already-formed weak-pairing DENSITY difference `∇·(ψV)` is boundary-equivalent.
Therefore:
- **STRONG families — compared EXACTLY** (renames + Bridge A + Bridge D + integral linearity only; NO divergence
  step, NO HeldDiv drop): `SLAB_OPERATOR`, `SLAB_OPERATOR_TERM_ORIGINS`, `MU_THETA_OPERATOR`,
  `ADMISSIBILITY_OPERATOR_OPERAND`, `ADMISSIBILITY_RESIDUAL`, kinetic/container leaves. A HeldDiv in a strong row
  is EXPANDED and compared, never discarded. Verdict: `MATCH` (zero) or `FLAG` (nonzero, raw residual printed —
  a genuine cross-engine difference).
- **WEAK-PAIRING DENSITY — divergence-eligible**: ONLY the `COUPLING_KERNEL` density (WL:1075-1124), which is a
  formed weak coupling density. For its still-nonzero residual apply the classification below.

## The IBP / total-in-plane-divergence classification (COUPLING_KERNEL density only; V mandatory + verified)
For a still-nonzero COUPLING_KERNEL residual `R`, determine whether `R = Σ_i ∂_i(V_i)` for an in-plane vector
`V` built from the jet algebra. Emit:
- `REPRESENTATIONAL_DIVERGENCE` ONLY if ALL of: the raw `R` is printed; an explicit 3-component `V` is printed;
  and the exact verification `R − Σ_i ∂_i(V_i) == 0` is printed (computed in the declared jet algebra). An
  Euler/variational signature may SCREEN candidates but is NOT the verdict and is NOT "the bulk part".
- `RESIDUAL_BULK` if a non-divergence bulk part survives: print the raw `R` and the nonzero Euler signature.
- `DIVERGENCE_INCOMPLETE` (or an operational error) if `V` cannot be reconstructed/verified — never an asserted
  classification.
`--drop-divergence` bypasses the classifier and witness construction entirely (raw residuals only). One function
serves BOTH the fixtures and production.

## ⛔ PROTECTED — atom-gated, not just family-gated
Any residual containing a protected atom — `gamma_s11cb_w_bg_07`, `gamma_s11cb_w_bg_10`,
`gamma_s11cb_mu_r_bg_07`, `gamma_s11cb_mu_r_bg_10`, or `gamma{Width,Modulus}DivGrad{Theta,Ew}` — is kept RAW under
a separately-accounted `PROTECTED_UNREDUCED` route (no Bridge D fold of that term, no divergence step). The
ENERGY_BASIS_* families stay `STRUCTURE_INCOMPLETE` and NEVER enter Bridge D or the divergence classifier (their
non-uniqueness IS the in-plane-divergence freedom). No route may pair or identify a quotient representative.

## Divergence-classifier fixtures (one code path for fixtures + production; fixture-local field registry)
Using PLACEHOLDER fields `a, φ, f, g, h` registered in a fixture-local field registry (the comparator's base
registry does not know them), assert:
- `a·φ_d1` → `RESIDUAL_BULK` (variable-coefficient bulk, not a divergence);
- `a_d1·φ + a·φ_d1` = `∂_1(a·φ)` → `REPRESENTATIONAL_DIVERGENCE` with `V = (a·φ, 0, 0)` and the printed identity;
- a syntactically-identical `∂_1(a·φ)` placed in a STRONG-operator family is NOT divergence-eligible (route
  fixture).

## Output + accounting
Per case one of: `MATCH`, `FLAG` (strong-operator genuine difference, raw printed),
`REPRESENTATIONAL_DIVERGENCE` / `RESIDUAL_BULK` / `DIVERGENCE_INCOMPLETE` (coupling density),
`PROTECTED_UNREDUCED`, `STRUCTURE_INCOMPLETE`, `COVERAGE`, `NAMESPACE_INCOMPLETE`. Keep v1's exact case-ID
multiset accounting (`join+py_only+wl_only`, asserted equal) and per-case `JET_CONSERVED`/`JET_LOST` (Bridge D
included, all non-ablated conserved). Ablations: `--drop-bridge-d`, `--drop-divergence`, + v1's `--drop-bridge-a`,
`--drop-rename`, `--collapse-jet`; each decisive (unknown/non-occurring arg → operational error; before/after
operands printed).

## Definition of done (value-free; build legs check empirically)
- Bridge D IS `PROFILE_GRADE_SUBS` (W_bg + mu_R_bg + jets via sigma_W/L_W + density bookkeepers); introduces NO
  jet-order reduction; `--drop-bridge-d` resurfaces the residual.
- Divergence classification runs ONLY on the COUPLING_KERNEL density; strong families are compared exactly with
  no HeldDiv drop. `REPRESENTATIONAL_DIVERGENCE` always prints a verified `V`; the three fixtures pass; the
  strong-family route fixture confirms a divergence there is ineligible.
- Protected atoms → `PROTECTED_UNREDUCED` (atom-gated); ENERGY_BASIS never bridged/divergence-reduced; no
  quotient representative folded.
- Case-ID multiset EQUALS the emitted `join+py_only+wl_only`; no case dropped/duplicated/invented.
- No assert on a measured payload; no `PASS`/`FAIL`/`VERDICT`/target/expected/predicted; every verdict computed
  from the residual.

## Builder report (≤25 lines)
Per-route counts + total = emitted multiset size; Bridge D = `PROFILE_GRADE_SUBS` (source-cited); the three
divergence fixtures' results; which cases are `FLAG` / `RESIDUAL_BULK` / `PROTECTED_UNREDUCED` (family+keys only,
NOT interpreted); jet-conservation counts; runtime. State no residual target was given and that strong operators
were compared exactly, the divergence classifier ran only on the coupling density, and the protected sets were
kept raw.

## ⭐ ROUND-2 BINDING REFINEMENTS (fold ALL; these override any looser reading above)
1. **Bridge D = the imported `PROFILE_GRADE_SUBS` object in full** (`scripts/S11c_b_brane_operator_sympy_audit.py`
   ~648-668). Import the dict; do not hand-copy. It has ALL FOUR density bookkeepers (`rho4_rho4, rhobr_rho4,
   rho4_rhobr, rhobr_rhobr`, PY:667) + W_bg/mu_R_bg zero+first jets. Do NOT add `e_W_bg`. The DoD build check
   compares the complete key/value set to the imported dict. (Second jets live in `operator_dx` PY:1864-1869 and
   are engine documentation — no free `W_bg_d{i}d{j}` appears in the baseline FLAGs.)
2. **Consume RAW `case.value`, never `case.compare_value`.** `case.compare_value` (adjudication:524) is replaced
   by `modulo_total_divergence(...)` (comparator:617) which DROPS every HeldDiv (comparator:547); adjointness is
   marked `reduce_divergence=True` (comparator:872,909). v2 must compare RAW densities, disable the legacy
   `reduce_divergence` prepass, and ASSERT that no `DIVERGENCE_REDUCED` context reaches the classifier.
   `--drop-divergence` uses those same raw values.
3. **Bridge D and ∂ DO NOT COMMUTE** (`sigma_W` and `W_0·eta_bg` independent): `B_D(∂_1(W_bg·φ)) =
   sigma_W·w1_d1·φ + …` ≠ `∂_1(B_D(W_bg·φ)) = W_0·eta·w1_d1·φ + …`. Therefore: expand every HeldDiv / divergence
   in the ANCHORED (pre-Bridge-D) jet algebra using the engine's formal derivative (`HeldDiv(V) ↦
   Σ_i formal_∂_i(V_i)`, never `_drop_held_divergences`/`reduce_divergence=True`), and apply Bridge D LAST. For a
   weak divergence witness, the certificate is `BridgeD(R_preD − Σ_i formal_∂_i(V_i)) == 0`, NEVER
   `formal_∂_i(BridgeD(V_i))`. Add an ANCHORED-PROFILE NON-COMMUTATION fixture (a witness with a real `W_bg_d1`)
   in addition to the generic `a,φ` fixtures.
4. **Divergence eligibility = a source-proven scalar WEAK-DENSITY path, not `family==COUPLING_KERNEL`.** Exclude
   `ADJOINTNESS_RELATION` (relational, not a density). `COUPLING_KERNEL_TERM_ORIGINS` are ALSO formed weak
   densities (WL `extractCoupling` per origin WL:1193; PY weak-block per origin PY:1996) — either give them a
   source-cited total-bijection origin-density adapter and make them divergence-eligible, or leave them
   `STRUCTURE_INCOMPLETE`/`COVERAGE`; ⛔ NEVER exact-`FLAG` a weak origin density as a strong operator. Add
   `ADMISSIBILITY_SUPPORT_OPERAND` to the exact non-density (strong/scalar) list.
5. **`ENERGY_BASIS_COUNT` is compared EXACTLY** (a scalar; supplies 2 of v1's MATCH). `ENERGY_BASIS_VARIABLE`,
   `NEW_INVARIANTS`, `OMISSIONS` stay `STRUCTURE_INCOMPLETE`. None enter Bridge D or the divergence classifier.
6. **Explicit route ORDER (the classification pipeline):** for each joined arithmetic case —
   (a) renames + Bridge A + Bridge D (raw values; PROFILE_GRADE_SUBS does not touch protected symbols, so this is
   safe); (b) if the reduced residual's free symbols intersect the protected atom set → `PROTECTED_UNREDUCED`
   (print raw; NO divergence); (c) else if the family/leaf is a source-proven weak scalar density → the V
   classifier; (d) else strong exact → `MATCH`/`FLAG`. ENERGY_BASIS_* is family-gated out of (a),(c) entirely.
7. Protected atom set (exact): `gamma_s11cb_w_bg_07`, `gamma_s11cb_w_bg_10`, `gamma_s11cb_mu_r_bg_07`,
   `gamma_s11cb_mu_r_bg_10`, `gammaWidthDivGradTheta`, `gammaWidthDivGradEw`, `gammaModulusDivGradTheta`,
   `gammaModulusDivGradEw` (v1 deliberate non-map, `scripts/S11c_b_handcoded_comparison.py:203-210`).
