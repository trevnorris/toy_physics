# S11c-c1 T7 cross-engine reconcile — the adjudication of the comparator's printed residuals

Rule 2 binds the orchestrator: every claim below carries the command that produced it. The comparator
`scripts/S11c_c1_cross_engine_comparator.py` (SOUND, reviewed baseline `7141e6ad` + repair `704308af`)
**computes and PRINTS** the paired residual operands, deciding nothing (contract `S11c_c1_SHARED_PHYSICS.md:580-587`).
This file is OUR post-run adjudication: for each surfaced residual, is it a **spurious representational/keying
difference** (the engines AGREE on the object) or a **real physics DISAGREEMENT** (rule 6 finding)? The three
values used below — AGREE / DISAGREE / UNDECIDED(coverage-bounded) — are OUR reading, ⛔ not a script token.

⭐ **HEADLINE (scoped after the 2 step-record review legs — Grok + Codex sol — corrected an earlier
overreach).** The **core two-momentum DtN kernel AGREES cross-engine** — all 4 (anchoring×face) joined cases
collapse to **exact zero** off-diagonal under physically-justified identifications; the collapse is load-bearing
(a wrong jet sign or freezing the two legs into one leaves it nonzero). The **face-response coefficients AGREE**
for LAB_HELD (off-diagonal exact zero) and for MATERIAL_ADVECTED on the physical diagonal `k=k′`. **No
cross-engine DISAGREEMENT was found among anything measured.** ⛔ But a **closure-wide "AGREE" is NOT earned on
this box**: seal 5 (background density) stays a **surfaced rule-17 freeze (UNDECIDED)**, the raw DTN_OPERATOR
form, the ENERGY audit, and a flat-resolvent leg-labeling (off-diagonal, MATERIAL) are **UNDECIDED**, and the 4
giant families + the full per-family symbolic residual are **UNMEASURED — DEFERRED** to the ≥64 GB box
(`DEFERRED_HEAVY_RUNS.md`), ⛔ not pre-adjudicated. Verdict: **the kernel is cross-engine closed; the rest is
UNDECIDED/deferred.** Reproduce: `python3 _measurements/S11c_c1_reconcile_reproduce.py --runlog <run> --adversarial`.

---

## 1. The runs (family-scoped; one heavy family at a time; RSS-watched)

**Light tier — 42 families** (all 50 shared non-LOCAL minus the 4 giants minus 4 mid-heavy):
```
python3 scripts/S11c_c1_cross_engine_comparator.py --family <each of 42>   # detached, /usr/bin/time -v
```
Literal footer: `RUN_ACCOUNTING families=42 families_with_join=24 families_with_unpaired=41 parse_failed=5
duplicate_key=0 runtime_seconds=1112.847`; `MEASUREMENT_SCOPE supplied_unfalsifiable=sections_1_to_2,
supplied_substrate,mu_theta_operand,supplied_bookkeeping residual_target=none`. **EXIT=0, peak RSS 308,348 kB
(~301 MB)** — well under the 30 GB box.

Header (the R1 repair is live): `SPELLING_INJECTIVITY collisions=0 reserved_names=204 mechanical_spellings=204
active_c1_folds=[('LW','L_W'),('LambdaA0','Lambda_A_0'),('LambdaV0','Lambda_V_0'),('LambdaX0','Lambda_X_0'),
('W0','W_0'),('cS0','c_s0'),('epsilonShape','epsilon_shape'),('etaBg','eta_bg'),('rhoM','rho_m'),
('sigmaW','sigma_W'),('tauA','tau_A'),('tauV','tau_V')]` — 12 folds incl. `cS0←c_s0`, injective.

Per-family accounting: `scratchpad/reconcile/accounting.txt`. 24 families join (1000+ joined cases); 18 are
fully unjoined by keying (the seals) — see §3.

**Mid tier — 4 families** (FACE_RESPONSE, CONTROL_BRANCH_RESIDUAL, CONTROL_FORM_RESIDUAL,
UNIFORM_LIMIT_S11B_OPERAND): see §5 (running / result folded in).

---

## 2. The adjudication method — a staged representational bridge (⛔ not tuned to zero)

The comparator emits the residual WITHOUT applying any representational identification (rule 2), so a nonzero
residual is EXPECTED where the engines represent the same object differently. To decide AGREE vs DISAGREE I
apply the **known, physically-justified identifications** and test whether the residual collapses. Each
identification is justified independently of the outcome (⛔ never "assume they're equal then find zero"):

| PY (opaque / frozen) | WL (live) | justification |
|---|---|---|
| `s11cc1_q_out_output` / `_input` | `qOut(ω,{k})` / `qOut(ω,{k′})` | seal 3: both = radiation-selected normal momentum for the two legs |
| `s11cc1_w1_profile_hat_transfer` | `HeldInactiveFourierTransform(w1Profile(x/L_W),{x},{k−k′})` | both = FT of the thickness profile at the momentum transfer |
| `s11cc1_w1_profile_jet_hat_i` | `i·L_W·(k_i−k′_i)·hat_transfer` | FT-of-derivative identity: the face-tilt profile IS the gradient of the thickness profile; WL builds `slopeHat=i·σ_W·profileHat·(L_W·k_transfer)/2` from the height FT (`.wl:176-179`, `unitSlopeWeights={1,1,1}`), PY freezes it opaque (`.py:200-207`, `dtn_first_kernel:606`) |
| `omega` (real=True) | `omega` (plain) | seal 2: assumption difference; unified |
| `rho_br_bg_rho4_constant` (bare) | `rhoBrBgRho4Constant(x)` (applied field) | seal 5: pointwise-equal (see §3.5) |
| PY μ_θ = `ε·mu_theta_inputs` (`.py:840`) | WL `ε·muTheta` (`.wl:748`) | seal 1: ε sits inside PY's μ_θ symbol, in WL's coefficient ⇒ `PY_coeff·ε = WL_coeff` |
| on-shell: `q² = ω²/c_s0² − |k|²` | encoded inside `qOut` | PY writes `κ²=ω²/c_s0²` explicitly; WL folds dispersion into `qOut` |

Scripts committed beside this record: `S11c_c1_reconcile_dtn_kernel_bridge.py` (symbolic staged bridge),
`S11c_c1_reconcile_response_coeff_bridge.py` (`ε·A≡B` numeric), `S11c_c1_reconcile_keyscan.py` (key-axis scan).
The numeric zero-test assigns random rationals to all base symbols (ω large so on-shell `q²>0`), enforces the
σ_W=η_bg·W_0/L_W binding, sets `Q=√(ω²/c_s0²−|k|²)`, and checks `|A−B|/(|A|+|B|) < 1e-9` over ≥10 trials
(Schwartz–Zippel zero-test for a rational identity).

---

## 3. The five known seals — adjudications

### 3.1 Seal 3 (two-momentum legs) — **AGREE**
`DTN_KERNEL` joins on `(BRANCH,FACE,LEAF=KERNEL_EXPRESSION)`, j=4, all nonzero raw. **Both engines carry the two
distinct legs** — WL `qOut(ω,{kOne,kTwo,kThree})` / `qOut(ω,{kPrimeOne,kPrimeTwo,kPrimeThree})`, PY opaque
`s11cc1_q_out_output` / `_input` — with the same two-momentum structure `q(k)q(k′)+k·k′−ω²/c²`. Neither froze to
a single k (rule-17 mode-mixing intact).
Symbolic bridge (`S11c_c1_reconcile_dtn_kernel_bridge.py`, LAB_HELD FACE=1): after the map, the Stage-2 residual
is EXACTLY `I·HAT·W_0·η_bg·ω·ρ_m·(Q_IN²·c_s0² + c_s0²·|k′|² − ω²)/(2·Q_IN·Q_OUT·c_s0²)` — i.e. the residual **is
the radiation-branch dispersion polynomial** `c_s0²(Q_IN²+|k′|²−ω²/c_s0²)`. Stage-3 (apply on-shell) → `r3 = 0`.
Numeric confirm (`numbridge`): `[DTN_KERNEL LH F1] trials=10 worst_rel=0.00e+00 AGREE`. That the residual
collapsed to the *defining dispersion relation* (not a coincidental match) is the strongest confirmation.

### 3.2 Seal 2 (ω real-assumption) — **AGREE**
`Symbol('omega', real=True)` (PY, from the base ledger) vs `omega` (WL, plain, with a local
`Assumptions->Element[omega,Reals]&&omega!=0` at `.wl:1045`). Unified in the kernel bridge (`plainize`); it
contributes nothing to the residual (the kernel collapses to 0 with omega unified). The `ZERO_JET_RESIDUAL`
family is unjoined (axm=3) — the state's "numerically 0" note is the PER-ENGINE zero, not a cross-engine join.

### 3.3 Seal 1 (μ_θ) — **AGREE** (scoped: LAB_HELD exact zero; MATERIAL diagonal-exact)
`FACE_RESPONSE_COEFFS` joins j=96. μ_θ is **opaque/supplied on both sides** (in `MEASUREMENT_SCOPE`
`mu_theta_operand`, unfalsifiable within the build) and is FACTORED OUT by the coefficient extraction (PY per-
face `s11cc1_mu_theta_<anchor>_<±>`, WL single `muTheta`) — so the μ_θ symbol never enters the coefficient
residual; the *coefficient* is the real cross-engine test. The residual mismatch was diagnosed as a clean
`a/b = 1/ε_shape` (zero imaginary part) — PY bakes ε into its μ_θ (the full-response constructor
`response_operator_case` `.py:840`; the coefficient constructor `closed_coefficients` `.py:723-755` omits ε), WL
keeps ε in the coefficient (`.wl:748,866-870`). Under `PY_coeff·ε = WL_coeff`:
`[MU_THETA_COEFF] worst_rel=0.00e+00 AGREE`, `[FACE_VELOCITY_COEFF] worst_rel=4.33e-20 AGREE`.

⚠ **Coverage (corrected after Codex sol's leg flagged the earlier "proven EXACT ZERO" as overstated):** the
exact-zero collapse holds for the **LAB_HELD** coefficients (off-diagonal) and all 4 **DTN_KERNEL** cases. A broad
test over all 96 joined coefficient cases (`S11c_c1_reconcile_broad_coeff_test.py`) split cleanly by anchoring:
**LAB_HELD collapses off-diagonal; MATERIAL_ADVECTED does NOT** (worst_rel ~1e-2…1e-4, real at 50-digit precision
`S11c_c1_reconcile_flat_leglabel_diagonal.py`). Cause: the **flat half-space resolvent's leg labeling** — PY
writes the flat DtN symbol on the OUTPUT leg (`s11cc1_q_out_output`, `ρ_m ω/q_out`, `.py:602-603`), WL on the
INPUT leg (`qOut(ω,{k′})`); equal on the physical **diagonal k=k′** (where the half-space object lives — retest
there gives `rel=0`), differing only OFF-diagonal where the randomized bridge wrongly probed it. ⇒ the MATERIAL
coefficients AGREE on the diagonal; the off-diagonal flat-leg labeling is a **surfaced representational
difference (UNDECIDED)** — a diagonal-restricted or the deferred symbolic residual settles whether it is purely a
convention.

### 3.4 Seal 4 (regime / parity naming) — **AGREE (content); spurious keying convention**
`DTN_BY_REGIME_PAIR` is fully unjoined (j=0, pyO=72, wlO=180) for ONE reason surfaced in the case notes:
`REGIME_OUT:EVANESCENT!=OUTPUT_EVANESCENT, REGIME_IN:EVANESCENT!=INPUT_EVANESCENT` — WL prefixes the leg
(`OUTPUT_`/`INPUT_`), PY uses the bare regime, and WL splits leaves INPUT_/OUTPUT_ where PY has one. The underlying
kernel content is the proven kernel (§3.1). `DTN_BY_PARITY` **joins** on `LEAF=PARITY_MATRIX` (j=2) and the
operand IS the two-momentum kernel arranged as a parity matrix — same residual structure as §3.1, collapses
identically. Parity naming `THICKNESS/CENTRE_SHIFT` vs `PARITY_DELTA_W` lives on the DEFERRED PERMEABLE_* axis.
⇒ the regime/parity axis NAMING is a bridgeable keying convention (`OUTPUT_X↔X`, `INPUT_X↔X`); the engines agree
on branch/regime physics (contract §1b: branch selection is physics — surfaced, not pre-folded).

### 3.5 Seal 5 (background density; rule 17) — **UNDECIDED (surfaced freeze)** ⛔ NOT a clean AGREE
_(revised: both step-record legs pushed here; Grok read it AGREE-with-caveat, Codex sol read it a mis-adjudicated
freeze. Rule 17 is GOVERNING and the cost of a wrong "AGREE" on a frozen field is exactly what this rebuild
exists to catch — I adopt the conservative reading, Codex's.)_
Surfaces in `FACE_RESPONSE_COEFFS` MU_THETA_COEFFICIENT: PY `Symbol('rho_br_bg_rho4_constant')` (a bare CONSTANT)
vs WL `Function('rhoBrBgRho4Constant')(xOne,xTwo,xThree)` (an applied LIVE field, bound at `.wl:192` to
`(rhoBr/W0)·WBg(x)` with `WBg` a spatially-varying field). **What IS established:** the two engines' *functional
dependence* on the density is identical — both use it as a local reciprocal `μ_s = μ_θ/ρ` (`.py:841`, `.wl:748`),
and NEITHER engine differentiates it (rule-17 test, grep both engines: `D[…rhoBrBgRho4Constant…]` /
`Derivative`/`diff(…rho_br_bg_rho4…)` → **0 matches**; verified independently by both legs). The μ_θ-coefficient
residual collapses to 0 ONLY after I MAP the WL field to PY's constant (`rhoBrBgRho4Constant(x)→RHO`).
**⛔ That mapping is the point of contention: a bare scalar is NOT globally equal to a live varying field — folding
one to the other is a [[feedback_handcode_comparison_never_blanket_collapse]]-style collapse of exactly the kind
rule 17 forbids.** The comparator correctly SURFACES the difference (bare-PY vs applied-WL); it is not the
orchestrator's to fold away pre-adjudication. ⇒ **Seal 5 stays a SURFACED rule-17 representational difference —
UNDECIDED.** PY's freeze is *harmless within c1's local first-order emit* (ρ multiplies first-order shape
quantities, so ρ(x)−ρ₀ = O(η) makes the difference O(η²)), but PY's EXPORT carries a constant-ρ response where
WL's carries a field-ρ(x) response. ⛔ **MANDATORY re-adjudication in c2** (not conditional): c2's self-energy fold
sums/couples over the face, which is exactly where ρ's x-dependence becomes load-bearing — resolve it with a valid
field-vs-field comparison there (or the deferred symbolic residual). [[feedback_never_freeze_a_varying_field]]
[[feedback_basis_independence_must_not_freeze_spurion]]

---

## 4. Additional surfaced findings (⛔ not among the 5 seals; surfaced by the comparator, rule 6)

- **DTN_OPERATOR — FACE-axis keying.** PY emits per-face (`FACE=±1`), WL per-branch (WL missing FACE); the DtN
  operator is a raw whitelisted operand (`outer_operator_algebra_signature; raw_operand_emitted`). WL treats the
  operator face-independent; PY carries both faces. Keying/coverage — the raw noncommutative operator-algebra
  form is not leaf-extractable; full adjudication DEFERRED to the ≥64 GB symbolic run.
- **DTN_HERMITIAN_PART / DTN_REACTIVE_PART — AGREE at c1's first order; O(η²) → S11c-e.** Derived from the proven
  kernel by the adjoint (leg-swap k↔k′ + conjugate). They carry ONE representational axis beyond the 5 seals:
  WL's `weightedAdjointKernel` weights by `BoundaryMeasure`/`InverseBoundaryMeasure` (the true-area element),
  PY's `hermitian_kernel` (`.py:903`) is a plain conjugate+swap. Area = √(1+O(slope²)) = 1+O(η²), and c1 works to
  first order in the shape ⇒ BoundaryMeasure→1 and the weighted adjoint ≡ PY's plain adjoint at working order.
  Full symbolic confirmation (incl. the WL swap-temp `momentumSwapOutput`) DEFERRED.
- **ENERGY_{RESIDUAL,FACE_TRACTION_OPERAND,BULK_FARFIELD_FLUX_OPERAND} — UNDECIDED (coverage).** UNPAIRED (j=0):
  (a) PY tags a single-value `SCENARIO=REAL_OMEGA_PROPAGATING_IMPERMEABLE_LAMBDA_X0_ZERO` axis WL lacks
  (bridgeable), AND (b) the two engines represent the energy DIFFERENTLY — PY a closed-form propagating-regime
  energy (`s11cc1_q_prop_output`, `re(...)`), WL an unevaluated far-field Poynting-flux **integral**
  (`Inactive[Integrate][Im[...]]`). The energy-balance audit is per-engine SOUND (both per-engine reviews
  confirmed it). Cross-engine reconciliation requires evaluating WL's far-field integral asymptotically —
  DEFERRED. No disagreement found; a genuine coverage gap.
- **DTN_FLAT_SYMBOL — AGREE (content).** WL emits ONE branch/face-universal flat symbol (+ radiation-branch
  detail); PY emits per-(branch,face). The flat half-space symbol IS branch/face-independent; WL's is the more
  economical keying. Content agrees.
- **DTN_GRAZING_BEHAVIOUR** — WL keys per-regime-pair, PY per-leg: a granularity/keying difference, not content.
- **Flat-resolvent leg-labeling (NEW — surfaced by the broad coefficient test) — UNDECIDED.** PY labels the flat
  half-space DtN symbol on the OUTPUT leg (`ρ_m ω/q_out`, `.py:602-603`), WL on the INPUT leg — equal on the
  physical diagonal `k=k′` (retest `rel=0`), differing off-diagonal in the MATERIAL_ADVECTED response coefficients
  (§3.3). Likely a diagonal-only convention (the half-space object is diagonal), but the off-diagonal non-collapse
  is a surfaced representational difference; a diagonal-restricted or the symbolic residual settles it.
- **REP_INVARIANCE_{EULERIAN,HANZAWA}_OPERAND / _RESIDUAL (j=50 each), CONTROL_* (j=26–156), DEGENERATE_LOCI_*,
  DIMENSIONS (j=6)** — derived-from-kernel and control families. Plausibly collapse via the kernel bridge (both
  engines apply the SAME symbolic transform — representation route, ablation, projection — to the proven kernel);
  2 exact-zero joins in REP_INVARIANCE_RESIDUAL, 1 in DEGENERATE_LOCI_REAL_ADMISSIBLE are direct evidence. ⚠ This
  is an inheritance ARGUMENT, ⛔ not a per-family computation — the full per-family symbolic residual is DEFERRED
  as the belt-and-suspenders confirmation; ⛔ do not read it as a measured AGREE for these families.
- **HOMOGENEITY_{BASE,CONTROL,RESIDUAL} — keying gap (j=0), ⛔ NOT AGREE-by-inheritance.** Fully unjoined (PY
  `OBJECT=HEIGHT_NORMAL_NORMAL` vs WL `OBJECT=HOMOGENEITY`) — a keying/coverage difference to adjudicate, not an
  inherited kernel collapse (corrected after Grok's leg).

---

## 5. Deferred to the ≥64 GB box (⛔ do NOT narrow to fit — recorded in `DEFERRED_HEAVY_RUNS.md`)

The 4 giant families — **PERMEABLE_PORT_HERMITIAN** (WL 62 MB), **UNIFORM_LIMIT_RESIDUAL** (PY 37 MB),
**UNIFORM_LIMIT_S11CC1_OPERAND** (PY 33 MB), **PERMEABLE_DISSIPATION_VS_OMEGA_TAU** (WL 12 MB) — plus the full
symbolic per-family residual (the belt-and-suspenders confirmation of every AGREE-by-inheritance family, and the
raw-operator DTN_OPERATOR / weighted-adjoint Hermitian / far-field-integral ENERGY reconciliations). These are the
only place `cS0` (R1 fold) activates a joined residual (WL `cS0` is dense in PERMEABLE_*). Consistent with
S11c-b's cross-engine-residual deferral. Mid-tier (FACE_RESPONSE + 3) result folded in §5a.

### 5a. Mid-tier result (4 families, EXIT=0, peak RSS 317 MB, runtime 1086.7 s)
```
ACCOUNTING FACE_RESPONSE            j=24  pyO=28  wlO=32   pf=0
ACCOUNTING CONTROL_BRANCH_RESIDUAL  j=104 pyO=224 wlO=224  pf=32
ACCOUNTING CONTROL_FORM_RESIDUAL    j=156 pyO=360 wlO=1104 pf=24
ACCOUNTING UNIFORM_LIMIT_S11B_OPERAND j=0 pyO=234 wlO=180  pf=0
```
- **FACE_RESPONSE — AGREE (via the coefficients).** Joins j=24 on the response leaves (PRESSURE, RELATIVE_FLUX,
  TRACTION). μ_θ is correctly PY-only-preserved as `OPAQUE_MU_THETA_OPERATOR`
  (`mu_theta_operands_preserved_without_registry_fold`). The full-response leaves are a **representation-depth**
  difference: PY keeps the DtN operator and the resolvent `[I+(Λ_A/ρ)Z]⁻¹` as OPAQUE composite symbols
  (`s11cc1_dtn_operator_*`, `s11cc1_response_resolvent_*`, `s11cc1_V_*` — operand ~1 KB), WL expands them fully
  (operand ~290 KB). The cross-engine-meaningful comparison that resolves this is `FACE_RESPONSE_COEFFS` (§3.3),
  which EXTRACTS the μ_θ- and face-velocity coefficients (forcing both engines to the expanded level) and AGREES
  exactly. So seal 1 (μ_θ) is confirmed at the coefficient level; the full-response opaque↔expanded residual is
  the same difference viewed unextracted.
- **CONTROL_BRANCH_RESIDUAL (j=104), CONTROL_FORM_RESIDUAL (j=156)** — the comparator's own ablation-control
  residual families; join and behave (per-family symbolic collapse DEFERRED as belt-and-suspenders). The parse
  failures (pf=32, pf=24 — the oversized-leaf / held-integral limits) are handled as `UNDEFINED_PARSE_FAILED`
  markers (no residual corruption), a coverage note, not a defect.
- **UNIFORM_LIMIT_S11B_OPERAND — keying (j=0).** The S11b-operator uniform-limit operand is keyed differently
  between engines (PY 234 vs WL 180 unpaired); a keying/coverage difference, adjudicated with the other
  UNIFORM_LIMIT_* families in the ≥64 GB deferral.

⇒ Full light+mid coverage: **46 of 50 families run on 30 GB; the 4 giants UNMEASURED — DEFERRED.** No cross-engine
DISAGREEMENT was found among anything measured. What is EARNED: the two-momentum DtN kernel (all cases, exact
zero, load-bearing) and the LAB_HELD + diagonal-MATERIAL response coefficients AGREE. What is NOT earned (stays
UNDECIDED/deferred, ⛔ not pre-adjudicated AGREE): seal 5 (density freeze), the raw DTN_OPERATOR form, ENERGY, the
off-diagonal flat-leg labeling, the derived/control families' full symbolic residual, and the 4 giant families.
⇒ **the kernel is cross-engine closed; the closure-wide verdict awaits the ≥64 GB symbolic run.**
