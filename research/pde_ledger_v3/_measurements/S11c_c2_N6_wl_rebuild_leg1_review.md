# Independent physics review — blind Wolfram N6 engine (S11c-c2) — REPORT

Artifact: `research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl` (working tree, 1029 lines,
Codex-astra). Reviewer: fresh independent leg. All kernel runs wrapped in `timeout 600`, ONE kernel at a
time. Ablations run on `/tmp` COPIES only; the working tree was never modified.

**VERDICT: BUILD CLEAR.** Every load-bearing object is reached by computation (no hand-typed/ tautological
payload); the affine split is an exact identity; all three amplitude/geometry controls bite one-sidedly; the
prediction is non-circular; the PIT is a genuine WL-derived finite-field probe with the off-diagonal
(k_out=k_in) blind spot explicitly excluded; the engine is import-free and byte-identically blind. Two
non-blocking parity observations are listed at the end.

---

## 1 · Independent derivation (done BEFORE reading the engine)

Script: `/tmp/n6_wl_review/independent_derivation.py`  · stdout: `/tmp/n6_wl_review/independent_derivation.stdout.txt`

I did not reconstruct the full nonlocal DtN/material fold (intractable by hand). From the physics named in
`S11c_c2_SHARED_PHYSICS.md` §5c + `_measurements/S11c_c2_N6_route2_spec_astra.md`, I derived the STRUCTURAL
identities the two checks must obey and the SIGN of each control response, in an independent sympy toy:

- **Affine split is an exact algebraic identity.** With `I(C,S)=closed(C,S)+bare(C)` (closed bilinear, bare
  linear in C only): `E−M = I(dc,ms) + closed(C_M,ds) + closed(dc,ds)` ⇒ **SPLIT_CHECK ≡ 0**. The carrier
  channel MUST use `ms` (not `es`) and MUST keep the bare `−C·p`; source/cross MUST be closed-response-only.
  My toy confirmed: using `es` in the carrier leaves a residual `−closed(dc,ds)` (SPLIT_CHECK ≠ 0); re-adding
  bare to source/cross leaves `−(bare(C_M)+bare(dc))`.
- **Carrier reconcile `C_E−C_M`** is a difference of two INDEPENDENT routes; no prescribed zero. A FORM change
  to the material covector/normal feeding `C_M` MUST move CARRIER_BRIDGE_RESIDUAL/CARRIER_CHANNEL while
  `C_E`, `es`, SOURCE_BRIDGE_RESIDUAL stay fixed.
- **Source naturality `R_cov = ms − ms_pred`** (commuting square), `ms_pred = source(μ_E∘Φ, V_E)` built ONLY
  from μ_E and Φ. Predicted responses I derived: advection-coeff knife (a_ρ 1→2) moves R_cov by `A·a_ρ+B·Da_ρ`;
  junk knife moves R_cov by `J_μ·e_W`; a circular prediction collapses R_cov to structural 0; the shipped
  ms_pred must be INVARIANT under any material-route corruption.
- **PIT/δ**: δ = min(1,D/(N−E)) with D,E read off the engine's OWN circuit; k_out≠k_in independent legs; joint
  singular rejection; real branch before modular reduction; no residual-zero exit.

Every prediction above was then confirmed against the engine (below).

---

## 2 · Ablations — each on a /tmp COPY, timeout 600, one kernel at a time

Generator (mechanical single-case + reduced-PIT edits): `/tmp/n6_wl_review/gen_copy.py`.
Comparator: `/tmp/n6_wl_review/cmp.py`. Copies: `A_baseline.wl A_carrierknife.wl A_phiknife.wl B_baseline.wl
B_junkknife.wl` (all `/tmp/n6_wl_review/`). Stdout: same names `.out`.

**Speed note (does NOT weaken any conclusion):** copies restrict to ONE case and reduce the PIT to 1 branch
cell / 1 prime / 8 draws (for review speed only — the physics construction is untouched). The engine emits
each object's full symbolic `ARITHMETIC` circuit, which is **cell- and prime-independent**; I compared that
circuit directly, so "MOVED / identical" verdicts hold across ALL cells/primes, not just the sampled one
(explicitly re-verified for the carrier knife — see §2.1). Same-case+same-seed also makes PROBE_NUMERATORS
directly comparable. Baseline A ran 172 s at 111 MB RSS, exit 0, all 46 tags emitted; no run hit the 600 s cap.

### 2.1 Carrier FORM knife — `materialNormalKnife: 0→1` (`invT[[1,4]] += knife·jet["WBg",{2}]`, line 289), case LAB_HELD/RHOBR
Literal `cmp.py A_baseline.out A_carrierknife.out` (full-line sha + PROBE_NUMERATORS):
```
RC_CARRIER_MATERIAL          =MOVED      PN_MOVED
RC_CARRIER_BRIDGE_RESIDUAL   =MOVED      PN_MOVED
RC_CARRIER_CHANNEL           =MOVED      PN_MOVED
RC_MATERIAL_OPERAND          =MOVED      PN_MOVED
RC_R_N6                      =MOVED      PN_MOVED
RC_CARRIER_EULERIAN          =identical  PN_identical
RC_SOURCE_BRIDGE_RESIDUAL    =identical  PN_identical
RC_EULERIAN_OPERAND          =identical  PN_identical
RC_SOURCE_EULERIAN           =identical  PN_identical
COV_R_COV / SOURCE_PREDICTED / SOURCE_ACTUAL = identical
```
Symbolic `ARITHMETIC`-only comparison (cell-independent): CARRIER_MATERIAL/BRIDGE_RESIDUAL/CHANNEL/R_N6 MOVED;
CARRIER_EULERIAN, SOURCE_BRIDGE_RESIDUAL, EULERIAN_OPERAND identical. Decisive extra finding: at baseline the
`ARITHMETIC` hash of CARRIER_MATERIAL EQUALS that of CARRIER_EULERIAN (`390d2117…`) — the two independent
routes genuinely COINCIDE at this case (that is why CARRIER_BRIDGE_RESIDUAL=0 in the baseline), and the knife
separates them (`e1f8c04a…` vs unchanged `390d2117…`). ⇒ `C_M` is a genuinely separate construction (material
covector via `geometriesM`), NOT cloned from `C_E`; the reconcile is able-to-fail; the corruption is strictly
one-sided (source velocity uses the knife-free `geometriesSourceM`, knife=0, so `ms`/COV are untouched).

### 2.2 Φ-coefficient (advection) knife — `actualAdvectionCoefficient: 1→2`, case LAB_HELD/RHOBR (a_ρ live)
Literal `cmp.py A_baseline.out A_phiknife.out`:
```
COV_R_COV                =MOVED     PN_MOVED
COV_SOURCE_ACTUAL        =MOVED     PN_MOVED
COV_R_COV_CONTROL_DELTA  =MOVED     PN_MOVED
COV_R_COV_INCREMENT      =MOVED     PN_MOVED
COV_SOURCE_CONTROL_DELTA =MOVED     PN_MOVED
COV_SOURCE_PREDICTED     =identical PN_identical
COV_SOURCE_BASELINE      =identical PN_identical
RC_CARRIER_EULERIAN / RC_SOURCE_EULERIAN = identical
```
R_cov moves even though the prediction's declared Φ (a_ρ+h_α truth table) is unchanged (SOURCE_PREDICTED
byte-identical). Baseline R_cov=0 (the commuting square commutes at the true coefficient — the desired
physics), and the knife makes it nonzero ⇒ the naturality check is able-to-fail. **SOURCE_PREDICTED identical
under a material-route corruption is direct proof of NON-CIRCULARITY** (ms_pred does not depend on μ_M).

### 2.3 θ-junk knife — `actualJunkCoefficient: 0→1` (`+ jk·junkMu·eW` into μ_M only, line 247), case MATERIAL_ADVECTED/RHO4
Literal `cmp.py B_baseline.out B_junkknife.out`:
```
COV_R_COV                =MOVED     PN_MOVED
COV_SOURCE_ACTUAL        =MOVED     PN_MOVED
COV_R_COV_INCREMENT      =MOVED     PN_MOVED
COV_R_COV_CONTROL_DELTA  =MOVED     PN_MOVED
COV_SOURCE_CONTROL_DELTA =MOVED     PN_MOVED
COV_SOURCE_PREDICTED     =identical PN_identical
COV_SOURCE_BASELINE      =identical PN_identical
RC_CARRIER_BRIDGE_RESIDUAL / RC_CARRIER_EULERIAN = identical
RC_R_N6                  =MOVED     PN_MOVED
```
R_cov moves (junk in μ_M reaches the source), prediction fixed. Coherent channel structure confirmed:
CARRIER_BRIDGE_RESIDUAL is IDENTICAL under a μ corruption because μ enters only the closed pressure SOURCE,
never the pressure-slot CARRIER coefficients (route-2 spec §4: "open pressure coefficients carry neither θ nor
μ"). So the three knives cleanly separate channels: the FORM/geometry knife moves the carrier; the two
μ-amplitude knives move the source. Also confirms computed structural absence for RHO4_CONSTANT
(`ADVECTION_ABSENCE`: a_ρ=0, DENSITY_GRADIENT={0,0,0}, MATERIAL_MU_TAG_DERIVATIVE=0) — emitted as the actual
computed 0, NOT an A−A control.

### 2.4 Non-circularity (ablation 4) — established without a separate run
Code (line 853): `muPred = muE /. predictedMap["MAP"]`, and (871) `predicted = sourceBind[sourceSolve, muPred,
velocitiesE[s], density3]` — ms_pred is built ONLY from μ_E and Φ and Eulerian velocity; never from μ_M/ms/the
pullback. This is proven operationally by §2.1–2.3: SOURCE_PREDICTED is byte-identical under EVERY material-route
corruption (carrier normal, advection coeff, junk). Had ms_pred been circular (built from the material pullback),
those corruptions would have moved it. They do not ⇒ genuinely independent. (A circular-prediction copy would
merely reconfirm this; the invariance evidence is stronger and already in hand.)

### 2.5 PIT probe (ablation 5) — WL-derived δ and off-diagonal blind spot excluded
From `A_baseline.out` LOCAL_PROBE (single reduced cell/prime):
```
"D" -> 45202, "E" -> 60342, "N" -> 1000000008, "PER_VALID_DRAW" -> 22601/499969833
"EXCLUDED_LOCI" -> {Equal[qOut,0], Equal[clearedDenominator,0], Equal[kOut,kIn]}
"ALL_ZERO_TABLE_SEMANTICS" -> ConditionalExpression[noNonzeroFound, conditionalFalseNegativeBound]
```
`PER_VALID_DRAW = 22601/499969833 = min(1, D/(N−E))` with D=45202 (numerator degree from `nodeDegree` over the
circuit) and E=60342 (excluded denominator degree from leaf/chart/qLeg circuits) — large non-round
circuit-derived integers, so δ is genuinely **WL-derived, not a transferred constant**. The k_out=k_in locus is
explicitly EXCLUDED and the sampler rejects any draw with `kOut===kIn` (line 783), so an off-diagonal-only
discrepancy CANNOT be hidden by diagonal collapse — the engine is structurally immune to the ablation-5 failure
mode. Both legs are sampled from independent chart coordinates per leg; singular rejection is JOINT across all
compared operands (lines 779–787); the real branch cell is chosen (rational stereographic/hyperboloid charts)
BEFORE modular reduction; bad-prime handled by worst-prime bound (Min[primeList]) + a conditional
`BAD_PRIME_CONDITION = Exists goodPrime …`.

---

## 3 · Structural / faithfulness checks (WHICH LINE computes each object)

- **Affine split exact.** `SPLIT_CHECK` = 288 entries, ALL zero (SparseArray numerators empty). First entry's
  ARITHMETIC is literally `A + (−1)·A`. `R_N6` is live (18 stored-nonzero). Carrier channel uses `ms`
  (line 891, affine=True → keeps bare `−dc·p`); source/cross use es−ms with affine=False → closed-response-only
  (lines 892–893); `R_COV_INCREMENT` = closed-response contraction, affine=False (line 896). Matches derivation.
- **Load-bearing objects all computed** (no hand-typed answer payload): energy `constructEnergy[]` (211, tensor
  contractions + divergence-quotient RowReduce, symbolic coeffs bRho/cCoupling/energyCoefficientN); μ_E=`el[…]`
  (843); μ_M=`materialAmplitude…["MU"]` pull-back-then-vary (846); μ_pred=`muE/.Φ` (853); C_E=`carrier[rowsE]`
  (868, EULERIAN_SLAB_ROWS); C_M=`carrier[rowsM]` (868, MATERIAL_FACE_FOLD); es/ms=`sourceBind` (869–870);
  channels=`buildContraction` (888–896); c1 kernel/response re-derived by `constructKernel`/`responseFamilies`
  (388/508) — genuine blind re-derivation, not an import. PROVENANCE tags the six INDEPENDENT builders
  (MU_E=EL, MU_M=PULLBACK_EL, V_E=EULERIAN_LEVEL_SET, V_M=MATERIAL_FLATTENING, C_E=EULERIAN_SLAB_ROWS,
  C_M=MATERIAL_FACE_FOLD).
- **Guards genuine, not tautological.** SLOT_GUARD_RESIDUAL first component = projectedRow − Σ_p C_p·p
  (slot-linearity, able-to-fail on any quadratic-in-slot term) + G_cross = D[row,p,q] + denominator pressure
  derivatives; CLOSURE guards compare closed-image vs carrier reconstruction. Baseline: SLOT_GUARD_RESIDUAL and
  CLOSURE_GUARD_RESIDUAL = ALL zero (rows are exactly slot-linear, reconstruction consistent). Emitted as
  measurements, never asserted.
- **Dimensional able-to-fail fired**: `CONTROL_ADDITION_CONSISTENCY -> Equal[2,1]` from the deliberate extra-W0
  pressure summand (deltaPPlus {−2,−2,1} + W0·deltaPMinus {−1,−2,1} → 2 distinct dims) — the incompatible-sum
  witness is genuinely produced.
- **No physics disposition / no VERDICT.** No verdict/assert/pass/fail tokens; all `Quit[…]` codes are
  structural/domain guards (bad name, dup emission, uncovered Φ domain, unsupported head, non-prime, sampler
  exhausted, pivot failure) — none is a residual-zero exit. Sampler loop condition is validity/count only
  (line 774). Operands are emitted before residuals (rcNames/covNames order, line 1017–1022). ALL_ZERO table =
  "no nonzero found", not a pass.

## 4 · Operational / blindness

- **Import-free**: grep for `Get|Import|Needs|ReadString|ReadList|OpenRead|<<|URLFetch|CloudGet|…` → NONE. The
  only path strings are in the header COMMENT (construction authorities), not live reads.
- **Isolated-vs-in-repo byte-identical**: same file `A_baseline.wl` run with CWD=/tmp vs CWD=repo-mathematica-dir
  → sha256 `38a4db20d958f59bd2cc4d58a75a95f4b5bfb6b99991dbff0c88e63a0e5db568` IDENTICAL. Combined with
  import-freedom ⇒ genuine blindness (output independent of location; nothing read from the repo).
- **No freeze of varying fields (M3/rule-17):** WBg kept live (`profileRules` scale jets to σ_W AFTER `td`
  forms live derivatives); `pruneBackgroundProducts` drops only σ_W^{≥2} products while keeping zero-jet WBg
  live; the two density representatives are the field-vs-constant pair; a_ρ uses the live `td[density4,i]`.

## 5 · Non-blocking observations (parity, not defects)

- **O1** — WL `phiMap` (236) records `MAX_RANK` but, unlike the SymPy sibling's `prolonged_phi` (which raises
  if `max_rank+1 > background_depth`), has no explicit rank-vs-depth halt. Harmless: `tdMany` prolongs to the
  FULL jet rank present and `profileRules` covers WBg/muRBg jets of any order, so nothing is silently
  truncated/frozen; the shipped μ_E is rank-2 within depth-3. A defensive-assert parity note only.
- **O2** — The shipped carrier FORM knife is a single off-diagonal covector contamination. It bites decisively
  (proves C_M independence + able-to-fail), so it is sufficient; it is less aggressive than a combined
  sign+off-diagonal+symbol-collapse knife. Observation on knife strength, not a defect.

## Artifacts (all absolute)
- Independent derivation: `/tmp/n6_wl_review/independent_derivation.py`, `…/independent_derivation.stdout.txt`
- Generator / comparator: `/tmp/n6_wl_review/gen_copy.py`, `/tmp/n6_wl_review/cmp.py`
- Ablation copies + stdout: `/tmp/n6_wl_review/{A_baseline,A_carrierknife,A_phiknife,B_baseline,B_junkknife}.{wl,out}`
- Blindness run: `/tmp/n6_wl_review/A_baseline_inrepo.out`

**BUILD CLEAR** — no blocking defects; O1/O2 are non-blocking parity notes.
