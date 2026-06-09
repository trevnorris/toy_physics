---
unit_id: 181
batch: V.2
verifier_model: Claude Opus 4.8 (1M context)
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 181

## Per-finding outcomes

### F1 — mathematica_transliteration (USER-AUTHORIZED full re-author)

**Classification:** resolved

**What changed:**
The `.wl` was fully re-authored (diff: 111→219 lines).
- The hand-typed literals `eps1Expected` (old wl:75-78) and `theta1Expected`
  (old wl:114-118) are GONE. Grep over the file finds no `*Expected` token at all.
- The `.py`-mirroring `s`-perturbation route `D[Log[…],s]/.s->0` (old wl:164-186)
  is GONE — no `Log[…s]` / `/. s -> 0` constructs remain.
- New derivation infrastructure: `directional[expr,pairs]` (wl:34-36) and
  `logDirectional[expr,pairs]` (wl:38) — first-order directional differentials in
  *named drift variables* (`overlapDrift, frequencyDrift, …`), a mechanism the
  SymPy script does not use (SymPy uses `D[eps,epsW]*epsW1` + `D[Log[…],s]/.s->0`).

**Assessment — INDEPENDENCE confirmed (genuine derivation, not a relabel):**
- **ε₁** is derived *two* independent ways and cross-checked:
  `splitDriftFromFullExpression = directional[splitBlock, splitVariableDrifts]`
  (wl:145, full directional derivative) vs `splitDriftFromProductLedger` =
  hand product-rule expansion `bareBlockDrift*splitMultiplier +
  wallBlock*directional[splitMultiplier,{{uSplit,uSplitDrift}}]` (wl:146-149),
  compared at wl:151-153. `epsilonOne` (wl:156) is the derived object, no literal.
- **Θ₁** is derived two ways: `thetaFromQuotient = logDirectional[trackingFactorDirect,…]`
  (wl:213) vs the factor-by-factor log-ledger `logDirectional[1+joint+uSplit]
  - logDirectional[1+joint] - logDirectional[1+uSplit]` (wl:214-218), compared at
  wl:220-223. `thetaOne` (wl:225) is derived, no `theta1Expected` literal.
- **Ξ₁** is derived as `xiOneFromShape = logDirectional[shapeFromDirectPort,
  placementDrifts]` (wl:159), cross-checked against the law (wl:165-168) AND against
  the support-loaded route `xiOneFromLoadedRoute` (wl:185-189).
This is a real second route (directional log-derivatives + product/factor ledgers),
not a re-typing with a different API. A transcription error in a single literal can
no longer pass both engines because the `.wl` has no shared literal target.

### F2 — tautological_check (both engines)

**Classification:** resolved

**What changed:**
The round-trip "support-loaded branch product law" check
(`R_target_loaded*Mtr - product_loaded`) is removed from BOTH engines.
- `.py`: diff shows exactly one removed line (old py:57); current file lines 54-58
  contain no such check.
- `.wl`: folded into the re-author; the old `expectZero["support-loaded branch
  product law", rTargetLoaded*loadMassFromSupport - productLoaded]` is absent.
Neither output transcript contains "branch product law".

**Assessment:** Correct. No other check was lost — the genuine support-blindness
content (`d/dζ` checks, reconstruction identities, spoiled-packet counter-check) is
intact, and the `.wl` re-author actually ADDED cross-checks (quotient-form,
loaded-Ξ₁). The `.py` is otherwise byte-identical (single-line deletion only).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `direct-selected transfer-shape identity = 0` (D1)
- `d/dzeta ln T^2 ... = 0`, `d/dzeta ln R_target ... = 0` (D2 blindness)
- `spoiled d/dzeta ln R_target = (…nonzero rational…)` — negative control fires
- `split-blocking drift eps_1 = 0`, `tracking-factor drift = 0`,
  `Xi_1 derived from T^2 matches defect law = 0`, `selected-branch identity = 0`

**Mathematica:** exit=0. Notable lines:
- `PASS: direct-selected transfer-shape identity`
- `PASS: d/dzeta ln T^2 (support-loaded route)`, `PASS: d/dzeta ln R_target …`
- `spoiled d/dzeta ln R_target = (…nonzero rational…)` then `PASS: spoiled …`
  (via `expectNonZero`, which FAILS the script if it were 0 — fail-on-purpose
  control survived)
- `PASS: split-blocking drift eps_1`, `PASS: tracking-factor drift`,
  `PASS: Xi_1 derived from transfer differential matches defect law`,
  `PASS: selected-branch identity`, `Stage 181 Mathematica audit passed.`

**Output freshness:** confirmed. Both `.txt` mtimes = 1780986153, newer than the
`.py` (1780979711) and `.wl` (1780979761) sources. Regenerated post-fix.

## Material-change assessment

`material_change`: false. The re-author and the F2 deletion change neither a derived
constant, sign, nor any emitted closed form. The printed Θ₁, Ξ₁, R₁, ε₁ forms match
the auditor's pre-fix values up to the `.wl`'s renamed symbols (e.g.
`uSplit↔deltaU`, `joint↔chi0`); the SymPy forms are unchanged. Downstream Stage 182
consumes the same closed forms. No unit > 181 is substantively affected.

## Side observations (non-blocking)

- The `.wl` deliberately renamed all physics symbols (e.g. `chi0→joint`,
  `deltaU→uSplit`, `zeta→supportShare`) — this is appropriate for an independent
  re-author and helps demonstrate the engines are not sharing a namespace.
- The `.wl` adds an extra `xiSupportRigid === 0` guard (wl:236) mirroring the
  SymPy support-rigid counter-check (py:142) — strengthens, not weakens, coverage.

## Verdict justification

Both findings are `resolved`. F1's 173-lesson independence requirement is met: the
hand literals and the `D[Log[…],s]/.s->0` py-mirror are gone, and ε₁/Θ₁/Ξ₁ are now
derived by genuine, can-fail directional/log-derivative routes cross-checked against
independent product- and factor-ledger expansions. F2's tautological round-trip is
removed from both engines with no collateral loss; the `.py` changed by exactly that
one line. All deliverables (D1, D2 blindness + surviving spoiled-packet negative
control, D3 ε₁/Ξ₁/selected-identity/Θ₁) verify with residuals → 0, both engines exit
0, outputs are fresh. Verdict: `verified`.
