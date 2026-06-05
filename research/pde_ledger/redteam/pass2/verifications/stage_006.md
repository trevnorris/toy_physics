---
unit_id: 006
batch: I.1
verifier_model: claude-opus-4-8
verify_date: 2026-06-04T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 006

Note: stage 006's fix was a PAPER-CARD edit (`paper/stages/stage_006.tex`), not a
script edit. The captured script diff (`exec_logs/stage_006_diff.patch`) is EMPTY
by design (0 bytes), and there are no exec logs — both engines were already
correct and were intentionally not touched. This verification is about
card↔engine ALIGNMENT, not a script change.

## Per-finding outcomes

### F1 — paper_misalignment (projected Ampère sign, target_mismatch)

**Classification:** resolved

**What changed:**
`paper/stages/stage_006.tex:38-41` (eq:stage006-ampere). The git diff of the card
shows exactly the directed edit and nothing else:
```
-\nabla\times\mathbf H_{\rm flux}-\partial_t\mathbf D_{\rm flux}
-=\mu_0\mathbf J_{\rm proj}+\mathbf L_{\rm mix}.
+-\nabla\times\mathbf H_{\rm flux}-\partial_t\mathbf D_{\rm flux}
++\mathbf L_{\rm mix}=\mu_0\mathbf J_{\rm proj}.
```
The card now reads `−∇×H_flux − ∂_t D_flux + L_mix = μ₀ J_proj`: the curl-H term
sign was flipped to `−`, `L_mix` was moved to the LHS, and `μ₀ J_proj` is alone
on the RHS. No other line of the card was modified. The appendix stage-006 row
(`paper/appendices/stage_appendix_part01.tex`) is a status/summary line that does
not restate the Ampère sign and was correctly left unchanged (per directive
`## Applied: F1`).

**Assessment:**
The corrected card matches the load-bearing assertion in BOTH engines:

- SymPy (`...stage006...sympy_audit.py:146-148`): `amp1_target =
  (∂_z H2 − ∂_y H3) − ∂_t D1 + L1`. Since `(∇×H)_1 = ∂_y H3 − ∂_z H2`, the leading
  group `= −(∇×H)_1`, so the verified law is `−∇×H − ∂_t D + L = μ₀ J`
  (printed at lines 153-155; residues asserted zero at 150-159).
- Mathematica (`...stage006...mathematica_audit.wl:65-68, 119-124`):
  `ampereCurl3[v][i] = Σ eps3[k,j,i] ∂_j v[k]`. Because `eps3[k,j,i] = −eps3[i,j,k]`
  (1↔3 index swap, odd permutation), `ampereCurl3 = −curl3`, and
  `ampereRearranged` (119-122) asserts the tensor divergence equals
  `ampereCurl3[Hflux] − ∂_t Dflux + leak`, i.e. `−∇×H − ∂_t D + L = μ₀ J`.

Both engines independently verify the same sign-flipped form; the card now agrees
with both. The card's Gauss law `+∇·D = +μ₀ρ` is consistent with the SAME
component map: the SymPy Gauss block (`...sympy_audit.py:125,131-135`) sets
`G10=D1` and `lhs0 = ∂_x D1 + ∂_y D2 + ∂_z D3 + L0 = μ₀ρ`, i.e. `+∇·D + Leak0 =
+μ₀ρ`. The same map (`G^{i0}=D_i`, `G^{ij}=ε^{ijk}H_k`) yields BOTH the corrected
Ampère sign and the unchanged Gauss sign, so the card is now internally
consistent. The Faraday/divB block and the field-split definitions were not
touched and continue to match. No collateral edit; no script change.

## Exec log assessment

**SymPy:** exit=n/a. No exec log was captured (paper-only edit; scripts not run).
This is correct by design — the scripts were already correct and the directive
explicitly forbids running/modifying them.

**Mathematica:** exit=n/a. Same — no exec log; scripts untouched.

**Output freshness:** the committed `.txt` outputs are unchanged and remain newer
than their scripts (sympy output 2026-05-25 17:24 > script 02:13; mathematica
output 17:29 > script 02:13). Since no script changed, no regeneration was needed.
`git status` confirms only `paper/stages/stage_006.tex` is modified; both
`.py` and `.wl` are clean.

## Material-change assessment

`material_change`: false. No script and no derived numeric/symbolic engine result
changed — only the paper card's displayed Ampère equation was corrected to match
what the engines already verified. Nothing downstream consumes stage 006's
vector Ampère sign (per the directive's cross-stage investigation: 005/008 stay
covariant in `F^{μν}`; 010/011/012/023 are sign-agnostic P2/bundle moments; the
EM-extension endgame 243/244/247 rides the scalar continuity-leakage / Poynting
branch). No downstream re-stale is induced by this fix.

## Side observations (non-blocking)

The original report's D5 note (card writes `F^{i0}=E_i` upper-index vs script
`F_{i0}=E_i` lower-index) remains a pure notational up/down-index difference whose
observable Faraday/divB signs still match the card; it was correctly not raised
as a finding and is out of scope for this verification.

## Verdict justification

The sole finding F1 (paper-vs-engine Ampère sign mismatch) is resolved: the card's
eq:stage006-ampere now reads `−∇×H_flux − ∂_t D_flux + L_mix = μ₀ J_proj`, which
matches the load-bearing assertions of both independent engines and is consistent
with the card's own (unchanged) Gauss law under the same component map. The edit is
exactly the directed change with no collateral edits; both scripts are clean and
their committed outputs are unchanged. The empty diff patch and absent exec logs
are expected by design for a paper-only fix. No regressions; `material_change:
false`. Verdict: verified.
