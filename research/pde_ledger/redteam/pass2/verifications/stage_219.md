---
unit_id: 219
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T19:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 219

## Per-finding outcomes

### F1 — paper_misalignment (paper_missing_script_claim; engine-coverage statement, not a value)

**Classification:** resolved (by deferral — outside scripts-only scope; no script change required or possible)

**What changed:**
Nothing in any script. The directive (`redteam/pass2/directives/stage_219.md`) correctly
applies NOTHING: the only finding is a stale `\stagefield{Verification}` annotation in the
published card (`paper/stages/stage_219.tex:11` reads "Mathematica audit: none yet") and the
parallel notes §9, while a passing dual-engine `.wl`
(`mathematica/.../stage219..._mathematica_audit.wl`, commit `1dfc3fe`) now exists and passes
M1–M7. The git diff patch `redteam/pass2/exec_logs/stage_219_diff.patch` is legitimately empty
(0 content lines) — this is the FIRST pass-2 batch needing zero script corrections, not a
"Codex failed to apply" condition. Both `.py` and `.wl` are byte-unchanged from HEAD.

**Assessment:**
This is a documentation-staleness in the direction of *understating* coverage (the card claims
less verification than exists), so it creates no false verification claim. It is a paper-side
edit — the user's call, not Codex's, and out of the scripts-only verifier scope. Per the
standing user decision it is deferred to PAPER_CLEANUP P4-51 (the VII.1 card `Verification`-field
doc-sync lag batched across the VII.1 cards). The deferred paper finding alone does NOT roll the
overall verdict to needs_rework. No script defect, no regression, no script edit prescribed.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L43: "Verified: chi_qW = Lambda / D0 and chi_qW^2 = P0 / D0"
- L59: "All Stage 219 symbolic checks passed."
- Asserts at py L46/59/85/110/130-131/164/199-202 all compare an engine-computed object
  (`K_red.det()`, Schur via `K_int.inv()`, `K_red.inv()` entries, `J^T·Kinv·J` product,
  the C6/C4/C2-assembled primitive kernel) against an independently-written paper closed form.

**Mathematica:** exit=0. Notable lines:
- L8/10/12/14/16/21/26: "PASS M1" … "PASS M7" (all residuals `{0,...}`).
- L19: "M6 x-degrees by {1, Exp[2 kappa x], Exp[4 kappa x]} = {4, 2, 0}"; L20: coefficient-list
  lengths `{5,3,1}`; L25: "M7 PositiveDefiniteMatrixQ = True"; L27: "Stage 219 Mathematica audit passed."

**Output freshness:** the orchestrator re-ran both engines directly; both logs are freshly
timestamped (2026-06-09T19:21:18-06:00), exit 0, deterministic (219's committed `.txt` is
byte-identical to the fresh run). Committed-`.txt` mtime is not a failure condition here.

## Material-change assessment

`material_change`: false. No script was edited; no derived result changed. No downstream unit
can be affected. (No `upstream_stale` propagation is warranted on a zero-edit unit.)

## Independence / non-tautology spot-check (scripts-only)

The `.wl` is a genuine independent recomputation, not a transliteration of the `.py`:
- M1–M5 each engine computes its own `Inverse[Kred]`/`Det[Kred]` and compares to the paper's
  claimed closed forms (RHS from the card, not from copying the other engine's algebra).
- M6 (load-bearing square-law suppression theorem) is *more* independent than the `.py`: the
  `.py` (L145–164) compares the matrix-product `delta_V_primitive` against a hand-assembled
  C6/C4/C2 closed form (two different construction paths → non-tautological). The `.wl`
  (L86–102) does NOT compare to a pre-baked answer; it scales `-2·dVp·x^6·Exp[4κx]`, substitutes
  `Exp[2κx]→z`, and *structurally proves* the family content via
  `PolynomialQ`/`Exponent`/`Coefficient`/`CoefficientList`: `Exponent[familyZ,z]==2`,
  `Exponent[familyZ,x]<=4`, `zDegrees==={4,2,0}`, `zeroQ[x5Coefficients]`, and all three
  z-coefficients nonzero. This degree-exclusion attack genuinely forbids the long-range
  `-1/x` (→ `x^5 z^2`) and `e^{-2κx}/x` (→ `x^5 z`) families, so the Output's suppression
  clause is actually exercised — a structurally distinct route the `.py` lacks.
- M7 adds `PositiveDefiniteMatrixQ[Knum]` plus the numeric-slice residuals (`Delta=140`,
  `D0=74/7`, `det=1480`, minors `{11,98,1480}`) — an independent positive-definiteness proof.

No tautological load-bearing assert found; no transliteration the auditor missed.

## Side observations (non-blocking)

The notes §9 "Script-backed status" sharing the same understated coverage is part of the same
deferred paper-side doc-sync (P4-51); it is outside scripts-only scope and raises no script
concern. The residual stale stage-number labels from pass 1 were not re-encountered as new
blockers here and remain the known deferred numbering-cleanup class.

## Verdict justification

Both engines pass (sympy_exit=0, mathematica_exit=0) on byte-unchanged scripts; the diff patch
is legitimately empty (zero-correction batch). Every load-bearing assertion compares an
engine-computed object against an independently-written paper closed form — none is
true-by-construction — and the `.wl`'s M6 structural-family/degree-exclusion proof plus M7
`PositiveDefiniteMatrixQ` make it a genuine independent recomputation, not a transliteration.
The single finding F1 is a stale card `Verification` annotation (paper-side, understating
coverage), correctly applying nothing and deferred to PAPER_CLEANUP P4-51 — outside the
scripts-only scope and explicitly non-blocking. Overall verdict: `verified`.
