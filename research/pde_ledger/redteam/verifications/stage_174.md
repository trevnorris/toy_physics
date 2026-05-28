---
unit_id: 174
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T16:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 174

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Codex made exactly four edits to
`mathematica/moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl`
(git diff: 7 insertions, 16 deletions, no other files touched):

- L26 — banner corrected `STAGE 157` → `STAGE 174`.
- L67-68 — the hand-typed four-term `b01Two` (`2*c1*dc1/w1^2 - 2*c1^2*dw1/w1^3 + ...`)
  replaced by `b0Eps = (c1 + eps*dc1)^2/(w1 + eps*dw1)^2 + (c2 + eps*dc2)^2/(w2 + eps*dw2)^2;`
  then `b01Two = FullSimplify[D[b0Eps, eps] /. eps -> 0, ...]`.
- L86-87 — hand-typed quotient-rule `z01Two` replaced by
  `z0Eps = (q1 + eps*q1p)/(delta1 + eps*delta1p) + (q2 + eps*q2p)/(delta2 + eps*delta2p);`
  then `z01Two = FullSimplify[D[z0Eps, eps] /. eps -> 0, ...]`.
- L105-106 — hand-typed four-term `n01Two` replaced by
  `n0Eps = (p1 + eps*p1p)^2/(delta1 + eps*delta1p)^2 + (p2 + eps*p2p)^2/(delta2 + eps*delta2p)^2;`
  then `n01Two = FullSimplify[D[n0Eps, eps] /. eps -> 0, ...]`.

**Assessment:**
The edits match the directive's prescribed text verbatim in all four spots. The three
load-bearing differentials are now reconstructed by genuine symbolic differentiation of
the *perturbed base quantities* rather than re-typed from the closed form:
`b0Eps`/`z0Eps`/`n0Eps` are `b0Two`/`z0Two`/`n0Two` with each primitive replaced by
`(base + eps*perturbation)`, and `D[..., eps] /. eps -> 0` reads off the first-order
coefficient. This is structurally independent of the SymPy script, which still carries the
hand-typed closed forms (`B01_two = 2*c1*dc1/w1**2 - ...` etc., sympy L78/L97/L116) — so a
transcription slip in either engine's differential would now diverge instead of agreeing,
which is exactly the cross-engine independence the finding demanded.

The check lines are unchanged in content: `expectZero["BdG weighted-average formula", b01Two/b0Two - deltaBWeighted]` (now L76), the Z analogue (L95), and the N analogue (L114).
Their line numbers shifted up because each hand-typed assignment collapsed from 4 lines to
2, but the asserted expression is byte-identical to the pre-fix version. The assertions
remain non-tautological: each compares the independently-derived `*01Two/*0Two`
log-slope against the weight-averaged candidate (`deltaBWeighted` etc.) built separately
from the `rho` weights.

`eps` hygiene is correct: it appears only inside the three `*0Eps` definitions and the
`D[..., eps] /. eps -> 0` calls (grep confirms 6 occurrences, all in L67/68/86/87/105/106),
is never assigned a value, is guaranteed fresh by the file-opening `ClearAll["Global`*"]`
(L1), and is not added to any `$Assumptions` list — matching both the directive instruction
and the `## Applied: F1` deviation note.

The exec transcript confirms the post-fix run: `BdG weighted-average formula = 0` / `PASS`,
`Z weighted-average formula = 0` / `PASS`, `N weighted-average formula = 0` / `PASS`, banner
reads `STAGE 174 — STATIC SELF-SIMILARITY DECOMPOSITION`, and the script reaches its final
`Stage 174 Mathematica audit passed.` line (the `fail` helper would have `Exit[1]` before
that). No collateral edits.

## Exec log assessment

**SymPy:** exit=0. The sympy script was not in scope for F1 (mathematica-only finding) and
was not modified. Its log is complete and error-free, e.g.:
- `weight identity = 0`
- `BdG weighted-average formula = 0`
- `Xi_load = (N0*(-B01 + K1 - Z01) + N01*(B0 - K + Z0))/(N0*(B0 - K + Z0))`
The log shows the older "STAGE 157" banner and no `PASS` lines because the `.py` predates
the PASS-printing convention; this is unrelated to F1 and the directive did not touch it.
No traceback; output ends cleanly with the carry-forward block.

**Mathematica:** exit=0. Notable lines:
- `STAGE 174 — STATIC SELF-SIMILARITY DECOMPOSITION`
- `BdG weighted-average formula = 0` followed by `PASS: BdG weighted-average formula`
- `Z weighted-average formula = 0` / `PASS`; `N weighted-average formula = 0` / `PASS`
- `Stage 174 Mathematica audit passed.`
All eight checks print `= 0` and `PASS`; the closing line confirms `Exit[0]` was reached.

**Output freshness:** confirmed. Mathematica script mtime `2026-05-28 15:55:34`; its output
log mtime `2026-05-28 16:13:48` (newer → regenerated post-fix). SymPy output mtime
`2026-05-28 16:10:21`, newer than the unchanged sympy script (`2026-05-11`); both logs are
fresh from the orchestrator's refresh.

## Material-change assessment

`material_change`: false.

The fix changes only *how* the Mathematica engine arrives at `b01Two`/`z01Two`/`n01Two` —
the derived values are provably identical to the prior hand-typed forms (the `expectZero`
residuals remain exactly 0), and no assertion, constant, or carry-forward result changed.
No downstream unit depends on any altered numeric or symbolic output. The auditor likewise
recorded `stop_cold: null` and noted the fix "changes no assertion or constant."

## Side observations (non-blocking)

- The SymPy script (`...sympy_audit.py`) still hand-types the three differentials. That is
  by design here — the cross-engine independence is now provided by the Mathematica side
  deriving them via `D[..., eps]`, which is exactly what the finding required. Not a defect.
- The SymPy transcript uses a different banner/print convention (no `PASS` lines, residual
  values printed without an explicit pass marker). Cosmetic only; out of scope for F1.

## Verdict justification

The single finding is fully resolved. Codex applied exactly the four edits the directive
prescribed, with no collateral changes; the three load-bearing differentials are now
independently reconstructed by symbolic perturbation differentiation in Mathematica while
the SymPy engine retains the hand-typed closed forms, restoring genuine cross-engine
confirmation; the banner is corrected; the three `expectZero` checks are unchanged and
still pass with zero residual per the freshly regenerated transcript; both engines exit 0;
and the change is non-material (no derived result altered). Verdict: verified.
