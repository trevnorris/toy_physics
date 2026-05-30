---
unit_id: 150
batch: IV.5
verifier_model: claude-opus-4-8-1m
verify_date: 2026-05-29T22:55:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 150

## Per-finding outcomes

### F1 — insufficient_verification (DISPLAY corroboration of the compact slope)

**Classification:** resolved

**What changed:**
- SymPy `scripts/...stage150...sympy_audit.py:41-47`: the slope assignment was
  restructured to build the printed form from FREE coefficient placeholders, then
  substitute the concrete definitions:
  ```python
  Aq_s, Cq_s = sp.symbols("Aq Cq")
  Sq_symbolic = Aq_s*k - Cq_s*Pi
  Sq = Sq_symbolic.subs({Aq_s: Aq, Cq_s: Cq})
  print("S_q(Pi) =", Sq_symbolic)        # compact: Aq*k - Cq*Pi  (k=pi/2)
  print("S_q(Pi) [expanded] =")
  sp.pprint(Sq)
  ```
- Mathematica `mathematica/...stage150...mathematica_audit.wl:41-45`: mirrored with free
  placeholders `aqS,cqS`:
  ```mathematica
  sQsymbolic = aqS*k - cqS*p;
  sQ = sQsymbolic /. {aqS -> aq, cqS -> cq};
  Print["S_q(Pi) = ", fmt[sQsymbolic]];
  Print["S_q(Pi) [expanded] = ", fmt[sQ]];
  ```
The diff (`stage_150_diff.patch`) shows ONLY these two display/assignment blocks changed;
no collateral edits to BCs, the residual, or the curvature logic. Matches the directive's
consult-Q6 approach (b) exactly. Directive carries an `## Applied: F1` block
(deviation: none).

**Assessment:**
Correct and addresses the finding. The printed compact form is **provably the asserted
slope**, not a decoupled display string:
- The SAME symbolic object is reused — `Sq = Sq_symbolic.subs(...)` (py:43) /
  `sQ = sQsymbolic /. {...}` (wl:42) — so the printed `Sq_symbolic`/`sQsymbolic` and the
  asserted `Sq`/`sQ` differ only by substituting the concrete `Aq/Cq` (`aq/cq`)
  definitions. The fabricated-display anti-pattern (hardcoded text) is avoided.
- The directive's can-fail self-test holds: deleting the `.subs`/`/.` step would leave
  `Sq`/`sQ` as the free-symbol form, making `T_q'(0)-S_q == 0` FAIL — i.e. the printed
  form and the asserted object are linked, not decoupled.
- Refreshed output line 5: SymPy prints `S_q(Pi) = pi*Aq/2 - Cq*Pi` (k=pi/2 substituted,
  `Aq`/`Cq` preserved as symbols); Mathematica prints
  `S_q(Pi) = -(cqS*p) + (aqS*Pi)/2` (`aqS`/`cqS` preserved as symbols). Both are followed
  by an expanded-value line. Neither shows the bare fully-substituted rational at line 5.
  Both forms are explicitly listed as acceptable in the directive verification (b).
- Load-bearing assertion `T_q'(0)-S_q == 0` (py:52 / wl:50) is UNCHANGED and passes in
  both engines. Curvature checks R(0), R'(0) (py:56-57 / wl:53-54) and R''(0)-target
  (py:60-64 / wl:56-59) are untouched and pass. SymPy reports
  `R''(0) = 3·Π·Σ·exp(Π)/(1-exp(Π))` and Mathematica `(-3*E^p*p*sigmaM)/(-1+E^p)` —
  equivalent forms of `-3 Σ Π/(1-e^{-Π})`, both with `R''(0) - target = 0`.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
```
S_q(Pi) = pi*Aq/2 - Cq*Pi
T_q'(0)-S_q = 0
R''(0) - target = 0
```

**Mathematica:** exit=0. Notable lines:
```
S_q(Pi) = -(cqS*p) + (aqS*Pi)/2
T_q'(0)-S_q = 0   /  PASS: T_q'(0)-S_q
R''(0) - target = 0   /  PASS: R''(0) - target
```
All in-file checks pass; both logs report `# exit_code: 0`.

**Output freshness:** confirmed. `stat` shows both `.txt` outputs at mtime
2026-05-29 22:49:24 vs both scripts at 22:46:09 — outputs are ~3 min newer, regenerated
post-fix. Committed `.txt` contents match the exec logs on the relevant lines, including
the compact line-5 slope.

## Material-change assessment

`material_change`: false.

This is a display/corroboration-only change. The slope VALUE (`Aq*k - Cq*Pi`), the
load-bearing `T_q'(0)-S_q` assertion, the residual identity, and all three curvature
checks are mathematically identical to the pre-fix committed state (the source slope was
already correct at 3e2b5c0). No derived result that a downstream unit could depend on has
changed. No specific downstream re-audit concern.

## Side observations (non-blocking)

- The directive's own "Orchestrator note" records that
  `notes/stages/review/stage_150_review.md` is mislabeled (it contains Stage 031
  content). This is a prose/notes issue outside the scripts-only loop and is correctly
  excluded from this fix; flagged here only for orchestrator awareness, not blocking.
- The original auditor's documented (non-finding) observation stands: neither engine
  independently re-derives `T_s`/`T_q` from the lane ODEs — both plug in the notes'
  closed forms. This was intentionally not raised as a directive finding (needs upstream
  ODE/BCs the auditor could not read) and is unaffected by F1. Not blocking.

## Verdict justification

The single finding F1 is fully resolved. Both transcripts now display the compact
coefficient-symbol slope (`pi*Aq/2 - Cq*Pi` / `-(cqS*p) + (aqS*Pi)/2`) at line 5,
followed by the expanded value, and that printed form is provably the same symbolic
object entering the load-bearing `T_q'(0)-S_q == 0` assertion (built from free
placeholders, then `.subs`/`/.` to the concrete coefficients — confirmed by the
can-fail self-test). The load-bearing slope assertion and all curvature checks are
unchanged and still pass; both engines exit 0; saved outputs are fresh and match the
logs. No regressions in the diff and no fabricated-display anti-pattern. Verdict:
**verified**.
