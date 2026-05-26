---
unit_id: 049
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 049

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
In `mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:52-53` the two-`Integrate` definitions were replaced with one integrator call (LHS, with explicit `Assumptions -> Element[n, Integers] && n >= 0 && l > 0`) and one bare independent closed-form expression `Sqrt[2 l]/((n + 1/2) Pi)` (RHS). The diff capture at `redteam/exec_logs/stage_049_diff.patch` shows exactly the before/after pair the directive specified — no other lines touched. The `expectZero["uniform overlap integral", overlapFromIntegral - overlapFormula]` call at line 56 is unchanged but now substantive.

**Assessment:**
Edit matches the directive byte-for-byte; no collateral changes. The new `overlapFormula` is an independently stated algebraic target, not a second integrator output, so the assertion can now fail if the integrator ever returns a non-equivalent form. No re-introduction of the `uniformDnOverlap` helper, so the structural mirror with SymPy stays broken (the v1 F2 intent is preserved). The saved Mathematica transcript prints `I_n from direct integral = (2*Sqrt[2]*Sqrt[l])/(Pi + 2*n*Pi)` and `I_n closed form = (Sqrt[2]*Sqrt[l])/((1/2 + n)*Pi)` — two algebraically equivalent but textually distinct forms — confirming the LHS/RHS are no longer identical expressions. `uniform overlap integral = 0` with `PASS: uniform overlap integral` still prints.

## Exec log assessment

**SymPy:** exit=n/a. No `redteam/exec_logs/stage_049_sympy.log` was captured for this iteration. Falling back to the canonical transcript `scripts/output/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.txt` (mtime 2026-05-26 02:55, newer than the SymPy script mtime 2026-05-22 16:55). All seven assertions print `= 0` and the carry-forward block emits cleanly, so the run completed normally.

Notable lines:
- `k_n satisfies D/N Neumann boundary = 0`
- `uniform overlap integral = 0`
- `exact twin-lane support ratio = 0`
- `lowest twin half-wave = 0`

**Mathematica:** exit=n/a. No `redteam/exec_logs/stage_049_mathematica.log` was captured. Falling back to `mathematica/output/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.txt` (mtime 2026-05-26 02:55, newer than the script mtime 2026-05-26 02:54 — i.e. regenerated after Codex's edit).

Notable lines:
- `I_n from direct integral = (2*Sqrt[2]*Sqrt[l])/(Pi + 2*n*Pi)`
- `I_n closed form = (Sqrt[2]*Sqrt[l])/((1/2 + n)*Pi)`
- `uniform overlap integral = 0` / `PASS: uniform overlap integral`
- `Stage 32 Mathematica audit passed.`

The two printed `I_n` forms differ textually (integrator output vs. independently stated closed form), confirming the RHS is not a copy of the LHS. Their FullSimplify difference is zero, which is the genuine substantive check the directive asked to restore.

**Output freshness:** confirmed. The Mathematica `.txt` (02:55) post-dates the `.wl` (02:54), and the SymPy `.txt` (02:55) post-dates its `.py` (2026-05-22). Both transcripts are post-fix.

## Material-change assessment

`material_change`: false.

The edit only restored an independent closed-form RHS to one Mathematica assertion. No symbolic intermediate, no carry-forward formula, and no printed result downstream of this stage changes. Both engines' `I_n`, `zeta_n^{phys}`, `zeta_n^{twin}`, `x`, and `zeta_0^{twin}` outputs remain identical to the pre-fix transcripts in substance (only the rendering of `I_n closed form` differs — now `(Sqrt[2]*Sqrt[l])/((1/2 + n)*Pi)` rather than `(2*Sqrt[2]*Sqrt[l])/(Pi + 2*n*Pi)`). The `PASS` lines are unchanged. Downstream units consuming the boxed `zeta` / `I_n` carry-forwards from stage 049 are unaffected.

## Side observations (non-blocking)

- No `redteam/exec_logs/stage_049_sympy.log` or `stage_049_mathematica.log` were produced for this iteration; only the diff capture is present. Verification fell back to the canonical `scripts/output/` and `mathematica/output/` transcripts as the prompt allows. Worth confirming the orchestrator is wiring exec logs into `redteam/exec_logs/` for v2+ iterations so future verifiers do not need the fallback.
- The SymPy script (`scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py`) mtime is 2026-05-22, predating this audit unit's v2; the saved SymPy transcript was nonetheless regenerated on 2026-05-26 02:55. This is consistent (SymPy never needed an edit for v2) but worth noting.
- The script banner still reads "STAGE 32" while the file name is `stage049`. Pre-existing; flagged in the v1 verification as well; not in scope for this directive.

## Verdict justification

The single finding F1 (a regression introduced by v1's F2 fix) was applied exactly as directed: the integrator self-comparison was replaced with an integrator-vs-independent-closed-form comparison. The diff shows no collateral edits, the Mathematica transcript confirms the two sides now print distinct algebraic forms while still reducing to zero under `FullSimplify`, and every `PASS` line is intact. SymPy was untouched and still passes all seven assertions. No exec logs were captured in `redteam/exec_logs/`, but the canonical output transcripts are fresh (post-edit mtimes) and show clean `PASS` on every assertion. No regressions, no new findings inside scope.
