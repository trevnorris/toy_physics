---
unit_id: 094
batch: IV.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 094

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
The static-limit block was de-tautologized symmetrically in both engines.

SymPy (`...stage094..._sympy_audit.py`): two accumulators are initialized before the
generic-cross loop (py:47-48) and summed inside it (py:53-54):
`Kg2_overlap += I_mass` where `I_mass = domega(Y00*Y)` (py:50), and
`Kg4_overlap += I_lap` where `I_lap = domega(Y00*(-lap_s2(Y)))` (py:52). The
static block (py:67-79) now sets `K_g2 = sp.simplify(Kg2_overlap)`,
`K_g4 = sp.simplify(Kg4_overlap)` — no bare `sp.Integer(0)` assignment — with
can-fail `assert K_g2 == 0` / `assert K_g4 == 0`, then derives `eps_2`/`eps_4`
from them, and adds individual `assert c_pole == sp.Rational(1,4)` /
`assert c_geom == sp.Rational(3,4)` (py:77-78) alongside the retained sum check.

Mathematica (`...stage094..._mathematica_audit.wl`): mirror accumulators
`kg2Overlap`/`kg4Overlap` initialized at wl:48-49 and summed in the `Do` loop at
wl:62-63 (`+ overlap` = `dOmega[y00*y]`, `+ lapCross` = `dOmega[y00*(-lapS2[y])]`).
Static block (wl:79-93): `Kg2 = FullSimplify[kg2Overlap,...]`, `Kg4 = FullSimplify[kg4Overlap,...]`,
with `expectZero["K_(g,2) overlap moment", Kg2]` / `expectZero["K_(g,4)...", Kg4]`,
derived `eps2`/`eps4`, plus new `expectZero["c_pole - 1/4", ...]` and
`expectZero["c_geom - 3/4", ...]`.

**Assessment:**
Correct and non-tautological. `K_g2`/`K_g4` are now the accumulated sums of the
same `domega`/`dOmega` overlap integrals proven zero mode-by-mode above — they are
0 *because* those integrals vanished, not by assignment. If any per-mode overlap
returned nonzero, the accumulator would be nonzero, `assert K_g2 == 0` /
`expectZero["K_(g,2)..."]` would fire, and it would propagate to `eps_2 != 0`.
The two named constants are now pinned individually (`c_pole == 1/4`, `c_geom == 3/4`)
rather than only via the complement-agnostic sum. Both engines were updated
symmetrically (same two integrals feed each accumulator). Edits match the
directive's required-change blocks exactly; the diff shows no collateral changes
beyond the four named hunks. Emitted values unchanged.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
`eps_2 = 0 ; eps_4 = 0 ; c_pole = 1/4 ; c_geom = 3/4` (L37);
all five `<Y00|Y2A> = 0` and `(-Delta)Y2A - 6Y2A = 0` residuals zero (L8-31);
`Stage 094 theorem verified...` (L39).

**Mathematica:** exit=0. Notable lines:
`K_(g,2) overlap moment = 0` / `PASS` (L78-79), `K_(g,4) overlap moment = 0` / `PASS`
(L80-81); new `PASS: c_pole - 1/4` (L87) and `PASS: c_geom - 3/4` (L89);
`eps_2 = 0; eps_4 = 0; c_pole = 1/4; c_geom = 3/4` (L92).

**Output freshness:** confirmed. Script mtimes 1780722767 (.py) / 1780722784 (.wl);
output mtimes 1780723321 for both .txt — outputs newer than scripts, regenerated
post-fix. The new derived-overlap PASS lines appear in the committed Mathematica
transcript, confirming the refresh.

## Material-change assessment

`material_change`: false. No emitted/derived value changed: `eps_2=0`, `eps_4=0`,
`c_pole=1/4`, `c_geom=3/4`, and `K_{g,2}=K_{g,4}=0` are identical pre- and post-fix.
The edit only strengthens how those values are established (consume vs. hardcode).
Downstream units depending on stage 094's outputs see no value change.

## Side observations (non-blocking)

- `Kg4_overlap` accumulates the `<Y00|(-Delta)Y2A>` Laplacian-weighted moment
  (= 6× the mass overlap by the eigenvalue identity), which is the natural
  Omega_Q^4-channel analog the directive specified; consistent with paper D2'.
  Not a finding.
- `I_grad`/`I_pot` remain computed in the loop but unused for the accumulators
  (as in the original); harmless, matches the prior structure.

## Verdict justification

The sole finding (tautological static-limit block) is fully resolved in both
engines symmetrically: `K_g2`/`K_g4` now derive from the accumulated, genuinely
computed S^2 overlap integrals with can-fail asserts, and `c_pole=1/4`/`c_geom=3/4`
are pinned individually. Both scripts exit 0, transcripts are fresh and print the
unchanged deliverable values, and the diff shows no regressions or out-of-scope
edits. Verdict: verified.
