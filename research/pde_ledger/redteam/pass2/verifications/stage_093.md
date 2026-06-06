---
unit_id: 093
batch: IV.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T23:30:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 093

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
Codex added a symbolic-eps anchor block after the four existing deliverable checks in
`mathematica/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.wl:45-51`:
- `ClearAll[e2, e4]; cPoleGen = (1 + e4)/(4*(1 + e2)^2);` (fresh free symbols, line 45-46)
- `expectZero["c_pole static-limit (symbolic)", (cPoleGen /. {e2 -> 0, e4 -> 0}) - 1/4];` (line 47)
- can-fail directional guard `If[TrueQ[FullSimplify[cPoleGen - 1/4] === 0], fail[...], pass["c_pole eps-dependent"]];` (line 48-50)
- off-static probe `expectZero["c_pole off-static probe (e2=1,e4=0)", (cPoleGen /. {e2 -> 1, e4 -> 0}) - 1/16];` (line 51)
The original lines 28-43 (four deliverable checks: 1/4, 3/4, 4/3, 1/3) are byte-for-byte unchanged; diff is purely additive (+8 lines). No collateral edits.

**Assessment:**
The edit matches the directive's required change exactly (fresh symbols `e2/e4` so the
`eps2=0; eps4=0;` bindings don't leak, the three new checks in the prescribed order). The
new checks are genuinely non-tautological:
- Static-limit symbolic: `(1+0)/(4*(1+0)^2) = 1/4` — a real substitution into the structured
  expression, so the `^2` and numerator are exercised (not the pre-collapsed literal).
- Off-static probe `e2=1,e4=0` targets `1/(4*4) = 1/16 ≠ 1/4`; it fails if the squared-denominator
  exponent or numerator/denominator structure is altered. The auditor's chosen probe correctly
  avoids the degenerate `e2=1,e4=3` point that re-hits 1/4.
- The directional guard takes the `pass` branch because `FullSimplify[cPoleGen - 1/4]` is not
  identically 0; it can fail iff the formula were spuriously eps-independent — a real can-fail.
Exec log shows all three new PASS lines plus the four unchanged deliverable PASS lines; exit 0.

## Exec log assessment

**SymPy:** exit=n/a. Stage 093 is Mathematica-only by design (card "SymPy audit: none yet";
manifest `is_status_only_candidate: True`); no SymPy script exists and no `stage_093_sympy.log`
is present (confirmed absent). Correctly classified n/a.

**Mathematica:** exit=0. Notable lines:
- `PASS: c_pole static-limit (symbolic)` (new)
- `PASS: c_pole eps-dependent` (new)
- `c_pole off-static probe (e2=1,e4=0) = 0` → `PASS: c_pole off-static probe (e2=1,e4=0)` (new)
- four original PASS lines (c_pole−1/4, c_geom−3/4, rho_alpha−4/3, zeta_req−1/3) unchanged
- `Stage 093 Mathematica audit passed.` / `# exit_code: 0`

**Output freshness:** confirmed. `.txt` mtime 2026-06-05 23:22:01 is newer than `.wl` mtime
2026-06-05 23:10:22. Committed `.txt` content matches the exec log line-for-line (includes all
three new PASS lines).

## Material-change assessment

`material_change`: false. The edit is additive verification only — it adds symbolic-anchor and
probe checks and leaves the four carry-forward deliverable values (c_pole=1/4, c_geom=3/4,
rho_alpha=4/3, zeta_req=1/3, cited by downstream stages 100–106) untouched. No derived result
changed; no downstream unit is affected.

## Side observations (non-blocking)

None. The probe and guard are self-contained and do not interact with the static-point block.

## Verdict justification

The single low-severity `insufficient_verification` finding is fully resolved. Codex's additive
edit anchors the obstruction formula's symbolic structure (static-limit reduction, eps-dependence
guard, and a 1/16≠1/4 off-static probe) without altering any deliverable value, the new checks
are genuinely can-fail (not tautological), and the Mathematica run exits 0 with all seven PASS
lines and a fresh committed output. Verdict: verified.
