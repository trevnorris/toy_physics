---
unit_id: 223
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T19:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 223

## Per-finding outcomes

The pass-2 auditor returned `verdict: clean` with **0 findings**. There is therefore
no directive file for stage 223 (confirmed absent at
`redteam/pass2/directives/stage_223.md` — expected), and Codex applied **no edits**
(the diff patch `redteam/pass2/exec_logs/stage_223_diff.patch` is 0 bytes — legitimately
empty, not a failed capture). Both `.py` and `.wl` are byte-identical to HEAD (empty
`git diff` for all four script/output files).

With no findings to re-check, this verification CONFIRMS that the `clean` verdict still
holds by (a) re-reading both current scripts, (b) confirming both fresh exec logs exit 0
with load-bearing PASS lines, (c) spot-checking that the load-bearing assertions are
non-tautological and the `.wl` is a genuine independent recomputation (not a
transliteration), and (d) confirming the 5PN R_Q values `138.814136942081` /
`137.502546600713` are emitted with no surviving `206/205`.

### Confirmation checklist

**Classification:** resolved (vacuously — 0 findings; clean verdict re-confirmed)

1. **Exec logs pass.** SymPy log ends `Stage 223 audit completed successfully.` then
   `# exit_code: 0`. Mathematica log ends `All Stage 223 Mathematica checks M1-M9
   completed successfully.` then `# exit_code: 0`. Every M1-M9 block prints `PASS:`;
   no `FAIL`; all reported residuals ≤ ~1e-13.

2. **Diff/edits.** Empty diff patch + clean `git diff` confirm zero Codex edits — the
   expected zero-correction VII.1 state. Nothing to regress.

3. **Load-bearing R_Q values present, no 206/205.** SymPy `.py` L372 pins the
   λ_W=0.2 scan row's wall residues as `138.814136942081, 137.502546600713`; the SymPy
   output L60 emits `(0.2, ..., 138.81413694208146, 137.50254660071312)`. A targeted
   grep for `206.814136942081` / `205.502546600713` (and `206.8141` / `205.5025`)
   across both scripts and both `.txt` outputs returns **no hits** — no surviving typo
   in any deliverable-carrying script/output. (The `.wl` scan starts at λ_W=2/5, L301,
   so it does not emit a λ_W=0.2 row; the corrected values are corroborated cross-engine
   via the shared `rqOmega` formula at the overlapping λ_W=0.4 row, which matches.)

4. **Compatibility surface (M3) non-tautological.** `.py` L154-156 builds `K_pole =
   3(M+B2+Z2)^2/(B4+Z4)+B0+Z0` (one-pole side) and `K_norm = N0/P0_target+B0+Z0`
   (normalization side) from disjoint constructions, then `sp.solve(sp.Eq(K_norm,K_pole),
   P0_target)[0]` solves the equality and L160-161 assert it equals the independently
   built boxed target `N0(B4+Z4)/(3(M+B2+Z2)^2)`. This is a genuine consistency theorem,
   not `x==x`. The `.wl` M3 (L157) reproduces with `First[Solve[kNormIso==kPoleIso,
   p0Target]]` and ADDS the independent M3b check (L161-164) `nZero/boxedTarget -
   3(...)^2/(b4+z4)` that SymPy lacks.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `kappa = 2*sqrt(2)/pi` (overlap constant emitted)
- `(0.2, 0.0005769708798430451, 29.31584648723137, 138.81413694208146, 137.50254660071312)` (corrected λ_W=0.2 R_Q row)
- `Stage 223 audit completed successfully.`

**Mathematica:** exit=0. Notable lines:
- `M1 kappa integral residual = 0` → `PASS` (native `Integrate`, not a transcribed closed form)
- `M3 solved target minus boxed target residual = 0` → `PASS`; `M3 equivalent N0/P0 compatibility equation residual = 0` → `PASS` (the extra M3b check)
- `M5 leading quartic coefficient is mass residual = 0` → `PASS` (check SymPy does not perform)
- `M7 positive omega roots residuals = {1.5e-15, 2.7e-15, 3.9e-16, 6.7e-15}` → `PASS` (NSolve @ WorkingPrecision 60)
- `All Stage 223 Mathematica checks M1-M9 completed successfully.`

**Output freshness:** confirmed. The committed `.txt` outputs (mtime 2026-06-02 16:37:47)
are newer than their scripts (`.py` 16:23, `.wl` 16:26). The orchestrator re-ran both
engines directly (fresh logs dated 2026-06-09T19:21:18); both exit 0 and the fresh log
bodies match the committed `.txt` deterministically. `git diff` for the two `.txt`
outputs is empty (byte-identical to HEAD), so committed-output mtime is not a concern.

## Independent-derivation spot-check (Mathematica)

The `.wl` mirrors the SymPy claim sequence (M1-M9 = the nine deliverables, which is
expected) but the derivation route uses native Mathematica primitives throughout and in
several places does strictly more:
- M1: `FullSimplify[Integrate[uConst fProfile, {s,0,len}]]` — genuine symbolic integral.
- M3: `First[Solve[kNormIso==kPoleIso, p0Target]]` — independently solves the stiffness
  equality; adds redundant M3b form absent in SymPy.
- M5: `CoefficientList[fyExpanded, y]` with `Exponent==4` AND `Last[coeffs]-mass==0`
  (leading coefficient equals mass) — a different decomposition plus an extra check.
- M7: `NSolve[poly==0, y, WorkingPrecision->60]` with independent positive-real
  selection — a different numeric kernel than SymPy's NumPy `np.roots`.

The only shared "operation" is `Solve`/`solve` on the same stiffness equality built from
the same physical moment definitions (`B0..Z4, N0`), which are the stage's physical
premises both engines must start from. No SymPy-specific choreography (substitution order,
helper naming, lambdify scaffolding) is echoed. **Not a transliteration — independent.**

## Material-change assessment

`material_change`: false. Codex made no edits; both scripts are byte-identical to HEAD.
No derived result changed, so no downstream unit can be affected by this verification.

## Side observations (non-blocking)

- The auditor's NOTE re. a deferred numbering class (formerly "Stage 240" in the notes
  title) is out of scripts-only scope and falls under the project-wide stem-keyed
  renumber pass (P4-53 residual class), not a pass-2 blocker. Not re-examined here.
- Symbol domains in both engines match the physical setup, and the `.wl` `$Assumptions`
  (L111-128) mirror the SymPy `positive=True` declarations and add the non-degeneracy
  guards `kWall-b0-z0 != 0`, `b4+z4 != 0`, `mass+b2+z2 != 0` that protect the divisions.
  No concern.

## Verdict justification

`verified`. Stage 223 is the zero-correction case: the auditor's `clean`/0-findings
verdict has no findings to re-resolve, so verification reduces to confirming the clean
state still holds. Both fresh exec logs exit 0 with all M1-M9 / sample PASS lines and
≤1e-13 residuals; the diff patch and `git diff` confirm zero Codex edits; the load-bearing
M3 compatibility surface is independently constructed on both sides and solved (not
tautological); the `.wl` is a genuine native recomputation (Integrate / Solve / NSolve@60
/ CoefficientList) with two checks SymPy lacks, not a transliteration; and the corrected
5PN R_Q values `138.814136942081` / `137.502546600713` are emitted with no surviving
`206/205` anywhere in the scripts or outputs. findings_total=0, findings_resolved=0,
material_change=false.
