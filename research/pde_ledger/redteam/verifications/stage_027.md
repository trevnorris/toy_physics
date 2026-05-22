---
unit_id: 027
batch: II.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 027

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl:91-96` the `wallStiffness[]` module was edited exactly as the directive specified. Per the captured diff (`redteam/exec_logs/stage_027_diff.patch`):

- The `Module` local-variable list now includes `gEta`.
- A new line builds the wall operator: `gEta = -tW*D[chi, {s, 2}] + (kEta + 6*tOmega)*chi;`.
- `kGeo` is now derived: `kGeo = FullSimplify[Integrate[chi*gEta, {s, 0, l}], Assumptions -> $Assumptions];` (was a literal closed-form assignment).
- `kGeoExpected` is now a distinct literal closed form: `kEta + 6*tOmega + tW*Pi^2*Sin[theta]^2/l^2;` (was previously `kGeoExpected = kGeo`).

No other lines, modules, or files were touched (the diff is exactly 5 changed lines, all inside `wallStiffness[]`). The assertion strings, the `K_geo(theta_max)` block, and the `branchSubstitution[]` `theta=0` check were left untouched as required.

**Assessment:**

The edit faithfully implements the directive. The two sides of the `expectZero["K_geo - expected", kGeo - kGeoExpected]` assertion are now structurally distinct: `kGeo` is the FullSimplify of a symbolic integral over `[0, l]` of `chi * (-tW * chi_ss + (kEta + 6*tOmega) * chi)`, while `kGeoExpected` is a hand-written closed form. The assertion can now fail (e.g., if `chi_ss` were sign-flipped, or if the basis prefactor `sqrt(2/l)` were wrong, or if the eigenvalue `(pi/l)^2` were misread).

The saved Mathematica output at `mathematica/output/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.txt` (mtime 2026-05-21 17:04:15, ~1 min 23 s after the script mtime of 17:02:52, so post-fix) confirms substance:

- Line 62: `K_geo(theta) = kEta + 6*tOmega + (Pi^2*tW)/(2*l^2) - (Pi^2*tW*Cos[2*theta])/(2*l^2)`. This is Mathematica's canonical FullSimplify rendering of `kEta + 6*tOmega + Pi^2*tW*Sin[theta]^2/l^2` (using `Sin[theta]^2 = (1 - Cos[2*theta])/2`). Crucially this is *not* the literal form of `kGeoExpected`, which is what one would see if the assignment had still been hard-coded. So the integral was genuinely performed and FullSimplify chose its own canonical form.
- Line 64: `PASS: K_geo - expected` — confirms FullSimplify of `(integral result) - (closed-form RHS)` simplifies to 0 under the assumptions. This is now a real check that the integral identity holds, not a self-comparison.
- Line 67: `PASS: K_geo(theta_max) - [K_eta + 6 T_Omega + 2 T_w pi^2/(11 L^2)]` — substitution of `Sin[theta] -> -Sqrt[2]/Sqrt[11]`, `Cos[theta] -> 3/Sqrt[11]` into the *derived* `kGeo` correctly reduces to the target. With the new `Cos[2*theta]` form, `TrigExpand[kGeo] /. maxSubs` is the genuine reduction path; the prior tautology defect is gone.
- Line 98: `PASS: K_geo(theta=0) - (K_eta + 6 T_Omega)` — `theta -> 0` substitution into the derived `kGeo` correctly reduces to `kEta + 6*tOmega`.

The new assertions are non-tautological: the integral and closed form are independently constructed expressions whose equality is checked by FullSimplify-of-difference. No collateral edits, no renamed symbols, no changes to assertion labels.

## Exec log assessment

**SymPy:** exit log was not captured by the orchestrator at `redteam/exec_logs/stage_027_sympy.log` (file not present). Per directive scope, SymPy script was not edited, so re-execution would not have been required; however, the saved sympy output at `scripts/output/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.txt` is fresh (mtime 2026-05-21 17:04:00, newer than the unchanged `.py` from 2026-04-01) and shows all assertions reaching `= 0` print lines (e.g. line 89: `K_geo(theta) = K_eta + 6*T_Omega + pi**2*T_w*sin(theta)**2/L**2`, line 90: `K_geo - expected = 0`).

**Mathematica:** exit log was not captured by the orchestrator at `redteam/exec_logs/stage_027_mathematica.log` (file not present). The saved `.txt` output is fresh (mtime 2026-05-21 17:04:15, newer than the post-edit `.wl` at 17:02:52) and ends with `Stage 10 Mathematica audit passed.` All five required PASS lines from the directive's verification list appear:

1. `PASS: K_geo - expected` (line 64)
2. `PASS: K_geo(theta_max) - [K_eta + 6 T_Omega + 2 T_w pi^2/(11 L^2)]` (line 67)
3. `PASS: K_geo(theta=0) - (K_eta + 6 T_Omega)` (line 98)
4. All Section I, II, IV, V PASS lines unchanged.
5. Final banner `Stage 10 Mathematica audit passed.`

The derived `K_geo(theta)` printed at line 62 has a different surface form than the prior tautological assignment (now contains `Cos[2*theta]` from FullSimplify canonicalisation), which is the load-bearing evidence that the integral was actually evaluated post-edit.

**Output freshness:** Confirmed. Both saved outputs have mtimes after the script mtimes — sympy output (17:04:00) > script (2026-04-01, unchanged); mathematica output (17:04:15) > script (17:02:52). No stale-output concern.

## Material-change assessment

`material_change`: false.

The edit only converts a self-referential assignment into a derived one; the symbolic value of `kGeo` (and therefore everything downstream — `d0`, `K_req(theta)`, the IV.2 `theta=0` reduction, the V.1 blind-angle no-go) is mathematically identical to what the prior hard-coded form already produced. Both engines printed the same `K_geo(theta_max) = K_eta + 6*T_Omega + 2*Pi^2*T_w/(11*l^2)` reduction before and after, and the SymPy reference value `K_eta + 6*T_Omega + pi^2*T_w*sin(theta)^2/L^2` is unchanged. No unit > 027 depends on a numerically different output of unit 027.

## Side observations (non-blocking)

- The orchestrator did not write `stage_027_sympy.log` or `stage_027_mathematica.log` into `redteam/exec_logs/` for this unit (only the diff was captured). Verification relied on the saved `.txt` outputs in `mathematica/output/` and `scripts/output/`, whose freshness mtimes are consistent with a post-edit run. This is not a finding against Codex; just a logging gap to flag for the orchestrator.
- The Mathematica `wallStiffness[]` module returns `{kappa, kGeo}`, where `kGeo` is the now-FullSimplify'd integral. `branchSubstitution[]` then uses `(TrigExpand[kGeo] /. theta0)` to recover the `theta=0` reduction. Since the new `kGeo` carries a `Cos[2*theta]` term, `TrigExpand` + substitution still reduces correctly (and the saved output confirms the PASS), so this is fine — just worth noting that the downstream call site is robust to FullSimplify's representation choice.

## Verdict justification

The single finding (F1) was applied exactly as specified — the directive's "After" code block matches the current file's lines 91–96, the diff shows only those four mechanical changes, and the post-edit Mathematica output shows the wall-stiffness assertion now passes through a genuine `Integrate`+`FullSimplify` evaluation (visible from the `Cos[2*theta]` canonical form on the `K_geo(theta) = ...` print line) rather than a self-comparison of a hard-coded expression. All downstream `theta_max` and `theta=0` reductions still PASS, confirming no regression. Material change is false because the symbolic equivalence is preserved.
