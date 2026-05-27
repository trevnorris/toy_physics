---
unit_id: 094
batch: IV.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 094

## Per-finding outcomes

### F1 — engine_disagreement

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage094_isotropic_geometry_decoupling_mathematica_audit.wl:60` now reads
```
cCross = FullSimplify[mu*overlap - tw*overlap - tOm*lapCross - kPot*overlap, Assumptions -> $Assumptions];
```
The `tw*gradCross` term was replaced with `tw*overlap`, matching the SymPy form `mu*I_mass - Tw*I_mass - TOm*I_lap - K*I_pot` at sympy line 52. The independent `gradCross` orthogonality assertion at .wl line 65 is retained.

**Assessment:**
The edit exactly matches the directive's required "Before/After" block. Both engines now compute the same symbolic identity for the bilinear `Y00`↔`Y2A` cross coefficient, faithful to the `T_w(∂_w η)²` angular factorization in the wall action. The Mathematica transcript still shows `PASS: Generic isotropic cross coefficient C_0,<label>` for all five Y2A labels (lines 20, 32, 45, 58, 71 of output). The new assertion is non-tautological — it still encodes the action structure and merely happens to evaluate to 0 because `overlap = 0` and `lapCross = 0` are independently asserted upstream. No collateral edits beyond the one-line fix.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy: lines 56–69 of `scripts/moving_throat_pde_stage094_isotropic_geometry_decoupling_sympy_audit.py` introduce `Omega_Q`, `K_pole` (positive symbols), set `K_g2 = K_g4 = 0` from orthogonality, compute `eps_2 = Omega_Q**2 * K_g2 / K_pole` and `eps_4 = Omega_Q**4 * K_g4 / K_pole`, assert each is 0, define `c_pole = 1/4`, `c_geom = 3/4`, and assert `c_pole + c_geom == 1` with a labeled print.
- Mathematica: lines 71–85 of the `.wl` insert the mirror block: `Kg2 = Kg4 = 0`, define `OmegaQ`, `Kpole`, compute `eps2`, `eps4`, `cPoleStatic = 1/4`, `cGeomStatic = 3/4`, with three new `expectZero` calls (`eps_2 (static limit)`, `eps_4 (static limit)`, `c_pole + c_geom - 1`) and a labeled `Print`.

**Assessment:**
The blocks match the directive's required-change patch sketches structurally and semantically. SymPy transcript line 32 prints `eps_2 = 0 ; eps_4 = 0 ; c_pole = 1/4 ; c_geom = 3/4` as specified. Mathematica transcript adds three new PASS lines at output lines 74, 76, 78, with the named summary line 79 `eps_2 = 0; eps_4 = 0; c_pole = 1/4; c_geom = 3/4`. The assertions are mildly tautological by construction (Kg2 and Kg4 are *assigned* zero, so eps_2 and eps_4 trivially vanish, and `1/4 + 3/4 == 1` is arithmetic) — but the directive explicitly authorized this as a *naming* exercise to register `K_g2`, `K_g4`, `eps_2`, `eps_4`, `c_pole`, `c_geom` as named symbolic objects in the script so the card's Check #1 has a script-side anchor. The orthogonality theorem (the genuine load-bearing claim) is verified upstream via the non-tautological angular integrals (A1–A4, A6–A10). The static-limit block is documentary; given the directive scope, that is acceptable. The Mathematica script uses `cPoleStatic`/`cGeomStatic` as local symbol names instead of `cPole`/`cGeom` from the patch sketch — minor cosmetic deviation, semantically identical, does not affect output or assertions.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- Line 27–31: `Generic isotropic cross coefficient C_0,<label> = 0` for all five Y2A.
- Line 32: `eps_2 = 0 ; eps_4 = 0 ; c_pole = 1/4 ; c_geom = 3/4`.
- Line 34: `Stage 094 theorem verified: isotropic l=0 <-> l=2 cross terms vanish exactly.`

**Mathematica:** exit=0. Notable lines:
- 34 total `PASS:` lines (counted): 1 Y00 norm + 5 Y2A norms + 5 overlap + 5 Laplacian-eigenvalue + 5 gradient-cross + 5 Laplacian-cross + 5 cCross + 3 static-limit = 34. Matches the orchestrator's expected count exactly.
- Line 73–78: three new `PASS:` lines for `eps_2 (static limit)`, `eps_4 (static limit)`, `c_pole + c_geom - 1`.
- Line 79: `eps_2 = 0; eps_4 = 0; c_pole = 1/4; c_geom = 3/4`.
- Line 81: `Stage 094 Mathematica audit passed.`

**Output freshness:** Both outputs are newer than their scripts.
- SymPy script mtime 2026-05-27 11:13:56; output mtime 2026-05-27 14:28:35 (output is newer — post-fix).
- Mathematica script mtime 2026-05-27 11:14:05; output mtime 2026-05-27 14:29:24 (output is newer — post-fix).
Both transcripts are post-edit.

## Material-change assessment

`material_change`: false.

The F1 fix changes which symbolic expression the Mathematica `cCross` represents, but the final asserted value (0) is unchanged; no downstream stage depends on the symbolic form, only the orthogonality result. The F2 additions are new *named* symbolic objects (K_g2, K_g4, eps_2, eps_4, c_pole, c_geom) with hard-coded values that match notes §3; they introduce no new derived quantities a downstream unit could depend on. No re-audit of units > 094 is needed on the basis of stage 094's edits.

## Side observations (non-blocking)

- The F2 SymPy block uses `Omega_Q` and `K_pole` as positive symbols but they never appear non-trivially (numerator is 0 either way). This is by design per the directive ("assigns `K_g2 = 0`, `K_g4 = 0` as the orthogonality result") and so the static-limit assertions are tautological-by-construction. A stronger version would have computed `eps_2`, `eps_4` from a downstream-passed formula and shown the orthogonality forces zero, but the directive explicitly accepted the weaker naming-only variant.
- The Mathematica script names the local constants `cPoleStatic`/`cGeomStatic` instead of the patch sketch's `cPole`/`cGeom`. No effect on transcript or assertions; record only for parity with directive text.
- Banner sweep STAGE 77/077 → STAGE 094 noted in the Applied block is reflected in the .wl line 26 banner string and the .py final print on line 71 — consistent with the renumbering pass and not a new issue.

## Verdict justification

Both findings are resolved with edits that match the directive exactly (F1 one-line replacement; F2 multi-line append at the specified insertion points in both engines). The SymPy transcript exits 0 with the expected new print line. The Mathematica transcript exits 0 with exactly 34 `PASS:` lines, matching the orchestrator-expected count of 30 orthogonality + 1 Y00 norm + 3 static-limit. Both outputs are fresher than their scripts (post-fix re-runs). No regressions; no collateral edits beyond the banner number sweep. Verdict: `verified`.
