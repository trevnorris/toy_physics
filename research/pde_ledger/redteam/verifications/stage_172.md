---
unit_id: 172
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 172

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
In `mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl`, the two `Normal[Series[...,{eps,0,1}]]/(eps*lam)` slope-extraction blocks were replaced by an implicit-differentiation route:

- L42-46 (was L42-43): `deltaU2`/`deltaP0` now obtained via `Solve[D[(u2 + t*deltaU2)*(d0 + t*dD0) + (d2 + t*dD2), t] == 0 /. t -> 0, deltaU2]` and the analogous relation `(p0 + t*deltaP0)*(d0 + t*dD0) - (n0 + t*dN0)` for `deltaP0`.
- L79-84 (was L76-77): `deltaU2Star`/`deltaU4Star` now obtained via `Solve`/`D` on the canonical-branch defining relations, including `(u4Star + t*deltaU4Star)*(d0 + t*dD0)^2 - ((d2Star + t*dD2)^2 - (d0 + t*dD0)*(d4Star + t*dD4))`.

Both blocks match the directive's "after" code verbatim, including the prefatory comment line.

**Assessment:**
Correct and complete. The new route is structurally independent of the SymPy script: SymPy truncates the explicit quotient `-dA2/dA0` via `series`, whereas the `.wl` now differentiates the implicit defining relations `u_2 D_0 + D_2 = 0`, `P_0 D_0 - N_0 = 0`, and `u_4 D_0^2 - (D_2^2 - D_0 D_4) = 0` at the base point. This is exactly the independence-restoring path F1 prescribed. The fresh dummy variable `t` is set to 0 inside each `Solve`/`D` and never leaks downstream, so no symbol pollution.

The seven `expectZero[...]` residual checks (call sites at wl L55, L59, L62, L90, L99, L103, L109) do NOT appear in the git diff hunks, confirming they are byte-unchanged; likewise the `kA`/`gA` definitions (L51-52), the canonical-branch setup, the print statements, and `Exit[0]`. The exec log confirms the printed slope forms are identical to the pre-fix forms recorded in the report: `delta u_2^(A) = -((dD2 + dD0*u2)/d0)`, `delta P_0^(A) = (d0*dN0 - dD0*n0)/d0^2`, `delta u_4^(A) on canonical branch = -1/81*(5*dD0 + 18*dD2 + 81*dD4)/d0`. All seven `PASS:` lines remain and the script exits 0. No collateral edit beyond the banner (F2, below). The SymPy script was correctly left untouched for F1 (its only diff is the banner). The new checks remain non-tautological — each still couples a slope (now derived by the alternative route) to the independently defined obstruction/coefficient, so a sign or factor error would survive `FullSimplify`.

### F2 — stale_banner_label

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.py:31`: `STAGE 155` -> `STAGE 172`.
- `mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:26`: `STAGE 155` -> `STAGE 172`.

**Assessment:**
Correct and minimal. The git diff shows only the single-line banner literal change in each file (the SymPy diff contains nothing else). Both `.txt` transcripts now open with `STAGE 172 — PHYSICAL SLOPE COLLAPSE OF THE LINEAR GROUPED OUTLET PROBLEM`. Print-only; no assertion or result affected.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 172 — PHYSICAL SLOPE COLLAPSE OF THE LINEAR GROUPED OUTLET PROBLEM` (banner fixed)
- `delta u_2^(A) = (-dD0*u2 - dD2)/D0`; `delta P_0^(A) = (D0*dN0 - N0*dD0)/D0**2`
- All seven residual lines print `= 0` (e.g. `K_A + D0*delta_u2 - (1/9-u2)*dD0 = 0`, `Delta_Q^(A) - delta P_0^(A)/P0 = 0`). No `AssertionError`.

**Mathematica:** exit=0. Notable lines:
- `STAGE 172 — ...` (banner fixed)
- `delta u_2^(A) = -((dD2 + dD0*u2)/d0)`; `delta u_4^(A) on canonical branch = -1/81*(5*dD0 + 18*dD2 + 81*dD4)/d0` (slope forms unchanged from report)
- Seven `PASS:` lines, then `Stage 172 Mathematica audit passed.`

**Output freshness:** confirmed. Script mtimes are 15:55 (both); output mtimes are 16:10 (sympy `.txt`) and 16:13 (mathematica `.txt`), both newer than their scripts. Logs are post-fix.

## Material-change assessment

`material_change`: false.

Both edits preserve the verified result set exactly. F1 changes only the derivation path of the slopes in the `.wl`; the printed slope expressions and all seven residual zeros are unchanged (confirmed against the report's recorded values and the exec log). F2 is a print literal. No derived quantity that a downstream unit could depend on changed. No downstream re-audit needed on account of this unit.

## Side observations (non-blocking)

- The notes header numbering offset (notes file reading "Stage 240" vs unit/stage 172) noted in the original report F2 is a prose matter outside scripts-only scope; not actionable here.
- `deltaKappa`/`deltaGamma`/`deltaQ` (wl L96-108, py L112-125) are still constructed from `kA`/`gA` rather than re-derived independently, but F1 only required the load-bearing slope identities (`deltaU2`, `deltaP0`, `deltaU4Star`) to use an independent route, which is done. Not a defect against the directive.

## Verdict justification

Both findings are `resolved`. F1's transliteration is cured: the `.wl` now derives `deltaU2`, `deltaP0`, `deltaU2Star`, and `deltaU4Star` by implicit differentiation of the defining relations — a route structurally distinct from the SymPy series-truncation — while the seven `expectZero` residual targets are byte-unchanged (absent from the diff) and the printed slope forms match the report exactly. F2's banner is corrected in both scripts with no collateral edits. Both engines exit 0 with all seven checks passing on freshly regenerated logs. No regressions, no material change.
