---
unit_id: 226
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 226

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created the required second-engine audit at the exact target path
`mathematica/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_mathematica_audit.wl`
(315 lines). It defines `expectZero` (line 22, `Exit[1]` on failure at line 31) and `expectClose` (line 35, `Exit[1]` at line 41) and independently re-derives the M1–M7 manifest, printing 34 `OK ` lines and ending `All Stage 226 Mathematica checks passed.` + `Exit[0]`. The independent exec re-run (`stage_226_mathematica.log`) confirms exit 0, 34 OK, 0 FAIL. Codex also applied the USER-AUTHORIZED notes renumber (notes-only) and the orchestrator re-ran SymPy (output `.txt` refreshed).

**Assessment — correct and non-tautological.**

INDEPENDENCE / anti-transliteration (verified against the SymPy `.py`):
- **M1 bridge** — `.wl` uses `Normal[Series[p0A,{eps,0,1}]]` then `Coefficient[...,eps]/(lam p0)` (lines 83-84); SymPy uses `diff(P0A,eps).subs(eps,0)/lam` (py:42). Different primitive (Series vs symbolic derivative). PASS.
- **Coefficient extraction / nullspaces (M5, M6)** — `.wl` uses `CoefficientArrays[Expand[expr],vars][[2]]` (`coeffRow`, line 46) then native `MatrixRank`/`NullSpace` (lines 232-233, 243-244); SymPy builds the matrix from per-variable `expr.coeff(v)` and calls `.rank()`/`.nullspace()` (py:217-220, 270-273). Different extraction route. PASS.
- **Projector norms (M5, M6)** — `.wl` `projectedNorm` (lines 48-52) uses `Orthogonalize[basis]` (Gram–Schmidt) → `Transpose[Q].Q` projector → `Sqrt[coeff.proj.coeff]`; SymPy uses `QRdecomposition()` → `Q*Q.T` (py:251-252, 306-307). Different orthonormalization. PASS.
- **Compiler** — both differentiate an eps-dressed bundle (`D[expr,eps]/.eps->0`). This is the one shared technique, but the directive explicitly authorized it ("a derivative is fine — it is a different surface form than Series"); the bundle is reassembled from the physical primitives with native Mathematica symbol names (lines 119-139, 194-209), not transcribed. PASS.
This is NOT a line-by-line port.

CLAIM COVERAGE (all 7 manifest items, 34 OK checks accounted for):
M1=1, M2=3, M3=8 (incl. exact one-pole identity), M4=2, M5=3 (rank 2, nullity 3, σ_even), M6=9 (rank 3, nullity 2, 6× M6c residuals, σ_transfer), M7=8 (4 even + 4 transfer budgets). Sum = 34 = the OK count exactly, 0 FAIL. No manifest item silently uncovered. The corridor checks verify basis-INVARIANT facts (ranks 2/3, nullities 3/2, σ_even≈2.6739, σ_transfer≈2.3156) rather than echoing raw basis vectors — exactly as the directive demanded (basis vectors are returned only up to basis/scale).

NON-VACUITY (spot-checks of load-bearing checks):
- M1 `expectZero[xi1Series - xiLoad]` — xi1Series from Series; xiLoad is the independent symbolic target `n01/n0 - d01/d0`. Not a self-equality; a wrong O(eps) coefficient hard-fails. PASS.
- σ_even / σ_transfer `expectClose[sigma, 2.6739.../2.3156...]` — σ is computed from the independently-derived `xiCoeffExact` projected onto the natively-computed nullspace, compared to external notes literals; a wrong projector or coefficient row yields a different σ → `Exit[1]`. PASS.
- M6c `transferMatrixExact.transferNull[[i]]` — `transferNull = NullSpace[transferMatrixExact]`, so this residual is zero-by-construction (the same redundancy the auditor flagged for SymPy A14). It is HARMLESS: the load-bearing pure-transfer checks are the rank-3/nullity-2 corridor dimension and σ_transfer (external-literal comparison), both substantive. M6c is a redundant confirmation, not a vacuous pass that carries the claim.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `sigma_even = 2.67386816837172775`, `sigma_transfer = 2.31561904386055442`, `rank = 2 nullity = 3`, `rank = 3 nullity = 2`, `All Stage 226 checks passed.`

**Mathematica:** exit=0. Notable lines: `OK M1 Xi1 equals load-defect bridge`; `OK M5 strict even-gate rank 2` / `nullity 3`; `sigma_even = 2.67386816837172774819...` + `OK M5 sigma_even`; `OK M6 transfer matrix rank 3` / `nullity 2`; `sigma_transfer = 2.31561904386055441584...` + `OK M6 sigma_transfer`; `All Stage 226 Mathematica checks passed.` 34 OK, 0 FAIL. The two engines agree on both corridor norms to ~17 sig-figs and on all base data.

**Output freshness:** confirmed. `.wl` mtime 2026-06-02 16:47:44; its output `mathematica/output/...mathematica_audit.txt` mtime 17:06:23 (newer, post-fix re-run). SymPy `.txt` mtime 17:06:23, refreshed by the orchestrator re-run; its diff is header-format-only (every printed numeric value byte-identical). SymPy `.py` is untracked-clean (NOT in git porcelain) — unmodified, as required.

## Material-change assessment

`material_change`: false.

A new second engine was added and the SymPy output was re-generated, but NO derived numeric result changed: σ_even, σ_transfer, ranks/nullities, base data, and all twelve budget literals are identical across both engines and to the pre-fix SymPy output. The notes edit is a pure label renumber (240→223, 241→224, 242→225, 243→226; supporting-file stage243→stage226) with every formula, matrix entry, and value unchanged. Downstream Part-VII units that consume the σ values / budgets see no change. No `upstream_stale` concern arises from a value perspective.

## Side observations (non-blocking)

- The notes renumber matches the directive's authorization exactly and is consistent with the auditor's "provenance-label drift" informational note (script comments / card use canonical 223–226; notes had been on the pre-renumber 240–243). The notes labels now agree with the `.py` comments and card. Not a finding; noting that the prior label mismatch the auditor flagged is now closed.
- M6c being zero-by-construction mirrors the auditor's A14 observation for the SymPy script; harmless and not raised as a new finding.

## Verdict justification

The single finding (missing Mathematica engine) is fully resolved: Codex added a genuinely independent `.wl` at the exact target path that re-derives M1–M7 with native primitives (Series, CoefficientArrays, MatrixRank/NullSpace, Orthogonalize-based projector) through a different decomposition than the SymPy script — not a transliteration. All 34 checks pass with `Exit[1]` hard-fail paths, both engines agree on every load-bearing value to ~17 digits, the authorized notes renumber is value-preserving and notes-only, the SymPy `.py` and paper card are untouched, and outputs are fresh. Verdict: verified; material_change: false.
