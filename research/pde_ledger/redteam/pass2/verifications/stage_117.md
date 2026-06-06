---
unit_id: 117
batch: IV.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-06T18:05:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 117

## Per-finding outcomes

### F1 — mathematica_transliteration (FULL re-author)

**Classification:** resolved

**What changed:**
The `.wl` (`mathematica/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.wl`) was re-authored to break the line-by-line `Series[]`-based transliteration of the `.py`. The diff (`exec_logs/stage_117_diff.patch`) touches ONLY this `.wl`.

- **§1–§4 (wl:37-44, 47-113):** A new `jetFromBalance[tag, den, num]` helper (wl:37-44) posits `trial = Sum[jetVars[k+1] z^k, {k,0,5}]`, imposes the defining relation `den·trial - num == 0` order-by-order via `Table[Coefficient[Expand[...], z, k] == 0, {k,0,5}]`, and `Solve`s the linear system for the jet coefficients. It replaces:
  - §1: `Normal[Series[(-3 s)/(s lambdaOut(βz)), …]]` → `jetFromBalance["scale/argument", lambdaOut/.z->beta z, -3]` (wl:47).
  - §2: `Normal[Series[(-3+rho)/lambdaR, …]]` → `jetFromBalance["pure Robin", lambdaR, -3+rho]` (wl:58).
  - §3: the `Series` of `lambdaOut - sigma/poleDen` → clears the pole denominator (`mixDenCleared = poleDen·lambdaOut - sigma`, `mixNumCleared = (-3-sigma)·poleDen`) and `jetFromBalance`s (wl:68-77).
  - §4: same pole-clearing approach (`hybDenCleared`, `hybNumCleared`) → `jetFromBalance` (wl:86-91); the residual `Series` checks (cancellation-trivial, compensated-collapse) were replaced by exact `lambdaHybExact = lambdaOut + rho - sigma/poleDen` evaluation and a `compCollapseJet = Sum[Coefficient[compCollapseNumerator, z, k] z^k]` numerator construction (wl:98-113).
- **§5 (wl:127-168):** Re-typed reduced forms (`rC, rhoC, sigmaC, kappaC, gammaC` literals) are GONE. Instead: `dW = 1-κ₀z²-iγ₀z⁵`, `coreMatrix = {{ks,lam},{lam,-kq dW}}`, `coreSource = {gs,gq}`, and `Solve[Thread[coreMatrix.{sCore,qCore}==coreSource]]` (wl:128-130); the Schur complement `deltaCoreEliminated = gs·sCore + gq·qCore` (wl:131-134). The named reduced pieces are then DERIVED via a symbolic-`dSym` variant (wl:136-157): `rC`/`rhoC` extracted by `Coefficient[...,dSym,1]/deltaDLead`, `sigmaTilde = (rhoC-deltaCoreD)(dSym+rC)`, `kappaC`/`gammaC` from a `Solve` of `coreShape - (1-κ_slot z²-iγ_slot z⁵)==0`, and `sigmaC` from `Solve[sigmaSlot(1+rC)==sigmaTilde]`. A hard guard (`coreSchurResidual === 0`, wl:158-160) ties the named forms back to the raw elimination.

**Assessment:**
The route is genuinely independent, not a re-typed mirror.
- **§1–§4:** Verified the `(den, num)` arguments are mathematically correct in each section (§1: `Λ_out(βz)·Y=-3`; §2: `Λ_R·Y=-3+rho`; §3/§4: pole-cleared `poleDen·Λ·Y = poleDen·Λ - σ` so `Y = (norm)/Λ_mix` — a different normalization scheme than the `.py`'s raw-Λ `-L2/L0`, `L2²/L0²-L4/L0`). The coefficient primitive is now an undetermined-coefficient `Solve`, not `sp.series(...).coeff()`. The matching targets (`1/9`, `4/81`), the branch-root assertions, and the odd-norm `χ_Q=1` / `1-9σγ` checks are unchanged.
- **§5:** Traced the elimination algebra: the Schur complement of `{{ks,lam},{lam,-kq dSym}}` against `{gs,gq}` is `rhoC - A/(dSym+rC)` with `A=(ks gq-lam gs)²/(ks²kq)`, `rC=lam²/(ks kq)`, `rhoC=gs²/ks`. The `.wl` recovers each of `rC`, `rhoC`, `sigmaTilde`(=A), `kappaC=κ₀/(1+rC)`, `gammaC=γ₀/(1+rC)`, `sigmaC=A/(1+rC)` from that rational function by `Coefficient`/`Solve` — all DERIVED. `sigmaC` and the §5 collapse are now consequences of the 2×2 `Solve`, not the re-typed literal `(ks gq-lam gs)²/(ks²kq(1+rC))`. The final `deltaCore` (wl:180-181) is built from `deltaCoreEliminated` (the eliminated Schur form), not the re-typed `rho_c - sigma_c/(...)`. The `coreSchurResidual === 0` guard is a real independence check: a mis-derived reduced piece would make it nonzero and `fail`.
- **Variable names genuinely differ from the `.py`:** `jetFromBalance`, `jetVars`, `trial`, `poleDen`, `mixDenCleared`/`mixNumCleared`, `hybDenCleared`/`hybNumCleared`, `lambdaHybExact`, `compCollapseJet`, `dW`/`dSym`, `coreMatrix`/`coreMatrixD`, `sCore`/`qCore`, `deltaCoreEliminated`/`deltaCoreD`, `sigmaTilde`, `coreShape`, `kappaSlot`/`gammaSlot`/`sigmaSlot`, `coreSchurTarget`/`coreSchurResidual` have NO counterpart in the `.py`.
- **No collateral edits:** the diff touches only the `.wl`; no paper/notes/numbering labels; the `lambdaOut` fingerprint, `lWForward`/`kappa0FromTube` forward law, `sigmaStar`, `deltaCoreExpected`, every `expectZero` label string, the κ₀/γ₀ provenance `Print` lines, and the §6 boolean wiring are preserved (§6 `nontrivialCompensated` now reads `deltaCoreResidual === 0`, a trivial rename of the same residual).

**Independence acceptance (a)–(f):**
(a) No `Series`/`Normal` extraction of §1–§4 ratios survives — `grep` shows the ONLY `Series` is wl:185, the permitted final §5 residual `Normal[Series[deltaCore - deltaCoreExpected, {z,0,5}]]`. ✓
(b) §5 contains a `Solve[Thread[coreMatrix.{sCore,qCore}==coreSource]]` elimination in place of the re-typed `sigmaC`. ✓
(c) Intermediate names/steps differ from the `.py` (list above). ✓
(d) Every `expectZero` prints `= 0` / `PASS` (12 PASS lines in the log). ✓
(e) All printed deliverables byte-identical to the committed transcript (confirmed below). ✓
(f) Script exits 0. ✓

## Exec log assessment

**SymPy:** exit=0. The `.py` is the unchanged reference engine (not in the working-tree diff). Notable lines (`exec_logs/stage_117_sympy.log`): `standalone mixed-pole kappa match = -1/9`; `core-balance surface branches = [g_s*(2*lam - sqrt(...))/(2*K_s), ...]`; `concrete core collapses to the compensated hybrid class = 0`.

**Mathematica:** exit=0 (orchestrator's independent re-run; `exec_logs/stage_117_mathematica.log` `# exit_code: 0`). Notable lines:
- `hybrid canonical-even branches = {{rho -> sigma, kappa -> 0}, {rho -> 4*sigma, kappa -> 1/3}}`
- `standalone mixed-pole kappa match = -1/9` / `sigma match = 0`
- `core-balance surface branches = {{gq -> (gs*lam)/ks - Sqrt[gs^2*kq*ks + gs^2*lam^2]/(2*ks)}, {gq -> (gs*lam)/ks + Sqrt[...]/(2*ks)}}`
- `concrete core collapses to the compensated hybrid class = 0` / `PASS`
- Capstone: exactly one nontrivial survivor `{compensated Robin-mixed core realization, True, True, True, ...}`; ends `Stage 117 Mathematica audit passed.`
The implicit `coreSchurResidual === 0` guard (wl:160) did not `fail` — the elimination reproduces the named reduced forms.

**Output freshness:** The committed `.txt` (`mathematica/output/moving_throat_pde_stage117_outlet_core_status_mathematica_audit.txt`) is NOT in `git diff --name-only HEAD` — it is UNCHANGED in git, as the directive required. A byte-for-byte `diff` of the committed `.txt` against the fresh exec-log payload (header/exit-code lines stripped) is IDENTICAL. The re-authored `.wl` reproduces the exact transcript, so no `.txt` regeneration was needed.

## Material-change assessment

`material_change`: false.

This is a method-only re-author: every emitted value, branch root, σ_*, core-balance branch, classification row, and survivor set is byte-identical to the prior transcript (the eliminated route reaches the same residuals as the re-typed forms). No derived result changes, so no downstream unit that depends on stage-117 outputs is affected. The orchestrator may still mark units > 117 `upstream_stale` per policy, but there is no substantive staling concern here.

## Side observations (non-blocking)

- The card's `\stagefield{Purpose}` "Stage 134" self-label (117+17 pre-renumber drift) noted by the auditor is a cosmetic numbering-band artifact already scoped to the dedicated numbering pass — not a math finding and outside scripts-only scope.
- The two un-exercised card Checks (Schur-complement-sign positivity, parent overlap ratios) remain diagnostic-only reminders, not `\stagefield{Output}` deliverables; unchanged by this re-author and not part of F1.
- `γ₀ = (1+r_c)/9` remains a correctly-labeled postulated ANSATZ (ansatz-catalog item), unchanged by the re-author. The §5 elimination derives κ_c/γ_c structurally as κ₀/(1+r_c), γ₀/(1+r_c); the ansatz status of γ₀ itself is untouched.

## Verdict justification

F1 is fully resolved. The `.wl` no longer transliterates the `.py`: §1–§4 derive the fingerprint jet coefficients from an undetermined-coefficient `Solve` of `den·Y==num` (the only surviving `Series` is the permitted final §5 residual), and §5 derives the core correction by an explicit 2×2 `coreMatrix` `Solve`/elimination with a hard `coreSchurResidual === 0` guard, so `sigmaC` and the collapse are in-engine consequences rather than re-typed literals. Intermediate variable names genuinely differ from the reference. The SymPy reference engine and the committed `.txt` transcript are unchanged in git; all 12 `expectZero` checks print `= 0`/`PASS`; both engines exit 0; all deliverables are byte-identical. This is method-only, so `material_change` is false.
