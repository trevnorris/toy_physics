---
unit_id: 239
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 239 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_239.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (row line 90; theorem block lines 1101-1174)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage239_rigid_mouth_physical_normal_form_exact_physical_to_microscopic_correction_compiler_and_cartesian_orbit_lock_theorem_mathematica_audit.txt`

## What the paper claims

The stage card `\stagefield{Output}` reads: "Cartesian orbit-lock theorem: rigid-mouth orbit lock is exactly $U=V=0$, equivalently $\mathcal T^2=\mathcal T_{\rm ref}^2$ and $\epsilon_\eta=\epsilon_{\eta,\rm ref}$." The `\stagefield{Derivation ledger}` states the stage "diagonalizes the packet in $(U,V)$, derives the exact microscopic correction $(\Delta_T,\Delta_{K_\eta},\Delta_\mu)=(0,-V,U-V)$, and identifies the Cartesian orbit-lock point." The notes enumerate the full deliverable set: (1) the diagonal physical log chart $U=\ln(\mathcal T^2/\mathcal T_{\rm ref}^2)$, $V=\ln(\epsilon_\eta/\epsilon_{\eta,\rm ref})$ with $M_{\rm phys}=I_2$; (2) the target-ratio factorization $R/R_{\rm ref}=\frac{1-\epsilon_{\rm ref}e^V}{1-\epsilon_{\rm ref}}e^{-U}$; (3) complementary projectors $P_{\mathcal T},P_\eta$ and the commuting transfer/dressing legs; (4) the compiler $C^{\rm dep}_{\rm phys}=\begin{psmallmatrix}0&0\\0&-1\\1&-1\end{psmallmatrix}$ with left inverse $L^{\rm dep}_{\rm phys}$; (5) the pure transfer image $(0,0,U)$ and dressing image $-V(0,1,1)$; (6) static/post-static/full correction packets; (7) support-blindness in $\zeta,M_{\rm mix}$; (8) the Cartesian orbit-lock theorem in finite and first-order form. Appendix part07 row (line 90) reproduces the chart, compiler, and orbit-lock condition verbatim, all marked `\StatusExactClosure{}` (open only on whether the actual branch realizes $U=V=0$).

## What the script claims to verify

The docstring enumerates 8 items matching the notes 1:1. The SymPy assertions exercise: the chart exponentials and diagonal map; the target-ratio formula reconstructed two ways from the selected-branch product identity $R\,\mathcal T^2=\Lambda_0(1-\epsilon_\eta)$; projector idempotence/orthogonality/completeness and the leg factorization; the inherited compiler matrix and its left inverse (with an independent `sp.solve` of the dependent-plane inversion); the two axis images; the three correction packets and their split; support-blindness propagated through the chart and all correction packets via a chain-rule substitution of vanishing derivatives; and the orbit-lock equivalence via `sp.solve` plus first-order perturbation. The Mathematica script verifies the same deliverable set (labels M1-M9).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Diagonal chart $U,V$, $M_{\rm phys}=I_2$ | py L71-76 / wl M1 (L98-100) | match |
| Target-ratio factorization | py L84-98 / wl M2 (L114-115) | match |
| Projectors + commuting legs | py L103-128 / wl M3 (L126-131) | match |
| Compiler $(0,-V,U-V)$ + left inverse | py L174-215 / wl M4-M5 (L146-154) | match |
| Pure transfer/dressing images | py L220-225 / wl M6 (L183-185) | match |
| Static/post-static/full corrections | py L230-241 / wl M7 (L191-195) | match |
| Support-blindness ($\zeta,M_{\rm mix}$) | py L248-305 / wl M8 (L219-225) | match |
| Cartesian orbit-lock + first-order | py L310-352 / wl M9 (L229-253) | match |

All paper deliverables are covered; no `extra` script-side checks beyond paper scope (the Stage 238 branch-formula carry, py L133-191 / wl section III, corroborates the inherited inputs and is referenced in the notes §0/§3).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 76 | `assert_matrix_zero(q_phys - [U,V])` | chart diagonal | yes |
| A2 | sympy | 84-98 | `assert_zero(ratio reconstructions)` | target ratio | yes |
| A3 | sympy | 106-115 | projector algebra `== 0` | projectors/legs | yes |
| A4 | sympy | 177-181 | `y_dep - [0,-V,U-V] == 0` | compiler | yes |
| A5 | sympy | 197-212 | `solve` inversion + `L·C - I == 0` | left inverse | yes |
| A6 | sympy | 310-316 | `solve(y_dep=0) == {U:0,V:0}` | orbit lock | yes |
| A7 | sympy | 262-305 | chain-rule support-blind `== 0` | support-blindness | yes |
| A8 | sympy | 341-348 | first-order compiler `== 0` | first-order form | yes |
| M4 | mathematica | 140-146 | `compilerJacobian = Table[D[...]]`; `J·UV - delta == 0` | compiler (DERIVED via Jacobian) | yes |
| M5 | mathematica | 150-154 | `PseudoInverse[J]`; `L·J - I == 0` | left inverse (DERIVED via PseudoInverse) | yes |
| M9 | mathematica | 229-235 | `Reduce[...]`; `Equivalent[set, U==0&&V==0]` | orbit lock (Reduce/Equivalent) | yes |

## Findings

None. The script holds up against every attack tried (see Verdict justification and Self-test notes).

## Independent-derivation check (Mathematica)

This is the checkpoint's load-bearing question. Pass-1 re-authored the `.wl` away from a transliteration; I confirm the re-author produced **genuinely different information flows** for the three load-bearing operations:

1. **Compiler matrix.** SymPy POSITS the 3×2 matrix as literal entries — `S_rm_dep = sp.Matrix([[0, 0], [0, -1], [1, -1]])` (py L174), `C_phys_dep = S_rm_dep * M_phys` — then asserts `C_phys_dep - sp.Matrix([[0,0],[0,-1],[1,-1]]) == 0` (py L177, self-referential by construction). Mathematica instead writes the boxed *dependent vector* `dependentDelta = {0, -physicalCoordinates[[2]], physicalCoordinates[[1]] - physicalCoordinates[[2]]}` (wl L135-139) and DERIVES the matrix by `compilerJacobian = Table[D[dependentDelta[[row]], physicalCoordinates[[col]]], ...]` (wl L140-144), then checks `compilerJacobian . physicalCoordinates - dependentDelta == 0` (wl L146). One POSITS the operator and confirms it echoes itself; the other DIFFERENTIATES the closed-form correction to recover the operator and confirms the operator regenerates the (linear) vector. Different flow.

2. **Left inverse.** SymPy POSITS `L_phys_dep = sp.Matrix([[0, -1, 1], [0, -1, 0]])` (py L206) and checks `L·C - I == 0` (py L212). Mathematica DERIVES it natively: `leftNative = PseudoInverse[compilerJacobian]` (wl L150), then checks `leftNative . compilerJacobian - IdentityMatrix[2] == 0` (wl L152) AND `leftNative . deltaSymbols - {DeltaMu - DeltaKeta, -DeltaKeta} == 0` (wl L153) — i.e. the Moore-Penrose inverse independently reproduces the paper's left-inverse coordinate formulas $U=\Delta_\mu-\Delta_{K_\eta}$, $V=-\Delta_{K_\eta}$. The `.py` cross-checks the same formulas via an algebra-system `sp.solve` (py L197-202), a third independent route. Genuinely different.

3. **Orbit lock.** SymPy uses `sp.solve([y_dep[1]=0, y_dep[2]=0], [U,V])` and compares the returned dict to `[{U:0, V:0}]` (py L310-316). Mathematica uses `Reduce[dependentDelta[[2]]==0 && dependentDelta[[3]]==0, {U,V}, Reals]` then `Equivalent[orbitLockReduced, U==0 && V==0]` (wl L229-235). Different solver semantics (Reduce over Reals + logical Equivalent vs. solve + structural dict compare).

These are not the same autodiff/solve dressed twice: where the `.py` hardcodes a matrix and a left-inverse, the `.wl` derives both from the boxed correction vector. No `mathematica_transliteration` finding.

## Engine cross-check

Both saved outputs are present and report all-pass. SymPy `.txt` (52 `[ok]` lines + "All Stage 239 symbolic checks passed."); Mathematica `.txt` (every M1-M9 block "PASS", `M9 Reduce orbit-lock set = U == 0 && V == 0`, "All Stage 239 symbolic checks passed."). The two engines agree on every shared quantity: chart map zero, target-ratio residual zero, compiler/left-inverse identities, all correction packets, support gradients, the orbit-lock set `{U=0,V=0}`, and the first-order dependent vector `{0,-dlnEps,dlnT2-dlnEps}`. No `engine_disagreement`.

## Verdict justification

Clean. Both engines are present and substantive; the checkpoint higher bar is MET. Attacks tried that failed: (a) **tautology hunt** — the `.py` does have self-referential checks (L76 `q_phys` vs `[U,V]` with `M_phys=I`; L177 `C_phys_dep` vs the same literal), but each of those load-bearing quantities is independently re-derived on the `.wl` side (Jacobian) and/or via a second `.py` route (`sp.solve` at L197), so the *claim* is non-tautologically established across the unit, which is the checkpoint standard. (b) **variable-independence self-test trap** — the support-blindness checks differentiate `Log[T2sbFn[zeta,Mmix]/T2ref]` w.r.t. `zeta,Mmix`, and `T2sbFn`/`epsEtaSbFn` genuinely depend on those vars, so the derivatives are non-trivially formed before the `supportRules` substitution sends the inner derivatives to 0 — not an identically-zero `D` (py L248-257/wl L200-209). (c) **branch/assumption check** — the `.wl` `$Assumptions` (`0<epsEta<1`, `eps!=1`, all positives) match the physical setup and the `.py` symbol domains; nothing is over-assumed to force a pass. (d) **paper alignment** — every `\stagefield{Output}` and notes deliverable maps to a script check; the appendix row (line 90) and theorem (line 1169) match the compiler, chart, and orbit-lock condition verbatim. Banner is canonical: `.wl` L65 reads `STAGE 239 — RIGID-MOUTH PHYSICAL NORMAL FORM`, and the carried-branch labels correctly say `Stage238`/`Stage236` (no stale `STAGE 222`/`Stage221` residue). Outputs are fresh (both `.txt` mtimes 2026-06-03 08:42, newer than `.py` 2026-05-11 and `.wl` 2026-06-03 08:36).

## Value Reconciliation (pass-2 augmentation)

This stage emits only **symbolic** boxed results (no numeric constants — there are no figures-of-merit, no pinned decimals). Reconciling each deliverable symbol against the card/notes/appendix:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Chart $U=\ln(\mathcal T^2/\mathcal T_{\rm ref}^2)$, $V=\ln(\epsilon_\eta/\epsilon_{\eta,\rm ref})$ | py L68-69 / wl L92-93 | tex L9; notes L93-97; appx L1104-1106 | MATCH |
| $M_{\rm phys}=I_2$ (diagonal map) | py L74 / wl L94 | notes L123 | MATCH |
| Target ratio $\frac{1-\epsilon_{\rm ref}e^V}{1-\epsilon_{\rm ref}}e^{-U}$ | py L83 / wl L108 | notes L137-141; appx L1118-1120 | MATCH |
| Compiler $(\Delta_T,\Delta_{K_\eta},\Delta_\mu)=(0,-V,U-V)$ | py L179 / wl L135-139; sy.txt L20, wl.txt L41 | tex L13; notes L239-253; appx L90,L1125-1127 | MATCH |
| Compiler matrix $C^{\rm dep}_{\rm phys}=[[0,0],[0,-1],[1,-1]]$ | py L174/L177 / wl L146 | notes L257-268 | MATCH |
| Left inverse $U=\Delta_\mu-\Delta_{K_\eta},\,V=-\Delta_{K_\eta}$ | py L202/L206 / wl L153 | notes L285-289; appx L1132-1134 | MATCH |
| Transfer image $(0,0,U)$; dressing image $-V(0,1,1)$ | py L223-224 / wl L183-184 | notes L307-333; appx L1141,L1148 | MATCH |
| Static/post-static/orbit packets $(0,0,-U)$, $(0,V,V)$, $(0,V,V-U)$ | py L234-237 / wl L191-193 | notes L349-398 | MATCH |
| Orbit-lock condition $U=V=0$ | py L315 / wl L234-235; wl.txt L126 | tex L15; notes L457-462; appx L90,L1153-1155,L1171 | MATCH |
| First-order vector $(0,-\delta\ln\epsilon_\eta,\delta\ln\mathcal T^2-\delta\ln\epsilon_\eta)$ | py L346 / wl L253 | notes L489-497; appx L1162-1166 | MATCH |

INTERNAL scaffolding (no finding): `assert_zero`/`assert_matrix_zero` flags, `expectZero`/`expectTrue`, `compilerJacobian . Mphys - compilerJacobian` identity-passthrough check, `Stage 238` branch-formula intermediates (`T2_stage238`, `Rtarget_stage238`, `q_nt_rigid_stage238`) which exist to corroborate the carried inputs, residual zeros, pass/fail prints.

reconciliation: complete; 10 deliverable values checked, 0 misaligned

## Self-test notes

Checked the three traps from the prompt. (1) Variable independence: the only `D[...]`/`sp.diff` over physical variables are the `compilerJacobian` (vector linear in U,V — Jacobian well-defined, non-degenerate columns) and the support-blind derivatives (integrand genuinely depends on `zeta,Mmix` via `T2sbFn`/`epsEtaSbFn` before the vanishing-derivative substitution), so no identically-zero derivative producing a silent pass. (2) No unbounded-domain integrals in this stage — parity trap N/A. (3) Trivial-case: substituting `U=V=0` collapses `dependentDelta` to `{0,0,0}` and the target ratio to 1, matching the asserted orbit-lock checks. Conclusion: no directive needed; verdict clean.
