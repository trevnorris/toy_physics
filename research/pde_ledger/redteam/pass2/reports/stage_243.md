---
unit_id: 243
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-10T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 243 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_243.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (rows at lines 84, 122-179 read)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage243_relaxed_constraint_branch_declaration_and_short_range_open_system_compiler_mathematica_audit.txt`

## What the paper claims

Stage 243 declares (does not yet assemble) the relaxed-constraint branch: a codimension-three lift of the Stage-242 front-end packet that legalizes three Session-I relaxations (transverse-current `J^w`, rigid-mouth `U=V=0`, positive-source closure) while preserving the same-charge short-range verdict. Quoted `\stagefield`-equivalent outputs / boxed deliverables: (1) leakage/work lane `S_{\rm leak}=-\frac{\sqrt2}{4}\ell_w j_0`, `\mathcal W_w=\frac{\sqrt{2\pi}}{8}\ell_w j_0 E_0`, both vanishing on `\ell_w=0`; (2) non-rigid response `U=\frac{k_Vf_U}{k_Uk_V-\chi_\lambda^2}`, `V=\frac{\chi_\lambda f_U}{k_Uk_V-\chi_\lambda^2}`, `V/U=\chi_\lambda/k_V`, admissibility `k_Uk_V-\chi_\lambda^2>0`, drain `\mathcal D_{UV}=\frac{\chi_\lambda^2k_Vf_U^2}{(k_Uk_V-\chi_\lambda^2)^2}\ge0`; (3) compensated source `\varsigma(z)=1+a\cos\pi z+b\cos2\pi z` with unit mean, quadratic rewrite `\varsigma(y)=1-b+ay+2by^2`, vertex `y_*=-a/4b`, `\varsigma(y_*)=1-b-a^2/8b`, boundary candidates `1\pm a+b`; (4) exact recovery slice `(\ell_w,f_U,a,b)=(0,0,0,0)`; (5) the load-bearing short-range theorem — static `\delta V_{\rm stat}` and dynamic `\Re\mathfrak V_{\rm dyn}` both lie in span`\{x^{-6}, e^{-2\kappa x}x^{-4}, e^{-4\kappa x}x^{-2}\}`, so `\lim_{x\to\infty}x\,\delta V_{\rm stat}=\lim x\,\Re\mathfrak V_{\rm dyn}=0`. Appendix row 84 + main theorem (lines 150-179) restate the same recovery map and short-range span. Checkpoint.

## What the script claims to verify

Both scripts verify, section by section: the exact Gaussian/odd-profile leakage and work integrals (with vanishing boundary term and the `\ell_w=0` reductions); the non-rigid `(U,V)` stationary solution, det-Hessian, ratio, and non-negative drain; the compensated source's exact unit mean, quadratic rewrite, interior vertex/value, and boundary values; the codimension-three recovery slice driving all lifted observables to their strict values; and the short-range kernel products plus their `x\to\infty` `1/x`-weighted limits vanishing for both the static and the linear-dynamic bundles. The `.wl` adds independent cross-routes (IBP closure, half-line work symmetry, `LinearSolve`, and `Series` asymptotics) absent from the `.py`.

## Paper ↔ script cross-check

| paper deliverable | script check | status |
|---|---|---|
| `S_leak = -√2/4 ℓ_w j0` | py:40,56 / wl:81,102 | match |
| `W_w = √(2π)/8 ℓ_w j0 E0` | py:41,57 / wl:85,104 | match |
| both vanish on ℓ_w=0 | py:58-59 / wl:105-106 | match |
| `U,V` solution | py:76-80,98-99 / wl:112-114,131-132 | match |
| `V/U = χ/k_V` | py:101 / wl:134 | match |
| admissibility det = k_U k_V − χ² | py:100 / wl:133 | match |
| drain `D_UV ≥ 0` form | py:102 / wl:135 | match |
| `∫ς dz = 1` | py:136 / wl:154 | match |
| quadratic rewrite `1−b+ay+2by²` | py:138 / wl:156-159 | match |
| vertex `y_*=−a/4b`, `ς(y_*)=1−b−a²/8b` | py:139 / wl:160-162 | match |
| boundary candidates `1±a+b` | wl:163-164 (py prints only) | match (wl covers) |
| recovery slice `(ℓ_w,f_U,a,b)=0` | py:145-172 / wl:168-181 | match |
| short-range span products | py:213-215 / wl:239-241 | match |
| `lim x δV_stat = lim x ReV_dyn = 0` | py:219-220 / wl:245-251 | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 55-57 | integral = closed form | leakage/work | yes |
| A2 | sympy | 58-59 | subs(ell_w,0)==0 | ℓ_w=0 reduction | yes |
| A3 | sympy | 98-106 | solve vs expected, det, ratio, drain, limits | non-rigid lane | yes |
| A4 | sympy | 136-139 | mean/quadratic/vertex | compensated source | yes |
| A5 | sympy | 167-172 | recovery subs ==0 | recovery slice | yes |
| A6 | sympy | 213-220 | products + x→∞ limits | short-range theorem | yes |
| B1 | math | 100-106 | boundary/IBP/leak/work/half-line | leakage/work (2 routes) | yes |
| B2 | math | 131-139 | Solve∧LinearSolve agree, det, ratio, drain | non-rigid lane (2 routes) | yes |
| B3 | math | 154-164 | mean/rewrite/vertex/boundary | compensated source | yes |
| B4 | math | 176-181 | recovery ==0 | recovery slice | yes |
| B5 | math | 239-251 | products + Limit ∧ Series limits | short-range theorem (2 routes) | yes |

All rows non-tautological and anchored. No D[]/diff over an absent variable.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is a genuine independent re-derivation, not a transliteration. Where the `.py` uses a single operation, the `.wl` adds a second, mechanically distinct route and cross-checks it:
- **Leakage:** `.py` computes `S_leak=∫W' j^w` directly (py:40). `.wl` computes the same integral (wl:81-84) but ALSO computes `ibpInterior=∫W j^{w'}` (wl:89-92) and asserts the integration-by-parts identity `Sleak+ibpInterior-boundary==0` (wl:101) — an information flow absent from the `.py`. It further cross-checks the work integral by half-line symmetry `2∫_0^∞ j^w E_w` (wl:93-96, asserted wl:103).
- **U/V solve:** `.py` uses `sp.solve` only (py:76). `.wl` solves both via `Solve` (wl:112) AND `LinearSolve[{{kU,-chiLam},{-chiLam,kV}},{fU,0}]` (wl:115-118), then asserts the two routes agree (wl:132). One derives the solution from gradient roots; the other from the linear stationarity system — different operations.
- **Short-range limits:** `.py` uses `sp.limit` (py:197-198). `.wl` uses `Limit` (wl:205-212) AND an independent `Series[x·expr,{x,Infinity,0}]` asymptotic expansion (wl:213-237, asserted wl:247-251). The series mechanism is a distinct algorithm, and prudently hedges Mathematica's known non-deterministic pole `Limit`.

No section has the `.wl` mirroring the `.py`'s core operation as the sole route. The pass-1 re-author is SUFFICIENT (unlike V.3-200): the load-bearing claims are each derived by a genuinely different mechanism in the second engine.

## Engine cross-check

Both saved outputs report all checks PASS. The `.py` prints the closed forms verbatim (e.g. `S_leak = -sqrt(2)*ell_w*j0/4`, `W_w = sqrt(2)*sqrt(pi)*E0*ell_w*j0/8`, `U = -f_U*k_V/(chi_lam**2 - k_U*k_V)` = paper form with denominator sign flipped both ways, `QQ=x**-6`, `lim x*deltaV_stat=0`). The `.wl` reports each `expectZero` residual as `0` with `PASS`. The two engines agree on every shared deliverable. The `.wl`'s `recovered varsigma = 0` (output line 81) is the residual `varsigmaRec-1`, equivalent to the `.py`'s `Recovered varsigma = 1` — convention difference, not disagreement.

## Verdict justification

Clean. I read the paper card, the full notes, and the appendix rows before the scripts. Every paper-side deliverable maps to a non-tautological, well-anchored assertion in BOTH engines, with values matching exactly. Attacks tried and failed: (1) variable-independence trap — every `subs`/`/.` reduction (`ell_w→0`, `f_U→0`, `chi→0`, `a=b→0`) acts on a symbol the expression genuinely depends on, so none is vacuous; (2) hardcoded-result trap — the `.wl`'s only literal form (`varsigmaY` at wl:148) is validated functionally against the trig original (wl:156-159), not asserted against itself; (3) transliteration trap — the `.wl` derives the three load-bearing results by mechanically distinct second routes; (4) symbol-domain — positivity (`k_U,k_V,χ,f_U,x,κ>0`) matches the paper's stability/short-range setup, and the boundary-term limit is taken under the safe local assumption block (wl:77-80); (5) sign/factor audit on `√2/4`, `√(2π)/8`, `−a/4b`, `−a²/8b`, kernel exponents `−2κ,−4κ` and powers `6,4,2` — all match paper. Checkpoint bar MET: both engines present and substantive, paper alignment exact, load-bearing short-range limits and leakage scalars re-derived (not trusted as literals). Banner canonical (STAGE 243 in both `.wl` source line 59 and output line 3). Outputs fresh (both `.txt` newer than their scripts).

## Self-test notes

Checked the three required traps: (1) Variable independence — listed dependencies for every `diff`/`subs`; `D[varsigmaY,y]/.y->yStar` (wl:161) and all recovery substitutions act on present variables, no identically-zero derivative. (2) Symmetry/parity — the leakage integrand `W'·j^w` is (odd·odd)=even on ℝ → nonzero `S_leak` is legitimate; the work integrand `j^w·E_w`=(odd·odd)=even → nonzero `W_w` legitimate and matches the half-line cross-check. (3) Trivial-case — recovery slice substitutions reduce each lifted observable to literal 0 / ς≡1 as asserted. No directive written (zero script-side findings).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 13 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `S_leak = -√2/4 ℓ_w j0` | py:43 out:9 / wl:97 | tex:70-72; md:206-213 | MATCH |
| `W_w = √(2π)/8 ℓ_w j0 E0` | py:44 out:11 / wl:98 | tex:70-73; md:226-233 | MATCH |
| `U = k_V f_U/(k_U k_V−χ²)` | py:79 out:19 / wl:113 | tex:94; md:285-290 | MATCH |
| `V = χ f_U/(k_U k_V−χ²)` | py:80 out:21 / wl:114 | tex:96; md:285-290 | MATCH |
| `V/U = χ/k_V` | py:83 out:24 / wl:124 | tex:98; md:301-304 | MATCH |
| det H = `k_U k_V − χ²` (admissibility >0) | py:82 out:23 / wl:123 | tex:103-104; md:317,323 | MATCH |
| `D_UV = χ²k_V f_U²/(k_U k_V−χ²)²` | py:85 out:25 / wl:126 | tex:108-110; md:329-334 | MATCH |
| `∫_0^1 ς dz = 1` | py:117 out:32 / wl:144 | tex:124; md:365-368 | MATCH |
| quadratic `ς(y)=1−b+ay+2by²` | py:118 out:33 / wl:148 | tex:128-129; md:388-391 | MATCH |
| `y_* = −a/(4b)` | py:126 out:34 / wl:149 | tex:135; md:398 | MATCH |
| `ς(y_*) = 1−b−a²/(8b)` | py:127 out:35 / wl:150 | tex:136; md:405 | MATCH |
| boundary values `1±a+b` | out:36 / wl:151-152 | tex:138; md:413-417 | MATCH |
| short-range span products `x^-6, e^{-2κx}/x^4, e^{-4κx}/x^2` + `lim x·δV=lim x·ReV_dyn=0` | py:188-198 out:54-63 / wl:191-237 | tex:170-187; md:545-595 | MATCH |

INTERNAL (scaffolding, no finding expected): boundary term =0, IBP-closure residual, half-line work residual, Solve-vs-LinearSolve residual, vertex-stationarity residual, recovery-slice residuals, series cross-route residuals, all PASS flags. The recovery map `(ℓ_w,f_U,a,b)=(0,0,0,0)` is itself a deliverable and MATCHES tex:157-162 / md:480-486 / appendix:84.

All emitted deliverable values reconcile against `.tex` and/or `.md`. No `paper_misalignment`.
