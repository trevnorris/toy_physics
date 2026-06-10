---
unit_id: 244
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-10T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 244 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_244.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_mathematica_audit.txt`

## What the paper claims

Stage 244 compiles the open-system leakage source `S_leak` and the scalar-photon work channel `J^w E_w` onto the Stage-242 selected-support packet, using an odd Gaussian one-mode closure `W_λ`, `φ_λ`. The card's body equations give: the exact leakage compiler `S_leak = √2 μ_w ρ₀/(2√π λ³) · E₀`; the bulk work `W_w^bulk = √2 μ_w q ρ₀/(2√π λ³) · E₀²` with `W_w^bulk = q E₀ S_leak`; the Session-I scalar `W_w^sess = 2μ_w q ρ₀/λ² · E₀² = 2√(2π)λ W_w^bulk = 4πqλ⁴/(μ_w ρ₀) · S_leak²`; the pullback `Π_tr = 32Λ(1−ε)/(3π²) = 16Λϱ/π²`, `E₀ = η_leak Π_tr`; the support/orbit split `∂_{R_tr}=∂_{R_target}=∂_{ε_η}=0`; and orientation parity (`S_leak` odd, both work scalars even in η_leak) plus recovery at η_leak=0. The appendix (lines 185–194) states the selected-branch results `S_leak = 8√2 η_leak μ_w ρ₀/(π^{5/2}λ³)·Λϱ` and `W_w^sess = 512 η_leak² μ_w q ρ₀/(π⁴λ²)·Λ²ϱ²`. No explicit `\stagefield{Output}` line; the card's `\stagefield{Verification}` line 4 states "Mathematica audit: none yet."

## What the script claims to verify

The SymPy script derives `S_leak` and `W_bulk` from scratch by integrating `W'(w)·j^w` and `J^w·E_w` over the real line (genuine integrals, not literals), checks the boundary term vanishes, confirms `W_bulk = q E₀ S_leak`, `W_sess = 2√(2π)λ W_bulk` and the quadratic law, derives `Π_tr` from `C_mix`, substitutes `E₀ = η_leak Π_tr`, then verifies every compiled `(1−ε)` and `ϱ` form against closed-form expected expressions, the support/orbit split via free-symbol disjointness with a positive support-coverage control, η-parity via sign substitution, and recovery at η=0. The `.wl` mirrors the same nine claim families (M1–M9) with independent Integrate / FullSimplify / FreeQ machinery.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `S_leak = √2 μ_w ρ₀/(2√π λ³)E₀` | py 44–56 / wl M2 | match |
| `W_bulk` closed form + `=qE₀S_leak` | py 62–77 / wl M3,M4 | match |
| `W_sess` (3 forms) | py 64–79 / wl M5 | match |
| `Π_tr` both forms, `E₀=η_leak Π_tr` | py 87–119 / wl M6 | match |
| Compiled `(1−ε)`+`ϱ` forms (incl. 128√2) | py 103–125 / wl M6 | match |
| Support/orbit split `∂=0` | py 135–143 / wl M7 (free-symbol form) | match |
| Orientation parity | py 161–172 / wl M8 | match |
| Recovery η=0 | py 145–155 / wl M9 | match |
| Card line 4 "Mathematica audit: none yet" | `.wl` present and passes | mismatch |

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 55 | `boundary == 0` | boundary vanishing | yes |
| A2 | sympy | 56 | `simplify(S_leak - S_expected)==0` | leakage compiler | yes |
| A3 | sympy | 76–79 | bulk/sess vs derived | work scalars | yes |
| A4 | sympy | 118–125 | pullback forms | Π_tr + compiled forms | yes |
| A5 | sympy | 141 | `orbit_syms.isdisjoint(free)` | support/orbit split | yes (w/ A6 control) |
| A6 | sympy | 143 | `support_syms.issubset(free)` | positive control (anti-vacuous) | yes |
| A7 | sympy | 170–172 | parity via subs | orientation parity | yes |
| A8 | sympy | 153–155 | recovery η=0 | recovery slice | yes |
| M1–M9 | mathematica | 98–236 | expectZero/expectTrue | same families | yes |

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** paper_missing_script_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_244.tex:4` quote: "SymPy audit: \StageFile{...sympy_audit.py}.  Mathematica audit: none yet."
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_mathematica_audit.wl` (exists, all M1–M9 pass per saved output)

**What's wrong:**
Pass-1 added an independent `.wl` (M1–M9 all PASS in the committed output), but the stage card's `\stagefield{Verification}` still says "Mathematica audit: none yet." The card is stale relative to the now-present second engine.

**Why this matters:**
The card under-reports verification coverage; a reader would believe the stage is single-engine. Cosmetic/prose only — no math is affected.

**Required change:**
Route to user (paper-side edit). Likely fix: update line 4 to cite the `.wl` path, mirroring how other dual-engine cards read.

**Verification:**
Card line 4 names the `mathematica/...stage244..._mathematica_audit.wl` file.

## Independent-derivation check (Mathematica)

Genuinely independent, not a transliteration. Both engines re-derive `S_leak`/`W_bulk` from the integrals `∫W'(w)j^w dw` and `∫J^w E_w dw` rather than echoing pre-baked algebra: SymPy uses `sp.integrate(sp.diff(W,w)*j_w, ...)` (py 44); Mathematica uses `Integrate[D[projector,w] currentW, ...]` (wl 75–83). The support/orbit split is implemented with different primitives — SymPy set algebra `orbit_syms.isdisjoint(free)` + `support_syms.issubset(free)` (py 141/143) vs Mathematica `FreeQ`/`Not[FreeQ]` (wl 187–215). Pullback substitutions and parity also use engine-native idioms. No line-by-line porting.

## Engine cross-check

Both outputs agree: `S_leak = √2 E₀ μ_w ρ₀/(2√π λ³)` (py txt line 10) ≡ wl M2=0; `W_bulk`/`W_sess` identical; compiled `(1−ε)` forms match; all M1–M9 PASS; SymPy all asserts pass. No residual disagreement.

## Verdict justification

The math holds under attack in both engines, with genuine integrals (no hardcoded results) and a real, load-bearing anti-vacuity control on the one place that could have been a self-test trap. The only finding is a stale prose Verification line on the card (a `.wl` now exists but the card says "none yet") — a low-severity `paper_missing_script_claim` routed to the user. Verdict: findings.

## Self-test notes

Self-test trap (F1 in pass-1): CONFIRMED fixed and load-bearing. The support/orbit split uses NO `sp.diff`/`D[]` w.r.t. orbit variables (which would be the vacuous trap); instead it asserts `orbit_syms.isdisjoint(free)` (py 141) guarded by the positive control `support_syms.issubset(free)` (py 143) — the output shows non-empty support coverage `[Lam, eta_leak, varrho]` (txt 39), so the disjointness is non-vacuous. The Mathematica M7 mirrors this with `Not[FreeQ[...,Lam]]` etc. The only `sp.diff`/`D[]` calls are `sp.diff(W,w)` (py 44) and `D[projector,w]` (wl 75): `w` genuinely appears in the Gaussian projector, so the derivative is non-zero and load-bearing. Parity uses sign substitution, not differentiation. Symmetry: odd integrand `W'(w)·j^w` is even×even (the Gaussian × derivative product) integrating to a nonzero, and `J^w E_w` even — both correctly nonzero; recovery at η=0 correctly zeroes the compiled forms.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `S_leak = √2μ_wρ₀/(2√π λ³)E₀` | py txt 10 / wl M2 | notes 193, tex 41 | MATCH |
| `W_bulk = √2μ_w qρ₀/(2√π λ³)E₀²` | py txt 16 / wl M3 | notes 222, tex 50 | MATCH |
| `W_sess = 2μ_w qρ₀/λ² E₀²` | py txt 19 / wl M5 | notes 245, tex 62 | MATCH |
| `Π_tr = 32Λ(1−ε)/(3π²)` | py txt 27 / wl M6 | notes 295, tex 84 | MATCH |
| `Π_tr = 16Λϱ/π²` | py txt 28 / wl M6 | notes 309, tex 87 | MATCH |
| `S_leak(ϱ) = 8√2 η μ_wρ₀/(π^{5/2}λ³)Λϱ` | py 104 / txt 53 | notes 342, appendix 187 | MATCH |
| `W_bulk(ε) = 512√2 η².../(9π^{9/2}λ³)` | py 105 / txt 31 | notes 355 | MATCH |
| `W_bulk(ϱ) = 128√2 η²μ_w qρ₀/(π^{9/2}λ³)Λ²ϱ²` | py 106 / txt 55 | notes 366 | MATCH (pass-1 196√2→128√2 fix present) |
| `W_sess(ε) = 2048 η².../(9π⁴λ²)` | py 107 / txt 32 | notes 379, tex 101 | MATCH |
| `W_sess(ϱ) = 512 η²μ_w qρ₀/(π⁴λ²)Λ²ϱ²` | py 108 / txt 56 | notes 390, appendix 194 | MATCH |

Internal scaffolding (no finding): boundary term =0, free-symbol overlap/coverage sets, recovery zeros, parity residuals, PASS flags.

reconciliation: complete; 10 deliverable values checked, 0 misaligned. The flagged `128√2` reconciles (notes line 366 now reads `128√2`, matching script line 106 and output line 55; the pass-1 `196√2→128√2` correction is present and consistent across the (1−ε) and ϱ forms).
