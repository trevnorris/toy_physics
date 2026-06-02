---
unit_id: 208
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy.md]
  paper_appendix: present
---

# Audit unit 208 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_208.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows at lines 47, 236, 833–876 read; the `\subsection{Pairwise mixed rays}` block 833–876 is the detailed narrative for this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The card's `\stagefield{Output}` states: "Pairwise mixed-ray cone, gradient and off-diagonal Hessian synergy laws, canonical two-ray screen, promotion/deferral rule, and minimal packet for the full ratio optimizer." Stage 208 opens the first genuinely mixed (two-coordinate) search sector. The notes enumerate six deliverables: (1) the exact one-parameter pairwise mixed-ray cone `s_ij(r) = (e_i + r e_j)/√(1+r²)`, r≥0; (2) the gradient-synergy law `k_ij(r) = (k_i + r k_j)/√(1+r²)` with derivative `(k_j − k_i r)/(1+r²)^{3/2}`, unique maximizer `r_grad = k_j/k_i`, and exact maximum `√(k_i²+k_j²) > max(k_i,k_j)`; (3) the off-diagonal curvature-synergy law `H_1,ij(r) = (h_ii + 2r h_ij + r² h_jj)/(1+r²)`, with cross weight `w_×(r) = 2r/(1+r²)` maximized to value 1 at r=1 (the equal-mix ray), plus the diagonal-neutrality law when `h_ij = 0`; (4) the mixed-ray curvature envelopes `κ_lo(r), κ_hi(r)` and the certified-bracket theorem using the carried root map `T(H0,k;c) = 2H0/(k+√(k²−2cH0))` with discriminants `Δ_lo, Δ_hi`; (5) the two canonical screen rays — gradient-optimal `(k_i e_i + k_j e_j)/√(k_i²+k_j²)` with slope `√(k_i²+k_j²)` and Rayleigh-quotient curvature, and equal-mix `(e_i+e_j)/√2` with slope `(k_i+k_j)/√2` and curvature `(h_ii+2h_ij+h_jj)/2` — which coincide iff `k_i = k_j`; (6) the promotion/deferral rule and the minimal data packet `P_ij^mix`. The appendix block (eqs. `app-part06-pairwise-cone` through `app-part06-pairwise-canonical-screens`) reproduces the cone, slope, gradient ratio, curvature, cross weight, and the two canonical screens verbatim, so the .tex, notes, and appendix agree.

## What the script claims to verify

The SymPy script builds `s = [1, r]ᵀ/√(1+r²)`, gradient `g = [−k_i, −k_j]ᵀ`, and forms the oriented slope `K_ij = gᵀs`, then asserts it equals `−(k_i + r k_j)/√(1+r²)` (line 55–58), the slope derivative law (63–66), gradient-optimal stationarity and value (74–76). Section II builds the quadratic form `H1ij = sᵀHs` for `H = [[h_ii,h_ij],[h_ij,h_jj]]` and asserts the weighted decomposition, diagonal neutrality, cross-weight derivative, stationarity at r=1, and max value 1 (90–104). Section III independently constructs the gradient-optimal direction `[k_i,k_j]ᵀ/√(k_i²+k_j²)` and equal-mix `[1,1]ᵀ/√2`, checking they match `s.subs(r,r_grad)` and `s.subs(r,1)`, and verifies their slope/curvature closed forms (112–142). Section IV builds the entrywise-envelope curvatures `κ_lo,κ_hi` and checks the weighted form (149–166). Section V forms discriminants `Δ_lo,Δ_hi` and bracket times `τ_lo,τ_hi` from the carried root map and verifies each `τ` satisfies the closure quadratic `H0 − k τ + ½ κ τ² = 0` (173–195). Section VI verifies the coincidence condition (`r_grad − 1 = (k_j−k_i)/k_i`) and the gradient-optimal cross weight `2 k_i k_j/(k_i²+k_j²)` (202–217).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) cone `s_ij(r)=(e_i+r e_j)/√(1+r²)` | line 44 construction; used throughout | match |
| (2a) slope law `k_ij=(k_i+r k_j)/√(1+r²)` | line 55–58 `mixed slope law` | match |
| (2b) slope derivative `(k_j−k_i r)/(1+r²)^{3/2}` | line 63–66 | match |
| (2c) maximizer `r_grad=k_j/k_i`, stationarity | line 68, 74 | match |
| (2d) max value `√(k_i²+k_j²)` | line 75–76 | match (value verified; the strict `> max(k_i,k_j)` inequality is qualitative, not asserted — acceptable) |
| (3a) curvature `H_1,ij=(h_ii+2r h_ij+r² h_jj)/(1+r²)` | line 96 | match |
| (3b) cross weight `w_×=2r/(1+r²)`, max 1 at r=1 | line 99–104 | match |
| (3c) diagonal-neutrality (h_ij=0) | line 97 | match |
| (4a) envelopes `κ_lo,κ_hi` | line 152–166 | match |
| (4b) root map / discriminants / brackets `τ_lo,τ_hi` | line 174–195 | match |
| (5a) gradient-optimal ray, slope, Rayleigh curvature | line 112–127 | match |
| (5b) equal-mix ray, slope, curvature | line 130–142 | match |
| (5c) coincide iff k_i=k_j | line 202–205 | match |
| (6) promotion/deferral rule, packet `P_ij^mix` | none (qualitative comparison rule / data-bundling statement) | not a symbolic identity; no script check expected |
| ALL of the above on a second engine | none — no `.wl` exists | missing (F1) |

`paper_alignment: aligned` — every symbolic deliverable the paper states is faithfully and non-tautologically exercised by the SymPy script. The sole gap is the absent second engine.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 55–58 | `simplify(K_ij + (k_i+r k_j)/√(1+r²)) == 0` | (2a) slope law | yes |
| A2 | sympy | 63–66 | `simplify(kij' − (k_j−k_i r)/(1+r²)^{3/2}) == 0` | (2b) | yes |
| A3 | sympy | 74 | `kij'.subs(r,r_grad) == 0` | (2c) stationarity | yes |
| A4 | sympy | 75 | `kij_grad − √(k_i²+k_j²) == 0` | (2d) | yes |
| A5 | sympy | 76 | `kij_grad² − k_i² − k_j² == 0` | (2d) value | yes (independent re-confirm of squared max) |
| A6 | sympy | 96 | `H1ij − (w_i h_ii + w_x h_ij + w_j h_jj) == 0` | (3a) | yes |
| A7 | sympy | 97 | `H1ij(h_ij→0) − (w_i h_ii + w_j h_jj) == 0` | (3c) | yes |
| A8 | sympy | 102 | `w_x' − 2(1−r²)/(1+r²)² == 0` | (3b) | yes |
| A9 | sympy | 103–104 | `w_x'(r=1)==0`, `w_x(1)−1==0` | (3b) max | yes |
| A10 | sympy | 122–127 | grad-ray = `s.subs(r,r_grad)`; `K_grad+√(...)==0`; Rayleigh curvature | (5a) | yes |
| A11 | sympy | 140–142 | eq-ray = `s.subs(r,1)`; slope; curvature | (5b) | yes |
| A12 | sympy | 159–166 | envelope weighted forms | (4a) | yes |
| A13 | sympy | 188–195 | `τ_lo,τ_hi` satisfy closure quadratic | (4b) | yes (genuine root-substitution, not tautology) |
| A14 | sympy | 202–205 | `r_grad−1 − (k_j−k_i)/k_i == 0` | (5c) | yes |
| A15 | sympy | 214–217 | `w_x(r_grad) − 2 k_i k_j/(k_i²+k_j²) == 0` | (3b)/(5) comparison | yes |
| — | mathematica | — | (none) | ALL | missing → F1 |

## Findings

### F1 — missing_verification_script

**Severity:** medium
**Files:**
- `(missing)` → target `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage208_pairwise_mixed_rays_and_off_diagonal_hessian_synergy_mathematica_audit.wl`

**What's wrong:**
Stage 208 is `is_status_only_candidate: False` and `is_checkpoint: False`, so the project dual-engine contract requires both a SymPy and a Mathematica audit. Only the SymPy script exists; there is no `.wl`. The paper card itself records the gap: `stage_208.tex:11` states "Mathematica audit: none yet." All of this stage's content is elementary symbolic linear algebra and single-variable calculus — quadratic forms `sᵀHs`, the derivative of `w_×(r)=2r/(1+r²)`, the slope law `(k_i+r k_j)/√(1+r²)` and its derivative, the maximizer `r_grad=k_j/k_i`, the Rayleigh-quotient curvatures, and the certified-bracket root relation. Mathematica can verify every one of these independently with native primitives (`Simplify`/`FullSimplify`, `D`, `Solve`/`Reduce`, matrix ops, `Series`+`Coefficient`). There is no genuine impossibility, so the gap is a finding rather than a justified single-engine carve-out.

**Why this matters:**
A single-engine verification has no cross-check: any transcription error in the SymPy closed forms (e.g., a wrong factor of 2 in the cross weight, or a sign in the slope) would pass undetected because the same engine authored both the construction and the comparison expression. The second engine is the only adversarial guard for these closed forms before Stage 209 consumes the slope/curvature/bracket packet.

**Required change:**
Add an independent Mathematica audit at the target path that re-derives and verifies the claim manifest (M1–M9 in the directive) using a DIFFERENT decomposition than the SymPy script — not a line-by-line port. See directive `stage_208.md`, F1.

**Verification:**
Verifier runs `redteam exec-mathematica 208`; a new `.wl` appears at the target path, contains substantive `expectZero`-style checks for M1–M9, and exits 0 with all checks passing. Spot-check 2–3 sections against the `.py` to confirm independent decomposition (transliteration is rejected).

## Independent-derivation check (Mathematica)

No `.wl` exists, so the transliteration check is moot at audit time. The directive's anti-transliteration guard pre-empts it: the new script must derive each claim via a different route (e.g., build the cross weight from `Coefficient[Expand[(sᵀHs)(1+r²)], h_ij]` rather than re-stating `2r/(1+r²)`; obtain `r_grad` via `Solve[D[k_ij,r]==0, r]` rather than substituting a pre-baked ratio; verify the bracket by `Reduce`-ing the closure quadratic rather than re-substituting a closed-form `τ`).

## Engine cross-check

Not applicable — only one engine present. This is the basis of F1.

## Verdict justification

The SymPy script is strong. I attacked every assertion for tautology and triviality and each survived: the slope, curvature, weight, envelope, and bracket expressions are independently constructed (matrix products `gᵀs`, `sᵀHs`, root-map substitutions) and then compared to the paper's closed forms, so a wrong closed form would make `expect_zero` raise. The derivative checks (A2, A3, A8, A9) operate on expressions that genuinely depend on `r`, so no derivative collapses to identically zero. The bracket check (A13) is a real root substitution into the closure quadratic, not a re-statement of the bracket formula. Paper, notes, and appendix all agree on the formulas, and the script matches them — `paper_alignment: aligned`. The output `.txt` (2026-05-11) is newer than the script (2026-04-14), so no `stale_output`. The only defect is the absent second engine, which the dual-engine contract requires and which Mathematica can genuinely supply. Verdict: `findings` (one finding, `missing_verification_script`/`missing_mathematica`). Not stop-cold: nothing is mathematically irreconcilable and the fix (adding a `.wl`) does not change any downstream-consumed value.

Note: the script's banners say "STAGE 191" (lines 35, 219) and the output header repeats it; this is a cosmetic label drift inside print strings only — it does not affect any assertion, file path, or the verified math, so it is recorded here but not raised as a finding (the audit contract corrects math defects, not print-string cosmetics).

## Self-test notes

I walked the four relevant traps. (1) Variable independence: every `diff(·, r)` in the script and every `D[·,r]` I prescribe in the manifest acts on an expression that actually contains `r` (`k_ij`, `w_×`), so no derivative is identically zero and no `assert_nonzero`/`assert_zero` passes trivially. (2) Parity/symmetry: no unbounded-domain integrals are involved (this stage is algebraic), so the parity trap does not apply. (3) Trivial-case pre-check: substituting `k_i=k_j` makes `r_grad=1` and the gradient-optimal and equal-mix rays coincide (M6), and substituting `h_ij=0` collapses the curvature to the diagonal interpolation (M4) — both reduce as the manifest claims. (4) Path spec: the directive names the full `.wl` target under `mathematica/` and the manifest is a requirement list, not an implementation route, so Codex designs the independent derivation.
