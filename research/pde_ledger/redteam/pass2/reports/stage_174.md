---
unit_id: 174
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage174_static_self_similarity.md]
  paper_appendix: present
---

# Audit unit 174 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_174.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage174_static_self_similarity.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 79, 194, 265, 560-609, 715, 1464, 1478, 1525)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage174_static_self_similarity_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage174_static_self_similarity_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage174_static_self_similarity_mathematica_audit.txt`

## What the paper claims

Stage 174 ("Static self-similarity") rewrites the Stage-173 remaining linear grouped 2.5PN loading defect `Xi_load := N01/N0 - D01/D0` as a *weighted failure of static self-similarity* relative to the wall baseline slope. The card's `\stagefield{Output}` reads: "Expresses \(\Xi_{\rm load}\) as weighted failure of static self-similarity between wall, BdG, conservative port, and outgoing port slopes." The notes/appendix make the deliverables explicit: (1) the exact static-slope decomposition `δ_D = ω_K δ_K − ω_B δ_B − ω_Z δ_Z` with the weight identity `ω_K − ω_B − ω_Z = 1` (from `D0 = K − B0 − Z0`); (2) the wall-referenced form `Xi_load = (δ_N−δ_K) + ω_B(δ_B−δ_K) + ω_Z(δ_Z−δ_K)`; (3) microscopic weighted-average forms of `δ_B`, `δ_Z`, `δ_N` (sums over modes/ports of logarithmic deformations, with weights summing to 1); (4) the static self-similarity theorem `Σ=0 ⇒ Xi_load = 0`; and the mismatch-field rewrite `Xi_load = Σ_r ρ^N Σ^N + ω_B Σ_α ρ^B Σ^B + ω_Z Σ_r ρ^Z Σ^Z`. This is an exact symbolic-closure stage; no numeric Family-1 constants appear.

## What the script claims to verify

Both scripts (docstring: "SymPy-backed audit for the static self-similarity decomposition") verify five things: (1) the weight identity `ω_K − ω_B − ω_Z − 1 = 0` and the `δ_D` weighted decomposition; (2) that `Xi_load = N01/N0 − D01/D0` equals the wall-referenced form `(δ_N−δ_K)+ω_B(δ_B−δ_K)+ω_Z(δ_Z−δ_K)`; (3) for explicit two-mode/two-port realizations, that the first-order perturbation of `B0 = Σ c²/ϖ²`, `Z0 = Σ Q/Δ`, `N0 = Σ P²/Δ²` equals the claimed normalized-weight average of logarithmic slopes; (4) that under the self-similar substitution `K1→δ*K, B01→δ*B0, Z01→δ*Z0, N01→δ*N0`, `Xi_load` collapses to 0; and (5) that the defect-field substitution `B01→(δ_K+Σ_B)B0` etc. reproduces `Σ_N + ω_B Σ_B + ω_Z Σ_Z`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `δ_D = ω_K δ_K − ω_B δ_B − ω_Z δ_Z` | `expect_zero("delta_D weighted decomposition", …)` py:57-60 / wl:46-49 | match |
| `ω_K − ω_B − ω_Z = 1` | `expect_zero("weight identity", omegaK-omegaB-omegaZ-1)` py:56 / wl:45 | match |
| `Xi_load = (δ_N−δ_K)+ω_B(δ_B−δ_K)+ω_Z(δ_Z−δ_K)` | `expect_zero("Xi_load wall-referenced form", …)` py:62-64 / wl:51-56 | match |
| `δ_B` = weighted avg of `2 δln(c_α/ϖ_α)` | `expect_zero("BdG weighted-average formula", …)` py:84 / wl:76 | match (2-mode realization) |
| `δ_Z` = weighted avg of `δln(Q_r/Δ_r)` | `expect_zero("Z weighted-average formula", …)` py:103 / wl:95 | match (2-port realization) |
| `δ_N` = weighted avg of `2 δln(P_r/Δ_r)` | `expect_zero("N weighted-average formula", …)` py:122 / wl:114 | match (2-port realization) |
| static self-similarity ⇒ `Xi_load = 0` | `expect_zero("Xi_load under static self-similarity", …)` py:130-131 / wl:121-125 | match |
| `Xi_load = Σ_N + ω_B Σ_B + ω_Z Σ_Z` (mismatch-field) | `expect_zero("Xi_load mismatch-field form", …)` py:135-136 / wl:127-131 | match |

`paper_alignment: aligned`. Every paper-side deliverable maps to a non-tautological script check; no script check is orphaned (the printed carry-forward formulas mirror the notes' boxed equations exactly).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 56 | `simplify(ω_K−ω_B−ω_Z−1)==0` | weight identity | yes |
| A2 | sympy | 57-60 | `simplify(δ_D − (ω_Kδ_K−ω_Bδ_B−ω_Zδ_Z))==0` | δ_D decomposition | yes |
| A3 | sympy | 62-64 | `simplify(Xi_load − Xi_wall_ref)==0` | wall-referenced form | yes |
| A4 | sympy | 84 | `simplify(B01_two/B0_two − δ_B_weighted)==0` | δ_B microscopic avg | yes |
| A5 | sympy | 103 | `simplify(Z01_two/Z0_two − δ_Z_weighted)==0` | δ_Z microscopic avg | yes |
| A6 | sympy | 122 | `simplify(N01_two/N0_two − δ_N_weighted)==0` | δ_N microscopic avg | yes |
| A7 | sympy | 130-131 | `simplify(Xi_self_similar)==0` | self-similarity theorem | yes |
| A8 | sympy | 135-136 | `simplify(Xi_sigma − (Σ_N+ω_BΣ_B+ω_ZΣ_Z))==0` | mismatch-field form | yes |
| B1-B8 | mathematica | 45,46-49,51-56,76,95,114,121-125,127-131 | `expectZero[…]` (FullSimplify→0, Exit[1] on fail) | same 8 claims | yes |

All assertions are subtraction-to-zero of two *independently-formed* expressions (the candidate identity vs. an explicit construction), so none is tautological. A2/A3 are the load-bearing ones; both have genuine content because `δ_D` and `Xi_wall_ref` are formed from the primitive symbols `K,B0,Z0,N0,…` and reduced, not assumed equal.

## Findings

None. See verdict justification.

## Independent-derivation check (Mathematica)

The `.wl` is **an independent re-derivation, not a transliteration** — confirmed strongly by the construction of the perturbation numerators. For all three microscopic sectors the two engines build the first-order slope by *different methods*:

- SymPy hand-codes the closed-form derivative, e.g. BdG (py:78):
  `B01_two = 2*c1*dc1/w1**2 - 2*c1**2*dw1/w1**3 + 2*c2*dc2/w2**2 - 2*c2**2*dw2/w2**3`
- Mathematica builds the perturbed expression and differentiates symbolically (wl:67-68):
  `b0Eps = (c1+eps*dc1)^2/(w1+eps*dw1)^2 + (c2+eps*dc2)^2/(w2+eps*dw2)^2;`
  `b01Two = FullSimplify[D[b0Eps, eps] /. eps -> 0, …]`

Same pattern for Z (py:97 explicit quotient-rule literal vs. wl:86-87 `D[z0Eps,eps]`) and N (py:116 literal vs. wl:105-106 `D[n0Eps,eps]`). That the two routes both reduce to the *same* weighted-average formula is exactly the cross-engine corroboration the second-engine policy wants. Sections 1 and 5 (the symbolic Xi_load substitutions) look structurally similar across engines, but that is the only natural way to perform the substitution on primitive symbols and is not evidence of a port; even there the engines emit *differently-presented* canonical forms (see Engine cross-check), which a line-by-line transliteration would not. Confidence the `.wl` is independent: high.

## Engine cross-check

Both engines pass every check and agree on the final `Xi_load` and mismatch-field prototype, modulo (independent) canonicalization:

- SymPy `Xi_load = (N0*(-B01 + K1 - Z01) + N01*(B0 - K + Z0))/(N0*(B0 - K + Z0))` (output:8)
- Mathematica `Xi_load = n01/n0 - (b01 - k1 + z01)/(b0 - k + z0)` (output:11)

These are algebraically identical: both equal `N01/N0 + (−B01+K1−Z01)/(B0−K+Z0) = δ_N − δ_D`. The mismatch-field prototype likewise matches (SymPy `(-B0*SigmaB + SigmaN*(B0-K+Z0) - SigmaZ*Z0)/(B0-K+Z0)`, output:30, equals Mathematica `sigmaN - (b0*sigmaB + sigmaZ*z0)/(b0-k+z0)`, output:38). The differing surface forms are an independent-canonicalization signature, reinforcing the non-transliteration call. `engines_agree: true`.

## Verdict justification

`clean`. I read the card, the notes (which carry the full boxed derivation), and the appendix block (rows 79, 560-609) before the scripts, and the scripts' eight checks map one-to-one onto the paper's eight deliverables with no mismatch, no extra, no missing. Attacks attempted that failed: (1) **tautology** — the load-bearing A2/A3 subtract a candidate identity from a *constructed* quantity (`δ_D` from primitives, `Xi_load` from `N01/N0−D01/D0`), so they can fail if the algebra is wrong; not tautological. (2) **symbol-domain abuse** — all symbols are `real` with the denominators (`K,B0,Z0,N0,c_i,w_i,Q_i,Δ_i,P_i ≠ 0`) declared `nonzero`; `D0=K−B0−Z0` is never forced nonzero but every division by `D0` keeps it symbolic (no false cancellation), and the self-similar substitution at A7 does not silently zero a denominator. (3) **self-similarity loophole** — the substitution `K1→δ*K,…,N01→δ*N0` makes every δ equal `δ*`, so `Xi_load→0` is a real consequence, not a hardcoded zero. (4) **stale Family-1 numerology** — no `4107`, `100π²`, `168π²`, or `168%` appears anywhere in the four files; this stage is purely symbolic and that watch class is inapplicable here. (5) **transliteration** — disproved by the divergent numerator-construction methods above. Outputs are fresh (sympy txt mtime 16:35:36 > py 16:35:28; mathematica txt 16:13:48 > wl 15:55:34) and their content matches the current scripts.

## Self-test notes

Checked variable-independence (the `D[…,eps]` derivatives in the `.wl` genuinely depend on `eps` via the perturbed numerators/denominators — no identically-zero derivative; the SymPy literals are real first-order forms), denominator-nonvanishing under the self-similar/defect-field substitutions (A7/A8 keep `D0` symbolic, no division-by-zero), and trivial-case substitution (setting all `δ`-numerators proportional → `Xi_load=0` follows; setting `Σ=0` → mismatch-field reduces to `Σ_N+ω_BΣ_B+ω_ZΣ_Z`, matching). No paper round-trip issue since no fix is prescribed. All traps clear.

## Value Reconciliation (pass-2 augmentation)

The scripts emit only **symbolic** closed-form deliverables (no numeric constants — confirmed: the sole integer literals are structural factors of 2, the weight-identity `1`, and exponents). Each emitted symbolic deliverable is checked against the `.tex` card / appendix and `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `ω_K − ω_B − ω_Z = 1` (weight identity) | py:56 / wl:45; out:5/6 | notes:110-113 (boxed); py/wl carry-print | MATCH |
| `δ_D = ω_K δ_K − ω_B δ_B − ω_Z δ_Z` | py:57-60 / wl:46-49; out:6/7 | notes:94-103 (boxed) | MATCH |
| `Xi_load = (δ_N−δ_K)+ω_B(δ_B−δ_K)+ω_Z(δ_Z−δ_K)` | py:62-64 / wl:51-56; out:7/9 | notes:127-135 & appendix:583-588 (eq:Xiload-selfsimilarity) | MATCH |
| `Xi_load` canonical form `δ_N − D01/D0` | out:8 (py) / out:11 (wl) | notes:6-14, 22-32 (Xi_load defn & boxed) | MATCH (algebraically; both engines) |
| `δ_B = Σ ρ^B · 2 δln(c_α/ϖ_α)` | py:84 / wl:76; out:13/16 | notes:162-186 (boxed) & appendix:592 | MATCH |
| `δ_Z = Σ ρ^Z · δln(Q_r/Δ_r)` | py:103 / wl:95; out:18/22 | notes:199-222 (boxed) & appendix:594 | MATCH |
| `δ_N = Σ ρ^N · 2 δln(P_r/Δ_r)` | py:122 / wl:114; out:23/28 | notes:235-258 (boxed) & appendix:596 | MATCH |
| static self-similarity ⇒ `Xi_load = 0` | py:130-131 / wl:121-125; out:28/34 | notes:354-356 (boxed) & card Output, appendix:609 | MATCH |
| mismatch-field `Xi_load = Σ_N + ω_B Σ_B + ω_Z Σ_Z` | py:135-137 / wl:127-132; out:29-30/36-38 | notes:300-309 (boxed) & appendix:599-607 (eq:Xiload-Sigma-fields) | MATCH |

INTERNAL scaffolding (accounted for, no finding): perturbation primitives `B01_two/B0_two`, `Z01_two/Z0_two`, `N01_two/N0_two`, `b0Eps/z0Eps/n0Eps`, the two-mode/two-port placeholders (`c1,c2,w1,w2,Q1,Q2,Δ1,Δ2,P1,P2` and their `_p`/`d_` perturbations), the substitution helper symbols `delta_star`, `SigmaB/SigmaZ/SigmaN`, and all PASS/`expect_zero`/`=0` residual flags.

reconciliation: complete; 9 deliverable values checked, 0 misaligned
