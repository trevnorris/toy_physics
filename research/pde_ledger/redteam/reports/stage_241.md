---
unit_id: 241
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 241 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_241.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (table row line 94; narrative lines 1237-1281)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The card's `\stagefield{Output}` states: "Primitive ranking theorem: the selected branch has a complete phase diagram with thresholds $\varrho_{W\Lambda}=2(1+\beta^2)/[3(2+\beta^2)]$ and $\varrho_{U\Lambda}=2(1+\beta^2)/[3(1+\beta+\beta^2)]$." The `\stagefield{Inputs}` import `\epsilon_*=1-3\varrho/2`, `\sigma=4/(3\varrho)-2`, `0<\varrho<2/3`, and the constructive coherent bound `0<\beta<2/11`. The `\stagefield{Derivation ledger}` enumerates four deliverables: (D1) the primitive coherent weights and `w_\chi=2w_Z`; (D2) the two thresholds `\varrho_{W\Lambda}`, `\varrho_{U\Lambda}`; (D3) their ordering; (D4) the three ranking regions. The notes add detail: the selected-branch reduction `\epsilon_*=1-3\varrho/2` is *derived* from the Stage-240 support law `\varrho=2(1-\epsilon_*)/3` (notes line 164); the strict twin-window inclusion `1/\varrho-2 < \sigma_{\rm sel} < 2/\varrho-2`; the branch-independent positivity identities for `w_Z-w_\Lambda`, `w_Z-w_W`, `w_W-|w_U|`; the factorized sign-transfer laws over a common denominator `D`; threshold ordering `0<\varrho_{W\Lambda}<\varrho_{U\Lambda}<2/3`; numeric windows from `0<\beta<2/11`; and representative region checks. The appendix narrative (lines 1237-1281) restates the same curve, threshold pair, ordering, and three region orderings verbatim.

## What the script claims to verify

The docstring (lines 2-25) and the eight assertion blocks claim to verify: (1) the selected-branch reduction `eps_* = 1-3varrho/2`, `sigma = 4/(3varrho)-2`; (2) strict twin-window inclusion via the two gap identities `1/(3varrho)` and `2/(3varrho)`; (3) the weight identities `w_chi=2w_Z` and the three positivity-form differences `w_Z-w_Lambda`, `w_Z-w_W`, `w_W-|w_U|`; (4) the two threshold-consistency identities relating `eps_sel` to the imported crossover conditions `1/(2+beta^2)` and `beta/(1+beta+beta^2)`; (5) the two factorized sign-transfer laws over `D`; (6) the threshold-ordering positivity forms; (7) numeric threshold endpoints at `beta=0` and `beta=2/11` plus the two derivative forms; (8) numeric region orderings at `beta=1/10` and three representative `varrho` samples. The output (mtime newer than script) shows all 31 checks PASS, exit 0.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Selected-branch reduction `eps_*=1-3varrho/2` (derived from `varrho=2(1-eps_*)/3`) | lines 58, 61 | **partial** — `eps_sel` is defined as the literal `1-3varrho/2` and the assertion is `eps_sel - (1-3varrho/2) == 0` (tautological; does not thread through the Stage-240 support law) |
| `sigma=4/(3varrho)-2` | lines 59, 62 | match — genuinely substitutes `2eps/(1-eps)` and reduces |
| Strict twin-window inclusion | lines 70-77 | match |
| Weight identities `w_chi=2w_Z`, `w_Z>w_Lambda`, `w_Z>w_W`, `w_W>|w_U|` | lines 90-106 | match (exact positive-form residuals) |
| Thresholds `varrho_WLambda`, `varrho_ULambda` (D2) | lines 118-132 + 139-152 | match — threshold-consistency + factorized sign laws place the zero exactly at the threshold |
| Threshold ordering `0<varrho_WLambda<varrho_ULambda<2/3` (D3) | lines 157-172 | match |
| Three ranking regions (D4) | lines 218-237 | match (numeric, single `beta=1/10`) |
| Numeric windows from `0<beta<2/11` | lines 179-205 | match (script `125/369` is the correct value; see F3 note re a notes typo) |
| Mathematica independent route | — | **missing** (card line 11: "Mathematica audit: none yet.") |

`paper_alignment: aligned` — every paper deliverable has a faithful script-side counterpart; the issues are a tautological first check, a wholly-absent second engine, and a notes-prose typo.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 61 | `eps_sel - (1 - 3varrho/2) == 0` | selected-branch reduction | **no** (tautological: `X-X`) |
| A2 | sympy | 62 | `sigma_sel - (4/(3varrho)-2) == 0` | sigma reduction | yes |
| A3 | sympy | 70-73 | `sigma_sel - (1/varrho-2) - 1/(3varrho) == 0` | twin-window lower gap | yes |
| A4 | sympy | 74-77 | `(2/varrho-2) - sigma_sel - 2/(3varrho) == 0` | twin-window upper gap | yes |
| A5 | sympy | 90 | `w_chi - 2 w_Z == 0` | `w_chi=2w_Z` | yes |
| A6 | sympy | 91-94 | `w_Z - w_Lambda - (1-eps)^2/N == 0` | `w_Z>w_Lambda` | yes |
| A7 | sympy | 95-102 | `w_Z - w_W - [beta^2 eps^2 + 3(eps-1/2)^2 + 1/4]/N == 0` | `w_Z>w_W` | yes |
| A8 | sympy | 103-106 | `w_W - |w_U| - eps(1-eps)(1-beta)/N == 0` | `w_W>|w_U|` | yes |
| A9 | sympy | 121-124 | `eps_sel - 1/(2+beta^2) - 3/2(varrho_WLambda-varrho) == 0` | WLambda threshold | yes (consistency w/ imported crossover) |
| A10 | sympy | 125-132 | `eps_sel - beta/(1+beta+beta^2) - 3/2(varrho_ULambda-varrho) == 0` | ULambda threshold | yes |
| A11 | sympy | 139-145 | `wL_sel - wW_sel - (2-3varrho)(2+beta^2)(varrho_WLambda-varrho)/D == 0` | WLambda sign law | yes (load-bearing for crossover) |
| A12 | sympy | 146-152 | `wL_sel - wU_sel - (2-3varrho)(1+beta+beta^2)(varrho_ULambda-varrho)/D == 0` | ULambda sign law | yes |
| A13 | sympy | 157-164 | `varrho_ULambda - varrho_WLambda - [pos form] == 0` | ordering gap 1 | yes |
| A14 | sympy | 165-172 | `2/3 - varrho_ULambda - [pos form] == 0` | ordering gap 2 | yes |
| A15-A18 | sympy | 179-194 | numeric threshold endpoints `beta=0, 2/11` | numeric windows | yes |
| A19-A20 | sympy | 196-205 | derivative forms `d varrho/d beta` | window monotonicity | yes |
| A21-A32 | sympy | 218-237 | numeric `assert_positive_sample` region orderings at `beta=1/10` | three regions | partial (single beta, float) |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.py:58`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.py:61`

**What's wrong:**
Line 58 defines `eps_sel = sp.simplify(1 - sp.Rational(3, 2) * varrho)`. Line 61 then asserts `assert_zero(eps_sel - (1 - sp.Rational(3, 2) * varrho), ...)`. The subtrahend is the byte-for-byte definition of `eps_sel`, so the residual is `X - X = 0` by construction. This assertion cannot fail for any physics; it merely confirms the assignment round-trips. Meanwhile the notes (line 158-173) state the selected-branch reduction is *derived* from the Stage-240 support law `\varrho = 2(1-\epsilon_*)/3`, i.e. the content of the claim is "solving the Stage-240 support relation for `eps_*` yields `1-3varrho/2`." The script never threads through that relation. (Note: the companion check A2 at line 62, which substitutes `2*eps_sel/(1-eps_sel)` and reduces to `4/(3varrho)-2`, is genuine and is *not* part of this finding.)

**Why this matters:**
The first listed deliverable in the docstring ("Selected-branch reduction `eps_* = 1-3varrho/2`") is reported as verified, but the check is vacuous. A wrong selected-branch law would still pass line 61. The PASS at output line 9 gives false assurance for the foundational reduction every later check builds on.

**Required change:**
Replace the tautological line-61 check with one that exercises the Stage-240 support law. Define the support law `varrho_from_eps = 2*(1 - eps_star)/3` (notes line 164) and assert that substituting `eps_star -> eps_sel` makes `varrho_from_eps - varrho` vanish (i.e. `assert_zero((2*(1 - eps_star)/3 - varrho).subs(eps_star, eps_sel), ...)`), so the residual passes through the support relation rather than re-stating `eps_sel`. Equivalently, `sp.solve(2*(1-eps_star)/3 - varrho, eps_star)` should equal `[1 - 3*varrho/2]`. Codex chooses the exact form; the acceptance criterion is that the residual is built from `varrho = 2(1-eps_*)/3` and not from a literal copy of `1-3varrho/2`.

**Verification:**
New check at/near line 61 references `2*(1 - eps_star)/3` (the support law), not a second literal `1 - 3*varrho/2`; output still shows `[ok]` for the selected-branch reduction; script exits 0.

### F2 — missing_verification_script (missing_mathematica)

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_241.tex:11` ("Mathematica audit: none yet.")
- no `.wl` exists under `mathematica/` for stage 241

**What's wrong:**
This unit is `is_checkpoint: False` but `is_status_only_candidate: False`, so both engines are required. Only the SymPy script exists. Every claim in this stage — rational-function identities (`w_chi=2w_Z`, the three positive-form weight differences), the factorized sign laws over the polynomial `D`, the threshold-ordering positivity forms, the exact numeric endpoints `125/369`, `250/441`, and the two `d varrho/d beta` derivative forms — is squarely within Mathematica's native capability (`Together`, `Factor`, `FullSimplify`, `D`, exact rational arithmetic, `Reduce`). There is no genuine impossibility, so the dual-engine rule requires an independent Mathematica route.

**Why this matters:**
A single-engine verification has no cross-check against a SymPy-specific `simplify`/`factor` quirk. The whole point of the second engine is an independent re-derivation from the physical premises (the weight formulas + the Stage-240 support law), not a transliteration.

**Required change:**
Codex writes a NEW independent-route Mathematica script (see directive claim manifest). It must be a genuinely different decomposition, not a line-by-line port of the `.py`.

**Verification:**
`mathematica/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_mathematica_audit.wl` exists, derives each manifest claim via native primitives, and exits 0 with all in-file checks passing. Independent-derivation check (transliteration scan) passes.

### F3 — paper_misalignment (notes_contradicts_script)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.md:577`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_sympy_audit.py:184`

**What's wrong:**
Notes line 577 boxes `\frac13 < \varrho_{W\Lambda} < \frac{193}{369}\approx 0.338753`. The script (line 184) asserts `varrho_WLambda(beta=2/11) = 125/369`. For `beta=2/11`, `beta^2=4/121`, `varrho_WLambda = 2(125/121)/[3(246/121)] = 250/738 = 125/369`, which equals `0.338753...`. The decimal printed in the SAME boxed line (`0.338753`) matches `125/369`, NOT `193/369` (`193/369 = 0.523...`). So the numerator `193` in the notes is a prose typo; the canonical value is `125/369` (the script value, which also matches the notes' own decimal). The paper card and appendix narrative do not carry this numeric window, so only the notes prose is wrong.

**Why this matters:**
A reader trusting the notes fraction `193/369` would get a value (`0.523`) that contradicts both the script and the decimal printed beside it, and that even exceeds `varrho_ULambda`'s lower bound. It is a localized typo, not a derivation error, but it should be corrected so the notes are internally consistent.

**Required change:**
None for Codex (notes are prose; routed to user). See `## Resolve before fix_loop` in the directive.

**Verification:**
After user resolution, notes line 577 reads `\frac{125}{369}\approx 0.338753` (or the user confirms the script is canonical and edits the notes accordingly). No script change.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot yet be assessed. F2 requires the new script to be an independent decomposition; the directive's claim manifest specifies a different route (build weights from `eps_*`, take the `w`-differences and factor, solve crossover conditions via `Reduce`/`Solve` rather than importing the threshold literals) so the verifier can confirm non-transliteration after Codex writes it.

## Engine cross-check

Only one engine present; `engines_agree: n/a`. SymPy output (mtime 2026-05-11T12:52, newer than script mtime 11:58) shows all 31 checks PASS, exit 0. No `stale_output`.

## Verdict justification

`verdict: findings`. The SymPy script is substantively sound and well-anchored on the bulk of its work: A2-A20 are genuine, non-tautological algebraic and rational-function identities that faithfully exercise every paper deliverable, and I confirmed by hand that the `w_Z-w_W` positive form, the threshold-consistency identities, and the sign-law factorizations reduce to zero through real algebra (not by construction). The defects are: (F1) the very first selected-branch-reduction check is a `X-X` tautology that does not thread through the Stage-240 support law it claims to verify; (F2) the second engine is entirely absent though Mathematica can independently verify every claim, violating the dual-engine rule; (F3) a notes-only numerator typo (`193/369` should be `125/369`, the script's correct value, matching the notes' own decimal) needing user resolution since notes are prose. Paper alignment is otherwise exact across card, notes, and appendix. No stop-cold condition: F1 is a local fix, F2 is additive, and F3 routes to the user without blocking the script fixes.

## Self-test notes

I hand-checked (no execution) the F1 replacement: `(2(1-eps_*)/3 - varrho)` with `eps_* -> 1-3varrho/2` gives `varrho - 2(3varrho/2)/3 = varrho - varrho = 0`, non-trivially threading the support law, so the new `assert_zero` is genuinely satisfiable-or-falsifiable (trap 3 clear). For F2's manifest I checked variable dependence: the `D[varrho_WLambda, beta]` claim depends on `beta` (nonzero derivative), and the crossover `Solve`/`Reduce` targets `eps_*`, both well-posed (trap 1 clear); the region orderings are at concrete `beta=1/10` rational points so no symmetric-domain parity issue arises (trap 2 n/a). I confirmed `125/369 = 0.338753`, `193/369 = 0.523`, `250/441 = 0.5669` by hand-reducing the rationals (trap: arithmetic), establishing the F3 typo direction. The `.wl` target path is in `mathematica/` per trap-4.
