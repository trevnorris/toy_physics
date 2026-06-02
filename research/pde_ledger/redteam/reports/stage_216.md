---
unit_id: 216
batch: VI.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-01T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: insufficient
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate.md]
  paper_appendix: present
---

# Audit unit 216 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_216.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part06.tex` (rows at lines 63, 211, 236; detailed narrative lines 1093-1154)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

The `\stagefield{Output}` states verbatim: "Unique five-coordinate simplex, exact face reduction to support-\(\le4\), support-cardinality-five gate, and proof that no higher-support local search exists." The card's `Derivation ledger` adds that the stage "derives the five-coordinate gradient-synergy and ten-way cross-leverage laws, states the five-way canonical screens, and proves the support-cardinality ceiling." The notes enumerate seven deliverables, of which the algebraically-verifiable ones are: (i) the unique positive spherical five-simplex `\(\Delta_5^+\)` with `\(\sum a_p^2 = 1\)` and its five codimension-one faces being the Stage-215/249 quadruple simplices; (ii) the gradient-synergy theorem with unique optimal point `\(\mathbf a_5^{\rm grad}=(k_\lambda,\dots,k_W)/\sqrt{\sum k_p^2}\)`, max slope `\(\sqrt{\sum k_p^2}\)`, and strict first-order improvement `\((k_5^{\rm grad})^2-(k_{Q_{\widehat p}}^{\rm grad})^2 = k_p^2>0\)` over each face; (iii) the ten-way cross-leverage law `\(w_\Sigma(\mathbf a)=2\sum_{p<q}a_pa_q=(\sum_p a_p)^2-1\)` with Cauchy slack `\(5\sum a_p^2-(\sum a_p)^2=\sum_{p<q}(a_p-a_q)^2\ge0\)`, giving `\(w_\Sigma\le4\)` with equality at the barycenter `\(\mathbf a_5^{\rm eq}=(1,1,1,1,1)/\sqrt5\)` where `\(w_\Sigma=4\)`; (iv) the fixed-simplex certified bracket `\(\tau_{5,\star}=2H_0/(k_5+\sqrt{k_5^2-2H_0\kappa_{5,\star}})\)`; (v) the canonical seven-row five-way screen packet (5 imported faces + grad + eq); and (vi) the support-cardinality ceiling (no support-`>5` rays, since the free space has exactly 5 axes). The appendix (lines 1093-1154) confirms the same `\(\Delta_5^+\)`, `\(k_5\)`, `\(\mathbf a_5^{\rm grad}\)`, `\(\mathbf a_5^{\rm eq}\)`, and `\(w_\Sigma\le4\)`. The card explicitly marks the Mathematica audit as "none yet."

## What the script claims to verify

The SymPy script computes the same family of quantities the paper states: it builds `a_grad = k/||k||` and computes `grad_norm_sq` (printed as 1) and `k_grad` (printed as `sqrt(sum k^2)`); it computes the five face gaps `k_norm_sq - face_max_sq` (each prints as `k_p**2`); it computes the four gradient ratios `r,s,t,u = k_c/k_lam, k_gam/k_lam, k_U/k_lam, k_W/k_lam`; it computes `identity_check = expand(w_sigma - (asum**2 - norm_sq))` (prints 0) and `slack_check = expand(slack - pair_sum)` (prints 0); it substitutes the barycenter into `w_sigma` to get `w_eq = 4`; and it forms the symbolic bracket `tau`. Crucially, the script contains NO assertion of any kind — no `assert`, no `sys.exit(1)`, no `raise`, no `expect_*`. Every result is merely appended to an output list and printed. The script's "PASS" status in the saved `.txt` reflects only that the script exits 0, which it always will regardless of whether any computed quantity is correct.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side coverage | Status |
|---|---|---|
| `\(\Delta_5^+\)` is unique support-5 simplex; 5 faces are quadruple simplices | Prints face count = 5 and labels as text strings only | partial (descriptive print, no check) |
| `\(\mathbf a_5^{\rm grad}=k/\|k\|\)`, `\(\|a_5^{\rm grad}\|^2=1\)` | `grad_norm_sq` computed, printed as 1; never asserted | partial (computed, not checked) |
| `\(k_5^{\rm grad}=\sqrt{\sum k_p^2}\)` | `k_grad` computed/printed; never asserted | partial |
| `\((k_5^{\rm grad})^2-(k_{Q_{\widehat p}}^{\rm grad})^2=k_p^2>0\)` for each face | five diffs computed/printed; never asserted | partial |
| Ratios `\(r,s,t,u=k_c/k_\lambda,\dots\)` | computed/printed; never asserted | partial |
| `\(w_\Sigma=(\sum a_p)^2-1\)` identity | `identity_check` computed, prints 0; never asserted | partial |
| Cauchy slack `\(5\sum a^2-(\sum a)^2=\sum(a_p-a_q)^2\)` | `slack_check` computed, prints 0; never asserted | partial |
| `\(w_\Sigma\le4\)`, equality at barycenter, `\(w_\Sigma(a_{\rm eq})=4\)` | `w_eq` computed, prints 4; bound itself never tested; never asserted | partial |
| Bracket `\(\tau_{5,\star}=2H_0/(k_5+\sqrt{k_5^2-2H_0\kappa})\)` | `tau` formed/printed; never asserted | partial |
| Lower-cardinality comparison `\(w_\Sigma=3,2,1\)` for quad/triple/pair | hardcoded as `Integer(3),(2),(1)` with provenance comments; never asserted | partial |
| Support-cardinality ceiling = 5 | printed as text only | partial (descriptive) |
| Independent second engine (Mathematica) | none | missing |

Every algebraic deliverable is *computed* and the computed value happens to match the paper, but NONE is asserted — so the script cannot fail if any of these were wrong. The paper alignment of the computed forms is exact (the legacy "Stage 250"/"Stage 249" numbering in the notes is internal renumbering, not a content mismatch), hence `paper_alignment: aligned`. The defect is verification strength, not paper↔script content disagreement.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| (none) | sympy | — | no `assert` / `sys.exit` / `raise` anywhere | n/a | no |

The script has zero assertions. `grad_norm_sq` (line 28), the face diffs (line 43), `identity_check` (line 69), `slack_check` (line 72), and `w_eq` (line 79) are all computed and printed (lines 32-33, 45, 75-76, 84) but never compared against an expected value with a failing check. Lines 80-82 (`w_quad_eq=3`, `w_triple_eq=2`, `w_pair_eq=1`) are hardcoded literals with provenance comments, printed only.

## Findings

### F1 — missing_verification_script

**Severity:** high
**Subtype:** script_doesnt_cover_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate_sympy_audit.py:5-115` (whole `main()`)

**What's wrong:**
The SymPy script contains no assertion of any kind. It computes the load-bearing quantities the paper requires — `grad_norm_sq` (line 28, should be `1`), `k_grad` (line 29, should be `sqrt(sum k^2)`), the five face gaps `diffs` (line 43, should each equal `k_p**2`), `identity_check` (line 69, should be `0`), `slack_check` (line 72, should be `0`), and `w_eq` (line 79, should be `4`) — and then merely appends each to the output list (`out.append(...)`, e.g. lines 32-33, 45, 75-76, 84) and prints them. There is no `assert`, `sys.exit(1)`, `raise`, or `expect_*`. The saved output's "Status: PASS" (output line 6) reflects only `Exit code: 0`, which is unconditional: if any of these computations were silently wrong (e.g. a sign flip in `w_sigma`, a wrong simplex constraint), the script would still exit 0 and still report PASS. This is the `script_doesnt_cover_claim` subtype: a script exists but contains no real assertion.

**Why this matters:**
Per the prompt, "Output transcripts that just say PASS" must be re-read for whether the assertion is even present. Here there is no assertion at all, so the entire stage is effectively unverified despite the green PASS line. None of the paper's six algebraic deliverables is actually defended against regression. A future edit could break any computed identity and the audit would not catch it.

**Required change:**
Add explicit failing assertions for each computed quantity that mirrors a paper deliverable. Each computed result that the script prints must be turned into a check that raises on mismatch. See the directive for the precise claim list (M1-M6).

**Verification:**
After the fix, the script must contain `assert`/`sys.exit(1)` (or an `expectZero`-style helper) covering each of M1-M6 below, and `python3 <path>` must exit 0 with all checks passing. The verifier confirms the new assertions appear in the file and that deliberately perturbing any one of them (mentally) would flip the exit code.

### F2 — missing_verification_script

**Severity:** high
**Subtype:** missing_mathematica
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` (no `.wl` for stage 216)

**What's wrong:**
Stage 216 is non-status-only (`is_status_only_candidate: False`, checkpoint: `False`) and computes genuinely verifiable symbolic algebra: the unit-norm of the gradient-optimal ray, the max-slope identity, the five per-face first-order gap `= k_p^2`, the cross-leverage identity `w_Sigma = (sum a)^2 - 1`, the Cauchy slack identity, the barycenter leverage value `4`, and the closed-form certified bracket. Every one of these is independently derivable in native Mathematica (`Solve`/`Reduce` for the Lagrange-multiplier optimum, `Simplify`/`Together`/`Expand` for the identities, `Series`/`Coefficient` or matrix quadratic forms for the leverage law, `Sqrt`/algebraic simplification for the bracket). The paper card itself records "Mathematica audit: none yet" (`stage_216.tex:11`). Per the dual-engine project contract and the prompt's checkpoint/non-status rule ("both scripts are required, and missing scripts are findings"), the absence of any `.wl` is a finding. There is no reason this stage genuinely *cannot* be verified by Mathematica — it computes concrete symbolic identities.

**Why this matters:**
A single-engine stage has no cross-check. The whole point of the second engine is to catch transcription/algebra errors the first engine's author would not see. With F1 also open (the SymPy side has no assertions), the stage currently has *zero* defended checks across *zero* independent engines.

**Required change:**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage216_unique_five_coordinate_positive_simplex_and_support_cardinality_5_gate.wl` that independently verifies M1-M6 (see directive). Codex designs the route; the directive states only the claim manifest and the anti-transliteration guard.

**Verification:**
After the fix, `math -script <path>` exits 0 with all in-file checks passing, and the `.wl` derives each claim via native Mathematica primitives through a different decomposition than the `.py` (not a line-by-line port).

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot be assessed yet. The directive's anti-transliteration guard requires the new `.wl` to re-derive the gradient optimum and the leverage law from the physical premises (Lagrange/Cauchy) using native Mathematica primitives via a different decomposition than the SymPy script — not a syntactic port of the `out.append` choreography.

## Engine cross-check

Only one engine present (SymPy), and it has no assertions; `engines_agree: n/a`.

## Verdict justification

The computed symbolic forms in the SymPy script match the paper and notes exactly — I tried to find a `target_mismatch` (wrong constant, wrong simplex constraint, wrong leverage bound) and found none: `grad_norm_sq=1`, `k_grad=sqrt(sum k^2)`, face gaps `=k_p^2`, the two identity residuals `=0`, and `w_eq=4` are all the right paper values, so there is no `paper_misalignment`. The defect is purely verification strength: the script contains no assertions whatsoever (F1, `script_doesnt_cover_claim`), so its PASS is unconditional, and there is no second engine (F2, `missing_mathematica`) despite the math being fully Mathematica-verifiable. Both are high-severity script-side findings, neither is `UNFIXABLE` (the fix is adding assertions and a second engine, not reconciling contradictory math) and neither is `CRITICAL_DOWNSTREAM` (adding checks changes no quoted-forward value). Verdict: `findings`.

## Self-test notes

(1) Variable independence: no `sp.diff`/`D[]` is proposed in the directive — the claims are pure algebraic identities and a Lagrange/Cauchy optimum, so the zero-derivative trap does not apply. (2) Symmetry/parity: no unbounded integrals; the `\(\Delta_5^+\)` constraint is a finite quadratic — no parity trap. (3) Trivial-case pre-check: I mentally substituted the barycenter `\(a_p=1/\sqrt5\)` into `w_Sigma` and got `2 * 10 * (1/5) = 4` (matches M5), and into the slack identity `5*1 - (5/\sqrt5)^2 = 5 - 5 = 0` (matches M4); the gradient optimum `k/||k||` has squared norm `(sum k^2)/(sum k^2)=1` (matches M1) and slope `sum k^2/||k|| = ||k|| = sqrt(sum k^2)` (matches M2); each face gap `||k||^2 - (||k||^2 - k_p^2) = k_p^2 > 0` (matches M3). All non-trivially true. (4) Path specs: F2 target is `mathematica/...wl`, F1 target is the existing `scripts/...py` — both named in full in the directive. (5) Paper round-trip: the assertions I prescribe assert exactly the paper's stated values (1, sqrt(sum k^2), k_p^2, 0, 0, 4, bracket form), introducing no new constant, so no new `paper_misalignment` is created.
