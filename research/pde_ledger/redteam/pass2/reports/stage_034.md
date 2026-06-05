---
unit_id: 034
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage034_softening_depth_normal_form.md]
  paper_appendix: present
---

# Audit unit 034 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_034.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage034_softening_depth_normal_form.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 58 references this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.txt`

## What the paper claims

Stage 034 (Part II, anchor MTDC-T5, claim status: exact closure) eliminates the selected eigenvector algebra by introducing a scalar softening-depth variable `x := A − λ_-` on the stable branch `0 ≤ x < A` (eq:app-stage034-x). The `\stagefield{Output}` reads verbatim: *"Stage~034 outputs the softening-depth variable \eqref{eq:app-stage034-x}, the loading law \eqref{eq:app-stage034-alpha-x}, the selected overlap \eqref{eq:app-stage034-s-x}, the normalization product \eqref{eq:app-stage034-N-x}, and the required support formula \eqref{eq:app-stage034-support-req}."* The five boxed deliverables are: (1) `x := A − λ_-`; (2) the loading law `α_0(x) = x(x+ΔK_ax) / [κ0²(x+ΔK_ax)+κ1²x]`; (3) the selected overlap `s_-(x) = [κ0²(x+ΔK_ax)+κ1²x]² / [κ0²(x+ΔK_ax)²+κ1²x²]`; (4) the invariant normalization product `N_-(x) = β0[κ0²(x+ΔK_ax)+κ1²x]⁴ / {κ0²(A−x)[κ0²(x+ΔK_ax)²+κ1²x²]²}`; and (5) the monotone loading derivative `dα0/dx = [κ0²(x+ΔK_ax)²+κ1²x²]/[κ0²(x+ΔK_ax)+κ1²x]² > 0`, used to obtain the required support formula `g_{B,req}²/ϖ² = x(x+ΔK_ax)/[κ0²(x+ΔK_ax)+κ1²x] − X²/(Ω_U²Δ0)`. The notes add the provenance of each form (secular law `1 = α0[κ0²/x + κ1²/(x+ΔK_ax)]`, the overlap identity `s_- = −dλ_-/dα0 = dx/dα0`, and `N_-(α0) = β0 s_-²/(κ0² λ_-)`).

## What the script claims to verify

Both engines build the secular quantities in the eigenvalue variable λ from the rank-1 secular sum `S1 = κ0²/(A−λ)+κ1²/(A+ΔK−λ)` (so `α_λ = 1/S1`, `s_λ = S1²/S1′`, `N_λ = β0 s_λ²/(κ0² λ)`), then independently posit the paper's candidate closed forms in `x` (`alpha_x`, `s_x`, `N_x`), and assert each λ-form equals the corresponding x-form under the substitution `λ = A − x`. They further verify the monotonicity form `dα0/dx` equals the manifestly-positive target form, the duality identity `s_x · dα0/dx = 1`, and that the required-support loading solved in x-form equals the same quantity solved in λ-form (`g_{B,req}²/ϖ²`). These are exactly the five paper deliverables plus the two supporting identities (overlap-as-inverse-derivative and the λ↔x support consistency).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `x := A − λ_-` | encoded as substitution `lam → A − x` throughout; not a standalone assertion (a definition) | match (definitional) |
| (2) `α_0(x)` loading law | `expect_zero("alpha(lambda=A-x) - alpha(x)", ...)` py:61 / wl:65 | match |
| (3) `s_-(x)` overlap | `expect_zero("s(lambda=A-x) - s(x)", ...)` py:62 / wl:66 | match |
| (4) `N_-(x)` normalization product | `expect_zero("N(lambda=A-x) - N(x)", ...)` py:63 / wl:67 | match |
| (5a) `dα0/dx` manifestly-positive form | `expect_zero("dalpha/dx - manifestly positive form", ...)` py:70 / wl:75 | match |
| (5b) `s_- = dx/dα0` duality (Checks item 2) | `expect_zero("s_x * d alpha / dx - 1", ...)` py:71 / wl:76 | match |
| (6) required support `g_{B,req}²/ϖ²` | x-form printed py:80-81 / wl:90; cross-checked vs λ-form `expect_zero("lambda-form vs x-form support loading", ...)` py:89 / wl:94 | match |

`paper_alignment: aligned` — every boxed deliverable maps to a substantive, non-tautological assertion in both engines; no extras, no mismatches.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 61 | `simplify(alpha_lam.subs(lam,A-x) - alpha_x) == 0` | (2) loading law | yes |
| A2 | sympy | 62 | `simplify(s_lam.subs(lam,A-x) - s_x) == 0` | (3) overlap | yes |
| A3 | sympy | 63 | `simplify(N_lam.subs(lam,A-x) - N_x) == 0` | (4) normalization product | yes |
| A4 | sympy | 70 | `simplify(dalpha_dx - dalpha_target) == 0` | (5a) dα0/dx form | yes |
| A5 | sympy | 71 | `simplify(s_x*dalpha_dx - 1) == 0` | (5b) `s_- = dx/dα0` | yes |
| A6 | sympy | 89-92 | `simplify(gBreq_lambda.subs(lam,A-x) - gBreq_x) == 0` | (6) support loading | yes |
| B1 | mathematica | 65 | `expectZero[(alphaLambda/.lambda->A-x) - alphaX]` | (2) loading law | yes |
| B2 | mathematica | 66 | `expectZero[(sLambda/.lambda->A-x) - sX]` | (3) overlap | yes |
| B3 | mathematica | 67 | `expectZero[(nLambda/.lambda->A-x) - nX]` | (4) normalization product | yes |
| B4 | mathematica | 75 | `expectZero[dAlphaDx - dAlphaTarget]` | (5a) dα0/dx form | yes |
| B5 | mathematica | 76 | `expectZero[sX*dAlphaDx - 1]` | (5b) `s_- = dx/dα0` | yes |
| B6 | mathematica | 91-96 | solve λ-form support, `expectZero[(gBReqLambda/.lambda->A-x) - gBReqSqOverVarpi2]` | (6) support loading | yes |

All twelve assertions are non-tautological: in each case the right-hand x-form is an independently posited candidate (the paper's claimed closed form) checked against the left-hand quantity *derived from the secular equation* (or from an independent solve, A6/B6). The substitution `λ = A − x` is the only bridge, so a wrong x-form would leave a nonzero residual. A5/B5 (`s_x · dα0/dx − 1`) is the strongest non-tautological check: `s_x` and `dalpha_dx` are constructed by completely separate algebra and the product collapsing to exactly 1 is the duality theorem, not an identity by construction.

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.txt` (mtime 2026-05-21 17:33:45)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.txt` (mtime 2026-05-21 17:33:50)
- vs scripts `...sympy_audit.py` (mtime 2026-06-03 15:59:11) and `...mathematica_audit.wl` (mtime 2026-06-03 15:59:11)

**What's wrong:**
Both saved `.txt` transcripts predate the scripts. The only difference is cosmetic banner stage-number tokens. Commit `e2a4780` ("numbering reconciliation Phase 1 (deterministic): 273 doc-only stage-label fixes") rewrote the self-banners only:
- sympy: `banner("STAGE 17.1 …")` → `banner("STAGE 34.1 …")` (and 17.2→34.2, 17.3→34.3), so the committed sympy output still prints `STAGE 17.1/17.2/17.3` (lines 3/31/47).
- mathematica: `banner["STAGE 017 …"]` → `banner["STAGE 034 …"]`, so the committed mathematica output still prints `STAGE 017` (line 3) and `Stage 034 Mathematica audit passed.` was already correct.
The commit diff confirms `no equation, value, variable, \label, or code logic changed`. Every result line in both transcripts (the alpha/s/N/dα0/support forms and all `= 0` / `PASS` residuals) matches what the current scripts produce. So the staleness affects only display labels, not any deliverable value.

**Why this matters:**
No mathematical impact. The captured residuals and closed forms are still authoritative; only the banner numbers lag the rename. A fresh run will regenerate identical math with corrected banner numbers.

**Required change:**
None for Codex — there is no script code change. The orchestrator's standard independent re-run (`exec-sympy 034`, `exec-mathematica 034`) refreshes both `.txt` transcripts; the banner tokens will then read 34.x / 034. No directive is written because there is no code edit to apply.

**Verification:**
After re-run, sympy output line 3 reads `STAGE 34.1 …` and mathematica output line 3 reads `STAGE 034 …`; all `= 0` / `PASS` lines remain unchanged.

## Independent-derivation check (Mathematica)

The `.wl` is a parallel re-derivation, not a verbatim transliteration. Both engines necessarily share the same minimal structure because there is essentially one route: build `S1`, differentiate, form `s = S1²/S1′` and `N = β0 s²/(κ0²λ)`, then test the candidate x-forms. The Mathematica script makes its own engine-native choices that distinguish it from a line-port:
- It uses `FullSimplify[…, Assumptions -> $Assumptions]` with an explicit physical assumption block (`0 < lambda < A`, `0 <= x < A`, all constants positive) at wl:29-32, where the SymPy script instead carries the domain on the symbol declarations (`positive=True`, `nonnegative=True`) at py:32-36.
- For the support loading it uses `Solve[…, Reals]` plus an explicit solution-count guard `If[Length[gBSolution] != 1, fail[...]]` (wl:85-86, 92), whereas SymPy uses `sp.solve(...)[0]` (py:78, 87) with no count check.
- `expectZero` reduces via `FullSimplify[Together[Expand[expr]]]` and tests `TrueQ[res === 0]` (wl:20-24); SymPy's `expect_zero` uses `simplify(expand(expr))` and `!= 0` (py:24-28).
These are genuinely independent normalizers and solvers acting on the same physics, which is the intent of the dual-engine policy. Not a `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on every emitted form and every residual:
- `alpha(x)`: sympy `x·(DeltaK+x)/(κ0_sq·(DeltaK+x)+κ1_sq·x)` (out line 9-12) = math `(x*(DeltaK+x))/(DeltaK*kappa0Sq+(kappa0Sq+kappa1Sq)*x)` (out line 6) — same up to factoring of the denominator.
- `s_-(x)`, `N_-(x)`, `dα0/dx`, `g_{B,req}²/ϖ²`: identical algebraic content across the two transcripts (sympy lines 13-25, 33-44, 50-55; math lines 7-8, 15, 21).
- All comparison checks: `alpha/s/N(lambda=A-x) − (x) = 0`, `dalpha/dx − target = 0`, `s_x·dα0/dx − 1 = 0`, `lambda-form vs x-form support loading = 0` — both engines report 0 / PASS (sympy lines 26-28,43-44,56-57; math lines 9-23). `engines_agree: true`.

## Verdict justification

Attacks tried and failed: (1) checked for tautology — the x-forms are independent candidates checked against secular-derived λ-forms, so a wrong x-form leaves a nonzero residual (not construction-guaranteed); A5/B5 `s_x·dα0/dx − 1` is a genuine duality theorem, not an identity. (2) Sign check on the overlap: `s_λ = S1²/S1′ = −dλ/dα0 = dx/dα0`, consistent with the notes' `s_- = −dλ_-/dα0` and the paper's `s_- = dx/dα0`; no sign error. (3) κ²-vs-κ symbol check: script's `ks0,ks1` are declared as `kappa0_sq,kappa1_sq` and used wherever the paper writes `κ0²,κ1²` — consistent. (4) Symbol domains (`0 < λ < A`, `0 ≤ x < A`, all positive constants) match the paper's stable-branch setup and do not over-assume. (5) Support-loading solve is linear in `gB²` so the unique-solution guard is sound. (6) Value reconciliation: all six symbolic deliverables present and matching in both `.tex` and `.md` (see below). The single finding is a cosmetic `stale_output` (banner number tokens only; no value impact), which the orchestrator's re-run resolves. Read paper, notes, and appendix; the script's claim matches the paper's claim exactly. Verdict `findings` solely to surface the informational stale-output; mathematically this stage is clean.

## Value Reconciliation (pass-2 augmentation)

All deliverables are symbolic (this stage emits no numeric constants — no `Pi_star`, `gamma_0`, benchmark numbers, etc.). Each boxed/labeled symbolic result emitted by the scripts is reconciled against the `.tex` card and `.md` notes below.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `x := A − λ_-` (softening-depth variable) | py:36/47 substitution, out via `lam→A−x`; wl:32 `0<=x<A` | `.tex:19` `\boxed{x:=A-\lambda_-}`; `.md:52` `x := A - lambda_-` | MATCH |
| `α_0(x) = x(x+ΔK)/[κ0²(x+ΔK)+κ1²x]` | py:48, sympy out:9-12; wl:46, math out:6 | `.tex:36-37` (eq alpha-x); `.md:78` | MATCH |
| `s_-(x) = [κ0²(x+ΔK)+κ1²x]²/[κ0²(x+ΔK)²+κ1²x²]` | py:49, sympy out:13-18; wl:50, math out:7 | `.tex:45-47` (eq s-x); `.md:98-99` | MATCH |
| `N_-(x) = β0[κ0²(x+ΔK)+κ1²x]⁴/{κ0²(A−x)[κ0²(x+ΔK)²+κ1²x²]²}` | py:50, sympy out:19-25; wl:54, math out:8 | `.tex:52-55` (eq N-x); `.md:119-120` | MATCH |
| `dα0/dx = [κ0²(x+ΔK)²+κ1²x²]/[κ0²(x+ΔK)+κ1²x]²  (>0)` | py:67, sympy out:33-43; wl:70, math out:15 | `.tex:71-73` (eq dalpha-dx); `.md:134-137` | MATCH |
| `s_- · dα0/dx = 1` (duality identity) | py:71; wl:76 | `.tex:97` (Checks item 2: `s_-=dx/d\alpha_0`); `.md:90-94` | MATCH |
| `α_mix = X²/(Ω_U²Δ0)` | py:76, sympy out:49; wl:84, math out:20 | `.tex:79` `\alpha_{\rm mix}=\frac{\mathcal X^2}{\Omega_U^2\Delta_0}`; `.md:157` | MATCH |
| `g_{B,req}²/ϖ² = x(x+ΔK)/[κ0²(x+ΔK)+κ1²x] − X²/(Ω_U²Δ0)` | py:78, sympy out:50-55; wl:87, math out:21 | `.tex:84-89` (eq support-req); `.md:169-171` | MATCH |

Internal scaffolding (accounted for, no finding expected in prose): `alpha_lam`/`alphaLambda`, `s_lam`/`sLambda`, `N_lam`/`nLambda` (λ-variable intermediate forms used only to bridge the substitution); `S1`/`s1`, `S1p`/`s1p` (secular sum and its derivative); `dalpha_target`/`dAlphaTarget` (the manifestly-positive comparison form, itself a deliverable but used as the assertion target); `gBreq_lambda`/`gBReqLambda` (λ-form support, cross-check intermediate); pass/fail flags, residual-zero values, the solution-count guards.

reconciliation: complete; 8 values checked, 0 misaligned

## Self-test notes

Traps checked: (1) Variable independence — `dα0/dx` (py:66/wl:69) differentiates `alpha_x`, which genuinely depends on `x`; not an identically-zero derivative. `S1′ = dS1/dλ` depends on λ; nonzero. (2) Symmetry/parity — N/A (no unbounded integrals; this stage is purely rational-function algebra). (3) Trivial-case — the load-bearing checks are exact symbolic-zero residuals confirmed in both transcripts (`= 0` / `PASS`), and `s_x·dα0/dx − 1 = 0` is the non-trivial duality collapse, confirmed nonzero-then-canceling rather than trivially true. (4) Path specs — N/A (no missing-script finding). (5) Paper round-trip — no fix prescribed (the only finding is a re-run, no code/derivation change), so no new paper_misalignment can be introduced. Conclusion: stage is mathematically clean; the lone finding is a cosmetic stale-output resolved by the orchestrator's standard re-run.
</content>
</invoke>
