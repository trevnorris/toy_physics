---
unit_id: 031
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md]
  paper_appendix: present
---

# Audit unit 031 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_031.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage031_selected_branch_reachability.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row at line 52; `\input{stages/stage_031}` at line 100)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.txt`

## What the paper claims

Stage 031 ("Reachability and Stable-Side Monotonicity", Part II, anchor MTDC-T5) proves that the selected lower-branch prefactor `P_{0,-} = β₀ s_-/λ_-` is strictly increasing on the stable loading interval `0 ≤ α₀ < α_crit`, turning normalization into a one-dimensional reachability theorem. The card's `\stagefield{Output}` reads verbatim: "Stage~031 outputs the derivative identities \eqref{eq:app-stage031-ds-dalpha}--\eqref{eq:app-stage031-dP-dalpha}, the refined threshold \eqref{eq:app-stage031-alpha-crit}, and the unique stable crossing theorem." The distinct deliverables are: (D1) the exact overlap derivative `ds_-/dα₀ = 2(ΔK_ax)²Π_κ/R_λ³ > 0` (boxed, eq stage031-ds-dalpha); (D2) the exact prefactor derivative `dP_{0,-}/dα₀ = β₀[(ds_-/dα₀)λ_- + s_-²]/λ_-² > 0` (boxed, eq stage031-dP-dalpha), using `dλ_-/dα₀ = -s_-`; (D3) the zero-loading start value `P_{0,-}(0) = β₀κ₀²/A` (eq stage031-P-zero); (D4) the refined softening threshold `α_crit = AB/(Bκ₀² + Aκ₁²)` (boxed, eq stage031-alpha-crit); (D5) the divergence `P_{0,-} → +∞` as `α₀ → α_crit⁻` (so the branch spans `P_{0,-}(0) ≤ P_{0,-} < +∞`); and (D6) the unique stable crossing theorem. The notes (`moving_throat_pde_stage031_*.md`) carry the same set with added derivational detail: `A = K_0 − Ξ_0`, `B = A + ΔK_ax`, `s_- = (v·e_-)²`, `λ_- = (A+B − α₀σ − R)/2`, `R = sqrt((ΔK_ax + α₀δ_κ)² + 4α₀²·KappaProd)`, the HF identity `dλ_-/dα₀ = −s_-`, the determinant identity `λ_-λ_+ = AB − α₀(Bκ₀² + Aκ₁²)`, and `λ_-(0)=A`, `s_-(0)=κ₀²`. The appendix row (line 52) summarizes: "Strict monotonicity of \(P_{0,-}\) and unique stable crossing," status `\StatusExactClosure{}`. The card states no free numeric constants — all results are exact symbolic identities (the constants `A, B, ΔK_ax, κ₀, κ₁, β₀` are upstream symbols, not pinned numbers).

## What the script claims to verify

The SymPy script (docstring lines 2-11) audits five things: the exact derivative of the selected overlap, the exact prefactor-derivative formula, initial values at α=0, the exact softening threshold + determinant factorization, and the stable-side divergence factorization. Concretely it defines `λ_∓ = (2A + DK − ασ ∓ R)/2` from the explicit `(σ, δ_κ)` parameterization (lines 39-45), DEFINES `s_- = −dλ_-/dα` via Hellmann-Feynman (line 47), then asserts: `ds_-/dα` equals the closed form `2 DK² KappaProd / R³` (line 50); the symbolic `dP_{0,-}/dα` equals `β₀[(ds)λ_- + s_-²]/λ_-²` (line 60); `λ_-(0)=A`, `s_-(0)=κ₀²(=x0)`, `P_{0,-}(0)=β₀κ₀²/A` (lines 67-69); the determinant factorization `λ_-λ_+ = A(A+DK) − α·T0` with `T0 = (A+DK)x₀ + Ax₁` (line 75); `λ_-(α_crit)=0` (lines 76-81); a threshold-radical-square identity (lines 84-88); the pole-form factorization `λ_-λ_+ = T0(α_crit − α)` (lines 97-100); and `P_{0,-} = β₀ s_- λ_+/[T0(α_crit − α)]` (lines 103-104). The Mathematica script is structurally independent: it builds the loaded 2×2 wall matrix `diag(a, a+dK) − α·v vᵀ` with `v={κ₀,κ₁}` (lines 37-38), extracts eigenvalues/eigenvectors via `Eigenvalues`/`Eigenvectors` (lines 42-53), DERIVES `s_- = (v·e_-)²/(e_-·e_-)` from the eigenvector overlap (lines 56-59), and VERIFIES the HF identity `D[λ_-,α] + s_- = 0` rather than assuming it (lines 63-67), then mirrors the same closed-form, initial-value, threshold (from `Det[M]=0`), and pole-factorization assertions (Parts II-VI).

## Paper ↔ script cross-check

| paper deliverable | script-side check (sympy / mathematica) | status |
|---|---|---|
| D1 `ds_-/dα = 2(ΔK_ax)²Π_κ/R³` | py:50 (`ds_exact − ds_expected==0`); wl:76 (`dsOverlap − dsExpected==0`) | match |
| D2 `dP_{0,-}/dα = β₀[(ds)λ_- + s_-²]/λ_-²` | py:60 (`dP_direct − dP_physical==0`); wl:87 (`dPDirect − dPClosed==0`) | match |
| D3 `P_{0,-}(0) = β₀κ₀²/A` | py:67-69 (`λ_-(0)=A`, `s_-(0)=κ₀²`, `P_{0,-}(0)=β₀κ₀²/A`); wl:45,92,93 | match |
| D4 `α_crit = AB/(Bκ₀² + Aκ₁²)` | py:73-75,81 (`T0`, det factorization, `λ_-(α_crit)=0`); wl:102-112 (`α_crit from Det[M]=0`, `λ_-(α_crit)=0`) | match |
| D5 divergence `P_{0,-}→+∞` at α_crit⁻ | py:97-104 (pole factorization `λ_-λ_+ = T0(α_crit−α)`, `P_{0,-}=β₀ s_- λ_+/[T0(α_crit−α)]`); wl:120-133 | match |
| D6 unique stable crossing (D1+D2 monotone, D3 start, D5 divergence ⇒ unique crossing) | corollary of D1/D2/D3/D5 above; no separate IVT assertion (prose theorem) | match (corollary) |

`dλ_-/dα = −s_-` (HF, used in D2) is verified explicitly in Mathematica (wl:63-67) and is the *definition* of `s_-` in SymPy (py:47) — covered between the two engines. Every script-side assertion maps to a paper deliverable; no orphan assertions, no uncovered deliverable. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 50 | `expect_zero(ds_exact − ds_expected)` | D1 | yes |
| A2 | sympy | 60 | `expect_zero(dP_direct − dP_physical)` | D2 | yes |
| A3 | sympy | 67 | `expect_zero(λ_-(0) − A)` | D3 | yes |
| A4 | sympy | 68 | `expect_zero(s_-(0) − x0)` | D3 | yes |
| A5 | sympy | 69 | `expect_zero(P_{0,-}(0) − β₀ x0/A)` | D3 | yes |
| A6 | sympy | 75 | `expect_zero(λ_-λ_+ − (A(A+DK) − α T0))` | D4 | yes |
| A7 | sympy | 81 | `expect_zero(λ_-(α_crit))` (post radical-replace) | D4 | yes |
| A8 | sympy | 85-88 | `expect_zero(T0²R²(α_crit) − root_crit²)` (radical guard) | D4 (justifies A7) | yes |
| A9 | sympy | 97-100 | `expect_zero(λ_-λ_+ − T0(α_crit − α))` | D5 | yes |
| A10 | sympy | 104 | `expect_zero(P_{0,-} − β₀ s_- λ_+/[T0(α_crit − α)])` | D5 | yes |
| M1 | mathematica | 45-46 | `λ_∓(0)` initial values from Eigenvalues | D3 | yes |
| M2 | mathematica | 67 | `expectZero(D[λ_-,α] + s_- )` (HF derived) | D2 prerequisite | yes |
| M3 | mathematica | 76 | `expectZero(dsOverlap − dsExpected)` | D1 | yes |
| M4 | mathematica | 87 | `expectZero(dPDirect − dPClosed)` | D2 | yes |
| M5 | mathematica | 92-93 | `s_-(0)=κ₀²`, `P_{0,-}(0)=β₀κ₀²/a` | D3 | yes |
| M6 | mathematica | 107 | `expectZero(α_crit(Det=0) − closed)` | D4 | yes |
| M7 | mathematica | 110-112 | `expectZero(λ_-(α_crit))` | D4 | yes |
| M8 | mathematica | 120-122 | `expectZero(Det[M] − T0(α_crit − α))` | D5 | yes |
| M9 | mathematica | 123-125 | `expectZero(λ_-λ_+ − Det[M])` | D5 | yes |
| M10 | mathematica | 131-133 | `expectZero(P_{0,-} − β₀ s_- λ_+/[T0(α_crit − α)])` | D5 | yes |

Every assertion is `simplify/FullSimplify(lhs − rhs) == 0` against an independently-formed target, not `x == x`. No tautology row; no orphan row.

## Findings

None. (One informational `stale_output` observation is recorded below in the verdict justification; it does not affect the substantive verdict and produces no directive.)

### Attacks attempted and survived

1. **Symbol-convention mismatch between engines (κ vs κ²).** The two scripts use different vocabularies: SymPy `x0, x1` are the *squared* kappas (`KappaProd = x0·x1 = κ₀²κ₁²`, `delta_kappa = x0 − x1 = κ₀²−κ₁²`, py:40-41), whereas Mathematica `kappa0, kappa1` are the *raw* kappas and the squaring happens in-formula (`rSquared = (dK + α(κ₀²−κ₁²))² + 4α²κ₀²κ₁²`, wl:72). I checked that both nonetheless land on the same paper object: paper notes line 36 define `R = sqrt((ΔK_ax + α₀δ_κ)² + 4α₀²·KappaProd)` with `KappaProd = κ₀²κ₁²`, and both `R²` expressions are identical under the mapping x0↔κ₀², x1↔κ₁². The `ds = 2 DK² κ₀²κ₁²/R³` form is identical in both (py:49, wl:74). No factor-of-κ error; the convention difference is consistent and is itself evidence of independent derivation.

2. **`s_-` defined via HF in SymPy — possible circular `ds` check.** SymPy DEFINES `s_- = −dλ_-/dα` (py:47), so `ds = diff(s_-,α) = −d²λ_-/dα²`; asserting it equals `2DK²·KappaProd/R³` is a real, non-tautological closed-form check on the second derivative of the explicitly-built `λ_-`. The potential gap (does the HF-defined `s_-` actually equal the eigenvector overlap `(v·e_-)²` the notes claim?) is closed by the Mathematica side, which computes `s_-` from `Eigenvectors` (wl:56-59) and then *verifies* `D[λ_-,α] + s_- = 0` as a derived identity (wl:67). So the HF relation is assumed on one engine and proved on the other — the claim is fully exercised between them.

3. **Branch-selection trap in `λ_-`.** SymPy picks the `−R` root explicitly (py:44). Mathematica picks the lower branch by `Sort[Eigenvalues]` (line 42) AND by identifying the root that reduces to `a` at α=0 (M1 checks `λ_-(0)=a`, `λ_+(0)=a+dK`, lines 45-46). Both select the same (stable, lower) branch; the initial-value checks pin the branch identity. No wrong-branch substitution.

4. **The manual radical replace in `λ_-(α_crit)` (py:77-80).** SymPy substitutes α=α_crit then replaces `sqrt(base)` with `root_crit = A²x₁ + (A+DK)²x₀` wherever `base − root_crit² == 0`. This is the one hand-guided step. It is justified two ways: (a) `root_crit > 0` (all symbols positive), so replacing the principal square root of a perfect square by its positive root is valid; (b) the separate "threshold radical square identity" (A8, py:84-88) independently confirms `T0²·R²(α_crit) = root_crit²` symbolically. Moreover the Mathematica side reaches `λ_-(α_crit)=0` (M7, wl:110-112) with NO manual replace — pure `FullSimplify` of the eigenvalue at the closed-form `α_crit`. The replace hides nothing.

5. **Threshold `T0`/`α_crit` factor check.** Paper `α_crit = AB/(Bκ₀² + Aκ₁²)`, `B = A+ΔK_ax`. SymPy `T0 = (A+DK)x₀ + Ax₁ = Bκ₀² + Aκ₁²`, `α_crit = A(A+DK)/T0 = AB/T0` (py:73-74) — exact match. Mathematica `alphaCritClosed = a(a+dK)/((a+dK)κ₀² + aκ₁²)` (wl:102-105) — exact match. Paper determinant identity `λ_-λ_+ = AB − α₀(Bκ₀² + Aκ₁²)` matches py:75 and wl:120-124. No factor-of-2 or sign slip.

6. **Symbol-domain / positivity assumptions.** SymPy: `A, DK, x0, x1, β₀` positive, `α` nonnegative, all real (py:34-37). Mathematica: identical (`a,dK,kappa0,kappa1,beta0 > 0`, `alpha ≥ 0`, all Reals, wl:29-30). These match the paper's stated stable-interval setup (`A = K_0 − Ξ_0 > 0`, `β₀ > 0`, `0 ≤ α₀ < α_crit`). The positivity is exactly what makes `R > 0`, `λ_- > 0` on-branch, and `root_crit > 0` — no over-strong assumption masks a branch; no assumption contradicts the physical setup.

7. **Variable-independence trap (self-test #1).** Every `diff` is taken w.r.t. `α`, and `λ_-`, `s_-`, `P_{0,-}`, `R` all genuinely depend on `α` (py:42-56; wl:42-81). No identically-zero derivative; the `expect_zero(... )` checks are residual-vanishing identities, not zero-derivative artifacts. There are no `assert_nonzero` controls in this stage, so the "trivially-passing assert_zero / failing assert_nonzero" failure mode does not arise.

8. **Monotonicity `> 0` not asserted as an inequality.** The card boxes D1/D2 as identities and derives `> 0` in the `\stagefield{Checks}` list by inspection of the closed form (nonneg numerator terms, positive denominator). The scripts verify the *identities* (D1, D2) exactly; positivity is a manifest corollary of the verified closed forms under the declared positive-symbol domain. This matches the paper's own treatment (the boxed equations are identities; `>0` is annotation), so it is not an `insufficient_verification` finding.

## Independent-derivation check (Mathematica)

The `.wl` is a genuinely independent second engine, not a transliteration. Evidence:
- **Different starting object.** SymPy hard-codes the explicit branch eigenvalues `λ_∓ = (2A + DK − ασ ∓ R)/2` (py:44-45). Mathematica builds the physical matrix `wallMatrix = {{a,0},{0,a+dK}} − α·Outer[Times,{κ₀,κ₁},{κ₀,κ₁}]` (wl:37-38) and extracts eigenvalues via `Sort[Eigenvalues[wallMatrix]]` (wl:42) — it never writes the closed-form `λ_∓` by hand.
- **`s_-` derived two different ways.** SymPy DEFINES `s_- = −dλ_-/dα` (py:47, i.e. assumes HF). Mathematica DERIVES `s_- = (v·e_-)²/(e_-·e_-)` from `Eigenvectors` (wl:56-59) and then *proves* HF as a check (`expectZero[D[λ_-,α] + s_-]`, wl:67). The HF identity is the script's own assertion M2, not a definition — this is the strongest non-transliteration signal.
- **`α_crit` derived two different ways.** SymPy pins `α_crit = A(A+DK)/T0` and back-checks the determinant factorization (py:74-75). Mathematica solves `Solve[Det[wallMatrix]==0, α]` and confirms the root equals the closed form (wl:98-107).
- **Distinct symbol convention** (raw κ vs squared x), noted in attack #1.

Conclusion: independent. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines pass every check (all residuals `= 0` in both transcripts). The shared closed forms agree once the κ vs κ² convention is mapped:
- `ds_-/dα`: py output line 5 `= 0` with displayed form `2·DK²·x₀·x₁/(4α²x₀x₁ + (DK+α(x₀−x₁))²)^{3/2}`; wl output line 15-17 `= 0` with `2 dK²κ₀²κ₁²/(4α²κ₀²κ₁² + (dK + α(κ₀−κ₁)(κ₀+κ₁))²)^{3/2}`. Identical under x₀=κ₀², x₁=κ₁² (note `(κ₀−κ₁)(κ₀+κ₁)=κ₀²−κ₁²`).
- `α_crit`: py output lines 72-75 `A(A+DK)/(Ax₁ + x₀(A+DK))`; wl output line 39 `a(a+dK)/((a+dK)κ₀² + aκ₁²)`. Identical.
- `λ_-(α_crit)=0`, det factorization `=0`, `λ_-λ_+ − Det = 0`, pole factorization `=0`, `P_{0,-}` factorization `=0`: all `0`/PASS in both transcripts (py lines 62-104; wl lines 5-52).
Engines agree.

## Verdict justification

`clean`. I read the paper card, the per-stage notes, and the appendix row (line 52, `\StatusExactClosure{}`), and built the six-deliverable model (D1 overlap derivative, D2 prefactor derivative, D3 start value, D4 threshold, D5 divergence factorization, D6 unique-crossing corollary). Each deliverable maps to a non-tautological, anchored assertion present in BOTH engines. I attacked: the κ/κ² convention difference (consistent, maps cleanly), the HF-as-definition circularity in SymPy (closed by the Mathematica eigenvector derivation + HF proof), the branch selection (pinned by initial values on both sides), the hand-guided radical replace (independently justified by A8 and corroborated replace-free by M7), the threshold factors/signs (exact match to paper), and the symbol-domain positivity (matches the stated stable interval, no masking). All survived. The card states no free numeric constants, so there is nothing to value-mismatch. **One informational note:** the SymPy `.py` mtime (2026-06-03 15:59) is newer than the SymPy `.txt` output (2026-05-26 00:43); git shows the only post-output change was a banner-string edit `"STAGE 14 AUDIT COMPLETE"` → `"STAGE 31 AUDIT COMPLETE"` (commit e2a4780, label-only numbering reconciliation; py:108). The saved transcript still prints `STAGE 14 AUDIT COMPLETE` (sympy.txt:120). No equation, value, assertion, or symbol changed — the captured numeric/symbolic content remains faithful to the current script. This is the standard `stale_output` signal (the orchestrator's independent re-run will refresh the banner line); it is informational, not blocking, and not a content finding. The Mathematica `.wl`/`.txt` pair is fresh (txt 2026-05-26 00:43:53 newer than wl 2026-05-25 23:42:57). Both engines pass independently and agree.

## Self-test notes

Checked: (#1) variable independence — every `diff`/`D` is w.r.t. `α`, on which `λ_-`, `s_-`, `P_{0,-}`, `R` genuinely depend; no zero-derivative trap, and there are no `assert_nonzero` controls to mis-fire. (#2) parity — N/A, no integrals; all checks are algebraic eigenvalue/rational identities. (#3) trivial-case pre-check — at α=0 the initial-value residuals (`λ_-(0)−A`, `s_-(0)−κ₀²`, `P_{0,-}(0)−β₀κ₀²/A`) reduce to 0 as asserted (confirmed in both transcripts). (#4 path spec) — N/A, both scripts present. (#5 paper round-trip) — no fix prescribed, no directive written, so no new `paper_misalignment` introduced. No directive written (zero findings).

## Value Reconciliation (pass-2 augmentation)

The scripts emit purely **symbolic** closed-form results; there are no named numeric constants, benchmarks, or figures of merit (`A, B, ΔK_ax, κ₀, κ₁, β₀, α₀` are upstream symbols, not pinned numbers, and the stage is `\StatusExactClosure{}` algebraic). Each emitted deliverable value is one of the card's boxed/displayed identities and is reproduced in BOTH the `.tex` card and the `.md` notes. Every deliverable reconciles (MATCH); no MISMATCH, no MISSING-DELIVERABLE.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `ds_-/dα = 2(ΔK_ax)²κ₀²κ₁²/R³` (D1) | py:49-50 / py.txt:5-12; wl:74-76 / wl.txt:15-17 | tex:17-22 (boxed eq stage031-ds-dalpha); md:38 | MATCH |
| `dP_{0,-}/dα = β₀[(ds)λ_- + s_-²]/λ_-²` (D2) | py:59-60 / py.txt:17; wl:83-87 / wl.txt:22-24 | tex:24-30 (boxed eq stage031-dP-dalpha); md:67-68 | MATCH |
| `dλ_-/dα = −s_-` (HF, supports D2) | py:47 (definition); wl:63-67 / wl.txt:9-10 (derived) | md:61-62 | MATCH |
| `λ_-(0) = A` (D3) | py:67 / py.txt:62; wl:45 / wl.txt:5-6 | md:92 ("λ_-(0) = A = K_0 − Ξ_0") | MATCH |
| `s_-(0) = κ₀²` (D3) | py:68 / py.txt:63; wl:92 / wl.txt:29-30 | md:94 ("s_-(0) = kappa_0^2") | MATCH |
| `P_{0,-}(0) = β₀κ₀²/A` (D3) | py:69 / py.txt:64; wl:93 / wl.txt:31-32 | tex:35-38 (eq stage031-P-zero); md:98 | MATCH |
| `α_crit = AB/(Bκ₀² + Aκ₁²)` (D4) | py:73-74 / py.txt:72-75; wl:102-107 / wl.txt:37-39 | tex:40-44 (boxed eq stage031-alpha-crit); md:112 | MATCH |
| `λ_-(α_crit) = 0` (D4) | py:77-81 / py.txt:70; wl:110-112 / wl.txt:40-41 | tex:45 (λ_-→0⁺); md:116 | MATCH |
| det identity `λ_-λ_+ = AB − α₀(Bκ₀²+Aκ₁²)` (D4/D5) | py:75 / py.txt:69; wl:120-124 / wl.txt:46-49 | md:107 ("λ_-λ_+ = AB − α₀(Bκ₀²+Aκ₁²)") | MATCH |
| pole form `P_{0,-} = β₀ s_- λ_+/[T0(α_crit−α)]`, `P_{0,-}→+∞` (D5) | py:103-104 / py.txt:93; wl:127-133 / wl.txt:50-52 | tex:46-53 (P_{0,-}→+∞, range); md:122-132 | MATCH |
| unique stable crossing (D6) | corollary of above (no separate assertion) | tex:55-65; md:136-158 | MATCH (corollary) |

Internal scaffolding (no prose expectation): pass/fail flags (`PASS:`/`FAIL:`), residual-zero check printouts, `T0`/`root_crit` intermediates, the `λ_+` branch and `λ_+(α_crit)`/`lam_+(0)` checkpoints, the threshold-radical-square guard identity (py:84-88), and the eigenvector normalization `(e_-·e_-)` in Mathematica.

reconciliation: complete; 11 deliverable values checked, 0 misaligned. (All emitted results are exact symbolic identities boxed/stated in the card and reproduced in the notes; no numeric constants exist to mismatch.)
