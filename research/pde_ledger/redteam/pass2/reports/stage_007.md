---
unit_id: 007
batch: I.1
auditor_model: claude-opus-4-8
audit_date: 2026-06-04T00:00:00Z
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 007 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_007.tex`
- notes: `(none)` — no `notes/stages/moving_throat_pde_stage007_*.md` exist (confirmed by directory listing)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 36 + `\input` line 93)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.txt`

## What the paper claims

Stage 007 ("Projection/reduction comparison", Part I anchor MTDC-T4) compares the
exact projection-first projected Maxwell law against the zero-mode reduction
channel, keeping observer dependence explicit. For a zero-mode ansatz
`A_mu(x,w)=a_mu(x)` the projected equation carries the weighted factors
`I_WZ = ∫WZ dw`, `I_WH = ∫WH dw`, `I_WS = ∫WS dw`
(eq:stage007-projected-integrals). The card's `\stagefield{Output}` is verbatim:
"Stage~007 exports the observer-dependent projected coefficients
\eqref{eq:stage007-effective-parameters} and marks the reduction-first brane law
as a matched channel rather than the controlling derivation." Those effective
parameters (eq:stage007-effective-parameters) are
`mu_eff^proj = mu0 · I_WS/I_WZ` and `1/xi_eff^proj = I_WH/(xi·I_WZ)` (i.e.
`xi_eff^proj = xi·I_WZ/I_WH`), declared as observer-kernel quantities until a
matching prescription is imposed. The "Mismatch guard" paragraph asserts that a
delta-like source and a smooth localized gauge profile do NOT automatically
reproduce the reduced brane coupling — source matching is a separate condition.
The appendix row 36 summarizes the deliverable as "Matched Gaussian zero-mode
channel and its coupling/gauge-parameter mismatch relative to the projection
law," status `\StatusExact{} / \StatusReduced{}`. The card states no numeric
constants of its own beyond the symbolic effective-parameter forms.

Distinct deliverables:
1. The projected-parameter formulas `mu_eff^proj = mu0·I_WS/I_WZ` and
   `xi_eff^proj = xi·I_WZ/I_WH` (the boxed Output).
2. The qualitative claim that these are observer-dependent (kernel-dependent)
   until a matching prescription is fixed.
3. The mismatch guard: matched-Gaussian-observer + delta-source projection does
   not equal the reduction-first coupling `mu0/Z_int`; there is a definite
   coupling/gauge-parameter mismatch factor.

## What the script claims to verify

The SymPy script (docstring lines 1-24) sets up the zero-mode comparison and
verifies, with concrete Gaussian profiles `Z=exp(-w²/λ²)`, `H=exp(-w²/ρ²)`, and
normalized smooth observer/source kernels, that (a) the projection-formula
overlaps `I_WZ`, `I_WH`, `I_WS` evaluate to the claimed closed forms; (b) the
projected `xi_eff` formula equals `xi·I_WZ/I_WH`; (c) a `w`-dependent field
mutation and a `w`-dependent source mutation each break the factorized
zero-mode projection by an explicit nonzero Gaussian-moment amount (mismatch
guard); (d) for the matched observer kernel `W=Z/Z_int` with a delta-localized
source, `mu0_eff^(proj)/mu0_eff^(red) = √2` and the matched gauge ratio is
`√(λ²+ρ²)/(√2·λ)`, reducing to 1 at `H=Z` (`ρ→λ`); (e) a regularized sharp
observer/source reproduces the sampling limits `I_WZ→1`, `I_WH→1`, the
projected `xi_eff→xi`, and shows `I_WS` diverges (no finite delta/delta
coupling). The Mathematica script independently evaluates the same Gaussian
integrals and asserts the same closed forms via numbered M1–M11c residual checks
with `Exit[1]` on any nonzero residual.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `mu_eff^proj = mu0·I_WS/I_WZ` (Output eq) | sympy `mu0_proj_match = mu0*I_WS_match/I_WZ_match` (line 148) + smooth `I_WS`/`I_WZ` overlaps (lines 90-97); wl M3/M4/M7 | match |
| `xi_eff^proj = xi·I_WZ/I_WH` (Output eq) | sympy "smooth xi_eff projection formula" (lines 102-105) + matched/regularized variants (lines 150, 187); wl M4c/M10c | match |
| Observer-dependent (kernel-dependent) | smooth/matched/eps kernels give distinct overlaps (sec 1-4); confirmed numerically distinct | match |
| Mismatch guard: matched+delta ≠ reduction coupling | sympy "delta-source ratio = √2" (line 165), `assert_nonzero(... ratio - 1)` (line 174); wl M8 | match |
| Gauge-parameter mismatch vs reduction | sympy "matched gauge ratio to reduction-first" (lines 166-169), `H=Z` collapse to 1 (lines 170-173); wl M8b/M8c | match |
| w-dependent mutations break factorization | sympy field/source mutation amplitudes + `assert_nonzero` (lines 111-128); wl M5/M6 | match |

`paper_alignment = aligned`. Every paper-side deliverable maps to a
non-tautological script-side check, and the script's effective-parameter forms
(`mu0·I_WS/I_WZ`, `xi·I_WZ/I_WH`) reproduce the card's eq:stage007-effective-parameters
exactly, including the `1/xi_eff = I_WH/(xi·I_WZ)` reciprocal convention.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 90-93 | `I_WZ_smooth - λ/√(λ²+σ²) == 0` | mu_eff overlap (I_WZ) | yes |
| A2 | sympy | 95-97 | `I_WS_smooth - 1/(√π√(σ²+τ²)) == 0` | mu_eff overlap (I_WS) | yes |
| A3 | sympy | 98-101 | `I_WH_smooth - ρ/√(ρ²+σ²) == 0` | xi_eff overlap (I_WH) | yes |
| A4 | sympy | 102-105 | `xi·I_WZ/I_WH - xi·λ√(ρ²+σ²)/(ρ√(λ²+σ²)) == 0` | xi_eff^proj formula | yes |
| A5 | sympy | 106-109 | projected-residual identity = `I_WZ·df - mu0·I_WS·j` | projected law structure | yes |
| A6 | sympy | 114-117 | field-mutation delta = `η λ³σ²/(2(λ²+σ²)^{3/2})` | mismatch guard | yes |
| A7 | sympy | 118 | `assert_nonzero(field mutation delta)` | mismatch guard (nonzero) | yes |
| A8 | sympy | 123-127 | source-mutation delta = `η mu0 x σ²τ²/(2√π(σ²+τ²)^{3/2})` | mismatch guard | yes |
| A9 | sympy | 128 | `assert_nonzero(source mutation delta)` | mismatch guard (nonzero) | yes |
| A10 | sympy | 138-140 | `Z_int=√πλ`, `Z2_int=√2√πλ/2`, `H_int=√πρ` | localization data | yes |
| A11 | sympy | 162-164 | `I_WZ_match = Z2_int/Z_int = √2/2`; `I_WH_match=ρ/√(λ²+ρ²)` | matched channel | yes |
| A12 | sympy | 165 | `mu0_proj_match/mu0_red - √2 == 0` | mismatch (coupling) | yes |
| A13 | sympy | 166-169 | matched gauge ratio = `√(λ²+ρ²)/(√2·λ)` | mismatch (gauge) | yes |
| A14 | sympy | 170-173 | `(xi ratio).subs(ρ,λ) - 1 == 0` | H=Z collapse | yes |
| A15 | sympy | 174 | `assert_nonzero(ratio - 1)` | mismatch guard (nonzero) | yes |
| A16 | sympy | 198-206 | regularized limits `I_WZ→1`,`I_WH→1`, `xi_eff→xi`, `I_WS=√2/(2√πε)` | sampling/divergence | yes |
| A17 | sympy | 207-210 | `assert_nonzero(I_WZ(0) - √2/2)` (mutation guard) | observer-dep guard | yes |
| M1–M11c | math | 21-243 | `FullSimplify[... target] === 0` / limit checks | mirror of A1–A16 | yes |

Every "Exercises" entry traces to a specific paper deliverable; no orphaned
assertions and no script-side checks that the paper does not mention.

## Findings

None. (All assertions hold up under attack; paper alignment exact.)

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration. It defines `Z`, `Hprofile`, `Wsmooth`,
`Ssmooth`, `Wmatch`, `Weps` and independently evaluates each overlap via
`Integrate[...,{w,-Infinity,Infinity}]` (e.g. `IWZsmooth = Integrate[Wsmooth[w]
Z[w], ...]`, line 47), then compares against the closed-form target with
`FullSimplify[... === 0]`. The shared targets (`λ/√(λ²+σ²)`, `√2/2`, etc.) are
the physics results, not borrowed algebra; both engines reach them by evaluating
the Gaussian integral from scratch rather than echoing each other's intermediate
steps. Representative correspondences:
- SymPy line 84 `I_WZ_smooth = sp.integrate(W_smooth*Z, (w,-oo,oo))` ↔ wl line 47
  `IWZsmooth = Integrate[Wsmooth[w] Z[w], {w,-Infinity,Infinity}]` — same
  integrand, independent evaluation; neither hardcodes the result.
- SymPy line 147 `I_WS_match = W_match.subs(w,0)` (delta-sampling) ↔ wl line 139
  `IWSmatch = Wmatch[0]` — both encode `∫δ(w)W=W(0)` directly; consistent
  convention, independently written.
- SymPy line 195 `sp.limit(I_WZ_eps, eps, 0, dir='+')` ↔ wl line 211
  `Limit[IWZeps, epsilon->0, Direction->"FromAbove"]` — same limit, native idiom
  each side.
The numbered M-series in the `.wl` is the project's standard mirror layout, not a
line-by-line port. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree. SymPy output reports `Z_int = sqrt(pi)*lam`,
`Z2_int = sqrt(2)*sqrt(pi)*lam/2`, `H_int = sqrt(pi)*rho`, `I_WZ(match)=sqrt(2)/2`,
`I_WH(match)=rho/sqrt(lam**2+rho**2)`, `ratio proj/reduction = sqrt(2)`,
`xi ratio = sqrt(2*lam**2+2*rho**2)/(2*lam)` (= `√(λ²+ρ²)/(√2·λ)`), regularized
`I_WZ(eps→0)=1`, `I_WH(eps→0)=1`, `STATUS: PASS`. Mathematica output shows every
residual M1–M11c = 0 and `STATUS: PASS`. The matched coupling ratio (√2), the
matched gauge ratio (`√(λ²+ρ²)/(√2 λ)`), the `H=Z` collapse to 1, and all sharp
limits coincide across engines. Note: SymPy line 185 evaluates
`integrate(W_eps*S_eps)` with `S_eps=W_eps`, while wl line 173 evaluates
`Integrate[Weps^2]` — algebraically the same integrand, both yielding
`√2/(2√π ε)`; this is a benign cosmetic difference, not an engine disagreement.

## Verdict justification

Clean. I attacked: (1) the Gaussian integrals — hand-verified `Z_int=√πλ`,
`Z2_int=√(π/2)λ=√2√πλ/2`, `I_WZ_smooth=λ/√(λ²+σ²)`, `I_WS_smooth=1/(√π√(σ²+τ²))`,
`I_WH_smooth=ρ/√(ρ²+σ²)`, `I_WZ_match=√2/2`, `I_WS_match=1/(√πλ)`,
`I_WS_eps=√2/(2√πε)` — all correct. (2) Tautology risk — each `assert_zero`
compares an integral SymPy evaluates symbolically against an independently
written closed form; none is `x==x`. (3) The delta-sampling shortcut
(`W_match.subs(w,0)`) — legitimate `∫δW=W(0)`, not a hidden assumption. (4)
Parameter coverage — smooth (σ,τ), matched (W=Z/Z_int), and regularized (ε)
channels each exercised; mutation guards use `assert_nonzero` so the mismatch
claim genuinely fails-loud if the moment vanished. (5) Paper round-trip — the
script's `mu_eff=mu0·I_WS/I_WZ` and `xi_eff=xi·I_WZ/I_WH` reproduce
eq:stage007-effective-parameters exactly, including the reciprocal-`xi`
convention; the matched √2 coupling mismatch and gauge mismatch realize the
appendix-row deliverable. Outputs are fresh (sympy/wl `.txt` mtimes 2026-05-25
17:24/17:29, both newer than the scripts at 02:14/02:15). I read the paper card,
confirmed no stage007 notes exist, and read the part-01 appendix row; the
script's claim matches the paper's claim. No `paper_misalignment`, no stale
output, no transliteration, no insufficient verification.

## Value Reconciliation (pass-2 augmentation)

Every RESULT/deliverable value emitted by the two scripts (per source + the
committed saved outputs) is enumerated below and located in the paper card
(`.tex`). No stage007 notes file exists, so the `.tex` card (which carries the
symbolic effective-parameter forms) is the natural carrier; per the augmentation
guards, a terse card legitimately omits intermediate Gaussian overlaps, so those
intermediates are classed INTERNAL (scaffolding that drives the deliverable
checks), not MISSING.

| value | source (py line / wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `mu_eff^proj = mu0·I_WS/I_WZ` | py docstring l.21, l.67, out l.16; wl mu0ProjMatch l.140 | stage_007.tex:25 (eq:stage007-effective-parameters) | MATCH |
| `xi_eff^proj = xi·I_WZ/I_WH` | py docstring l.22, l.68, out l.17; wl xiProjMatch l.142 | stage_007.tex:27 (`1/xi_eff = I_WH/(xi I_WZ)`) | MATCH |
| matched coupling ratio proj/red = `√2` | py l.165, out l.44 "ratio proj/reduction = sqrt(2)"; wl M8 l.146 | stage_007.tex:36 (appendix row 36, "coupling … mismatch") + card "Output" mismatch claim | MATCH (qualitative; magnitude in appendix-row prose) |
| matched gauge ratio = `√(λ²+ρ²)/(√2 λ)` | py l.166-169, out l.47; wl M8b l.154 | stage_007.tex:36 (appendix row, "gauge-parameter mismatch") | MATCH (qualitative; magnitude is scaffolding for the mismatch claim) |
| H=Z collapse: gauge ratio → 1 | py l.170-173, out (implicit); wl M8c l.163 | stage_007.tex (Mismatch-guard / matched-channel narrative) | MATCH (qualitative) |
| projection vs reduction: `mu0/Z2_int` not `mu0/Z_int` | py out l.49-50; wl M7/M8 | stage_007.tex:31-35 (Mismatch guard paragraph) | MATCH (qualitative) |
| `mu0_eff^(red) = mu0/Z_int` (reference) | py docstring l.23, out l.43; wl mu0Red l.141 | (reduction-first reference; card cites the brane law, not the value) | INTERNAL/MATCH (comparison reference) |
| `xi_eff^(red) = xi·Z_int/H_int = xi λ/ρ` (reference) | py l.72, out l.46; wl xiRed l.143 | (reduction-first reference) | INTERNAL (comparison reference) |

INTERNAL scaffolding values (Gaussian overlaps and limits that exist only to
drive the deliverable checks; not expected in the terse card): `Z_int=√πλ`,
`Z2_int=√2√πλ/2`, `H_int=√πρ`, `Z2_int/Z_int=√2/2`, `I_WZ_smooth=λ/√(λ²+σ²)`,
`I_WS_smooth=1/(√π√(σ²+τ²))`, `I_WH_smooth=ρ/√(ρ²+σ²)`,
`xi_eff^(proj,smooth)=xiλ√(ρ²+σ²)/(ρ√(λ²+σ²))`, `I_WZ_match=√2/2`,
`I_WH_match=ρ/√(λ²+ρ²)`, `I_WS_match=1/(√πλ)`, field-mutation moment
`ηλ³σ²/(2(λ²+σ²)^{3/2})`, source-mutation moment
`η mu0 x σ²τ²/(2√π(σ²+τ²)^{3/2})`, `I_WZ_eps=λ/√(ε²+λ²)`, `I_WH_eps=ρ/√(ε²+ρ²)`,
`I_WS_eps=√2/(2√πε)`, `xi_eff^(proj,eps)=xiλ√(ε²+ρ²)/(ρ√(ε²+λ²))`,
`I_WZ(ε→0)=1`, `I_WH(ε→0)=1`, `xi_eff(ε→0)=xi`.

The two stage-level deliverable forms (`mu_eff^proj`, `xi_eff^proj`) reconcile
exactly to the card's eq:stage007-effective-parameters; the matched coupling and
gauge mismatches reconcile to the appendix-row deliverable and the card's
Mismatch-guard / "matched channel" narrative (the card is intentionally symbolic
and states no numeric `√2`, which the augmentation guards permit since the
magnitude is scaffolding that substantiates the qualitative mismatch deliverable,
not a `\stagefield{Output}` quantity in its own right).

reconciliation: complete; 8 deliverable-level values checked, 0 misaligned

## Self-test notes

Variable-independence: the mutation amplitudes (A6/A8) differentiate
`F_mut = f_test + η x w²` w.r.t. `x`, which `f_test` (and the η-term via `x`)
genuinely depend on, so the derivative is non-trivial and the matching
`assert_nonzero` (A7/A9) is not vacuous. Symmetry/parity: all integrands are
even Gaussians (and even×even products) over a symmetric domain, so the nonzero
overlaps are legitimate; the mutation `w²` factor preserves evenness, so its
Gaussian second moment is the nonzero `λ³σ²/(2(λ²+σ²)^{3/2})` form, confirmed by
hand. Trivial-case: the matched coupling ratio is `√2≠1`, the regularized
`I_WS=√2/(2√πε)` diverges as ε→0, and the sharp-sampling limits give exactly 1,
all consistent with the asserted literals; no `assert_zero` reduces to a
trivially-true identity. No directive written (zero findings).
