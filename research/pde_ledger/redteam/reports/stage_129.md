---
unit_id: 129
batch: IV.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage129_mouth_boundary_layer.md]
  paper_appendix: present
---

# Audit unit 129 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_129.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage129_mouth_boundary_layer.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only an `\input{stages/stage_129}` line at L1292; no narrative row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.txt`

## What the paper claims

The paper card (`stage_129.tex` L15-17) states the headline result verbatim as a quote:
"Electrochemical zero-flux law yields the exponential source profile \(\sigma_\Pi\)."

The notes (`moving_throat_pde_stage129_mouth_boundary_layer.md`) supply the
concrete form. Starting from the GNLS+localized-Maxwell free energy
`F_mouth[σ] = ∫_0^L [Θ_σ σ(ln(σ/σ_*) − 1) + V_m(z) σ] dz` with `V_m(z) ≈ V_1 z`
and the linear Onsager current `J_σ = −M_σ σ ∂_z μ_σ^chem = −M_σ(Θ_σ σ' + V_1 σ)`,
imposing the stationary zero-flux branch `J_σ = 0` gives `Θ_σ σ'(z) + V_1 σ(z) = 0`.
The notes then derive the normalized profile (L98-120):
`σ_Π(z) = Π e^{−Π z/L} / (L(1 − e^{−Π}))` with `Π := V_1 L / Θ_σ > 0`,
normalized by `∫_0^L σ(z) dz = 1`. The single dimensionless bias is
`Π_m = V_1 L / Θ_σ` (L131-141).

Deliverables: (D1) the profile `σ_Π(z)` has the stated exponential form;
(D2) it normalizes to unity on `[0,L]`; (D3) it satisfies the stationary
zero-flux ODE `Θ_σ σ' + V_1 σ = 0`; (D4) the identification
`Π = V_1 L / Θ_σ` is the load-bearing definition tying the two parameter
groups together.

## What the script claims to verify

Both engines do the same four steps. The SymPy script
(`...stage129...sympy_audit.py:13-31`) defines the exponential profile
`sigma = Pi*exp(-Pi*z/L)/(L*(1-exp(-Pi)))`, prints it, computes
`∫_0^L sigma dz` and asserts it equals 1 (L26-27), substitutes
`V_1 → Pi*Theta/L` into `J_σ = −M(Θ σ' + V_1 σ)` and asserts the result
is zero (L19, L28-29), and finally asserts the ODE residual
`Θ σ' + V_1 σ` vanishes under the same substitution (L23, L30-31).
The Mathematica script (`...stage129...mathematica_audit.wl:33-46`) runs
the identical three checks via `expectZero`. Neither script back-solves
the ODE to derive the exponential form — they verify that the closed-form
σ written down in the paper satisfies (a) normalization, (b) zero flux,
and (c) the boundary-layer ODE with the stated bias `Π = V_1 L / Θ_σ`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1: exponential profile `σ_Π(z) = Π e^{−Πz/L} / (L(1−e^{−Π}))` | Both scripts define `sigma = Pi*exp(-Pi*z/L)/(L*(1-exp(-Pi)))` (py L13; wl L33) | match (verbatim form) |
| D2: normalization `∫_0^L σ dz = 1` | sp.integrate(sigma,(z,0,L)) == 1 (py L17, L26-27); Integrate[sigma,{z,0,lM}]-1 == 0 (wl L37-38) | match |
| D3: zero-flux ODE `Θ σ' + V_1 σ = 0` under `V_1 = Π Θ/L` | residual after sub `V1 → Pi*Theta/L` simplifies to 0 (py L23-24, L30-31; wl L44-46) | match |
| D4: bias `Π = V_1 L / Θ_σ` | substitution `V1 → Pi*Theta/L` IS the identification; final print restates `Pi_m = V1*L/Theta` (py L34; wl L50) | match |

Paper alignment is `aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | L26-27 | `if simplify(∫_0^L sigma − 1) != 0: raise` | D2 (normalization) | yes |
| A2 | sympy | L28-29 | `if simplify(J_sub) != 0: raise` (after `V1 → Pi*Theta/L`) | D3 + D4 (zero-flux at the stated bias) | yes |
| A3 | sympy | L30-31 | `if simplify(res.subs(V1: Pi*Theta/L)) != 0: raise` | D3 (ODE residual) | yes (this is algebraically the same identity as A2 up to the `−M` prefactor, so it is partly redundant; still non-tautological) |
| A4 | mathematica | L38 | `expectZero["profile normalization", Integrate[sigma,{z,0,lM}]-1]` | D2 | yes |
| A5 | mathematica | L40-42 | `expectZero["zero-flux current", jSub]` after `v1 → piM*thetaSigma/lM` | D3 + D4 | yes |
| A6 | mathematica | L44-46 | `expectZero["boundary-layer ODE residual", residual]` after `v1 → piM*thetaSigma/lM` | D3 | yes (parallel to A3) |

No assertion is tautological by construction. The exponential form is
declared up-front (not algebraically guaranteed to satisfy the ODE), and
the substitution `V_1 = Π Θ/L` is the load-bearing identification: had
the script used `V_1 = 2 Π Θ/L`, the residual would be `(V_1 − Π Θ/L) σ ≠ 0`
and the assertion would fail. The normalization check is a real integration:
`∫_0^L (Π/L)(1−e^{−Π})^{−1} e^{−Π z/L} dz = (1−e^{−Π})/(1−e^{−Π}) = 1`,
not algebraically pre-baked.

## Findings

None. The script's three checks correspond one-to-one with the paper's
three operational deliverables (D2, D3, D4), and D1 (the form of σ) is
the shared starting expression.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally parallel to the SymPy script
(same sigma, same `v1 → piM*thetaSigma/lM` substitution, same three
checks in the same order). However, "independent derivation" of a
substitution-and-integration identity on a single-line exponential
profile has very little room to differ across engines — there is no
meaningful alternative algebraic path to verify `∫_0^L (Π/L)
e^{−Πz/L}/(1−e^{−Π}) dz = 1` or to check `(V_1 − ΠΘ/L) σ = 0`. The
Mathematica engine uses its native `Integrate`, `D`, and
`FullSimplify[Together[Expand[...]]]` pipeline with $Assumptions on
positivity, which is a different simplification stack than SymPy's
`simplify`. Both reach 0 independently. I do not call this a
`mathematica_transliteration` finding because the algebra is genuinely
two-line and there is no first-principles re-derivation alternative the
.wl could have taken that the .py didn't.

## Engine cross-check

Both engines report the same three results: normalization = 1
(sympy output L10; mathematica output L15), zero-flux current = 0
(L11; L17-19), ODE residual = 0 (L12; L20-22). Final identification line
also agrees: `Pi_m = V1*L/Theta` (sympy) vs `Pi_m = V1*L/Theta_sigma`
(mathematica), modulo the conventional symbol name for Θ_σ. Engines
agree.

## Verdict justification

The unit cleanly proves what the paper card and notes claim: a
single-line exponential `σ_Π(z) = Π e^{−Πz/L} / (L(1−e^{−Π}))` is the
zero-flux stationary solution of the GNLS+localized-Maxwell mouth
boundary-layer ODE `Θ_σ σ' + V_1 σ = 0` with `Π = V_1 L / Θ_σ`, and it
normalizes to 1 on `[0,L]`. Both engines verify the three operational
identities, both exit 0, the outputs are newer than the scripts
(sympy: 11 May 12:45 vs script 11:56; mathematica: 11 May 13:10 vs
script 11:56), and the substituted-bias relation `V_1 = ΠΘ/L` is
load-bearing (it is not tautological — varying the constant would fail
the assertion). Attacks attempted: (1) is the substitution
`V1 → Pi*Theta/L` algebraically guaranteed? No — it is the exact value
that makes `Θ σ' + V_1 σ` vanish for the given σ, which is the paper's
claim. (2) is the normalization integral trivial? No — it requires a
real exponential integration. (3) Does the script verify a different
identity than the paper requires (e.g., a sign-flipped Onsager current,
a different normalization interval, wrong potential slope)? No — every
sign, interval, and constant matches the notes. (4) Are there
deliverables in the notes that the script misses? No — D1-D4 are all
covered. Verdict: clean.

One cosmetic, non-finding-grade note: the Mathematica banner at
`...stage129...mathematica_audit.wl:26` reads
`banner["STAGE 112 — EXPLICIT GNLS + LOCALIZED-MAXWELL MOUTH BOUNDARY LAYER"]`
(see also output L11). The final pass line at L53 correctly says
"Stage 129 Mathematica audit passed." and the file's own name is for
stage 129. This is a stale string-literal label, not a math defect; it
does not affect any assertion or paper claim. Mentioned for completeness;
not raised as a finding because it does not fit any of the ten
categories and the audit transcript already self-identifies as stage 129
elsewhere.

## Self-test notes

Checked: (1) the substitution `V_1 = Π Θ/L` is the load-bearing identification — varying it would make the residual `(V_1 − ΠΘ/L)σ ≠ 0`, so the asserts can fail; not tautological. (2) Normalization integral is a real exponential integration, not algebraic by construction. (3) Symbol positivity (`Pi, L, Theta, V1, M, sigma_star` all positive in both engines) is justified by the physical setup (Π > 0, L > 0, Θ_σ > 0, M_σ > 0 per notes L40-83); no over-strong assumption hides a branch. (4) Outputs are fresher than scripts (mtimes checked). No traps fired.
