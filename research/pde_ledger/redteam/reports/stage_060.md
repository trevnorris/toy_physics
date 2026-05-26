---
unit_id: 060
batch: III.2
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26T00:00:00Z
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage060_entropic_microclosure.md
  paper_appendix: present
---

# Audit unit 060 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_060.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage060_entropic_microclosure.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 98 + `\input` at line 238)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage060_entropic_microclosure_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage060_entropic_microclosure_mathematica_audit.txt`

## What the paper claims

The paper card (`stage_060.tex:20`) states the stage proves: "A positive-density microclosure for the transport source family and the definition of the microscopic gain variable." The card body (`stage_060.tex:14-18`) reproduces the source-transport equation of Stage~056 from an Onsager form `J=-M(sigma) partial_s (delta F / delta sigma)` with a strictly positive mobility and convex local entropy, and asserts that the stationary zero-flux state is the Stage-39 exponential branch. The notes file enumerates eight specific deliverables: (1) the explicit free-energy functional `F[sigma,phi]`, (2) its Euler-Lagrange variations (chemical potential, support bulk equation, and Robin/Neumann BCs `T_X phi_s(0)=K_m phi(0), phi_s(L)=0`), (3) the Onsager current and Einstein relation `D_sigma=M_sigma Theta_sigma`, (4) recovery of the exponential family with `Pe=(Lambda_phi/Theta_sigma) Delta phi` and `Sigma_Pe(x)=Pe e^{Pe x}/(e^{Pe}-1)`, (5) the microscopic coupling `Xi_micro = Lambda_phi^2 L^2/(Theta_sigma T_X) = chi_sigma Lambda_phi^2 L^2/T_X = mu_sigma Lambda_phi^2 L^2/(D_sigma T_X)`, (6) passivity `dF/dt <= 0` under no-flux BCs. The appendix row (`stage_appendix_part03.tex:98`) summarises: "Positive source/free-energy closure reproducing the transport source family." All paper-side claims are anchored either in the .tex card or the notes.

## What the script claims to verify

The SymPy and Mathematica scripts together verify: (i) the chemical potential `mu = Theta log(sigma/sigma_*) - Lambda phi` from variation of the free-energy density (sympy:38-43, wl:43-48); (ii) the Onsager current decomposition `J = -M_sigma Theta sigma_s + M_sigma Lambda sigma phi_s` and the Einstein-relation substitution `M_sigma Theta -> D_sigma` (sympy:46-52, wl:52-62); (iii) the affine-drop ODE annihilation by `sigma=Cnorm e^{gamma s}` (sympy:62-63, wl:75-84) followed by an explicit normalisation `Cnorm = a/(e^{aL}-1)` and rescaling check `Sigma_from_rescale = Sigma_Pe` (sympy:65-78, wl:86-96); (iv) a derived-rate identification `gamma_derived = Lambda Delta_phi/(Theta L)` via `solve` on the ansatz, equivalent to `Pe = (Lambda/Theta) Delta phi` after `Pe=gamma L` (sympy:80-89, wl:97-107); (v) a constant-source K_X=0 BVP for `phi` with the Robin/Neumann BCs and a K_m -> infinity limit yielding `Delta phi = Lambda L^2 sigma_0/(2 T_X)` as a scale-confirming sanity check (sympy:91-116, wl:109-127); (vi) the gain `Xi_micro = Lambda^2 L^2/(Theta T_X)` and its susceptibility and phenomenological substitutions, anchored to the BVP-validated scale `Lambda L^2/T_X` (sympy:118-127, wl:129-142); (vii) dissipation-density positivity `M sigma mu_s^2 >= 0` under M, sigma > 0 (sympy:135-145, wl:144-156); (viii) the support bulk EL `-Lambda sigma + K_X phi - T_X phi'' = 0` from the full free-energy density (sympy:148-153, wl:158-165).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Free-energy → mu = Theta log(sigma/sigma_*) - Lambda phi | A1/B1 (sympy:42-43, wl:48) | match |
| Onsager current J = -M sigma partial_s mu | A2/B2 (sympy:46-49, wl:52-57) | match |
| Einstein relation D_sigma = M_sigma Theta | A3/B3 (sympy:51-52, wl:58-62) | match |
| Affine-drop exponential branch annihilates J=0 | A4/B4 (sympy:62-63, wl:75-84) | match |
| Sigma_Pe(x) = Pe e^{Pe x}/(e^{Pe}-1) | A5/B5 (sympy:76, wl:94) — derived via rescaling | match |
| Normalisation of Sigma_Pe to 1 on [0,1] | A6/B6 (sympy:77-78, wl:95-96) | match |
| Pe = (Lambda/Theta) Delta phi (i.e. gamma = Lambda dphi/(Theta L)) | A7/B7 (sympy:88-89, wl:106-107) via solve | match (Pe = gamma*L identification implicit, not explicitly named) |
| Support bulk EL -Lambda sigma + K_X phi - T_X phi'' = 0 | A8/B8 (sympy:152-153, wl:165) | match |
| Robin/Neumann BCs T_X phi'(0)=K_m phi(0), phi'(L)=0 | used in BVP (sympy:103-104, wl:116-117); not derived from variation | partial (asserted-not-derived; used downstream) |
| phi_from_Phi scale Lambda L^2/T_X validated via BVP | A9/B9 (sympy:113-116, wl:125-127), K_m->infty limit yields Lambda L^2 sigma_0/(2 T_X) | match |
| Xi_micro = Lambda^2 L^2/(Theta T_X) | A10/B10 (sympy:125, wl:140) | match |
| Xi_micro susceptibility / phenomenological forms | A11-A12/B11-B12 (sympy:126-127, wl:141-142) | match |
| Passivity dF/dt <= 0 (no-flux) | dissipation density M sigma mu_s^2 >= 0 verified (sympy:141-143, wl:151-156); the integration-by-parts → integral form is stated in print only | partial (density-level positivity confirmed; full integral identity narrated, not derived) |

`paper_alignment: aligned` — every paper-side deliverable maps to a substantive script-side check. Two partial entries (BCs not derived from variation; full Lyapunov integral identity stated in print rather than derived from continuity + integration by parts) are below the finding bar because (a) the BCs are exercised in the downstream BVP and the bulk EL is fully derived from the variational density, and (b) the dissipation-density positivity check captures the essential physics statement; the integration-by-parts step is a pure calculus identity given continuity, and continuity itself is a Stage-056 result not a Stage-060 claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42-43 | `mu_expr - (Theta log(sigma/sigma_*) - Lambda phi) == 0` | mu from free-energy variation | yes |
| A2 | sympy | 48-49 | Onsager current decomposition | J = -M sigma partial_s mu | yes (verifies chain-rule expansion) |
| A3 | sympy | 51-52 | D_sigma substitution | Einstein relation | yes (verifies M Theta = D in the current expression) |
| A4 | sympy | 62-63 | exponential ansatz annihilates affine-drop ODE | recovery of exp branch | yes (substantive) |
| A5 | sympy | 68-69 | `integrate(Csol*exp(a s), 0..L) - 1 == 0` | normalisation of trial family | yes |
| A6 | sympy | 76 | `Sigma_from_rescale - Sigma_x == 0` | Sigma_Pe(x) form via rescaling | yes (now derived, not hand-typed) |
| A7 | sympy | 77-78 | `integrate(Sigma_x, 0..1) - 1 == 0` | normalised Sigma_Pe family | yes |
| A8 | sympy | 88-89 | `gamma_derived - Lambda dphi/(Theta L) == 0` via solve | Pe identification | yes (gamma is solved for, not asserted) |
| A9 | sympy | 115-116 | BVP K_m->infty limit equals `Lambda L^2 sigma_0/(2 T_X)` | scale-confirming check for phi_from_Phi | yes (true BVP derivation) |
| A10 | sympy | 125 | `Xi_micro - Lambda^2 L^2/(Theta T_X) == 0` | Xi_micro headline | partial (definitional; load-bearing through BVP scale validation in A9) |
| A11 | sympy | 126 | `Xi_micro.subs(Theta,1/chi) - chi Lambda^2 L^2/T_X == 0` | susceptibility form | partial (substitution check) |
| A12 | sympy | 127 | `Xi_micro.subs(Theta,D/M) - M Lambda^2 L^2/(D T_X) == 0` | phenomenological form | partial (substitution check) |
| A13 | sympy | 141-143 | `sp.ask(Q.nonnegative(dissipation_density), ...)` | passivity (density level) | yes (genuine positivity check) |
| A14 | sympy | 152-153 | `EL_phi - (-Lambda sigma + K_X phi - T_X phi'') == 0` | support bulk EL | yes |
| B1-B14 | math | corresponding | Same identities | corresponding paper claims | mirror of A1-A14 — see independent-derivation check below |

Every assertion traces to a specific paper-side deliverable, addressing the v2 paper-alignment gate. The three partial rows (A10-A12) are not orphaned — they are anchored to A9's BVP scale validation; the susceptibility/phenomenological forms are simple symbol renames of an already-grounded quantity.

## Findings

(None — verdict is clean.)

## Independent-derivation check (Mathematica)

The Mathematica script is now meaningfully independent from the SymPy script in two of three previously flagged areas:

1. **Rescaling derivation**: both engines now compute `cNormSol = a/(Exp[a*ell]-1)` then rescale to verify `Sigma_Pe = Pe e^{Pe x}/(e^{Pe}-1)`. The wl uses `deltaDrop -> pe*theta/lambdaPhi` (wl:91) while sympy uses the same substitution `subs(dphi, Pe*Theta/Lam)` (sympy:73). This is structurally parallel rather than independent, but acceptable since both engines arrive at the rescaling via the explicit normalisation + change of variables (no hand-typed Sigma_Pe survives).

2. **Pe identification**: both engines call `Solve` / `sp.solve` on the ansatz against `J=0` (wl:101-104, sympy:84-86). Same algorithm; symbolic computer algebra paths are inherently similar here.

3. **Xi_micro construction**: the wl script (wl:132) builds `xiMicro = lambdaPhi^2*ell^2/(theta*tX)` from a dimensional combination and verifies the three forms via consistency substitutions (wl:133-138). The SymPy script (sympy:122-127) uses the older `phi_from_Phi` route and verifies the same identities. Different intermediate constructions, same endpoint — this is the F3-required independence delta from v1, and it landed.

4. **BVP**: both engines do `DSolve` / `sp.dsolve` on the constant-sigma K_X=0 EL with matching BCs (wl:111-127, sympy:93-116), arriving at the same `K_m -> infty` limit `Lambda L^2 sigma_0/(2 T_X)`. Same physical claim, parallel codepaths.

Calling this a transliteration would be unfair given the v1 corrections; calling it fully independent would be generous. The two engines do agree on every line of the .txt outputs, which is the relevant cross-check.

## Engine cross-check

Both engines pass all assertions with zero residual. The numeric/symbolic outputs match:
- `mu_sigma^(chem)` form: matches.
- `J expanded`: matches (sympy line 8, wl line 9).
- `gamma_derived = Delta_phi*Lambda_phi/(L*Theta)` (sympy) vs `gammaDerived = (deltaDrop*lambdaPhi)/(ell*theta)` (wl) — symbol-renamed identical.
- `Delta_derived = L**2*Lambda_phi*sigma_0/(2*T_X)` (sympy) vs `(ell^2*lambdaPhi*sigma0)/(2*tX)` (wl) — identical.
- `Xi_micro = L**2*Lambda_phi**2/(T_X*Theta)` (sympy) vs `(ell^2*lambdaPhi^2)/(theta*tX)` (wl) — identical.
- Support bulk EL form: identical structure.

No engine disagreement.

## Outputs freshness

- sympy script mtime: 2026-05-22 17:55; sympy output mtime: 2026-05-22 17:59 → output newer.
- wl script mtime: 2026-05-22 17:58; wl output mtime: 2026-05-22 17:59 → output newer.

Both outputs are fresh relative to their producing scripts.

## Verdict justification

Verdict is `clean`. Read the paper card and notes first, then the scripts. The script covers every numbered deliverable in the notes and the headline `Output` line of the paper card. The microscopic gain `Xi_micro = Lambda^2 L^2/(Theta T_X)` is now anchored to a BVP-derived scale check (sympy:115-116, wl:126-127), not a hand-typed `phi_from_Phi`. The exponential family is derived by explicit normalisation and rescaling rather than hand-typed (A5/A6). The rate identification goes through a genuine `solve` rather than substitution into its own definition (A8). The dissipation positivity is checked at the density level via `sp.ask` and `Reduce[ForAll[...]]` rather than via a tautological cancellation. The support bulk EL is recovered from the full free-energy density. 

Attacks tried and failed: (1) substituted `Csol = a/(e^{aL}-1)` into the normalisation integral and confirmed it integrates to 1 for generic non-zero `a`, ruling out the v1 piecewise-degenerate-branch bug; (2) traced `Sigma_from_rescale = L * Csol * exp(a s)` under `s -> xL, deltaDrop -> Pe*Theta/Lambda` and confirmed it equals `Pe exp(Pe x)/(exp(Pe)-1)`; (3) verified the BVP residual `phi(s) = -Lambda sigma_0 s^2/(2 T_X) + C1 s + C2` with BCs `T_X C1 = K_m C2`, `C1 = Lambda sigma_0 L/T_X` yields `Delta = Lambda sigma_0 L^2/(2 T_X)`, K_m-independent in this limit — matching the script output; (4) checked the `gamma_solved` branch selection (sympy:86, wl:104) — `solve` returns `gamma = 0` (degenerate) and `gamma = Lambda dphi/(Theta L)`, and the script's `[g for g in gamma_solved if g != 0][0]` correctly picks the nonzero branch; (5) verified the dissipation density `M sigma mu_s^2` is provably nonnegative under the declared assumptions. The only borderline items — BCs asserted rather than re-derived from boundary-term variation, and Lyapunov integral identity narrated rather than derived from continuity — are below the finding bar because the BCs are exercised in the downstream BVP and the integral identity is a calculus consequence of continuity (a Stage-056 result, not a Stage-060 claim).

Paper alignment is exact on every load-bearing identity. No `paper_misalignment` finding applies. No downstream propagation concern: `Xi_micro = Lambda^2 L^2/(Theta T_X)` is consistent across paper, notes, and both scripts; Stage 061 (which consumes this) inherits an unambiguous definition.

## Self-test notes

(1) Variable independence: the `gamma_solved` step solves a single-variable polynomial-in-the-exponential equation; gamma genuinely appears in the equation. The BVP `dsolve` step has both BCs depending on `K_m, T_X, L, Lambda, sigma_0` — none of the variables collapse trivially. (2) Symmetry/parity: no symmetric-domain integrals are involved; the normalisation integral over `[0,L]` (resp. `[0,1]`) is one-sided. (3) Trivial-case pre-check: setting `dphi -> 0` in `gamma_derived = Lambda dphi/(Theta L)` gives gamma = 0, matching the degenerate branch that the script's selector excludes — confirms the nonzero-branch selection logic is sound. Substituting `sigma_0 = 0` in the BVP collapses to `phi'' = 0`, giving constant `phi` and `Delta = 0`, matching `Lambda * 0 * L^2/(2 T_X) = 0`. (4) Path specifications: no missing-script findings. (5) Paper round-trip: no fixes prescribed, so no risk of introducing new misalignment.
