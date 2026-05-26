---
unit_id: 027
batch: II.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-25T00:00:00Z
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage027_nonconstant_axial_family.md
  paper_appendix: present
---

# Audit unit 027 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_027.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage027_nonconstant_axial_family.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 44)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage027_nonconstant_axial_family_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage027_nonconstant_axial_family_mathematica_audit.txt`

## What the paper claims

The stage card's `\stagefield{Output}` states:

> "Stage 027 outputs `kappa(theta)`, the blind-angle no-go (eq:app-stage027-blind-angle), the max-coupling angle (eq:app-stage027-theta-max), and the profile-dressed branch equation (eq:app-stage027-Kgate)."

Expanded across the card body and notes the distinct deliverables are:
- **D1** Exact overlap law `kappa(theta) = kappa_0 cos(theta) + kappa_1 sin(theta)` with `kappa_0 = 2*sqrt(2)/pi` and `kappa_1 = -4/(3*pi)` (eqs eq:app-stage027-kappa-theta, eq:app-stage027-kappa0-kappa1).
- **D2** Maximum overlap magnitude `|kappa|_max = sqrt(kappa_0^2 + kappa_1^2) = 2*sqrt(22)/(3*pi)` (eq:app-stage027-kappa-max).
- **D3** Blind angle `tan(theta_blind) = 3*sqrt(2)/2` with `C = G_W = R = P = N_0 = 0` and a positive-quadrupole no-go on that branch (eq:app-stage027-blind-angle plus subsequent block).
- **D4** Max-coupling angle `tan(theta_max) = kappa_1/kappa_0 = -sqrt(2)/3`, with `sin^2(theta_max) = 2/11` (eq:app-stage027-theta-max).
- **D5** Profile-dressed wall stiffness `K_geo(theta) = K_eta + 6*T_Omega + pi^2*T_w*sin^2(theta)/L^2` (eq:app-stage027-Kgeo).
- **D6** Profile-dressed normalization gate `K_geo(theta) = K_req(theta)` with `K_req` taking the explicit form in eq:app-stage027-Kreq, `Delta(theta) = Omega_U^2 Omega_W^2 - lambda_R^2 kappa(theta)^2`, and `Q(theta)` as written (eq:app-stage027-Kgate).
- **D7** Profile normalization `int_0^L chi_theta(s)^2 ds = 1` (paper Checks item 1; follows from u0/u1 orthonormality).

The notes file echoes the same deliverables and provides equivalent amplitude form `kappa(theta) = rho cos(theta - theta_max)`; no additional numerical claims that diverge from the stage card. The appendix row 44 is one-line and consistent: "Profile overlap kappa(theta), blind-angle no-go, max-coupling angle, and dressed stiffness equation."

## What the script claims to verify

The SymPy script verifies, in order: (i) orthonormality of `u_0 = 1/sqrt(L)`, `u_1 = sqrt(2/L) cos(pi s/L)` and `f_0 = sqrt(2/L) sin(pi s/(2L))` on `[0, L]` (Section I); (ii) the closed-form overlap constants `kappa_0 = 2*sqrt(2)/pi`, `kappa_1 = -4/(3*pi)`, the linear overlap law `kappa(theta) = kappa_0 cos + kappa_1 sin`, and the maximum magnitude `rho = 2*sqrt(22)/(3*pi)` (Section II.1); (iii) `kappa(blind) = 0`, `kappa(max) = rho`, and `sin^2(theta_max) = 2/11` evaluated at the exact trigonometric substitutions consistent with `tan(theta_blind) = 3*sqrt(2)/2` and `tan(theta_max) = -sqrt(2)/3` (Section II.2); (iv) the wall-stiffness expectation `K_geo(theta) = K_eta + 6*T_Omega + T_w*pi^2*sin^2(theta)/L^2` *derived* via `Integrate[chi*(-T_w chi_ss + (K_eta + 6 T_Omega) chi), {s,0,L}]` and its reductions at theta=0 and theta=theta_max (Section III); (v) the profile-dressed Stage-8/9 reduced quantities `Delta, Q, P, B0` and their `theta=0` recovery (Section IV); and (vi) the blind-angle no-go via `P_blind = 0`, `N0_blind = 0`, `lhs_blind = 0` against the strictly-positive GR target `54*G*c_s^5/(5*a^5*c^5)` (Section V).

The Mathematica script verifies the same six blocks. Two minor differences from the SymPy script: (a) Section V omits the explicit `lhs(blind) == 0` assertion, only checking `P(blind) = 0` and `N0(blind) = 0` (which together imply `lhs(blind) = 0` since the lhs numerator is `mhat^2 * N0` and the denominator is `K_geo - B0 - Q/Delta`, generically nonzero in the free parameters); (b) the `K_geo(theta)` expression is built from a real `Integrate[chi*gEta, {s, 0, l}]` and only collapses to `sin^2(theta)` form via `TrigExpand`, so the printed pre-expansion form is `kEta + 6 tOmega + Pi^2 tW/(2 l^2) - Pi^2 tW Cos[2 theta]/(2 l^2)` — confirmation that the integral is genuinely computed.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1 `kappa(theta)`, `kappa_0`, `kappa_1` | sympy 114-119 / mathematica 70-72 | match |
| D2 `|kappa|_max = 2*sqrt(22)/(3*pi)` | sympy 120 / mathematica 73 | match |
| D3 blind angle, no-go | sympy 141 (kappa(blind)=0) + 271-273 (P,N0,lhs); mathematica 85,172-173 | match (tan(theta_blind) value is implicit in the cos/sin substitution rule; the substitution is algebraically equivalent and kappa(blind)=0 is the operative no-go content) |
| D4 max-coupling angle, sin^2 = 2/11 | sympy 142-143 / mathematica 86-87 | match (tan(theta_max) implicit in cos/sin rule; sin^2 = 2/11 explicit) |
| D5 `K_geo(theta) = K_eta + 6 T_Omega + pi^2 T_w sin^2/L^2` | sympy 156-161 (derived from integral) / mathematica 94-99 (derived from integral) | match |
| D6 normalization gate `K_geo = K_req` with K_req form | sympy 247-256 / mathematica 160-166 (constructs K_req from previously-anchored B0, Q/Delta, P^2, target; prints) | partial (K_req is *built* from anchored components and printed; never asserted against a literal paper-side closed form. The components Delta, Q, P, B0 are all anchored at A5/B5, and `target = 54 G c_s^5/(5 a^5 c^5)` is the literal GR-side normalization constant. The construction is therefore traceable, but no explicit `assert K_req - K_req_paper == 0` exists.) |
| D7 `int chi^2 ds = 1` | sympy 87-89 / mathematica 49-51 (orthonormality of u0, u1, u0-u1) | implicit (follows from cos^2 + sin^2 = 1 and the three orthonormality checks; not asserted as a direct `int chi^2 - 1 == 0` line) |

`paper_alignment: aligned`. The two "partial/implicit" rows (D6, D7) are both cases where the script verifies the load-bearing pieces and assembles the final object correctly but does not write a final tautology-free assertion line. Neither rises to a paper_misalignment because the script does not assert a *different* identity than the paper — it just stops one assertion short, and the missing assertion is algebraically guaranteed by the anchored building blocks (D7) or by component anchoring of the K_req expression (D6). I considered raising D6 as `insufficient_verification` but the K_req formula in the paper has no independent existence beyond what the components define; checking `K_req - (paper-side K_req)` would compare the script's `B0 + Q/Delta + mhat^2 P^2/(target Delta^2)` against the paper's `lambda_B^2 kappa^2/varpi^2 + Q/Delta + mhat_rad^2 kappa^2 (Omega_U^2 lambda_W + lambda_R lambda_U)^2 / (N_Q^target Delta^2)`, and that reduces to checking `B0 = lambda_B^2 kappa^2/varpi^2` (already done) and `P^2 = kappa^2 (Omega_U^2 lambda_W + lambda_R lambda_U)^2` (already done via P-expected). So the substantive content is in fact anchored.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 87-90 | `int u_i u_j - delta_ij == 0`, `int f0^2 - 1 == 0` | D7 (component) | yes |
| A2 | sympy | 114-120 | `kappa0 - 2*sqrt(2)/pi`, `kappa1 + 4/(3*pi)`, `kappa - (k0 c + k1 s)`, `rho - 2*sqrt(22)/(3*pi)` | D1, D2 | yes |
| A3 | sympy | 141-143 | `kappa(blind) == 0`, `kappa(max) - rho == 0`, `sin^2(theta_max) - 2/11 == 0` | D3 (kappa side), D4 | yes |
| A4 | sympy | 161, 170-173 | `K_geo - (K_eta + 6 T_Omega + T_w pi^2 sin^2/L^2) == 0` (post-integration), `K_geo(theta_max) - ...` | D5 | yes |
| A5 | sympy | 221-224 | `Delta/Q/P/B0 - expected == 0` | D6 (component) | yes |
| A6 | sympy | 230-238 | `kappa(0) - kappa0`, `K_geo(0) - (K_eta+6T_Omega)`, `Delta(0) - ...` | D5 (limit) | yes |
| A7 | sympy | 271-273 | `P(blind) == 0`, `N0(blind) == 0`, `lhs(blind) == 0` | D3 (no-go) | yes |
| B1 | mathematica | 49-52 | `int u_i u_j` and `int f0^2 - 1` checks | D7 (component) | yes |
| B2 | mathematica | 70-73 | overlap constants and law | D1, D2 | yes |
| B3 | mathematica | 85-87 | blind/max trigonometric values | D3 (kappa side), D4 | yes |
| B4 | mathematica | 94-99, 102-104 | `kGeo - kGeoExpected == 0` (now genuinely derived via `Integrate[chi*gEta]`), `K_geo(theta_max) - ...` | D5 | yes |
| B5 | mathematica | 142-145 | `delta/q/p/b0 - expected == 0` | D6 (component) | yes |
| B6 | mathematica | 150-155 | recovery at theta=0 | D5 (limit) | yes |
| B7 | mathematica | 172-173 | `P(blind) == 0`, `N0(blind) == 0` | D3 (no-go) | yes |

All assertions are anchored to a paper deliverable. No tautological assertions remain in the current Mathematica script — the prior round's F1 on `wallStiffness[]` (a hard-coded `kGeo = kGeoExpected` self-comparison) has been fixed: lines 94-95 now construct `gEta = -tW*D[chi, {s, 2}] + (kEta + 6*tOmega)*chi` and compute `kGeo = FullSimplify[Integrate[chi*gEta, {s, 0, l}], Assumptions -> $Assumptions]`. The saved output line 62 confirms genuine integration: `K_geo(theta) = kEta + 6*tOmega + (Pi^2*tW)/(2*l^2) - (Pi^2*tW*Cos[2*theta])/(2*l^2)`, which is the `Integrate` result *before* the trig identity `Sin^2(theta) = (1 - Cos[2 theta])/2` collapses it — a signature only a real integral leaves behind.

## Findings

None.

## Independent-derivation check (Mathematica)

The Mathematica script remains structurally parallel to the SymPy script: identical section banners, 1:1 function name correspondence (`finiteThroatBasis` ↔ `finite_throat_basis`, `overlapLaw` ↔ `overlap_law`, etc.), identical intermediate variables, identical substitution dictionaries for blind/max angles. This pattern would normally invite a `mathematica_transliteration` finding, but:

1. The Section III wall-stiffness integral — the one place where genuinely independent computation can disambiguate the engines — is now an actual `Integrate[chi*gEta, {s, 0, l}]` in Mathematica (lines 94-95), not a hard-coded re-statement. The saved-output evidence (`(1 - Cos[2 theta])/2` pre-expansion form on line 62) confirms the integral is really being computed by Mathematica's integration engine, independently of SymPy's. The two engines arrive at the same `K_geo` from independent quadrature.
2. The remaining "parallel" content (orthonormality of `u_0, u_1, f_0`, the linear overlap law, the trivial `Delta/Q/P/B0` component identities) admits essentially one symbolic form. Forcing different surface algebra would be stylistic refactor outside the auditor's scope.

I therefore decline to file a `mathematica_transliteration` finding. The substantive independent-derivation requirement is met by Section III post-F1; the rest is unavoidable structural mirroring of single-form algebra.

## Engine cross-check

The two engines run to completion (exit 0) and their printed closed forms agree:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `kappa0` | `2*sqrt(2)/pi` | `(2*Sqrt[2])/Pi` |
| `kappa1` | `-4/(3*pi)` | `-4/(3*Pi)` |
| `kappa(theta)` | `2*(-2*sin + 3*sqrt(2)*cos)/(3*pi)` | `(6*Sqrt[2]*Cos - 4*Sin)/(3*Pi)` (algebraically identical) |
| `rho` | `2*sqrt(22)/(3*pi)` | `(2*Sqrt[22])/(3*Pi)` |
| `sin^2(theta_max)` | `2/11` | `2/11` |
| `kappa(blind)` | `0` | `0` |
| `kappa(max)` | `2*sqrt(22)/(3*pi)` | `(2*Sqrt[22])/(3*Pi)` |
| `K_geo(theta)` | `K_eta + 6 T_Omega + pi^2 T_w sin^2/L^2` | `kEta + 6 tOmega + Pi^2 tW/(2 l^2) - Pi^2 tW Cos[2 theta]/(2 l^2)` (= same after `Sin^2 = (1 - Cos[2 theta])/2`) |
| `K_geo(theta_max)` | `K_eta + 6 T_Omega + 2 pi^2 T_w/(11 L^2)` | `kEta + 6 tOmega + 2 Pi^2 tW/(11 l^2)` |
| `P(blind)` | `0` | `0` |
| `N0(blind)` | `0` | `0` |
| Section IV `Delta/Q/P/B0/Z0/N0/D0` | (printed; expected residuals = 0) | (printed; expected residuals = 0) |

The Section IV polynomial-in-`kappa` expressions differ on the surface — SymPy prints in the `(2 sin - 3 sqrt(2) cos)` representation, Mathematica in the `(6 Sqrt[2] Cos - 4 Sin)` representation — but `(6 sqrt(2) cos - 4 sin) = -2 (2 sin - 3 sqrt(2) cos)`, so the squared forms are identical and the linear forms differ only by an overall sign that cancels in the `expectZero` residuals. No engine disagreement.

## Verdict justification

The two scripts hold up against attack and they jointly cover every distinct deliverable on the paper card. The orthonormality integrals are real integrals on `[0, L]`; the overlap-law and rho identities are non-tautological (LHS by integration, RHS literal closed form); the blind/max-angle algebraic substitutions correctly satisfy `cos^2 + sin^2 = 1` and reproduce the paper-stated tangent values; the wall-stiffness integral is now genuinely derived in *both* engines (the prior round's tautological Mathematica `wallStiffness[]` block has been replaced by a real `Integrate[chi*gEta]`); the Section IV reduced-branch quantities collapse to their expected closed forms; and the blind-angle no-go is driven by `kappa(blind) = 0` propagating to `P, N_0, lhs` against the strictly-positive GR target. The two minor "partial" rows in the paper-cross-check (D6's missing top-level `K_req` literal assertion, D7's missing direct `int chi^2 = 1` assertion) are anchored at the component level and reduce to identities already in the assertion list; they do not constitute paper_misalignment or insufficient_verification once the underlying anchors are recognized.

Outputs are fresh (script mtimes Apr 1 and May 21; output mtimes both May 21 17:04, after the most recent script touch). The cosmetic docstring drift in the SymPy module (the file's internal docstring still names "Stage 10" rather than "Stage 027") and the Mathematica final-line "Stage 10 Mathematica audit passed." are pre-renumbering carryovers that do not affect the verified math; they are common across the renumbered Part II ledger and the user has previously declined to file findings on docstring-only stage-number drift.

Attacks tried that the scripts survived:

- Sign of `chi_ss`: `u1_ss = -(pi/L)^2 u1`, so `-T_w*<chi, chi_ss> = +T_w sin^2(theta) pi^2/L^2`. Matches the asserted positive sign.
- Orthogonality `int_0^L cos(pi s/L) ds = 0` (so `int u0 u1 = 0`); the saved output confirms.
- Blind-substitution consistency: `cos^2 + sin^2 = 2/11 + 9/11 = 1`. Max-substitution: `9/11 + 2/11 = 1`. Both unit-normed.
- `tan(theta_blind) = sin/cos = (3/sqrt(11)) / (sqrt(2)/sqrt(11)) = 3/sqrt(2) = 3*sqrt(2)/2`. Matches paper eq:app-stage027-blind-angle.
- `tan(theta_max) = sin/cos = (-sqrt(2)/sqrt(11)) / (3/sqrt(11)) = -sqrt(2)/3`. Matches paper eq:app-stage027-theta-max.
- `kappa(blind)` arithmetic: `(2 sqrt(2)/pi)(sqrt(2)/sqrt(11)) + (-4/(3 pi))(3/sqrt(11)) = 4/(pi sqrt(11)) - 4/(pi sqrt(11)) = 0`. Holds.
- `kappa(max) - rho`: `(2 sqrt(2)/pi)(3/sqrt(11)) + (-4/(3 pi))(-sqrt(2)/sqrt(11)) = 6 sqrt(2)/(pi sqrt(11)) + 4 sqrt(2)/(3 pi sqrt(11)) = (18 + 4) sqrt(2)/(3 pi sqrt(11)) = 22 sqrt(2)/(3 pi sqrt(11)) = 2 sqrt(22)/(3 pi) = rho`. Holds.
- `K_geo(theta_max)` arithmetic: `pi^2 T_w/L^2 * (2/11) = 2 pi^2 T_w/(11 L^2)`. Matches.
- Mathematica integral signature: line 62 of saved output shows `(1 - Cos[2 theta])/2` form, confirming `Integrate` is actually being invoked (a hard-coded RHS would print as `Sin[theta]^2` directly).
- Blind-angle propagation to `C = lambda_B kappa`, `G_W = lambda_W kappa`, `R = lambda_R kappa`: all vanish trivially when `kappa(blind) = 0`. Paper's enumeration `C=G_W=R=P=N_0=0` is therefore covered (even though only P and N_0 carry their own explicit assertions).

Verdict: `clean`. The unit's paper claims are aligned with the script claims; both engines verify the load-bearing physics; no findings raised. `paper_alignment: aligned`. No directive file is written.

## Self-test notes

I checked three traps before finalizing the verdict-of-clean:

1. **Tautological "expected" pattern**: I re-confirmed that the current Mathematica `wallStiffness[]` no longer self-compares a hard-coded `kGeo` against itself — line 94 builds `gEta` from `D[chi, {s, 2}]` and line 95 computes `kGeo` by `Integrate`, with `kGeoExpected` (line 96) a *separate* literal symbol. The saved output's pre-`TrigExpand` form `(1 - Cos[2 theta])/2` is the smoking-gun that an actual integral was evaluated. The prior-round F1 has held.
2. **Blind/max substitutions hide the angle identities**: I verified `tan(theta_blind) = 3/sqrt(2) = 3 sqrt(2)/2` and `tan(theta_max) = -sqrt(2)/3` are algebraically implied by the cos/sin rules the scripts use, so the paper's tangent-form statements are anchored even though the scripts never write `assert tan(theta_blind) - 3*sqrt(2)/2 == 0`. The substitution rules are unit-normed (`cos^2 + sin^2 = 1` for both), so they truly correspond to the paper's angles, not to arbitrary numerical fits.
3. **Missing `int chi^2 = 1` assertion**: I confirmed this is algebraically forced by the three orthonormality checks plus `cos^2 + sin^2 = 1`, so flagging it as `insufficient_verification` would be pedantic. Likewise the `K_req` printed expression is built from B0, Q/Delta, and P^2 — all of which carry their own anchored `- expected == 0` assertions — so the absence of a final `K_req - K_req_paper == 0` line does not create an unverified gap.

No findings, no directive.
