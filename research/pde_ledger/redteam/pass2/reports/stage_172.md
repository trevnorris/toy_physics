---
unit_id: 172
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage172_physical_slope_collapse.md]
  paper_appendix: present
---

# Audit unit 172 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_172.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage172_physical_slope_collapse.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows/blocks at lines 75, 461-487; anchor descriptions 1463)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.txt`

## What the paper claims

Stage 172 ("Physical slope collapse") rewrites the Stage 171 microscopic outlet-obstruction
pair in terms of the physical grouped-response slopes. The `\stagefield{Output}` states
verbatim: "Identifies \(\mathfrak K_1=-D_0u_2^{(1)}\) and \(\mathfrak G_1=D_0P_1\), reducing
the linear grouped outlet problem to physical slopes." The notes flesh this out into six
boxed deliverables: (1) the exact first-order variations \(\delta u_2^{(A)}=-(\delta D_2+u_2\delta D_0)/D_0\)
and \(\delta P_0^{(A)}=(\delta N_0-P_0\delta D_0)/D_0\); (2) the obstruction-pair rewrite
\(\mathcal K_A=-D_0\delta u_2^{(A)}+(1/9-u_2)\delta D_0\), \(\mathcal G_A=D_0\delta P_0^{(A)}\);
(3) the canonical-even collapse \(\mathcal K_A=-D_0\delta u_2^{(A)}\) at \(u_2=1/9\); (4) the
hidden-even consistency relation \(\delta u_4^{(A)}=(8/9)\delta u_2^{(A)}\iff \delta D_4=(2/3)\delta D_2+(1/27)\delta D_0\),
together with the slope \(\delta u_4=-(5\delta D_0+18\delta D_2+81\delta D_4)/(81D_0)\) on
\((u_2,u_4)=(1/9,4/81)\); (5) the physical direct-outlet maps
\(\delta\kappa_W^{(A)}=-3(1-\sigma_*)/\sigma_*\,\delta u_2^{(A)}\) and
\(\delta\gamma_W^{(A)}=-(1-\sigma_*)/(9\sigma_*)\,\delta P_0^{(A)}/P_0\); and (6) the
even-preserving residual \(\Delta_Q^{(A)}=\delta P_0^{(A)}/P_0\). The appendix block (lines
461-487) and row (line 75) restate items 1, 2 and 5. The card is an exact-closure original-stage
audit card; no novel numeric constants are introduced beyond the canonical \(u_2=1/9,u_4=4/81\).

## What the script claims to verify

Both scripts (identical six-check structure, declared in the SymPy docstring lines 7-13) compute
the first-order variations \(\delta u_2^{(A)}\), \(\delta P_0^{(A)}\) (SymPy by `sp.series`,
Mathematica by implicitly differentiating the defining relations), then assert, via `expect_zero`
/ `expectZero` residuals: (1) `K_A + D0*delta_u2 - (1/9-u2)*dD0 == 0`; (2) `G_A - D0*delta_P0 == 0`;
(3) the same minus the `(1/9-u2)*dD0` term after `u2 -> 1/9`; (4) `delta_u4 - (8/9)delta_u2 == 0`
after substituting the even-consistency relation `dD4 = (2/3)dD2 + (1/27)dD0`; (5) the two
direct-outlet residuals `delta_kappa + 3(1-sigma)/sigma*delta_u2 == 0` (at u2=1/9) and
`delta_gamma + (1-sigma)/(9 sigma P0)*delta_P0 == 0`; and (6) `Delta_Q - delta_P0/P0 == 0`.
`K_A`, `G_A` carry the literal Stage-171 definitions; `delta_kappa`, `delta_gamma`, `Delta_Q`
carry the Stage-170 direct-outlet/normalization maps. Every assertion combines an
independently-derived slope with one of the obstruction/outlet definitions, so none is tautological.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| \(\delta u_2^{(A)}=-(\delta D_2+u_2\delta D_0)/D_0\) | `delta_u2` printed = `(-dD0*u2 - dD2)/D0` (sympy out L5; wl out L5) | match |
| \(\delta P_0^{(A)}=(\delta N_0-P_0\delta D_0)/D_0\) | `delta_P0` printed = `(D0*dN0 - N0*dD0)/D0**2` (out L6) | match |
| \(\mathcal K_A=-D_0\delta u_2+(1/9-u_2)\delta D_0\) | A1 (py L62-65 / wl L55-58) | match |
| \(\mathcal G_A=D_0\delta P_0^{(A)}\) | A2 (py L66-69 / wl L59) | match |
| Canonical-even collapse at \(u_2=1/9\) | A3 (py L72-75 / wl L62-65) | match |
| \(\delta u_4=(8/9)\delta u_2\iff\delta D_4=(2/3)\delta D_2+(1/27)\delta D_0\) | A4 (py L101-104 / wl L90-93) | match |
| \(\delta\kappa_W^{(A)}=-3(1-\sigma_*)/\sigma_*\,\delta u_2^{(A)}\) | A5 (py L115-118 / wl L99-102) | match |
| \(\delta\gamma_W^{(A)}=-(1-\sigma_*)/(9\sigma_*)\,\delta P_0^{(A)}/P_0\) | A6 (py L120-123 / wl L103-106) | match |
| \(\Delta_Q^{(A)}=\delta P_0^{(A)}/P_0\) | A7 (py L126-129 / wl L108-109) | match |

`paper_alignment: aligned` — every boxed deliverable in the notes and every Output/appendix claim
has a faithful, non-tautological script-side residual check. No extra/orphaned assertions.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy/wl | py62-65 / wl55-58 | `simplify(K_A + D0*du2 - (1/9-u2)dD0) == 0` | obstruction-pair K (notes §2) | yes |
| A2 | sympy/wl | py66-69 / wl59 | `simplify(G_A - D0*dP0) == 0` | obstruction-pair G (notes §2) | yes |
| A3 | sympy/wl | py72-75 / wl62-65 | `(K_A + D0*du2)|_{u2=1/9} == 0` | canonical-even collapse (notes §2) | yes |
| A4 | sympy/wl | py101-104 / wl90-93 | `(du4 - 8/9 du2)|_{dD4=(2/3)dD2+(1/27)dD0} == 0` | hidden-even relation (notes §4) | yes |
| A5 | sympy/wl | py115-118 / wl99-102 | `(dkappa + 3(1-σ)/σ du2)|_{u2=1/9} == 0` | direct κ map (notes §5) | yes |
| A6 | sympy/wl | py120-123 / wl103-106 | `dgamma + (1-σ)/(9σ P0) dP0 == 0` | direct γ map (notes §5) | yes |
| A7 | sympy/wl | py126-129 / wl108-109 | `Delta_Q - dP0/P0 == 0` | even-preserving residual (notes §6) | yes |

All rows "yes": each residual feeds an independently-computed `delta_u2`/`delta_P0`/`delta_u4`
(from series/Solve) into the literal Stage-170/171 obstruction-and-outlet definitions, so a sign
or factor error anywhere in the derivation would produce a nonzero residual and trip the assertion.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent re-derivation**, not a transliteration. Decisive evidence: the
SymPy script derives the slopes by Taylor-expanding the *explicit* quotient
`u2A = -DA2/DA0` (py L52: `sp.series(u2A, eps, 0, 2).removeO()`), whereas the Mathematica
script deliberately differentiates the *implicit* defining relation and carries an explicit
comment to that effect:

- wl L42: `(* Independent route: differentiate the implicit defining relations, not the explicit quotient. *)`
- wl L43-44: `du2sol = Solve[D[(u2 + t*deltaU2)*(d0 + t*dD0) + (d2 + t*dD2), t] == 0 /. t -> 0, deltaU2]; deltaU2 = FullSimplify[deltaU2 /. First[du2sol], ...]`
- (SymPy counterpart, py L52) `delta_u2 = sp.simplify((sp.series(u2A, eps, 0, 2).removeO() - u2) / (eps * lam))`

The same divergence holds for `delta_P0` (wl L45-46 `Solve[D[(p0 + t*deltaP0)(d0+t*dD0) - (n0+t*dN0), t]==0, deltaP0]` vs py L53 series) and for `delta_u4` (wl L81-84 implicit
`Solve` on `u4·D0² = D2² - D0 D4` vs py L94-95 series). The Mathematica also uses different
symbol names (lowercase `d0,n0` vs `D0,N0`) and a different verification harness (`Solve`-based
slopes + `FullSimplify[Together[Expand[...]]]`). The two engines converge on the identical slope
expressions, but by genuinely distinct algebraic routes. Confidence: **high** (independent).

## Engine cross-check

Both engines emit the identical slopes and all seven residuals = 0:

| quantity | SymPy out | Mathematica out |
|---|---|---|
| `delta u_2^(A)` | `(-dD0*u2 - dD2)/D0` (L5) | `-((dD2 + dD0*u2)/d0)` (L5) |
| `delta P_0^(A)` | `(D0*dN0 - N0*dD0)/D0**2` (L6) | `(d0*dN0 - dD0*n0)/d0^2` (L6) |
| `delta u_4` canonical | `(-5*dD0 - 18*dD2 - 81*dD4)/(81*D0)` (L23) | `-1/81*(5*dD0 + 18*dD2 + 81*dD4)/d0` (L26) |
| all 7 residuals | `= 0` | `= 0` + `PASS` |

These are the same expressions modulo display normalization (case + sign factoring). Engines agree.

## Verdict justification

`clean`. I read the paper card, notes, and appendix block first and built the six-deliverable
model, then attacked the scripts: I checked that each `expect_zero` residual feeds an
independently-derived slope (not a pre-baked constant) into the literal Stage-170/171
obstruction/outlet definitions — so none is tautological; I checked the literal `1/9` in `K_A`,
the `(2/3,1/27)` even-consistency coefficients, the `8/9` relation, the `3(1-σ)/σ` and
`(1-σ)/(9σ)` outlet prefactors, and the `-9σ/(1-σ)` `Delta_Q` factor all match the notes and
appendix verbatim; I checked symbol domains (`D0,N0,sigma` declared `nonzero`, which is required
for the divisions and is consistent with the physical setup). I tried to find a stale-output
mismatch (both `.txt` are newer than their scripts and reproduce the current source exactly) and
a transliteration (the engines use materially different derivation routes — explicit-quotient
series vs implicit-relation `Solve`). Every attack failed. Paper claim and script claim match
exactly across all seven deliverables. This is the last unit of the 105–175 first-pass
orchestrator-direct band, and the transliteration watch ends clean here.

## Value Reconciliation (pass-2 augmentation)

All deliverable values the scripts emit are symbolic closed forms (this stage pins no new
numeric figures-of-merit; the only numeric constants are the canonical `u_2=1/9`, `u_4=4/81`,
carried in from upstream and used as branch substitutions). Watch-item constants
(168π²/100π², 4107, Family-1 radius √(4107−100π²)/(10π)) do **not** appear in this stage —
confirmed by grep over both scripts; this is a purely symbolic slope-collapse stage, so that
class is inapplicable here.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| \(\delta u_2^{(A)}=-(\delta D_2+u_2\delta D_0)/D_0\) | py L52-55 / wl L43-48; out L5 | notes L129-133; appendix eq (slope) | MATCH |
| \(\delta P_0^{(A)}=(\delta N_0-P_0\delta D_0)/D_0\) | py L53-56 / wl L45-49; out L6 | notes L134-137 | MATCH |
| \(\mathcal K_A=-D_0\delta u_2+(1/9-u_2)\delta D_0\) | py L58,62-65 / wl L51,55-58 | notes L141-145 (box) | MATCH |
| \(\mathcal G_A=D_0\delta P_0^{(A)}\) | py L59,66-69 / wl L52,59 | notes L148-153 (box) | MATCH |
| canonical collapse \(\mathcal K_A=-D_0\delta u_2\) at \(u_2=1/9\) | py L72-75 / wl L62-65 | notes L161-164 (box); tex L15 Output | MATCH |
| \(\mathfrak K_1=-D_0u_2^{(1)},\ \mathfrak G_1=D_0P_1\) | carry-forward print py L132-133 / wl L113-114 | tex L15 Output; notes L31-35,197-200; appendix L75, eq:app-part05-KG-physical | MATCH |
| \(\delta u_4=-(5\delta D_0+18\delta D_2+81\delta D_4)/(81D_0)\) | py L95-98 / wl L81-87; out L23/L26 | notes L240-243; appendix eq:app-part05-u4-slope-stage173 | MATCH |
| \(\delta u_4=(8/9)\delta u_2\iff\delta D_4=(2/3)\delta D_2+(1/27)\delta D_0\) | py L100-104 / wl L89-93 | notes L247-249 (box), L228-232; appendix eq:app-part05-hidden-even-D41-general | MATCH |
| \(\delta\kappa_W^{(A)}=-3(1-\sigma_*)/\sigma_*\,\delta u_2^{(A)}\) | py L112,115-118 / wl L96,99-102 | notes L278-283 (box); appendix eq:app-part05-dk-dg-physical | MATCH |
| \(\delta\gamma_W^{(A)}=-(1-\sigma_*)/(9\sigma_*)\,\delta P_0^{(A)}/P_0\) | py L113,120-123 / wl L97,103-106 | notes L285-292 (box); appendix eq:app-part05-dk-dg-physical | MATCH |
| \(\Delta_Q^{(A)}=\delta P_0^{(A)}/P_0\) | py L125-129 / wl L108-109 | notes L335-339, L344-350 (box) | MATCH |
| canonical \(u_2=1/9,\ u_4=4/81\) | py L83-84 / wl L68-69 | notes L74-80; appendix L513 | MATCH |

INTERNAL (scaffolding, no finding expected): `eps`, `lam` (perturbation bookkeeping symbols),
`t` (Mathematica differentiation parameter), `D2 = -u2*D0`, `D4_star = D0*(u2_star^2 - u4_star)`
(intermediate definitions), `K_A`/`G_A` literal Stage-171 carry-ins (verified upstream),
`delta_kappa`/`delta_gamma` literal Stage-170 carry-ins, `Solve`/`series`/`removeO` machinery,
PASS flags and banners.

reconciliation: complete; 12 deliverable values checked, 0 misaligned

## Self-test notes

Variable-independence: each `delta_*` genuinely depends on the differentiation/expansion
variable (`eps` for SymPy series, `t` for Mathematica `Solve`), and on `dD0/dD2/dD4/dN0`; no
identically-zero derivative trap — the printed nonzero slopes confirm dependence. No integrals,
so the parity trap is N/A. Trivial-case spot check: substituting the even-consistency relation
`dD4=(2/3)dD2+(1/27)dD0` into `delta_u4 - (8/9)delta_u2` collapses the
`-(5dD0+18dD2+81dD4)/(81D0)` numerator to `(8/9)·(-(dD0/9+dD2)/D0)`, i.e. residual 0, matching
output L24/L27. Symbol domains (`D0,N0,sigma` nonzero) are required for and consistent with the
quotient definitions. No directive written (zero findings).
