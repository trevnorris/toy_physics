---
unit_id: 173
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage173_axisymmetric_loading_mismatch.md]
  paper_appendix: present
---

# Audit unit 173 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_173.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage173_axisymmetric_loading_mismatch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 77, 265, 490-552)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage173_axisymmetric_loading_mismatch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage173_axisymmetric_loading_mismatch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.txt`

## What the paper claims

The stage computes the weak-axisymmetric physical grouped slopes directly from the actual moving-throat grouped operator/transfer moment expansions, then collapses the remaining first-order defect to a single scalar. The card's `\stagefield{Output}` states verbatim: "On the even-preserving branch, the remaining grouped defect is the scalar \(\Xi_{\rm load}=N_{01}/N_0-D_{01}/D_0\)." The notes and appendix enumerate the full deliverable set: (D1) the general slope law \(u_2^{(1)}=-(D_{21}+u_2 D_{01})/D_0\); (D2) the canonical \(u_4^{(1)}=-(5D_{01}+18D_{21}+81D_{41})/(81D_0)\) on \((u_2,u_4)=(1/9,4/81)\); (D3) the prefactor slope \(P_1/P_0=N_{01}/N_0-D_{01}/D_0\); (D4) the hidden-even operator law \(D_{41}=\tfrac23 D_{21}+\tfrac1{27}D_{01}\); (D5) the even-preserving collapse \(D_{21}=-D_{01}/9\), \(D_{41}=-D_{01}/27\); (D6) the scalar \(\Xi_{\rm load}=N_{01}/N_0-D_{01}/D_0\) with lane form \(\Delta_Q^{(20)}=\epsilon\Xi,\ \Delta_Q^{(21)}=\tfrac\epsilon2\Xi,\ \Delta_Q^{(22)}=-\epsilon\Xi\). The Stage-024 axisymmetric signature is \((\lambda_{20},\lambda_{21},\lambda_{22})=(1,1/2,-1)\).

## What the script claims to verify

Both scripts build the lane expansions \(D_{A,0}=D_0+\epsilon\lambda D_{01}\), etc., series-expand the exact response definitions \(u_2=-D_2/D_0\), \(u_4=(D_2^2-D_0 D_4)/D_0^2\), \(P_0=N_0/D_0\) to first order in \(\epsilon\), and extract the slope coefficients. They then assert: the series-derived \(u_2^{(1)}\) equals the boxed law (with \(u_2\to -D_2/D_0\)); the canonical \(u_4^{(1)}\) equals \(-(5D_{01}+18D_{21}+81D_{41})/(81D_0)\); \(P_1/P_0=N_{01}/N_0-D_{01}/D_0\); the hidden-even residual vanishes; \(D_{21}=-D_{01}/9\) and \(D_{41}=-D_{01}/27\) follow from \(u_2^{(1)}=0\) plus the hidden-even relation solved via `solve`/`Solve`; and \(\Xi_{\rm load}\) plus its three lane multiples. This is exactly the paper's deliverable set.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1 `u2^(1) = -(D21+u2 D01)/D0` | py L54-57 / wl L45 `u2 slope identity` | match |
| D2 canonical `u4^(1)` | py L68-71 / wl L56 `u4 canonical formula` | match |
| D3 `P1/P0 = N01/N0 - D01/D0` | py L72-75 / wl L57 `P1/P0 formula` | match |
| D4 hidden-even `D41 = 2D21/3 + D01/27` | py L78-85 / wl L60-61 residual + solve | match |
| D5 `D21=-D01/9`, `D41=-D01/27` | py L88-95 / wl L64-73 | match |
| D6 `Xi_load` + lanes (1, 1/2, -1) | py L98-111 / wl L76-84 | match |

`paper_alignment: aligned`. Every paper-side deliverable has a faithful, non-tautological script-side check, and constants (1/9, 4/81, 1/27, the (1,1/2,-1) signature) all match. The only finding is engine-policy, not math.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 54-57 | `expect_zero(u21 - (-(D21+u2*D01)/D0).subs(u2:-D2/D0))` | D1 | yes |
| A2 | sympy | 68-71 | `expect_zero(u41_can + (5D01+18D21+81D41)/(81D0))` | D2 | yes |
| A3 | sympy | 72-75 | `expect_zero(P1_ratio - (N01/N0 - D01/D0))` | D3 | yes |
| A4 | sympy | 78-82 | `expect_zero(hidden_even_residual)` | D4 | yes |
| A5 | sympy | 94 | `expect_zero(u21_zero_D21 + D01/9)` | D5 | yes |
| A6 | sympy | 95 | `expect_zero(D41_even + D01/27)` | D5 | yes |
| A7 | math | 45 | `expectZero["u2 slope identity", ...]` | D1 | yes |
| A8 | math | 56 | `expectZero["u4 canonical formula", ...]` | D2 | yes |
| A9 | math | 57 | `expectZero["P1/P0 formula", ...]` | D3 | yes |
| A10 | math | 61 | `expectZero["hidden-even residual", ...]` | D4 | yes |
| A11 | math | 72-73 | `expectZero["D21+D01/9"/"D41+D01/27"]` | D5 | yes |

All assertions trace to a specific deliverable and are non-tautological: the LHS of each `expect_zero` is independently produced by a series expansion / `solve`, and the RHS is the paper's target formula. Substituting a wrong constant (e.g. 18→17 in A2, or 1/9→1/8 in the canonical substitution) would make the residual nonzero, so the checks can fail.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl:1-99`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage173_axisymmetric_loading_mismatch_sympy_audit.py:1-120`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. The variable choreography, intermediate objects, helper functions, substitution order, solve sequence, print order, and the trailing carry-forward text block are identical; only the syntax is translated to Mathematica.

Corresponding excerpts:

(1) Lane expansion + slope extraction — SymPy L37-48:
```python
D0A = D0 + eps*lam*D01
...
u2A = sp.expand(sp.series(-D2A/D0A, eps, 0, 2).removeO())
...
u21 = sp.simplify(sp.diff(u2A, eps).subs(eps, 0)/lam)
```
Mathematica L32-39:
```mathematica
d0A = d0 + eps*lam*d01;
...
u21 = FullSimplify[Coefficient[Series[-d2A/d0A, {eps, 0, 1}] // Normal, eps, 1]/lam, ...];
```
Same `-d2A/d0A` ratio, same eps-series, same "order-1 coefficient over lam". Identical variable names (`d0A`, `u21`, `u41`, `p1`).

(2) Even-preserving collapse — SymPy L88-92:
```python
u21_zero_D21 = sp.simplify(sp.solve(sp.Eq(u21_can, 0), D21)[0])
...
D41_even = sp.simplify(D41_hidden.subs(D21, u21_zero_D21))
```
Mathematica L64-69:
```mathematica
u21ZeroD21 = FullSimplify[d21 /. First[Solve[u21Can == 0, d21]], ...];
...
d41Even = FullSimplify[d41 /. First[Solve[(u41Can == 8 u21Can/9) /. d21 -> u21ZeroD21, d41]], ...];
```
Same `Solve[u21==0, d21]`/`First` choreography, same substitution into the hidden-even relation.

(3) The carry-forward `Print` block (SymPy L113-119, Mathematica L87-93) is byte-for-byte identical text.

**Why this matters:**
The second-engine policy requires Mathematica to derive the result independently from the physical premises so the two engines provide genuinely independent confirmation. A transliteration only re-checks that the same algebra was typed correctly in two syntaxes; a shared conceptual error (e.g. an incorrect lane signature or an incorrect canonical \(D_2/D_0=-1/9\) relation) would pass both engines identically. This is the recurring transliteration class flagged across the 105-175 band.

**Required change:**
Re-author the `.wl` to derive the slope laws by an independent route rather than echoing the `.py` choreography. Concretely: instead of porting the `Series`/`Coefficient` pipeline that mirrors SymPy's `series`/`diff`, derive each first-order slope via direct implicit differentiation of the defining relations at \(\epsilon=0\) — e.g. for \(u_2^{(A)}=-D_{A,2}/D_{A,0}\), compute \(\partial_\epsilon u_2^{(A)}|_{0}/\lambda = -(D_{21} D_0 - D_2 D_{01})/D_0^2\) using `D[..., eps]` on the ratio directly (no intermediate `d2A/d0A` named exactly as in the .py), and confirm it equals the boxed \(-(D_{21}+u_2 D_{01})/D_0\) with \(u_2=-D_2/D_0\). Use distinct intermediate naming and remove the verbatim carry-forward text block (or replace it with a Mathematica-native summary). The acceptance criterion is that the `.wl` no longer reproduces the `.py`'s exact variable choreography and that all six `expectZero` checks still pass.

**Verification:**
After re-authoring, the verifier confirms (a) the `.wl` exits 0 with all `expectZero` checks PASS, and (b) the structural correspondence is broken — the `.wl` no longer shares the `.py`'s named intermediate `d2A/d0A`-style choreography and verbatim print block. The output residuals must still be 0 for all six identities.

## Independent-derivation check (Mathematica)

Not independent. See F1. The `.wl` mirrors the `.py` step-for-step: same lane-expansion variables, same `Series`→order-1-coefficient extraction (`Coefficient[Series[...], eps, 1]` vs `diff(...).subs(eps,0)`), same canonical substitutions `d2 -> -(1/9) d0`, `d4 -> -(1/27) d0`, same `Solve`/`First` for the even-preserving collapse, and a byte-identical carry-forward print block. Confidence: high that this is a transliteration.

## Engine cross-check

Both engines present and agree exactly. Side-by-side of the load-bearing outputs:

| quantity | sympy output | mathematica output |
|---|---|---|
| u2^(1) general | `(-D0*D21 + D01*D2)/D0**2` | `(d01*d2 - d0*d21)/d0^2` |
| u4^(1) general | `(-D0**2*D41 + D0*(D01*D4 + 2*D2*D21) - 2*D01*D2**2)/D0**3` | `(-2*d01*d2^2 + 2*d0*d2*d21 + d0*d01*d4 - d0^2*d41)/d0^3` |
| P1 general | `(D0*N01 - D01*N0)/D0**2` | `(-(d01*n0) + d0*n01)/d0^2` |
| u4 canonical | `(-5*D01 - 18*D21 - 81*D41)/(81*D0)` | `-1/81*(5*d01 + 18*d21 + 81*d41)/d0` |
| Xi_load | `N01/N0 - D01/D0` | `-(d01/d0) + n01/n0` |

All identical up to ordering. Every `expect_zero`/`expectZero` returned 0 / PASS in both transcripts. No `engine_disagreement`.

## Verdict justification

The math is correct and fully aligned with the paper. I attacked each assertion: the `u2 slope identity` is not tautological because `u21` is produced by an independent series expansion and is then compared to the boxed formula re-expressed in operator variables; the canonical and even-preserving checks would fail under any perturbed constant (18, 81, 1/9, 1/27); the `solve`-derived `D21`/`D41` are confirmed against the paper's stated values; the (1,1/2,-1) signature matches Stage 024. I read the paper card, notes, and appendix and confirmed the six deliverables map one-to-one to the six load-bearing checks. The single defect is that the Mathematica engine is a transliteration of the SymPy script rather than an independent re-derivation (F1), which weakens the two-engine guarantee. Verdict: `findings` (one `mathematica_transliteration`, medium). No `paper_misalignment`, no stale outputs, no stop-cold.

## Value Reconciliation (pass-2 augmentation)

Procedure: enumerated every labeled RESULT/deliverable value emitted by the scripts (from `.py`/`.wl` source plus the two saved `.txt` transcripts) and located each in the `.tex` card and `.md` notes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `u2^(1) = -(D21+u2 D01)/D0` | py L50/56, wl L41/45; sympy.txt L5, math.txt L5 | card `\stagefield{Derivation ledger}` (implied) + notes L40-42, L156-158; appendix L505-507 | MATCH |
| `u4^(1) = -(5D01+18D21+81D41)/(81D0)` | py L65/70, wl L53/56; sympy.txt L14, math.txt L15 | notes L185-188; appendix L508-514 | MATCH |
| `P1/P0 = N01/N0 - D01/D0` | py L66/74, wl L54/57; sympy.txt L15, math.txt L19 | card Output L15; notes L53-59, L161-167; appendix L516-520 | MATCH |
| `D41 = 2 D21/3 + D01/27` (hidden-even) | py L85, wl carry L91; sympy.txt L23 | notes L197-200; appendix L528-530 | MATCH |
| `D21 = -D01/9` (even-preserving) | py L89, wl L65; sympy.txt L28, math.txt L31 | notes L72-74, L214-216; appendix L523-525 | MATCH |
| `D41 = -D01/27` (even-preserving) | py L92, wl L70; sympy.txt L29, math.txt L32 | notes L82-84, L221-223; appendix L533-535 | MATCH |
| `Xi_load = N01/N0 - D01/D0` | py L98, wl L76; sympy.txt L36, math.txt L41 | card Output L15; notes L89-92, L277-279; appendix L538-542 | MATCH |
| `Delta_Q^(20)/eps = Xi_load` | py L105/109, wl L82; sympy.txt L37, math.txt L42 | notes L97-103, L289; appendix L547 | MATCH |
| `Delta_Q^(21)/eps = Xi_load/2` | py L106/110, wl L83; sympy.txt L38, math.txt L43 | notes L100/291; appendix L549 | MATCH |
| `Delta_Q^(22)/eps = -Xi_load` | py L107/111, wl L84; sympy.txt L39, math.txt L44 | notes L102/293; appendix L551 | MATCH |
| lane signature `(lam20,lam21,lam22)=(1,1/2,-1)` | py L101-103, wl L79-81 | notes L22-24, L139-141; card Checks L21 `(1,1/2,-1)` | MATCH |
| canonical `u2=1/9, u4=4/81` (substitutions D2=-D0/9, D4=-D0/27) | py L60-61, wl L48-49 | notes L177-181, L240-243; appendix L513 | MATCH |

INTERNAL (scaffolding, no finding expected): `u2^(1) general`, `u4^(1) general`, `P1 general` (intermediate operator-variable forms before canonical substitution); residual print values (all 0); PASS flags; `D41 from hidden-even relation` intermediate; banner strings.

No stale-constant class (`168π²`/`100π²`) or Family-1 radius `√(4107−100π²)/(10π)` appears in this stage — those constants are not part of stage 173 (a symbolic operator-slope stage). Searched the watch terms; none present.

reconciliation: complete; 12 values checked, 0 misaligned.

## Self-test notes

(1) Variable independence: the only derivative-style operations are `series`/`Coefficient` order-1 extraction; the integrands (`-D2A/D0A`, `(D2A^2-D0A*D4A)/D0A^2`, `N0A/D0A`) genuinely depend on `eps` through the `eps*lam*D01` terms, so the order-1 coefficients are non-trivial (confirmed by the nonzero general-form outputs). (2) No unbounded integrals, so no parity trap. (3) Trivial-case spot check: setting `D01=0, D21=0` makes `u21_zero_D21=0` consistent and `Xi_load=N01/N0` nonzero, matching expectation; the asserted residuals reduce to literal 0. F1 is engine-policy only and prescribes a re-derivation route (direct `D[ratio, eps]`) that does not introduce any new constant, so the paper round-trip is preserved.
