---
unit_id: 131
batch: IV.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md"]
  paper_appendix: present
---

# Audit unit 131 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_131.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage131_parent_mouth_threshold.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_131}` at line 1296 + a `\clearpage`; no separate appendix narrative row for this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.txt`

## What the paper claims

The card (status `\StatusExactClosure`) states the canonical compensated branch is the parent micro-threshold quoted verbatim as: "Canonical branch is the parent threshold \(T_m-q_*A'_0=\Pi_*\Theta_\sigma/L\)." (`stage_131.tex:16`). The notes expand this: translate the Family-1 compensation point \(\Pi_*\approx 1.50882951349316\) into a direct parent-level threshold on the localized-Maxwell + confinement mouth data. Deliverables: (D1) the canonical branch condition \(\Pi_m=\Pi_*\), equivalently \(V_1=\Pi_*\Theta_\sigma/L\) (notes Sec.1, boxed); (D2) the localized-Maxwell/mechanical split \(V_1=T_m-q_*A_0'\) giving \(T_m-q_*A_0'=\Pi_*\Theta_\sigma/L\) (notes Sec.2, boxed = card line 16); (D3) the linearized deviation law with slope \(\mathfrak g'(\Pi_*)\approx 0.0714453558083195\) (notes Sec.3, boxed). The card's Checks demand positivity of the mouth source, zero-flux/boundary-layer normalizations, and that the Family-1 compensation point sits on the LOWER branch, not the singular equal-normalized branch (card lines 21-25).

## What the script claims to verify

The SymPy script (a) confirms the carried-forward F1 lower-branch value \(g_-^{F1}\) closed form \((2\sqrt{4107-100\pi^2}-37\sqrt3)/(20\pi)\) equals its literal `0.758035078944663`; (b) solves the monotone bias function \(g_\Pi(\Pi)=2\Pi(2\Pi e^\Pi+\pi)/((4\Pi^2+\pi^2)(e^\Pi-1))\) against \(g_\Pi=g_-^{F1}\) for the root \(\Pi_*\) and asserts it matches notes Sec.1 `1.50882951349316`; (c) computes the slope \(g'(\Pi_*)\) and asserts it matches notes Sec.3 `0.0714453558083195`; (d) prints the parent bias mismatch formula \((T_m-q_*A_0')-\Pi\Theta_\sigma/L\); and (e) performs branch discrimination per Checks item 3 — lower-branch membership (residual `<1e-30`), exclusion of the singular equal-normalized branch \(g_{\rm nat}=1\) (separation matches notes \(\Delta g_-=0.241964921055337\)), and exclusion of the upper branch \(g_+^{F1}\) (separation `>1`). The Mathematica script mirrors the same physical claims but reaches \(\Pi_*\) by an independent denominator-cleared polynomial `FindRoot` with a bracketing seed pair.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1: canonical branch \(\Pi_m=\Pi_*\), \(V_1=\Pi_*\Theta_\sigma/L\) | `Pi_star` solved from `gPi - g_minus`; asserted == `1.50882951349316`; `V1 = Pi*Theta_sigma/L`, `V1_*` printed (py:24-29,42-45) | match |
| D2: split \(T_m-q_*A_0'=\Pi_*\Theta_\sigma/L\) (card line 16) | `threshold_residual = (Tm - qstar*A0p) - Pi*Theta_sigma/L` printed (py:34-35); identical symbolic form | match |
| D3: slope \(g'(\Pi_*)\approx 0.0714453558083195\) | `gprime_star = diff(gPi,Pi)` at `Pi_star`; asserted == `0.0714453558083195` (py:31-32,47-53) | match |
| Checks item 3: lower (not singular) branch | (4a) membership `<1e-30`; (4b) singular `g_nat=1` excluded, separation == \(\Delta g_-\); (4c) upper branch excluded (py:55-88) | match |
| Checks items 1-2 (positivity, zero-flux/boundary-layer normalizations) | not directly exercised — these are imported-input qualitative gates from upstream stages 125/129, not algebraic deliverables of this stage | n/a (carry-forward gate, not a numeric claim) |

Dominant pattern: aligned. The card's only displayed equation (line 16) and all three boxed notes deliverables are exercised by non-tautological assertions. Checks items 1-2 are upstream-import qualitative conditions, not algebraic targets verifiable within this stage's scope.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 17-19 | `abs(N(g_minus_exact) - literal) < 1e-14` | upstream input to D1/D3 (g_-^F1) | yes |
| A2 | sympy | 42-44 | `abs(Pi_star_num - 1.50882951349316) < 1e-14` | D1 | yes |
| A3 | sympy | 50-52 | `abs(gprime_num - 0.0714453558083195) < 1e-14` | D3 | yes |
| A4 | sympy | 66-68 | `\|g_Pi(Pi_*) - g_-^F1\| < 1e-30` | Checks item 3 (lower membership) | yes |
| A5 | sympy | 74-77 | `\|sing_sep - Delta_g_-\| < 1e-12` | Checks item 3 (singular excluded) | yes |
| A6 | sympy | 78-80 | `sing_sep > 1e-3` | Checks item 3 (singular excluded) | yes |
| A7 | sympy | 85-87 | `upper_sep > 1` | Checks item 3 (upper excluded) | yes |
| M1 | mathematica | 37-38 | `expectApprox[g_minus, literal, 1e-14]` | upstream input | yes |
| M2 | mathematica | 61-62 | `expectApprox[piStar, 1.50882951349316, 1e-14]` | D1 | yes |
| M3 | mathematica | 65-66 | `expectApprox[slope, 0.0714453558083195, 1e-14]` | D3 | yes |
| M4 | mathematica | 76-80 | `lowerResidual < 1e-30` | Checks item 3 (lower) | yes |
| M5 | mathematica | 83-88 | `\|singSep - Delta_g_-\| < 1e-12 && singSep > 1e-3` | Checks item 3 (singular) | yes |
| M6 | mathematica | 91-95 | `upperSep > 1` | Checks item 3 (upper) | yes |

Notes on anchoring. A2/M2 are genuine: `Pi_star` is found by a numerical root solve of the monotone `gPi == g_minus`, NOT pinned to the literal; the assertion then independently confirms the solved root reproduces the documented value. A3/M3 differentiate the closed-form `gPi` and evaluate at the solved root — neither hardcoded nor tautological. A4/M4 (membership) could superficially look tautological because `Pi_star` was solved from `gPi == g_minus`, but the membership check uses a much tighter tolerance (`1e-30`) than the solver's anchor tolerance and confirms the root quality; combined with A5-A7, it is the discriminator that distinguishes lower from singular/upper branches as the card's Checks item 3 demands. A1/M1 verify the carried-forward `g_-^F1` closed form against its decimal literal, defending against a back-fit constant.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the SymPy script. Three corresponding sections justify this:

1. Root finding. SymPy: `Pi_star = sp.nsolve(gPi - g_minus, 1.5, tol=1e-30, ...)` — a single-seed Newton solve on the rational equation `gPi == g_minus` (py:24). Mathematica deliberately takes a DIFFERENT route (wl:41-47): it clears denominators into a polynomial-in-`(piM, Exp[piM])` residual `gThresholdResidual[p] := 40*Pi*p*(2*p*Exp[p]+Pi) - 20*Pi*gMinus*(4*p^2+Pi^2)*(Exp[p]-1)` and applies `FindRoot[... == 0, {piM, 1.4, 1.6}]` with a bracketing seed PAIR. The denominator-cleared polynomial form and the two-point bracket are an independent numerical strategy, not a syntactic port of `nsolve`. The inline comment (wl:41) explicitly states this design intent.

2. Branch closed forms. Both engines build `g_-^F1` and `g_+^F1` from the same upstream stage122 closed form `(2√(4107-100π²) ∓ 37√3)/(20π)`; this is shared input data (a carried-forward constant), not algebra echoed between engines, which is the legitimate way both engines anchor to the same physical premise.

3. Branch discrimination. Both compute lower/singular/upper separations, but the choreography is independent: SymPy uses `sp.N(..., 40)`/`sp.Float` comparisons; Mathematica uses `N[..., 40]`/`TrueQ[Abs[...] < 10^-30]`. Same physical targets, separately coded.

Verdict on independence: PASS. No `mathematica_transliteration` finding.

## Engine cross-check

The two engines agree to the precision they claim:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `Pi_*` | `1.50882951349315558300555075595` (sympy.txt:2) | `1.5088295134931555830055507559542749...` (math.txt:7) |
| `g'(Pi_*)` | `0.0714453558083195211894603019881` (sympy.txt:4) | `0.07144535580831952118946030198813575...` (math.txt:9) |
| Parent bias formula | `-A0p*q_* + T_m - Pi*Theta_sigma/L` (sympy.txt:5) | `-(a0Prime*qStar) - (piM*thetaSigma)/lM + tM` (math.txt:10) |
| Lower-branch residual | `5.5e-42` (sympy.txt:8) | PASS, `<1e-30` (math.txt:15) |
| Singular separation | `0.241964921055337173080319109586` (sympy.txt:9) | PASS, == `Delta g_-` (math.txt:16) |
| Upper separation | `2.03991691306063058319190804376` (sympy.txt:10) | PASS, `>1` (math.txt:17) |

`Pi_*` agrees to ~29 digits; `g'(Pi_*)` agrees to ~29 digits; the symbolic bias formula is identical up to term ordering and symbol-name choice. No `engine_disagreement`.

## Verdict justification

`clean`. I read the card, the single notes file, and the part04 appendix context (which only `\input`s the stage). The card's sole displayed equation `T_m-q_*A'_0=\Pi_*\Theta_\sigma/L` (line 16) is exercised verbatim by both scripts' bias-mismatch formula and the `Pi_*` anchor; both boxed notes constants (`Pi_*=1.50882951349316`, slope `0.0714453558083195`) are reproduced by genuine root-solve + differentiation, not hardcoded; the carried-forward `g_-^F1` closed form is checked against its literal so no back-fit constant slips through. Attacks I tried and which failed: (i) tautology on A2/A4 — defeated because `Pi_*` is solved from the monotone `gPi == g_minus` and the membership check uses a strictly tighter tolerance than the anchor; (ii) hardcoded-result on the closed form — defeated by A1/M1 checking the closed form against the decimal and the closed form being anchored upstream at stage122 line 56; (iii) transliteration — defeated, the Mathematica root route is denominator-cleared polynomial `FindRoot` with a bracket, structurally distinct from SymPy `nsolve`; (iv) wrong-branch — the script explicitly excludes both the singular `g_nat=1` supremum (separation matches notes `Delta g_- = 0.241964921055337`) and the upper branch, satisfying Checks item 3; (v) symbol-domain — `Pi`/`piM` positive-real for the physical bias and `pi`/`Pi` (constant) are correctly distinguished, with positivity justified by the card's reduced branch. Checks items 1-2 are upstream-import qualitative gates (positivity, zero-flux/boundary-layer normalizations from stages 125/129), not algebraic deliverables verifiable inside this stage, so their absence from the script is not a finding.

## Self-test notes

Checked traps: (1) variable independence — `sp.diff(gPi, Pi)`/`D[gPi, piM]` differentiate w.r.t. the bias parameter that `gPi` genuinely depends on (it appears in numerator, denominator and exponent), so the slope is genuinely nonzero (`0.0714…`), not an identically-zero derivative. (2) No unbounded-domain integrals in this stage, so the parity trap does not apply. (3) Trivial-case: the singular and upper separations are positive nonzero literals (`0.2419…`, `2.0399…`) and the lower-membership residual collapses to `~1e-42`, all consistent with the monotone `g_Π: 2/π → 1` law from stage130. (4) No missing-script finding, so path specs n/a. (5) Paper round-trip: no fix prescribed, so no new misalignment risk. No directive written.

## Value Reconciliation (pass-2 augmentation)

Enumeration of every RESULT/deliverable value the scripts emit, with doc location and status. The committed outputs (`sympy_audit.txt`, `mathematica_audit.txt`, both mtime 2026-05-29 15:39:50, NEWER than both scripts at 15:31:13 / 15:31:33 — outputs are fresh) are the authoritative record.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `g_-^F1 = 0.758035078944663` (closed form `(2√(4107-100π²)-37√3)/(20π)`) | py:15-16, wl:35-36; checked vs literal in both outputs | notes stage130:98, stage122:63/56 (carried-forward upstream); referenced for this stage in notes Sec.3 context | MATCH |
| `Pi_* = 1.50882951349316` (printed `1.5088295134931555830…`) | py:24,28 / wl:46,54; sympy.txt:2, math.txt:7 | notes Sec.1 line 8, line 41, line 64; (stage130:105) | MATCH |
| `V1_* = Pi_* · Theta_sigma/L` (symbolic) | py:27,29 / wl:50,55; sympy.txt:3, math.txt:8 | notes Sec.1 boxed line 41 `V_1 = Π_*Θ_σ/L ≈ 1.50882951349316 Θ_σ/L` | MATCH |
| `g'(Pi_*) = 0.0714453558083195` | py:31-32 / wl:51,56; sympy.txt:4, math.txt:9 | notes Sec.3 line 92, line 100 | MATCH |
| Parent bias mismatch formula `(T_m - q_*A_0') - Π Θ_σ/L` | py:34-35 / wl:52,57; sympy.txt:5, math.txt:10 | card line 16 `T_m-q_*A'_0=Π_*Θ_σ/L`; notes Sec.2 boxed line 61-65 | MATCH |
| `Δg_- = 0.241964921055337` (singular-branch separation `g_nat - g_Π(Π_*)`) | py:72-73 / wl:83-84; sympy.txt:9 | notes stage122:104 `Δg_- = 1 - g_-^F1 ≈ 0.241964921055337` (cross-stage carry, anchored upstream) | MATCH |

INTERNAL items (verification scaffolding, no prose deliverable expected; no finding):
- upper-branch separation `2.03991691306063…` (sympy.txt:10 / math.txt:17) — a discrimination residual proving the upper branch `g_+^F1` is excluded; not a stated deliverable.
- `g_+^F1` closed form / its numeric value — discriminator only (the stage's deliverable is the LOWER branch); the closed form is anchored upstream (stage122:56) but the upper value is not a stage-131 deliverable.
- `g_nat = 1` — the supremum constant used as the singular-branch reference (notes Sec.1 / stage130:87 `2/π → 1`); used as a comparison anchor, not a reported result.
- lower-branch membership residual `~5.5e-42` — near-zero verification residual.
- closed-form-vs-literal residuals (`~1.7e-16`, etc.) — tolerance-check residuals.

reconciliation: complete; 6 deliverable values checked, 0 misaligned
