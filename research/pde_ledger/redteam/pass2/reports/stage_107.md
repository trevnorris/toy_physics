---
unit_id: 107
batch: IV.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage107_general_dtn_deformation.md]
  paper_appendix: present
---

# Audit unit 107 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_107.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage107_general_dtn_deformation.md` (only file present)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (subsection "General isotropic deformation algebra", lines 348–397; band summary line 86; input row line 1248)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage107_general_dtn_deformation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage107_general_dtn_deformation_mathematica_audit.txt`

## What the paper claims

Stage 107 introduces the **first explicit general isotropic `l=2` DtN deformation family** and derives the exact map from deformation data to the retarded quadrupole normalization `chi_Q`. The card quotes: *"Deformed branch \(S\Lambda_2^{\rm out}(\beta z)+\Sigma_0+\Sigma_2z^2+\Sigma_4z^4+i\Sigma_5z^5\) gives explicit \(\chi_Q\)."* The notes and appendix (eqs `app-part04-general-dtn-deformation`, `app-part04-Ydef-expanded`, `app-part04-chiQ-general`) make the deliverables precise:
1. The deformed branch expands as `Lambda_def = L0 + L2 z² + L4 z⁴ + i L5 z⁵` with `L0=-3S+Σ0, L2=Sβ²/3+Σ2, L4=Sβ⁴/9+Σ4, L5=Sβ⁵/9+Σ5`.
2. The normalized response `Ŷ_2^def = 1 - (L2/L0)z² + (L2²/L0² - L4/L0)z⁴ - i(L5/L0)z⁵`.
3. The canonical-even matching conditions `-L2/L0 = 1/9`, `L2²/L0² - L4/L0 = 4/81` give `Σ2 = -(3Sβ²-3S+Σ0)/9`, `Σ4 = -(3Sβ⁴-3S+Σ0)/27`.
4. With those imposed, the exact odd normalization is **`chi_Q = 3(Sβ⁵+9Σ5)/(3S-Σ0)`** (equivalently `chi_Q-1 = (3S(β⁵-1)+Σ0+27Σ5)/(3S-Σ0)`).

The appendix additionally records a linearized corollary `Δ_Q = 5b + a0/3 + 9a5 + O(ε²)` (eq `app-part04-linear-DeltaQ-triple`); this is an appendix-only consequence, not boxed in the notes' deliverable set or the stage-card quote.

## What the script claims to verify

Both engines build `Lambda_def` from `S·Lambda_out(βz)` plus the additive even/odd throat-core slots, extract `L0,L2,L4,L5`, and verify (a) the closed-form normalized expansion `Ŷ_2^def` matches the engine's own Taylor series of `L0/Lambda_def` (`expect_zero('normalized expansion direct-formula', ...)`); (b) the even-matching `Solve` for `(Σ2,Σ4)` returns a unique solution equal to the boxed formulas (`expect_zero('Sigma2 exact formula', ...)`, `'Sigma4 exact formula'`); and (c) substituting that solution into `chi_Q = (-L5/L0)/(1/27)` yields `3(Sβ⁵+9Σ5)/(3S-Σ0)` (`expect_zero('chi_Q - 3(...)/(3S - Sigma0)', ...)`). The docstring/banner topic is "general isotropic DtN deformation algebra."

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| Deformed-branch coefficients `L0,L2,L4,L5` | py L33–40 / wl L36–39 extract + print; consumed by all later checks | match |
| Normalized expansion `Ŷ_2^def` (eq `Ydef-expanded`) | `expect_zero('normalized expansion direct-formula', Y_direct − Y_formula)` (py L42–44 / wl L41–52) | match |
| Even-matching `Σ2,Σ4` (boxed, notes 95–97) | `Solve(m2=1/9, m4=4/81)` then `expect_zero('Sigma2/Sigma4 exact formula')` (py L54–67 / wl L62–69) | match |
| `chi_Q = 3(Sβ⁵+9Σ5)/(3S-Σ0)` (eq `chiQ-general`) | `expect_zero('chi_Q - 3(Sβ⁵+9Σ5)/(3S-Σ0)')` (py L69–74 / wl L71–76) | match |
| `chi_Q-1` form (eq `chiQ-minus-one-general`) | not separately asserted; algebraically `chi_Q − 1`, derivable from the verified `chi_Q` | partial (non-load-bearing) |
| Linearized `Δ_Q = 5b + a0/3 + 9a5` (appendix eq `linear-DeltaQ-triple`) | not exercised | missing (appendix-only corollary, outside card/notes deliverable set) |

`paper_alignment = aligned`: every boxed/load-bearing deliverable in the notes and card is faithfully verified by both engines with matching constants. The two unverified items are (i) an algebraically-equivalent restatement and (ii) an appendix-only linearization not in the stage card's stated scope; neither rises to a `paper_misalignment`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 44 | `expect_zero(Y_direct − Y_formula)` (engine series vs hand closed form) | claim 2 (`Ŷ_2^def`) | yes |
| A2 | sympy | 54–56 | `Solve(m2=1/9, m4=4/81)` unique-solution guard | claim 3 (even match) | yes |
| A3 | sympy | 60–67 | `expect_zero(Σ2 − boxed)`, `expect_zero(Σ4 − boxed)` | claim 3 | yes |
| A4 | sympy | 71–74 | `expect_zero(chi_even − 3(Sβ⁵+9Σ5)/(3S−Σ0))` | claim 4 (`chi_Q`) | yes |
| A5 | mathematica | 52 | `expectZero[yDirect − yFormula]` | claim 2 | yes |
| A6 | mathematica | 62–63 | `Solve[m2==1/9, m4==4/81]` unique-solution guard | claim 3 | yes |
| A7 | mathematica | 68–69 | `expectZero[Σ2 + (...)/9]`, `expectZero[Σ4 + (...)/27]` | claim 3 | yes |
| A8 | mathematica | 73–76 | `expectZero[chiEven − 3(sNorm β⁵+9 σ5)/(3 sNorm−σ0)]` | claim 4 | yes |

All eight rows are non-tautological. A1/A5 compare an **independently engine-computed Taylor series** of `L0/Lambda_def` against a hand-typed closed form (not a re-substitution of the same expression), so they can fail. A2/A6 use `Solve` (not back-substitution of a pre-baked answer). A4/A8 substitute the solved even-match values into the `chi_Q` definition and compare to the boxed odd-normalization — also genuinely falsifiable. No hardcoded answer is asserted against itself.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage107_general_dtn_deformation_mathematica_audit.wl:33–76`
- (parallel) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py:27–74`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Every step corresponds one-to-one in the same order, on the same expressions, with the same routine choices:

| step | sympy | mathematica |
|---|---|---|
| outgoing branch | L27 `Lambda_out = -3 + z**2/3 + z**4/9 + I*z**5/9` | L33 `lambdaOut = -3 + z^2/3 + z^4/9 + I*z^5/9` |
| build deformation | L30 `Lambda_def = expand(S*Lambda_out.subs(z,beta*z) + Sigma0 + Sigma2*z**2 + …)` | L34 `lambdaDef = Expand[sNorm*(lambdaOut /. z->beta*z) + sigma0 + sigma2*z^2 + …]` |
| extract coeffs | L33–36 `L0=…subs(z,0)`, `L2=coeff(z,2)`, `L4=coeff(z,4)`, `L5=im(coeff(z,5))` | L36–39 `l0=…/z->0`, `l2=Coefficient[…,z,2]`, `l4=Coefficient[…,z,4]`, `l5=Coefficient[…,z,5]/I` |
| series | L42 `Y_direct = sp.series(L0/Lambda_def, z, 0, 6)` | L41 `yDirect = Series[l0/lambdaDef, {z,0,5}]` |
| hand closed form | L43 `Y_formula = 1 - (L2/L0)*z**2 + (L2**2/L0**2 - L4/L0)*z**4 - I*(L5/L0)*z**5` | L42 `yFormula = 1 - (l2/l0)*z^2 + (l2^2/l0^2 - l4/l0)*z^4 - I*(l5/l0)*z^5` |
| even match | L54 `solve([Eq(m2,1/9), Eq(m4,4/81)], [Sigma2,Sigma4])` | L62 `Solve[{m2==1/9, m4==4/81}, {sigma2,sigma4}, Reals]` |
| assertion order | normalized-expansion → Σ2 → Σ4 → chi_Q | normalized-expansion → Σ2 → Σ4 → chi_Q (identical) |

Both engines call the SAME `Series[]`/`series()` on the SAME expression `L0/Lambda_def` and extract coefficients the same way — the canonical transliteration signature flagged in the 105–175 first-pass orchestrator-direct band watch. The second-engine policy requires an independent route (e.g., derive `chi_Q` from `−L5/L0` via the geometric-series obstruction algebra without a `Series[]` on the identical ratio, or invert the normalized law symbolically rather than re-typing the same coefficient extraction).

**Why this matters:**
A transliterated second engine does not independently corroborate the algebra; it re-runs the same choreography in a different syntax. If the SymPy derivation embedded a subtle convention error (e.g., wrong parity convention on `L5`, sign on the odd slot), the `.wl` would reproduce it rather than catch it. The math here is in fact correct and both engines agree, but the independence guarantee the dual-engine policy exists to provide is absent.

**Required change:**
This is a re-author-vs-accept decision at the user level (per the IV.1 lesson: dual-engine re-authoring is a USER-LEVEL call, not orchestrator/Codex discretion). If re-authoring is authorized, rewrite the `.wl` so its `chi_Q` and normalized-law derivation does NOT reuse `Series[l0/lambdaDef, …]` on the identical ratio and does NOT mirror the `.py` coefficient-extraction order — e.g., obtain `Ŷ_2^def`'s odd coefficient by directly inverting `L0/Lambda_def = 1/(1 + (L2 z² + L4 z⁴ + i L5 z⁵)/L0)` via an explicit geometric expansion (`Sum`/manual `O(z^6)` truncation), then read off `−L5/L0` analytically. Keep all emitted values byte-identical to the current output.

**Verification:**
The rewritten `.wl` must not contain `Series[l0/lambdaDef, …]`; its derivation order should differ from the `.py`; and `math -script` must still print `chi_Q under canonical-even matching = (-3*(9*sigma5 + beta^5*sNorm))/(sigma0 - 3*sNorm)` and all four PASS lines, exiting 0. The verifier confirms the route is independent (no `Series[]` on the shared ratio) and the output is byte-identical to the committed `.txt`.

## Independent-derivation check (Mathematica)

Not independent — transliteration. See F1's step-by-step correspondence table. The decisive evidence is `Series[l0/lambdaDef, {z,0,5}]` (wl L41) vs `sp.series(L0 / Lambda_def, z, 0, 6)` (py L42): both engines apply the same series operator to the same rational expression, then compare against an identically-typed hand closed form, in the same assertion order. There is no second route to `chi_Q`.

## Engine cross-check

Both engines produce identical results. SymPy: `chi_Q = 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0)`; Mathematica: `chi_Q = (-3*(9*sigma5 + beta^5*sNorm))/(sigma0 - 3*sNorm)` — algebraically the same expression (numerator/denominator both negated). `Σ2_evenmatch = -S*beta**2/3 + S/3 - Sigma0/9` agrees on both sides. All `expect_zero`/`expectZero` checks emit `0`/`PASS`. `engines_agree = true`. The agreement is expected and unremarkable precisely *because* the second engine is a port (F1).

## Verdict justification

All eight load-bearing assertions are non-tautological, correctly anchored to the paper's deliverables, and the emitted constants (`L0…L5`, `Σ2`, `Σ4`, `chi_Q`) match the notes (lines 50–54, 95–97, 104–106) and the appendix (eqs `Ydef-expanded`, `chiQ-general`) exactly. Attacks tried and failed: (1) is `Y_direct − Y_formula` a round-trip? No — `Y_direct` is the engine's own series of `L0/Lambda_def`, `Y_formula` is independently typed; a wrong coefficient sign would break it. (2) Is the even match a hardcoded back-substitution? No — `Solve` derives `(Σ2,Σ4)` and a unique-solution guard precedes the formula comparison. (3) Symbol-domain traps? `S,β` are `nonzero` and the `3S−Σ0 ≠ 0` non-vanishing-denominator assumption is declared on both sides; no positivity over-assumption. (4) Does the script verify the paper's chi_Q? Yes, exactly. The sole defect is that the `.wl` is a line-by-line transliteration of the `.py` (F1), defeating the dual-engine independence requirement. Verdict: **findings** (one medium `mathematica_transliteration`); resolution direction (re-author vs accept) is a user-level call.

## Self-test notes

Checked: (a) variable-independence — no `diff/D` here; series is on `L0/Lambda_def`, which depends on all of `S,β,Σ0,Σ2,Σ4,Σ5`, so the `z`-series is non-degenerate. (b) The even-match `Solve` target values `1/9` and `4/81` match the canonical even fingerprint `z²/9, 4z⁴/81` (notes 84–85). (c) Trivial-case: setting `Σ5=0, β=1, S=1, Σ0=0` gives `chi_Q = 3(1)/(3) = 1` (canonical), consistent with the card's `chi_Q=1` baseline — the formula degenerates correctly. (d) Confirmed the proposed F1 re-author keeps emitted values identical and does not introduce a different constant than the paper states (paper round-trip clean).

## Value Reconciliation (pass-2 augmentation)

Deliverable-level table (symbolic results; this stage emits no numeric figures-of-merit):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Lambda_out = -3 + z²/3 + z⁴/9 + i z⁵/9` | py L27 / wl L33; out L5 | notes L13; appendix eq `app-part04-DtN-out` (L308) | MATCH |
| `L0 = -3S + Σ0` | py L33 / wl L36; out L7 | notes L50 | MATCH |
| `L2 = Sβ²/3 + Σ2` | py L34 / wl L37; out L8 | notes L51 | MATCH |
| `L4 = Sβ⁴/9 + Σ4` | py L35 / wl L38; out L9 | notes L52 | MATCH |
| `L5 = Sβ⁵/9 + Σ5` | py L36 / wl L39; out L10 | notes L53 | MATCH |
| `Ŷ_2^def = 1 - (L2/L0)z² + (L2²/L0² - L4/L0)z⁴ - i(L5/L0)z⁵` | py L43–45 / wl L42,50; out L12 | notes L66–72; appendix eq `app-part04-Ydef-expanded` (L362) | MATCH |
| `Σ2_evenmatch = -Sβ²/3 + S/3 - Σ0/9` | py L58 / wl L66; out L16 | notes L95 (`Σ2=-(3Sβ²-3S+Σ0)/9`) | MATCH |
| `Σ4_evenmatch = -Sβ⁴/9 + S/9 - Σ0/27` | py L59 / wl L67; out L17 | notes L96 (`Σ4=-(3Sβ⁴-3S+Σ0)/27`) | MATCH |
| `chi_Q = 3(Sβ⁵+9Σ5)/(3S-Σ0)` | py L70,79 / wl L72; out L20,26 | notes L104–106; appendix eq `app-part04-chiQ-general` (L373) | MATCH |

INTERNAL (scaffolding, no finding expected): `m2`, `m4` (intermediate even/odd ratios used to set the `Solve` targets), the unique-solution length guard, and all `expect_zero`/`PASS` residual flags.

Note: the `chi_Q-1` restatement (appendix eq `app-part04-chiQ-minus-one-general`, L378; notes L110–113) and the linearized `Δ_Q = 5b + a0/3 + 9a5` (appendix eq `app-part04-linear-DeltaQ-triple`, L395) are paper-side corollaries the script does not separately emit; the former is `chi_Q − 1` (algebraically covered), the latter is an appendix-only linearization outside the stage card's stated scope — neither is a missing *emitted* value, so no reconciliation finding.

`reconciliation: complete; 9 deliverable values checked, 0 misaligned`
