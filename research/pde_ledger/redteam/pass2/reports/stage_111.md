---
unit_id: 111
batch: IV.2
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
  notes_stage_files: [moving_throat_pde_stage111_mixed_sidechannel_pole.md]
  paper_appendix: present
---

# Audit unit 111 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_111.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage111_mixed_sidechannel_pole.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (subsec "Robin and mixed-pole tests", eq:app-part04-mixed-pole-outlet at L414-423; the appendix `\input`s `stages/stage_111` at L1256, no separate tabular row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage111_mixed_sidechannel_pole_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage111_mixed_sidechannel_pole_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.txt`

## What the paper claims

The card's quoted body claim is: "A standalone hidden mixed pole cannot preserve the canonical even branch unless it is absent." The notes flesh this out: take the canonical outgoing branch `Lambda_2^out(z) = -3 + z^2/3 + z^4/9 + i z^5/9` and add one isotropic Schur-complement-style hidden mixed pole `- sigma_W/(1 - kappa_W z^2 - i gamma_W z^5)`, expand to `O(z^6)`, giving coefficients `L0 = -(3+sigma_W)`, `L2 = 1/3 - sigma_W kappa_W`, `L4 = 1/9 - sigma_W kappa_W^2`, `L5 = 1/9 - sigma_W gamma_W`. The two canonical-even fingerprint conditions `-L2/L0 = 1/9` and `L2^2/L0^2 - L4/L0 = 4/81` force the two boxed no-go results `kappa_W = -1/9` (incompatible with passive `kappa_W>0`) and then `sigma_W = 0` (the pole must vanish). The notes additionally box the outgoing-normalization shift `chi_Q^mix = 3(1 - 9 sigma_W gamma_W)/(3 + sigma_W)`, its linearization `chi_Q^mix = 1 - sigma_W(1/3 + 9 gamma_W) + O(sigma_W^2)`, and the equivalent branch-selection triple `(b, a_0, a_5) = (0, -sigma_W, -sigma_W gamma_W)`. The appendix (L414-423) restates the pole model and the `kappa_W=-1/9 ⟹ sigma_W=0` no-go verbatim.

Five distinct deliverables: (D1) `kappa_W = -1/9`; (D2) `sigma_W = 0`; (D3) `chi_Q^mix = 3(1-9 sigma gamma)/(3+sigma)`; (D4) linearization `1 - sigma(1/3+9 gamma)`; (D5) the coefficient expansion `L0..L5` (carries D1–D4).

## What the script claims to verify

Both engines build `Lambda_mix(z) = Lambda_out(z) - sigma/(1 - kappa z^2 - i gamma z^5)`, series-expand through `z^5`, and extract `L0, L2, L4, L5` (with `L5 := coeff(z^5)/I`). They then solve `-L2/L0 = 1/9` for `kappa` (→ `-1/9`), substitute that into `L2^2/L0^2 - L4/L0 = 4/81` and solve for `sigma` (→ `0`), and form `chi_mix = (-L5/L0)/(1/27)` plus its `sigma`-linearization. The four load-bearing assertions check `kappa_match = -1/9`, `sigma_match = 0`, `chi_mix = 3(1-9 sigma gamma)/(3+sigma)`, and `chi_mix_lin = 1 - sigma(1/3+9 gamma)`. The Mathematica script adds a FIFTH assertion (`chi_Q^mix routes agree`, L75) cross-checking `chiMix` against an independently-routed `chiMixAlt` computed directly from the pole's geometric series (L54-62), bypassing the `L0/L5` extraction.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| D1 `kappa_W = -1/9` (z^2 even match) | `assert kappa_match + 1/9 == 0` (py L32; wl L71) | match |
| D2 `sigma_W = 0` (z^4 even match) | `assert sigma_match == 0` (py L33; wl L72) | match |
| D3 `chi_Q^mix = 3(1-9 sigma gamma)/(3+sigma)` | `assert chi_mix - 3(1-9 sigma gamma)/(3+sigma) == 0` (py L34; wl L73) | match |
| D4 linearization `1 - sigma(1/3+9 gamma)` | `assert chi_mix_lin - (1 - sigma(1/3+9 gamma)) == 0` (py L35; wl L74) | match |
| D5 `L0..L5` expansion | printed L0/L2/L4/L5 (py L16-20; wl L41-45) feeding D1–D4 | match |

Dominant pattern: aligned. `paper_alignment = aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 32 | `assert simplify(kappa_match + 1/9) == 0` | D1 | yes |
| A2 | sympy | 33 | `assert simplify(sigma_match) == 0` | D2 | yes |
| A3 | sympy | 34 | `assert simplify(chi_mix - 3(1-9 sigma gamma)/(3+sigma)) == 0` | D3 | yes |
| A4 | sympy | 35 | `assert simplify(chi_mix_lin - (1 - sigma(1/3+9 gamma))) == 0` | D4 | yes |
| A5 | mathematica | 71 | `expectZero["kappa_match + 1/9", kappaMatch + 1/9]` | D1 | yes |
| A6 | mathematica | 72 | `expectZero["sigma_match", sigmaMatch]` | D2 | yes |
| A7 | mathematica | 73 | `expectZero[chi closed form]` | D3 | yes |
| A8 | mathematica | 74 | `expectZero[chi linearized]` | D4 | yes |
| A9 | mathematica | 75 | `expectZero["chi_Q^mix routes agree", chiMix - chiMixAlt]` | D3 (independent cross-route) | yes |

All nine assertions are anchored to paper-side deliverables and non-tautological — each is the residual of a quantity computed by the coefficient-matching pipeline against an independently-stated closed-form target (the paper's boxed result), not a definition compared to itself. The `1/27` scale in `chi_mix` is the canonical normalization (`chi_Q=1` on the unperturbed branch: with `sigma=0`, `-L5/L0 = (1/9)/3 = 1/27`), so it is anchored to the appendix's `chi_Q=1` canonical claim, not a magic constant.

## Findings

None. (First-pass `mathematica_transliteration` F1 has been remediated — see Independent-derivation check.)

## Independent-derivation check (Mathematica)

The first-pass audit (`redteam/reports/stage_111.md`, F1, low severity) flagged the `.wl` as a line-by-line port of the `.py`: identical `Lambda_out` literal, identical `Series` on the same subtracted-pole expression, identical `L0/L2/L4/L5` extraction, identical `Solve` order, identical four assertions. That directive explicitly authorized re-deriving exactly ONE of the four identities ("Codex picks one; do not rewrite all four"), suggesting the `chi_Q` route specifically.

The current `.wl` carries that remediation (L54-62):
- `poleSeries = Series[sigma/(1 - kappa*z^2 - I*gamma*z^5), {z, 0, 5}]` — expands the POLE alone, not `lambdaMix`.
- `imagPart5 = Coefficient[Normal[poleSeries], z, 5]/I` — yields `sigma*gamma`, the pole's own `z^5` imaginary coefficient.
- `chiMixAlt = FullSimplify[27*(1/9 - imagPart5)/(3 + sigma), ...]` — reconstructs `chi_Q^mix` from the canonical `Lambda_out` z^5 coefficient `1/9` and `-L0 = 3+sigma` directly, WITHOUT going through the extracted `L5`/`L0` of `lambdaMix`. Output line 10 confirms `(3 - 27*gamma*sigma)/(3 + sigma)`, byte-equal to the primary `chiMix` route (output line 13), and L75 asserts the two routes agree.

This is a genuinely different algebraic path for D3 (different intermediate object, different coefficient source, no reuse of the `lambdaMix` pipeline), satisfying the first-pass directive's "pick one" scope. The TWO no-go core deliverables D1/D2 (`kappaMatch`, `sigmaMatch`) remain a transliteration of the `.py` (same `Solve[-l2/l0 == 1/9, kappa]` / `Solve[(l2^2/l0^2 - l4/l0) == 4/81, sigma]` choreography). Per the sanctioned remediation scope this partial independence is the intended state, so I classify the `.wl` as **partial** independence and do NOT re-raise `mathematica_transliteration`. The pole-check idiom is handled cleanly: there is no `Limit`-near-pole; both engines use `Series`/coefficient extraction, and `chiMixAlt` avoids any `=!= Infinity` style test.

## Engine cross-check

Both engine outputs agree at the level they claim:

| Quantity | SymPy (`...sympy_audit.txt`) | Mathematica (`...mathematica_audit.txt`) |
|---|---|---|
| `L0` | `-sigma - 3` | `-3 - sigma` |
| `L2` | `-kappa*sigma + 1/3` | `1/3 - kappa*sigma` |
| `L4` | `-kappa**2*sigma + 1/9` | `1/9 - kappa^2*sigma` |
| `L5` | `-gamma*sigma + 1/9` | `1/9 - gamma*sigma` |
| `kappa_match` | `-1/9` | `-1/9` |
| `sigma_match` | `0` | `0` |
| `chi_Q^mix` | `3*(-9*gamma*sigma + 1)/(sigma + 3)` | `(3 - 27*gamma*sigma)/(3 + sigma)` |
| `chi_Q^mix linearized` | `-9*gamma*sigma - sigma/3 + 1` | `1 - sigma/3 - 9*gamma*sigma` |
| `chi_Q^mix independent route` | (n/a — sympy has no alt route) | `(3 - 27*gamma*sigma)/(3 + sigma)` |

All shared pairs are algebraically identical (`3*(-9 g s + 1)/(s+3) = (3-27 g s)/(3+s)`). The Mathematica alt-route value equals its primary value, and its L75 `chi_Q^mix routes agree = 0` passes. No `engine_disagreement`.

## Verdict justification

The nine script-side assertions cleanly verify all five paper-side deliverables. The card's body sentence ("a standalone hidden mixed pole cannot preserve the canonical even branch unless it is absent") is captured by the conjunction of A1/A5 (`kappa_W=-1/9`, incompatible with passivity) and A2/A6 (`sigma_W=0`). I attacked: (a) tautology — each assertion is a residual of a pipeline-computed quantity against an independently-stated boxed target, so none can pass by construction; (b) the `1/27` scale and `4/81`/`1/9` RHS constants — all three are the canonical `Lambda_out` values (verified by substituting `sigma=0`), so they are anchored, not hardcoded answers; (c) symbol assumptions — `z, sigma, kappa, gamma` are declared `real` only with no positivity that could mask the `kappa=-1/9` no-go (the `.wl` adds only `sigma != -3` to keep `L0` nonzero, which is physically required, not a hidden simplification crutch); (d) the `L5 = coeff(z^5)/I` extraction — correct, the imaginary `z^5` term is the only imaginary piece in both `Lambda_out` and the pole; (e) transliteration — the first-pass finding was remediated within its authorized "pick one" scope via the independent `chiMixAlt` route, leaving partial (sanctioned) independence. Verdict: `clean`.

## Self-test notes

Variable independence: no `sp.diff`/`D` derivatives appear; the `Series` expansions are in `z` (and `chi_mix_lin` in `sigma`), and both `Lambda_mix` and `chi_mix` genuinely depend on those variables, so no identically-zero-derivative trap. Parity/symmetry: not applicable (finite Taylor coefficient matching, no symmetric-domain integral). Trivial-case pre-check: at `sigma=0` the canonical anchors reduce exactly (`-L2/L0 = 1/9`, `L2^2/L0^2 - L4/L0 = 4/81`, `chi_mix = 1`), confirming the RHS constants are the canonical values and the no-go solves return `kappa=-1/9`, `sigma=0` as shown in both transcripts. No directive written (zero findings).

## Value Reconciliation (pass-2 augmentation)

Deliverable values emitted by the scripts (from source + committed `.txt` outputs):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `L0 = -(3+sigma)` | py L17 / wl L42; out py L2, wl L6 | notes L30 region (stage111 expansion `-(3+sigma_W)`, L28-29) | MATCH |
| `L2 = 1/3 - kappa*sigma` | py L18 / wl L43; out py L3, wl L7 | notes L31 (`1/3 - sigma_W kappa_W`) | MATCH |
| `L4 = 1/9 - kappa^2*sigma` | py L19 / wl L44; out py L4, wl L8 | notes L33 (`1/9 - sigma_W kappa_W^2`) | MATCH |
| `L5 = 1/9 - gamma*sigma` | py L20 / wl L45; out py L5, wl L9 | notes L35 (`1/9 - sigma_W gamma_W`) | MATCH |
| `kappa_W = -1/9` | py L22/L27 / wl L47; out py L6, wl L11 | notes L46 (boxed); appendix L423 | MATCH |
| `sigma_W = 0` | py L23/L28 / wl L48; out py L7, wl L12 | notes L53 (boxed); appendix L423 | MATCH |
| `chi_Q^mix = 3(1-9 sigma gamma)/(3+sigma)` | py L24/L29 / wl L52; out py L8, wl L13 | notes L60-66 (boxed) | MATCH |
| `chi_Q^mix linearized = 1 - sigma(1/3+9 gamma)` | py L25/L30 / wl L64; out py L9, wl L14 | notes L67-74 (boxed series) | MATCH |
| `Lambda_mix(z)` full expansion | py L16 / wl L41; out py L1, wl L5 | notes L26-37 (boxed expansion) | MATCH |

INTERNAL (scaffolding, no finding expected): `kappaMatch+1/9` residual, `sigmaMatch` residual, `chiMix - chiMixAlt` route-agreement residual, `imagPart5` (= sigma*gamma, intermediate), `poleSeries` (intermediate), PASS/FAIL flags, `chi_Q^mix (independent route)` print (equals the boxed `chi_Q^mix`, already reconciled above).

`reconciliation: complete; 9 values checked, 0 misaligned`
