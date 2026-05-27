---
unit_id: 111
batch: IV.2
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage111_mixed_sidechannel_pole.md]
  paper_appendix: present
---

# Audit unit 111 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_111.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage111_mixed_sidechannel_pole.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (the appendix only `\input`s `stages/stage_111`; no separate row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage111_mixed_sidechannel_pole_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage111_mixed_sidechannel_pole_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.txt`

## What the paper claims

The stage card's quoted Output is: "A standalone hidden mixed pole cannot preserve the canonical even branch unless it is absent." The notes flesh this out: the canonical outgoing branch is augmented by a single isotropic Schur-complement-style hidden mixed pole `sigma_W / (1 - kappa_W z^2 - i gamma_W z^5)`, expanded to order `z^5`. Matching the even `l=2` fingerprint at `z^2` forces `kappa_W = -1/9` (boxed result), which is already incompatible with the assumed `kappa_W > 0`. Formally inserting that value into the `z^4` match then forces `sigma_W = 0` (boxed). The notes also state the outgoing-normalization shift `chi_Q^mix = 3(1 - 9 sigma_W gamma_W) / (3 + sigma_W)` and its linearization `chi_Q^mix = 1 - sigma_W(1/3 + 9 gamma_W) + O(sigma_W^2)`, which encode the branch-selection triple `(b, a_0, a_5) = (0, -sigma_W, -sigma_W gamma_W)`.

## What the script claims to verify

Both scripts construct `Lambda_mix(z) = Lambda_out(z) - sigma/(1 - kappa z^2 - i gamma z^5)`, expand to `z^5`, extract `L0, L2, L4, L5` (with `L5` defined as `coeff(z^5)/i`), solve `-L2/L0 = 1/9` for `kappa` and `(L2^2/L0^2 - L4/L0)|_{kappa=kappa_match} = 4/81` for `sigma`, then form `chi_mix = 27 (-L5/L0)` and its linearization in `sigma`. The four assertions verify (i) `kappa_match = -1/9`, (ii) `sigma_match = 0`, (iii) `chi_mix = 3(1 - 9 sigma gamma)/(3 + sigma)`, (iv) `chi_mix_lin = 1 - sigma(1/3 + 9 gamma)`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `kappa_W = -1/9` (z^2 match no-go) | `assert kappa_match + 1/9 == 0` (sympy L32; math L60) | match |
| `sigma_W = 0` (z^4 match forces absence) | `assert sigma_match == 0` (sympy L33; math L61) | match |
| `chi_Q^mix = 3(1 - 9 sigma gamma)/(3 + sigma)` | `assert chi_mix - 3(1 - 9 sigma gamma)/(3+sigma) == 0` (sympy L34; math L62) | match |
| Linearization `chi_Q^mix = 1 - sigma(1/3 + 9 gamma)` | `assert chi_mix_lin - (1 - sigma(1/3 + 9 gamma)) == 0` (sympy L35; math L63) | match |
| Triple `(b, a_0, a_5) = (0, -sigma_W, -sigma_W gamma_W)` | implicit in linearization check (no separate extraction) | partial (re-parameterization only; covered) |

Dominant pattern: aligned. `paper_alignment = aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 32 | `assert simplify(kappa_match + 1/9) == 0` | kappa_W = -1/9 no-go | yes |
| A2 | sympy | 33 | `assert simplify(sigma_match) == 0` | sigma_W = 0 no-go | yes |
| A3 | sympy | 34 | `assert simplify(chi_mix - 3(1-9 sigma gamma)/(3+sigma)) == 0` | chi_Q^mix closed form | yes |
| A4 | sympy | 35 | `assert simplify(chi_mix_lin - (1 - sigma(1/3+9 gamma))) == 0` | linearization | yes |
| A5 | mathematica | 60 | `expectZero["kappa_match + 1/9", kappaMatch + 1/9]` | kappa_W = -1/9 | yes |
| A6 | mathematica | 61 | `expectZero["sigma_match", sigmaMatch]` | sigma_W = 0 | yes |
| A7 | mathematica | 62 | `expectZero[chi closed form]` | chi_Q^mix closed form | yes |
| A8 | mathematica | 63 | `expectZero[chi linearized]` | linearization | yes |

All assertions are anchored to paper-side deliverables and are non-tautological (each is the residual of an independently-computed expression against the paper's quoted result).

## Findings

### F1 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage111_mixed_sidechannel_pole_sympy_audit.py:1-37`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.wl:26-66`

**What's wrong:**
The `.wl` script is a line-by-line port of the `.py` script: same variable choreography, same intermediate quantities, same algorithm order, same four assertions in the same order. Compare:

- SymPy L8: `Lambda_out = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + I*z**5/sp.Integer(9)`
- Math L33: `lambdaOut = -3 + z^2/3 + z^4/9 + I*z^5/9`

- SymPy L9: `Lambda_mix = sp.expand(sp.series(Lambda_out - sigma/(1 - kappa*z**2 - I*gamma*z**5), z, 0, 6).removeO())`
- Math L34: `lambdaMix = Expand[Normal[Series[lambdaOut - sigma/(1 - kappa*z^2 - I*gamma*z^5), {z, 0, 5}]]]`

- SymPy L22-25 vs Math L47-53: identical pipeline (`kappa_match` from `-L2/L0 = 1/9`, then `sigma_match` from `L2^2/L0^2 - L4/L0 = 4/81` with kappa substituted, then `chi_mix = (-L5/L0)/(1/27)`, then `chi_mix_lin = series(chi_mix, sigma, 0, 2)`).

The Mathematica engine does not derive the four paper-side identities by an independent path; it just re-runs SymPy's coefficient-matching pipeline in Wolfram syntax.

**Why this matters:**
The second-engine policy exists so that a bug in the choreography of one script (wrong series order, wrong coefficient extraction, wrong substitution order) does not silently propagate. Two engines walking the same algorithmic steps with the same variables only catch CAS-engine-specific bugs (simplification quirks), not author-introduced choreography errors. The math here happens to be correct, so this is informational, not blocking.

**Required change:**
Restructure the Mathematica script so it derives at least one of the four assertions by an independent path. Suggested independent re-derivations (Codex picks one; do not rewrite all four):
- Compute `chi_Q^mix` directly from the closed-form geometric series for the pole (without extracting `L5` separately): `chi_Q^mix = -27 Im[Series[sigma/(1 - kappa z^2 - I gamma z^5), {z, 0, 5}] /. {kappa -> -1/9}]/Coefficient[...,z,0]` or similar — produce the boxed formula without going through `L5`.
- Or: verify the `z^2` and `z^4` no-go by direct substitution into the boxed `Lambda_2^mix` expansion stated in the notes, rather than by `sp.series` + `Coefficient`.

The four assertions must remain at least as strong as they are now. Do not remove any existing `expectZero` call.

**Verification:**
The verifier inspects the Mathematica script and confirms that at least one of `chiMix`, `kappaMatch`, or `sigmaMatch` is computed by a route that does not literally mirror the SymPy pipeline (different intermediate names is not sufficient; the algebraic path must differ). All four `expectZero` lines still exit 0.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration of the `.py`. Both build `lambdaMix` by subtracting the same pole, both extract the same four coefficients, both solve the same equations in the same order, both form `chi_mix = (-L5/L0)/(1/27)`. The Mathematica `Solve[..., Reals]` and `FullSimplify` calls are the only Mathematica-specific touches; the algorithmic spine is identical. See F1 for the side-by-side quotes.

## Engine cross-check

Both engine outputs agree:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `L0` | `-sigma - 3` | `-3 - sigma` |
| `L2` | `-kappa*sigma + 1/3` | `1/3 - kappa*sigma` |
| `L4` | `-kappa**2*sigma + 1/9` | `1/9 - kappa^2*sigma` |
| `L5` | `-gamma*sigma + 1/9` | `1/9 - gamma*sigma` |
| `kappa_match` | `-1/9` | `-1/9` |
| `sigma_match` | `0` | `0` |
| `chi_Q^mix` | `3*(-9*gamma*sigma + 1)/(sigma + 3)` | `(3 - 27*gamma*sigma)/(3 + sigma)` |
| `chi_Q^mix linearized` | `-9*gamma*sigma - sigma/3 + 1` | `1 - sigma/3 - 9*gamma*sigma` |

All eight pairs are algebraically identical. No `engine_disagreement` finding.

## Verdict justification

The four script-side assertions cleanly verify all four paper-side identities quoted in the notes (the two no-go conditions `kappa_W = -1/9` and `sigma_W = 0`, the closed-form `chi_Q^mix`, and its linearization). The card's body sentence "a standalone hidden mixed pole cannot preserve the canonical even branch unless it is absent" is captured by the conjunction of the two no-go assertions. The Inputs/Outputs are consistent with the script's symbol setup. I attacked the assertions for tautology (each is the residual of an independently-computed expression against a closed-form target — not tautological), for hidden assumptions (`sigma`, `kappa`, `gamma` are declared real only, with no positivity that would mask the kappa=-1/9 result), and for symmetry (the L5 extraction via `coeff(z,5)/I` is correct because the i-coefficient comes from the same imaginary piece in both `Lambda_out` and the pole expansion). The remaining concern is `mathematica_transliteration`: the second engine does not exercise independent algebraic choreography. Verdict: `findings`, one low-severity informational item.

## Self-test notes

I checked: (a) the series order — SymPy `series(..., 0, 6)` and Math `Series[..., {z, 0, 5}]` both retain through `z^5` and discard `O(z^6)`, consistent; (b) the `L5 = coeff/I` definition — extracting the imaginary part of the `z^5` coefficient as a real symbol — both engines do this identically; (c) absence of positivity assumptions on `kappa`, `sigma`, `gamma` that could have masked the boxed no-go conclusions; (d) the boxed `chi_Q^mix = 3(1 - 9 sigma gamma)/(3 + sigma)` derived in the script matches the notes' boxed expression; (e) the algebraic round-trip — at `kappa = -1/9` and `sigma = 0`, both `kappa_match + 1/9` and `sigma_match` reduce to the literal `0` displayed in the saved transcripts. The directive's "independent re-derivation" instruction was checked against the paper: any acceptable rewrite must keep the four assertions and produce the same closed forms, so the fix cannot introduce a new `paper_misalignment`.
