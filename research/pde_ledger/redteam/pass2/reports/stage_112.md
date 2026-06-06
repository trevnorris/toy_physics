---
unit_id: 112
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
  notes_stage_files: [moving_throat_pde_stage112_hybrid_robin_mixed_compensation.md]
  paper_appendix: present
---

# Audit unit 112 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_112.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage112_hybrid_robin_mixed_compensation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows 425–456 "Compensated Robin–mixed outlet"; also 27, 78, 449, 454)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage112_hybrid_robin_mixed_compensation_mathematica_audit.txt`

## What the paper claims

Stage 112 ("Exact Robin–Mixed Compensation Law") forms the hybrid outlet
`Lambda_2^hyb(z) = Lambda_2^out(z) + rho_R - sigma_W/(1 - kappa_W z^2 - i gamma_W z^5) + O(z^6)`,
expands its low-frequency coefficients (`L0 = -3 + rho_R - sigma_W`, `L2 = 1/3 - sigma_W kappa_W`,
`L4 = 1/9 - sigma_W kappa_W^2`, `L5 = 1/9 - sigma_W gamma_W`), and imposes the canonical conservative-even
fingerprint (`-L2/L0 = 1/9`, `L2^2/L0^2 - L4/L0 = 4/81`). This yields exactly two branches:
the trivial cancellation `(rho_R, kappa_W) = (sigma_W, 0)` and the nontrivial compensated branch
`(rho_R, kappa_W) = (4 sigma_W, 1/3)`. On the nontrivial branch the card/appendix/notes give the
outgoing-normalization factor `chi_Q^hyb = (1 - 9 sigma_W gamma_W)/(1 - sigma_W)` (appendix
eq:app-part04-hybrid-chiQ line 449; notes lines 76–81), and canonical preservation holds **iff
`gamma_W = 1/9`** (\stagefield{Output}-equivalent quote in the card body, line 16:
"Nontrivial compensated branch has \(\rho_R=4\sigma_W\), \(\kappa_W=1/3\), and preserves odd normalization
iff \(\gamma_W=1/9\)"; appendix line 454; notes lines 84, 106). The notes additionally record the
Stage-92 deformation data `(b, a0, a5) = (0, 3 sigma_W, -sigma_W gamma_W)` (notes line 98) and the linearized
preservation condition `a0/3 + 9 a5 = sigma_W(1 - 9 gamma_W) = 0` (notes line 102), an independent route to
`gamma_W = 1/9`; and the scaled-collapse identity `Lambda_2^hyb = (1 - sigma_W) Lambda_2^out + O(z^6)` at
`gamma_W = 1/9` (notes line 89). Distinct deliverables: (D1) two-branch even-fingerprint solution; (D2)
`chi_Q^hyb` closed form on branch B; (D3) `gamma_W = 1/9` preservation condition (the checkpoint claim);
(D4) Stage-92 `(b, a0, a5)` data + linearized condition; (D5) scaled-collapse identity at gamma=1/9.

## What the script claims to verify

Both engines build `Lambda_hyb` by series-expanding the hybrid outlet, extract `L0/L2/L4/L5`, solve the two
even-fingerprint equations, and assert the solution set is exactly `[{kappa:0, rho:sigma}, {kappa:1/3, rho:4 sigma}]`.
They then form `chi_A` (branch A) and `chi_B` (branch B) as `(-L5/L0)/(1/27)` evaluated on each solution and
assert `chi_A = 1 - 9 sigma gamma` and `chi_B = (1 - 9 sigma gamma)/(1 - sigma)`. The checkpoint preservation
claim is closed two ways: (i) `chi_B(gamma=1/9) = 1`; (ii) the factorization `numer(chi_B - 1) = sigma(1 - 9 gamma)`
plus the degenerate `sigma=0` case. The Mathematica script additionally runs an **independent** Stage-92
`(b, a0, a5)` reconstruction directly from `lambdaHyb /. solB` coefficients, builds the linearized
`presCond = a0/3 + 9 a5`, asserts `presCond = sigma(1 - 9 gamma)`, and runs `Reduce[presCond == 0 && sigma != 0, gamma]`
to FORCE `gamma_W = 1/9` on the nontrivial branch. Both close with the scaled-collapse identity at gamma=1/9.

## Paper ↔ script cross-check

| Deliverable | Script-side check | Status |
|---|---|---|
| D1 two-branch even-fingerprint solve | py L36–43 / wl L41–53 `Solve[{-l2/l0==1/9, l2^2/l0^2-l4/l0==4/81}]` asserted == the two explicit branches | match |
| D2 `chi_Q^hyb = (1-9 sigma gamma)/(1-sigma)` on branch B | py L51–53 / wl L76,82 `chi_B - (1-9 sigma gamma)/(1-sigma) == 0` | match |
| D2′ `chi_A = 1-9 sigma gamma` (branch A) | py L46–48 / wl L75,81 | match |
| D3 `gamma_W = 1/9` (CHECKPOINT) | py L54,58–61 (gamma=1/9⟹chi=1 + factorization `sigma(1-9gamma)` + sigma=0 case); wl L69–71 `Reduce[...&&sigma!=0]`→gamma=1/9 (forces the hard ⟹ direction) | match |
| D4 `(b,a0,a5)=(0,3 sigma,-sigma gamma)` + linearized cond | wl L61–68 only (SymPy omits) | match (verified on MMA side; present in notes L98,102) |
| D5 scaled-collapse identity at gamma=1/9 | py L63–65 / wl L85–87 `(lambdaHyb/.solB) - (1-sigma)lambdaOut /. gamma->1/9 == 0` | match |

`paper_alignment: aligned` — every paper deliverable has a faithful, non-tautological script-side check.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 43 | `sols == [{kappa:0,rho:sigma},{kappa:1/3,rho:4 sigma}]` | D1 | yes |
| A2 | sympy | 48 | `chi_A - (1-9 sigma gamma) == 0` | D2′ | yes |
| A3 | sympy | 53 | `chi_B - (1-9 sigma gamma)/(1-sigma) == 0` | D2 | yes |
| A4 | sympy | 54 | `chi_B(gamma=1/9) - 1 == 0` | D3 (⟸) | yes |
| A5 | sympy | 58–59 | `numer(chi_B_general-1) - sigma(1-9 gamma) == 0` | D3 (factorization) | yes |
| A6 | sympy | 60–61 | `chi_B_general(sigma=0) - 1 == 0` | D3 (degenerate branch) | yes |
| A7 | sympy | 63–65 | `(Lambda_hyb/.solB) - (1-sigma)Lambda_out @ gamma=1/9 == 0` | D5 | yes |
| B1 | mathematica | 41–42 | `Solve[...]`, `Length==2`, sorted | D1 | yes |
| B2 | mathematica | 50–53 | `(rho/.solA)-sigma`, `kappa/.solA`, `(rho/.solB)-4 sigma`, `kappa/.solB - 1/3` all 0 | D1 | yes |
| B3 | mathematica | 64–66 | `b=0`, `a0-3 sigma`, `a5+sigma gamma` on solB == 0 | D4 | yes |
| B4 | mathematica | 67–68 | `presCond - sigma(1-9 gamma) == 0` | D3/D4 (linearized) | yes |
| B5 | mathematica | 69–71 | `Reduce[presCond==0 && sigma!=0,gamma]` → `gamma=1/9` asserted | D3 (⟹, hard direction) | yes |
| B6 | mathematica | 73 | `chiBgen(sigma=0) - 1 == 0` | D3 (degenerate) | yes |
| B7 | mathematica | 81–83 | `chiA-(1-9 sigma gamma)`, `chiB-(1-9 sigma gamma)/(1-sigma)`, `chiB(gamma=1/9)-1` == 0 | D2/D2′/D3 | yes |
| B8 | mathematica | 85–87 | scaled identity on branch B @ gamma=1/9 == 0 | D5 | yes |

## Findings

None. See per-check analysis below.

## Independent-derivation check (Mathematica)

**Partial-but-substantively-independent → NOT a transliteration finding.**

The shared spine — `lambdaOut`, `Series[..]`, coefficient extraction `l0..l5`, the two-branch `Solve`, and the
`chiA/chiB = (-l5/l0)/(1/27)` construction with mirrored assertions — IS parallel choreography between the two
engines (py L27–53 ↔ wl L33–53, 75–87). Taken alone that spine would be a candidate `mathematica_transliteration`.
However, the load-bearing checkpoint direction is reached by **two genuinely different derivations**:

- SymPy (py L55–61): forms `chi_B_general` and proves preservation by factoring the **rational function**
  `chi_B - 1`: `numer(chi_B-1) = sigma(1-9 gamma)`, plus the `sigma=0` degenerate case. The "⟹ gamma=1/9 on the
  nontrivial branch" step is left to the reader of the factor `sigma(1-9 gamma)`.
- Mathematica (wl L55–71): reconstructs the **Stage-92 linearized deformation data** `(b, a0, a5)` *directly from
  the branch-B coefficients of `lambdaHyb`* (a path the `.py` never computes), forms the linearized invariant
  `presCond = a0/3 + 9 a5`, proves `presCond = sigma(1-9 gamma)`, and then explicitly **forces** `gamma=1/9` via
  `Reduce[presCond==0 && sigma!=0, gamma]` followed by `Solve` extraction and an `expectZero` on `(gamma/.rule)-1/9`.

This `Reduce[...&&sigma!=0]` branch has **no `.py` analogue** and closes the hard ⟹ direction of the biconditional
rigorously (the `.py` only structurally exhibits it). Per the prompt's own criterion ("a genuine
`Reduce[...&&sigma!=0]` branch in `.wl` with no `.py` analogue is evidence of independence"), this is the stronger
side and clears the checkpoint independence bar. Verdict: **partial, leaning independent** — independent on the
checkpoint claim, parallel on the upstream scaffolding. No finding raised.

## Engine cross-check

Both engines emit identical results:

| Quantity | SymPy output | Mathematica output |
|---|---|---|
| `Lambda_hyb` | `-I*gamma*sigma*z**5 - kappa**2*sigma*z**4 - kappa*sigma*z**2 + rho - sigma + I*z**5/9 + z**4/9 + z**2/3 - 3` (L1) | `-3 + rho - sigma + z^2/3 - kappa*sigma*z^2 + z^4/9 - kappa^2*sigma*z^4 + (I/9)*z^5 - I*gamma*sigma*z^5` (L5) — algebraically identical |
| even solutions | `[{kappa:0,rho:sigma},{kappa:1/3,rho:4 sigma}]` (L2) | `{{rho->sigma,kappa->0},{rho->4 sigma,kappa->1/3}}` (L6) |
| `chi_Q branch A` | `-9*gamma*sigma + 1` (L3) | `1 - 9*gamma*sigma` (L27) |
| `chi_Q branch B` | `(9*gamma*sigma - 1)/(sigma - 1)` (L4) | `(-1 + 9*gamma*sigma)/(-1 + sigma)` (L28) — same rational |
| scaled identity | `0` (L7) | `0` (L36) |

All in agreement. `engines_agree: true`.

## Verdict justification

`clean`. I attacked the four standard vectors plus the checkpoint and transliteration vectors and all held.
**Re-derivation of the checkpoint:** independently expanding `1/(1-kappa z^2 - i gamma z^5) = 1 + kappa z^2 +
kappa^2 z^4 + i gamma z^5 + O(z^6)` gives `L0=-3+rho-sigma`, `L5=1/9-sigma gamma`; on branch B (`rho=4 sigma,
kappa=1/3`) `L0=-3(1-sigma)`, so `chi_B = 27(-L5/L0) = (1-9 sigma gamma)/(1-sigma)`. Setting `chi_B=1` gives
`1-9 sigma gamma = 1-sigma`, i.e. `sigma(1-9 gamma)=0`, so on the nontrivial branch (`sigma≠0`) `gamma=1/9` —
matching the printed literal, which I did NOT trust but reconstructed. The script proves BOTH directions: `⟸`
(gamma=1/9 ⟹ chi_B=1) at py L54 / wl L83, and the hard `⟹` (nontrivial branch forces gamma=1/9) via the
factorization py L58–59 and the explicit `Reduce[...&&sigma!=0]` wl L69–71 — not a re-substitution of the answer.
`gamma=1/9` is FORCED by the factorization/Reduce, not asserted as a bare literal. Tautology check: the two-branch
`Solve` and the `chi_B` factorization can genuinely fail (e.g. a wrong sign in `Lambda_out` would break A3/A5);
the `chi_B_general` literal at py L55 is independently anchored to the derived `chi_B` by the proven equality at
py L53, so it is not a hardcoded-in-a-vacuum value. Symbol domains are sound (`sigma!=1`, `rho-sigma!=3` keep the
relevant denominators nonzero; SymPy all-real, no positivity traps). Outputs are fresh (both `.txt` mtimes
post-date their scripts). Engines agree exactly. Checkpoint higher bar: **cleared**.

## Self-test notes

Checked: (1) variable-independence — no spurious `diff`/`D` derivatives in either script; all assertions are
algebraic identities. (2) Trivial-case pre-check — substituting `gamma=1/9` into `chi_B=(1-9 sigma gamma)/(1-sigma)`
gives `(1-sigma)/(1-sigma)=1` literally, confirming A4/B7. (3) Round-trip/hardcoded — the only literal candidate,
`chi_B_general` (py L55), is proven equal to the derived `chi_B` at py L53, so the downstream factorization check is
not self-referential. (4) Reduce branch — `Reduce[sigma(1-9 gamma)==0 && sigma!=0, gamma]` correctly yields the
single root `gamma=1/9` (the `sigma=0` root is excluded by the `sigma!=0` constraint), matching the degenerate-case
check that any gamma works at sigma=0. No directive written (findings_count = 0).

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation (script emits ↔ docs):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| even-branch solutions `(rho,kappa) = (sigma,0)` and `(4 sigma, 1/3)` | py out L2 / wl out L6 | card L16 ("\(\rho_R=4\sigma_W\), \(\kappa_W=1/3\)"); notes L43–48; appendix L440–442 | MATCH |
| `chi_Q branch A = 1 - 9 sigma gamma` | py out L3 / wl out L27 | notes L78–81 (chi_Q^hyb numerator `1-9 sigma_W gamma_W`); appendix L449 (branch-B form) | MATCH (branch-A intermediate; numerator form matches notes/appendix) |
| `chi_Q branch B = (1 - 9 sigma gamma)/(1 - sigma)` | py out L4 / wl out L28 | notes L76–81; appendix eq L449 | MATCH |
| `gamma_W = 1/9` (preservation condition; CHECKPOINT) | py L54/58 + out L5,8 / wl L71 + out L23 | card L16; notes L84,106; appendix L454 | MATCH |
| `(b, a0, a5) = (0, 3 sigma, -sigma gamma)` (Stage-92 data) | wl out L15–20 (b,a0,a5 checks) | notes L98 | MATCH (verified MMA-only; lives in notes) |
| linearized cond `a0/3 + 9 a5 = sigma(1 - 9 gamma)` | wl out L21–22 | notes L102 | MATCH |
| scaled-collapse identity `Lambda_hyb = (1-sigma)Lambda_out @ gamma=1/9` (residual 0) | py out L7 / wl out L35–36 | notes L89; appendix L456 | MATCH |

Internal scaffolding (no prose expectation, no finding): per-coefficient `L0/L2/L4/L5` intermediate forms;
the all-zero `expectZero`/`expect_zero` residuals and `PASS`/`stage112: PASS` flags; the `degenerate sigma=0`
self-check residual.

reconciliation: complete; 7 deliverable values checked, 0 misaligned.
