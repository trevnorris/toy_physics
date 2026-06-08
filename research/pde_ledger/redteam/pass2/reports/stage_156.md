---
unit_id: 156
batch: IV.6
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md]
  paper_appendix: present
---

# Audit unit 156 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_156.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage156_renormalized_canonical_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows at lines 941–971, the "Renormalized canonical point" subsection, and the `\input{stages/stage_156}` at line 1346)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage156_renormalized_canonical_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage156_renormalized_canonical_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage156_renormalized_canonical_branch_mathematica_audit.txt`

## What the paper claims

The stage proves that the lower compensated Family-1 branch *survives* full core–mouth co-evolution, but only after a unique upward renormalization of the mouth traction. The notes state the self-consistent restoration condition `g_fp(Sigma0) = g_*` has a unique positive root on the analyzed interval, and box the deliverables: `Sigma0_can ≈ 4.651033550168867`, `T_hat_m,can ≈ 1.4467083664567613` (via the closure `Sigma0 = (20/9) T_hat^2`), the restored fixed point `g_can = g_* ≈ 0.758035078944663` with `R_can = 1/4`, and the selected `S_can ≈ 0.6703621156734617`, `Pi_can ≈ 3.871564377479002`. The card's quoted block states "Exact compensation is restored at `Sigma_0^can ≈ 4.65103`, `T_hat_m,can ≈ 1.44671`," and its Checks require deviations about the renormalized point and `R = 1/4`. The notes additionally report relative shifts vs the original canonical point (`Sigma0` ratio − 1 ≈ 1.5754, `Pi` ratio − 1 ≈ 1.5659, `T_hat` ratio − 1 ≈ 0.6048) and characterize the result as a ~60.48% traction increase. The status tag is `\StatusNumerical{}` — explicitly a numerical placement inside the exact map, not a closed-form PDE theorem.

## What the script claims to verify

Both scripts discretize the exact self-consistent Family-1 map on an N=2401 grid (trapezoid weights), define the transport operators `Ts`/`Tq`, the moments `g`/`S`, the loading ratio `R = (g − rF1)^2 / (1 + rF1^2)`, and the update `sig_new = normalize(exp(−Sigma0·(Ts − R·Tq)))`. They solve the fixed point by Picard iteration to 1e-13, bracket-and-bisect the root of `g_fp(Sigma0) = g_*` on [3,6], then compute `Sigma0_can`, `g_can`, `S_can`, `R_can`, `Pi_can = Sigma0_can·(1 − R_can·S_can)`, and `T_hat_can = sqrt(9·Sigma0_can/20)`. The `.py` asserts (via `raise AssertionError`) that `g_can = g_*`, `R_can = 1/4`, and that `Sigma0_can`/`S_can`/`Pi_can`/`T_hat_can` match the boxed notes values within 1e-10/1e-11. The `.wl` additionally scans `g_fp` on a 0.5-step grid and asserts (`expectTrue`) strict monotone increase, then runs the same six numeric checks (`expectApprox`).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Unique root `g_fp(Sigma0)=g_*` (monotone map) | `.wl` scan + `expectTrue["monotone increase on scan grid"]` (wl:124); bisection in both | match |
| `Sigma0_can ≈ 4.651033550168867` | py:127–128 assert; wl:147 `expectApprox` | match |
| `T_hat_m,can ≈ 1.4467083664567613` via `Sigma0=(20/9)T_hat^2` | `T_hat = sqrt(9·Sigma0/20)` (py:109, wl:133) + check py:133, wl:152 | match |
| `g_can = g_* ≈ 0.758035078944663` | py:123, wl:148 | match |
| `R_can = 1/4` | py:125–126, wl:150 | match |
| `S_can ≈ 0.6703621156734617` | py:129–130, wl:149 | match |
| `Pi_can ≈ 3.871564377479002` | `Pi=Sigma0·(1−R·S)` (py:108, wl:132) + check py:131, wl:151 | match |
| Relative shifts (1.5754 / 1.5659 / 0.6048) | printed py:119–121, wl:143–145 | match (printed, not asserted) |

`paper_alignment: aligned`. Every boxed/quoted deliverable has a corresponding script-side check. The relative-shift values are printed (not asserted) but agree with the notes to ~13 digits, so they are reported, not load-bearing.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 123–124 | `abs(g_can − g_star) > 1e-10 → raise` | restored fixed point `g_can=g_*` | yes |
| A2 | sympy | 125–126 | `abs(R_can − 0.25) > 1e-10 → raise` | `R_can = 1/4` | yes |
| A3 | sympy | 127–128 | `abs(Sigma0_can − 4.651…867) > 1e-10 → raise` | `Sigma0_can` deliverable | yes |
| A4 | sympy | 129–130 | `abs(S_can − 0.670…617) > 1e-11 → raise` | `S_can` deliverable | yes |
| A5 | sympy | 131–132 | `abs(Pi_can − 3.871…002) > 1e-10 → raise` | `Pi_can` deliverable | yes |
| A6 | sympy | 133–134 | `abs(T_hat_can − 1.446…613) > 1e-10 → raise` | `T_hat_can` deliverable | yes |
| A7 | mathematica | 124 | `expectTrue["monotone increase", Min[scanDiffs] > 0]` | uniqueness of root (monotone map) | yes |
| A8 | mathematica | 147 | `expectApprox sigma0Can ~ 4.651…867` | `Sigma0_can` | yes |
| A9 | mathematica | 148 | `expectApprox gCan ~ 0.758…663` | `g_can=g_*` | yes |
| A10 | mathematica | 149 | `expectApprox sCan ~ 0.670…617` | `S_can` | yes |
| A11 | mathematica | 150 | `expectApprox rCan ~ 0.25` | `R_can=1/4` | yes |
| A12 | mathematica | 151 | `expectApprox piCan ~ 3.871…002` | `Pi_can` | yes |
| A13 | mathematica | 152 | `expectApprox tHatCan ~ 1.446…613` | `T_hat_can` | yes |

None of these is tautological: the checked quantities are produced by Picard iteration + bisection, then compared to the boxed targets; a wrong discretization or wrong root would change the computed value and trip the check. The comparison literals are the notes' published values used as regression anchors, which is legitimate for a numerical stage because the value is genuinely re-derived in-script.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration of the `.py`, despite structural parallelism (expected for a shared numerical scheme). Evidence of genuine independence:
- It is a numerical fixed-point solve, not symbolic algebra being re-typed; each engine evaluates the discretized map in its own native floating-point arithmetic. The committed outputs diverge in the last 1–2 digits (`Sigma0_can = 4.651033550168867` py vs `4.651033550168874` wl; `Pi_can = 3.871564377479002` vs `3.871564377479008`), the signature of two independent numerical runs rather than echoed algebra.
- The `.wl` adds a numerical-stabilization shift absent from the `.py`: `phShift = ph - Min[ph]; normalize[Exp[-phShift]]` (wl:70–72) vs the `.py`'s raw `np.exp(-ph)` (py:64). Since `normalize` divides out the constant factor, this is mathematically invariant but a deliberately different implementation.
- The `.wl` does **more** than the `.py`: an explicit `g_fp` monotonicity scan on [3,6] with `expectTrue["monotone increase on scan grid"]` (wl:119–124), which directly substantiates the notes' uniqueness claim and is absent from the `.py`. A line-by-line port would not add checks the source lacks.
- The `.wl` uses memoization (`fixedPointAt[s] = …`, `gFp[s] = …`, wl:89–90), an independent implementation choice.

## Engine cross-check

Both engines agree on every deliverable to ~13–14 significant digits; all differences are float noise:

| value | sympy out | mathematica out |
|---|---|---|
| Sigma0_can | 4.651033550168867 | 4.651033550168874 |
| g_can | 0.758035078944663 | 0.7580350789446629 |
| S_can | 0.6703621156734617 | 0.6703621156734617 |
| R_can | 0.2500000000000005 | 0.2500000000000005 |
| Pi_can | 3.871564377479002 | 3.871564377479008 |
| T_hat_can | 1.4467083664567613 | 1.4467083664567622 |

The `.wl`'s `expectApprox` diffs (output lines 19–30) are all ≤ 7.1e-15, comfortably inside the 1e-10/1e-11 tolerances. No `engine_disagreement`.

## Verdict justification

Clean. I attacked: (1) tautology — the checks compare iteration-derived values to notes anchors, not a value to itself; they can and would fail on a wrong discretization. (2) hardcoded result — the boxed numbers are regression targets, but the quantities are genuinely recomputed via the fixed-point solve, so no number is asserted against itself. (3) transliteration — the `.wl` diverges in stabilization, adds a monotonicity check, memoizes, and produces last-digit-distinct output, so it is an independent numerical confirmation, not a port. (4) the warned stale-constant patterns (`168π²`/`100π²`/`4107`) do not appear in this stage. (5) symbol/assumption errors — none; this is real-valued numpy/Mathematica numerics with no domain assumptions to violate. (6) the susceptibility closure `T_hat = sqrt(9·Sigma0/20)` matches the notes `Sigma0=(20/9)T_hat^2` exactly. Outputs are newer than scripts (no `stale_output`). Paper card, notes, and appendix all carry the same deliverables and the same `\StatusNumerical{}` status, and every script check traces to a paper-side deliverable. The only doc divergence is sub-tolerance last-digit roundoff between the appendix renderings and the notes/script (see reconciliation), which is not a value mismatch.

## Self-test notes

Checked: variable independence (no `diff`/`D` in this numerical script — N/A); parity/symmetry of integrals (the `Ts`/`Tq` quadratures are cumulative sums over [0,1], not symmetric-domain vanishing claims — N/A); trivial-case (the `R_can=1/4` and `g_can=g_*` targets are the defining compensation conditions, and the bisection drives `g_fp→g_*` by construction, so the checks reduce correctly). No directive is written (zero findings).

## Value Reconciliation (pass-2 augmentation)

Authoritative record: script source + committed `.txt` outputs (both fresh; outputs newer than scripts). Doc carriers: stage card `stage_156.tex`, notes `..._renormalized_canonical_branch.md`, and the part-04 appendix "Renormalized canonical point" subsection (lines 941–971).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Sigma0_can ≈ 4.651033550168867` | py:111 / wl:135 (sympy out L5 `4.651033550168867`; wl out L8 `…874`) | notes L31 `…867`; card L16 `≈4.65103`; appendix L956 `…876` | MATCH (notes/script exact; card to 6 sf; appendix last-digit roundoff ~1e-14) |
| `g_can = g_* ≈ 0.758035078944663` | py:112 / wl:136 (out L6/L9) | notes L60 `…663`; appendix L963 `g_can=g_*` | MATCH |
| `S_can ≈ 0.6703621156734617` | py:113 / wl:137 (out L7/L10) | notes L69 `…4617`; appendix L967 `…673462` | MATCH (appendix last-digit roundoff ~1e-16) |
| `R_can = 1/4` | py:114 / wl:138 (out L8/L11 `0.25000…05`) | notes L61 `R_can=1/4`; card Checks; appendix L965 | MATCH |
| `Pi_can ≈ 3.871564377479002` | py:115 / wl:139 (out L9/L12) | notes L72 `…002`; appendix L969 `…009` | MATCH (appendix last-digit roundoff ~7e-15) |
| `T_hat_can ≈ 1.4467083664567613` | py:116 / wl:140 (out L10/L13) | notes L43 `…613`; card L16 `≈1.44671`; appendix L958 `…62` | MATCH (notes/script exact; appendix last-digit roundoff ~1e-15) |
| `Sigma0 ratio − 1 ≈ 1.5754070949…` | py:119 / wl:143 (out L13/L17) | notes L114 `…223031` | MATCH (~1e-14) |
| `Pi ratio − 1 ≈ 1.5659389234…` | py:120 / wl:144 (out L14/L18) | notes L121 `…213572` | MATCH (~1e-14) |
| `T_hat ratio − 1 ≈ 0.6048074946…` | py:121 / wl:145 (out L15/L19) | notes L128 `…616844` | MATCH (~1e-15); notes L135 "60.48%" MATCH |

Internal scaffolding (accounted for, no finding expected in prose): input constants `rF1=1.77799353547498`, `g_star=0.758035078944663`, `Pi_star=1.50882951349316`, `Sigma0_star=1.80594111095636`, `T_hat_star=0.901484054174204` (carried-forward upstream "original canonical point" values, also reflected in notes §3 `(Pi_*,T_hat_*)≈(1.5088,0.9015)`); grid params N=2401, dx, kappa=π/2; iteration tolerance 1e-13; bisection iters 55/60; pass/fail flags and `expectApprox` diffs.

Note on the notes' "224.59% mouth-bias increase" (notes L136): this percentage is a notes-prose figure that does **not** correspond to any value the scripts emit (the script's Pi ratio − 1 is 156.59%, not 224.59%), so it is outside script→doc reconciliation scope; it is a notes-internal prose number, not a script deliverable, and is flagged here only as an observation (no script-side counterpart to reconcile against).

reconciliation: complete; 9 deliverable values checked, 0 misaligned (appendix last-digit divergences are sub-tolerance numerical roundoff between independent captures, not value mismatches).
