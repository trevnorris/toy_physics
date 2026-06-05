---
unit_id: 079
batch: III.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage079_family1_quadrupole_pe_map.md]
  paper_appendix: present
---

# Audit unit 079 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_079.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage079_family1_quadrupole_pe_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 136)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage079_family1_quadrupole_pe_map_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage079_family1_quadrupole_pe_map_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.txt`

## What the paper claims

Stage 079 converts a required quadrupole support-ratio demand `zeta_req` into the transport bias required on the explicit Family-1 branch. It fixes the branch factor `A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2)` with `kappa_F1 = 12321/5` (carried from Stages 073-075) and `y_F1` the Robin root of `y tan y = 37`, `0 < y < pi/2`, giving the boxed value `A_F1 ≈ 1.00005192880220`. The `\stagefield{Output}` line names two boxed deliverables: the demand map `zeta_F1(Pe) = A_F1 Omega_Pe^2` (eq. app-stage079-zetaF1) and the hard constructive ceiling `zeta_max^(F1) = A_F1 pi^2/4 ≈ 2.46752922945601` (eq. app-stage079-ceiling). The notes add the supporting endpoint facts `zeta_F1(0) = A_F1`, `lim_{Pe→∞} zeta_F1 = A_F1 pi^2/4`, the carried Stage-056 closed form for `Omega_Pe` and its limits (`Omega(0)=1`, `Omega(∞)=pi/2`), and the small-demand expansion `zeta_F1(Pe) = A_F1[1 + ((4-pi)/pi) Pe + O(Pe^2)]`. The appendix row summarizes the stage as "Demand-to-transport map … `zeta_F1(Pe)=A_F1 Omega_Pe^2` and hard ceiling."

## What the script claims to verify

Both scripts fix `kappa_F1 = 12321/5`, `eta_F1 = 37`, solve the Robin root `y_F1` from `y tan y = 37`, form `A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2)`, and carry the Stage-056 closed form `Omega = pi Pe (2 Pe e^Pe + pi)/[(4 Pe^2 + pi^2)(e^Pe - 1)]`. They then assert: `Omega(0+) = 1`; `Omega(∞) = pi/2`; `zeta_F1(0+) = A_F1`; `zeta_F1(∞) = A_F1 pi^2/4`; and that the small-Pe series of `zeta_F1 = A_F1 Omega^2` equals `A_F1(1 + ((4-pi)/pi) Pe)`. The Mathematica script additionally asserts the Robin residual is ~0, splits the small-Pe check into separate constant- and linear-coefficient checks, and adds an independent slope check `Omega'(0+) = (4-pi)/(2 pi)`. The labeled-result prints reproduce `A_F1`, `zeta_F1(0)`, and `zeta_max^(F1)`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `A_F1 = (kappa_F1+pi^2/4)/(kappa_F1+y_F1^2) ≈ 1.00005192880220` | py L37/40, wl L40/43; both print value; constructed from declared inputs | match |
| Demand map `zeta_F1(Pe) = A_F1 Omega_Pe^2` (boxed eq-zetaF1) | py L55 / wl L67-68 form, plus endpoint + series assertions exercising it | match |
| Ceiling `zeta_max^(F1) = A_F1 pi^2/4 ≈ 2.46752922945601` (boxed eq-ceiling) | py L63 `zeta_F1(inf) - A_F1 pi^2/4 == 0`, wl L70; printed L75/62 | match |
| Endpoint `zeta_F1(0) = A_F1` (notes) | py L62, wl L69 | match |
| Carried `Omega(0)=1`, `Omega(∞)=pi/2` (notes/Stage 056) | py L52-53, wl L58-59 | match |
| Small-demand expansion `zeta_F1 = A_F1[1+((4-pi)/pi)Pe+O(Pe^2)]` (notes §4) | py L66-69; wl L72-76 (split) + L77-79 slope | match |
| Robin root `y_F1 ≈ 1.52948248371470` (eq-yF1) | py L36 nsolve (printed L39); wl L39 FindRoot + asserted residual L46 | match |
| Inversion-existence statement `A_F1 ≤ zeta_req ≤ zeta_max` (notes §3) | py L76 / wl narrative — printed as a logical statement, follows from the two endpoint asserts + monotonicity carried from Stage 056 | match (no extra check required; it is a corollary of the asserted endpoints + carried monotonicity) |

`paper_alignment = aligned`. No `mismatch`, no `extra`, no `missing` deliverable.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 52 | `expect_small(Omega0 - 1)` | carried Omega(0)=1 (demand map endpoint) | yes |
| A2 | sympy | 53 | `expect_small(Omega_inf - pi/2)` | carried Omega(∞)=pi/2 (ceiling factor) | yes |
| A3 | sympy | 62 | `expect_small(zeta0 - A_F1)` | `zeta_F1(0)=A_F1` | yes |
| A4 | sympy | 63 | `expect_small(zeta_inf - A_F1 pi^2/4)` | ceiling `zeta_max^(F1)` | yes |
| A5 | sympy | 69 | `expect_small(series(zeta_F1) - A_F1(1+((4-pi)/pi)Pe))` | small-demand expansion | yes |
| B1 | mathematica | 46 | `expectApprox(yF1 tan yF1 - 37, 0)` | Robin root `y_F1` | yes |
| B2 | mathematica | 58 | `expectZero(omega0 - 1)` | Omega(0)=1 | yes |
| B3 | mathematica | 59 | `expectZero(omegaInf - Pi/2)` | Omega(∞)=pi/2 | yes |
| B4 | mathematica | 69 | `expectApprox(zeta0Sym - aF1, 0)` | `zeta_F1(0)=A_F1` | yes |
| B5 | mathematica | 70 | `expectApprox(zetaInfSym - aF1 Pi^2/4, 0)` | ceiling | yes |
| B6 | mathematica | 75 | `expectApprox(Coefficient(seriesDiff,pe,0), 0)` | expansion const coeff | yes |
| B7 | mathematica | 76 | `expectApprox(Coefficient(seriesDiff,pe,1), 0)` | expansion linear coeff | yes |
| B8 | mathematica | 79 | `expectApprox(omegaPrime0 - (4-Pi)/(2Pi), 0)` | Omega slope (drives expansion) | yes |

Every script-side check traces to a specific paper-side deliverable. No orphaned scaffolding.

## Findings

None. The audit is clean.

## Independent-derivation check (Mathematica)

The Mathematica script is an independent re-derivation, not a transliteration:
- The two engines share only the (legitimately carried) Stage-056 closed form for `Omega` and the Stage-057 form for `A_F1` — these are *inputs*, not the algebra under test, and the paper/notes require them to be carried verbatim.
- The limit and series machinery is independent: SymPy uses `sp.limit(...)` / `sp.series(zeta_F1, Pe, 0, 2).removeO()` (py L46-47, L56-57, L66); Mathematica uses `Limit[...]` / `Normal[Series[..., {pe,0,1}]]` (wl L52-53, L67-68, L72). These are distinct CAS pathways, not a line-by-line port.
- The Mathematica script adds checks the SymPy script does not have: it asserts the Robin residual (wl L46), splits the small-Pe residual into separate `Coefficient[...,pe,0]` and `Coefficient[...,pe,1]` checks (wl L75-76), and adds an independent slope check `Omega'(0+) = (4-Pi)/(2 Pi)` via `Limit[D[omega,pe],...]` (wl L77-79). That extra, non-mirrored structure is the opposite of transliteration.
No `mathematica_transliteration` finding.

## Engine cross-check

The engines agree to full reported precision; differences are confined to the last 1-2 digits, consistent with SymPy `nsolve` (root good to ~14 digits) vs Mathematica `FindRoot` at `WorkingPrecision -> 80`:

| Quantity | SymPy output | Mathematica output | Paper boxed |
|---|---|---|---|
| `y_F1` | 1.52948248371469963… | 1.5294824837146996499… | 1.52948248371470 |
| `A_F1` | 1.00005192880219520… | 1.00005192880219532… | 1.00005192880220 |
| `zeta_F1(0)` | 1.00005192880219492… | 1.00005192880219532… | (=A_F1) |
| `zeta_max^(F1)` | 2.46752922945601123… | 2.46752922945601223… | 2.46752922945601 |

All assertions pass in both transcripts (SymPy residuals at 1e-16, Mathematica residuals at 1e-40…1e-55 and exact `0` for the symbolic `Omega` limits). No `engine_disagreement`.

Note (not a finding): the SymPy transcript prints `Robin residual = -1.42e-14`, larger than the `nsolve(tol=1e-34)` convergence target. This is an evaluation-precision artifact of the steep `tan` near `y≈1.529` (the residual is `y tan y - 37` evaluated at finite working precision, not the root error), and it is only *printed*, not asserted, in SymPy. The root itself is accurate: Mathematica's `WorkingPrecision -> 80` FindRoot drives the same residual to 2.27e-55 and asserts it (wl L46, PASS), and every downstream value built from `y_F1` (A_F1, zeta endpoints) matches the paper to the full quoted precision. Coverage is adequate.

## Verdict justification

Clean. I attacked the small-Pe check for the classic factor-of-2 trap: because `zeta_F1 = A_F1 Omega^2` and `Omega = 1 + ((4-pi)/(2pi))Pe`, squaring doubles the slope to `((4-pi)/pi)`, and the script's `expected_series = A_F1(1+((4-pi)/pi)Pe)` correctly carries that doubled coefficient — the check would fail loudly if the factor of 2 were wrong. I checked whether SymPy's `expect_small` could silently pass a `Pe`-dependent residual: it cannot — a leftover `Pe` term reaches `complex(val)` and raises `TypeError` (fail-loud), and the transcript confirms the residual simplified to literal `0`. I confirmed every assertion is non-tautological: each compares an independently-Taylor-expanded / limit-evaluated quantity against a separately-written closed form, never `x == x`. I verified all carried inputs against upstream sources: `kappa_F1=12321/5` and `eta=37` (Stage 074 notes L95/L11; Stage 075 card L7), the `A_F1` form (Stage 057 card L19 / notes L102), the `Omega` closed form and its `1+((4-pi)/(2pi))Pe` expansion (Stage 056 notes L204/L176) — all match the script verbatim. Paper alignment is exact: both `\stagefield{Output}` boxed deliverables (demand map + ceiling) and the notes' endpoint/expansion/Robin-root facts each have a faithful, non-tautological script-side check, and the two boxed numerics (`A_F1`, `zeta_max^(F1)`) match the transcripts to full precision. Both engines are present, independent, and agree. Outputs are fresh (both `.txt` newer than their scripts). I read the paper card, notes, and appendix row before the scripts, and the scripts' claim matches the paper's claim.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 values checked, 0 misaligned

Deliverable-level table:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `y_F1 ≈ 1.52948248371470` | py out L5 `1.52948248371469963…`; wl out L5 `1.5294824837146996499…` | tex L20 `y_{\rm F1}\simeq1.52948248371470`; md L61 | MATCH |
| `A_F1 ≈ 1.00005192880220` | py out L6/L27 `1.00005192880219520…`; wl out L6 `1.00005192880219532…` | tex L33 (boxed) `A_{\rm F1}\simeq1.00005192880220`; md L63/L87 | MATCH |
| demand map `zeta_F1(Pe) = A_F1 Omega_Pe^2` | py out L15/L26; wl out L22; asserted (A3-A5/B4-B8) | tex L39 (boxed eq-zetaF1); md L81 | MATCH |
| `zeta_F1(0) = A_F1 ≈ 1.0000519288…` | py out L16/L28 `1.0000519288021949…`; wl out L23 `1.0000519288021953…` | tex (=A_F1 box, L33); md L87 `zeta_F1(0) = A_F1 ≈ 1.00005192880219520364…` | MATCH |
| `zeta_max^(F1) = A_F1 pi^2/4 ≈ 2.46752922945601` | py out L17/L29 `2.46752922945601123…`; wl out L24 `2.46752922945601223…` | tex L46 (boxed eq-ceiling) `\simeq2.46752922945601`; md L35/L90/L158 | MATCH |
| `Omega(0)=1`, `Omega(∞)=pi/2` (carried branch limits) | py out L11-12; wl out L17-18 | md L77 / Stage-056 notes L208/L212; tex via eq-zetaF1 factor | MATCH |

INTERNAL (verification scaffolding, no prose expected; raise no finding): `kappa_F1 = 12321/5` and `eta_F1 = 37` (declared inputs, do appear in tex L7/L27 and notes as inputs — reported, not deliverables of this stage); the printed Robin-residual diagnostic (`-1.42e-14` / `2.27e-55`); all `... - target = 0` residual lines; `Omega'(0+) = (4-Pi)/(2Pi)` (an intermediate driving the asserted expansion, present in Stage-056 notes); the `small-Pe` printed polynomial line (py out L20); pass/PASS flags and tolerances.

Every emitted deliverable value reconciles against the `.tex` card and/or the `.md` notes. No MISMATCH, no MISSING-DELIVERABLE.

## Self-test notes

Variable-independence: the only derivative in either script is `D[omega, pe]` (wl L77); `omega` genuinely depends on `pe`, so the slope is non-trivial (it evaluates to `-1/2 + 2/Pi ≠ 0`, wl out L37) — no identically-zero-derivative trap. Symmetry/parity: no unbounded symmetric-domain integrals in this unit, so the parity trap does not apply. Trivial-case/factor-of-2: I confirmed the squaring of `Omega`'s linear slope correctly doubles to `(4-pi)/pi` in the asserted `zeta_F1` expansion, and that `expect_small` fails loudly (TypeError) on any residual still carrying `Pe` — so the small-Pe check is substantive, not a silent pass.
