---
unit_id: 139
batch: IV.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-07T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage139_family1_actual_mouth_gains.md]
  paper_appendix: present
---

# Audit unit 139 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_139.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage139_family1_actual_mouth_gains.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (stage card embedded via `\input{stages/stage_139}` at line 1312; no separate summary row — the lines mentioning "139" at 893/951 belong to the later co-evolving-fixedpoint section, not Stage 139)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.txt`

## What the paper claims

The stage card (`stage_139.tex:7`) states its purpose: "Stage~139 is a coupled mouth fixed point and gain selection ledger step. Its audit target is the verification output quoted below," with the claim status `\StatusExactClosure / \StatusNumerical` (`:4`) and the block narrative "Natural and exact-compensated gain pairs are evaluated on the Family--1 branch" (`:16`). The card is deliberately terse and carries no numbers — it is "a derivation ledger entry, not an unconditional actual-branch theorem" (`:27`). The detail lives in the notes. The notes (`moving_throat_pde_stage139_*.md`) enumerate the deliverables: the Family-1 radius `r_F1 ≈ 1.77799353547498` (notes:22), the imported `Π_* ≈ 1.50882951349316` and `S_q(Π_*) ≈ 0.658075937605429` (notes:27-28), the natural equal-normalized branch (`g_c = 1`) with `R_q^nat ≈ 0.145454452260421` (notes:64), `M_s^nat,* ≈ 1.66854252965624` (notes:74), `M_q^nat,* ≈ -0.242696939724365` (notes:81); the exact-compensated branch with `g_c ≈ 0.758035078944663` (notes:98), `R_q = 1/4` (notes:100), `M_s^comp,* ≈ 1.80594111095636` (notes:110), `M_q^comp,* ≈ -0.451485277739090` (notes:116); and the comparison ratios `M_s^comp/M_s^nat − 1 ≈ 0.0823464663669` (notes:128) and `|M_q^comp|/|M_q^nat| ≈ 1.86028418097` (notes:132). The card's Checks (`:21-25`) require: the gain pair against outlet consistency, the self-matched susceptibility closure, and that fixed points are recorded as numerically located, not closed-form constants.

## What the script claims to verify

The SymPy script (and its line-for-line Mathematica mirror) computes `r_F1` from the Stage-121 closed form `√(12/π²·(37/20)²−1)` (py:6), imports `Π_*` and `S_q(Π_*)` as Stage-134 literals (py:8-9), then derives both gain pairs: natural (`R_q^nat = (1−r_F)²/(1+r_F²)`, py:13) and compensated (`g_- = r_F − √(1+r_F²)/2`, py:21; `R_q^comp = (g_−−r_F)²/(1+r_F²)`, py:22). It asserts each emitted deliverable against the notes literals at `tol = 1e-12` (py:44-60). Two checks are deliberately non-trivial: (R2) the `g_-^F1 = 0.758035…` value (py:52), documented as the falsifiable lower-branch discriminator (a branch/sign typo gives the upper branch `g_+ ≈ 2.79` and FAILS); and (R1) an independent reconstruction of `S_q(Π_*)` from the Stage-134 closed-form kernel at `κ=π/2` (py:66-71), confirming the imported literal via a route that does not reuse `M_s = Π_*/(1−R_q S_q)`. The outlet-form residual (py:76-79) and `R_q^comp = 1/4` (py:84) are explicitly labelled structural/definitional sanity, not load-bearing literal checks.

## Paper ↔ script cross-check

| Paper-side deliverable (notes line) | Script-side check | Status |
|---|---|---|
| `r_F1 ≈ 1.77799353547498` (notes:22) | py:6 closed form, asserted py:44 | match |
| `S_q(Π_*) ≈ 0.658075937605429` (notes:28) | py:9 import, independently reconstructed py:66-71 | match |
| `R_q^nat ≈ 0.145454452260421` (notes:64) | py:13, asserted py:47 | match |
| `M_s^nat,* ≈ 1.66854252965624` (notes:74) | py:14, asserted py:55 | match |
| `M_q^nat,* ≈ -0.242696939724365` (notes:81) | py:15, asserted py:56 | match |
| `g_c^comp ≈ 0.758035078944663` (notes:98) | py:21, asserted py:52 (branch discriminator) | match |
| `R_q^comp = 1/4` (notes:100) | py:22, asserted py:84 (definitional) | match |
| `M_s^comp,* ≈ 1.80594111095636` (notes:110) | py:23, asserted py:59 | match |
| `M_q^comp,* ≈ -0.451485277739090` (notes:116) | py:24, asserted py:60 | match |
| `M_s^comp/M_s^nat − 1 ≈ 0.0823464663669` (notes:128) | py:35 printed | match (printed deliverable) |
| `|M_q^comp|/|M_q^nat| ≈ 1.86028418097` (notes:132) | py:36 printed | match (printed deliverable) |
| outlet consistency `Π_* = M_s + M_q S_q` (notes:33; card check :22) | py:76-79 structural residual = 0 | match (true by construction, correctly labelled non-load-bearing) |

`paper_alignment: aligned`. Every notes deliverable has a faithful script-side check; no script-side value lacks a notes counterpart.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 44 | `abs(rF − 1.77799353547498) < 1e-12` | r_F1 (notes:22) | yes |
| A2 | sympy | 47 | `abs(Rq_nat − 0.145454452260421) < 1e-12` | R_q^nat (notes:64) | yes |
| A3 | sympy | 52 | `abs(g_minus − 0.758035078944662…) < 1e-12` | g_c^comp / lower-branch discriminator (notes:98) | yes |
| A4 | sympy | 55 | `abs(Ms_nat − 1.66854252965624) < 1e-12` | M_s^nat,* (notes:74) | yes |
| A5 | sympy | 56 | `abs(Mq_nat − (−0.242696939724365)) < 1e-12` | M_q^nat,* (notes:81) | yes |
| A6 | sympy | 59 | `abs(Ms_comp − 1.80594111095636) < 1e-12` | M_s^comp,* (notes:110) | yes |
| A7 | sympy | 60 | `abs(Mq_comp − (−0.451485277739090)) < 1e-12` | M_q^comp,* (notes:116) | yes |
| A8 | sympy | 71 | `abs(Sq_recon − Sq_star) < 1e-12` (kernel route) | S_q(Π_*) (notes:28); card check :23 | yes (independent route) |
| A9 | sympy | 84 | `abs(Rq_comp − 1/4) < 1e-25` | R_q^comp definitional (notes:100) | partial (definitional; self-labelled non-load-bearing) |
| A10 | mathematica | 67 | `expectApprox r_F1` | r_F1 | yes |
| A11 | mathematica | 68 | `expectApprox R_q^nat` | R_q^nat | yes |
| A12 | mathematica | 72 | `expectApprox g_-^F1 value` | g_c^comp / discriminator | yes |
| A13 | mathematica | 73-76 | `expectApprox M_s/M_q nat/comp` | the four gains | yes |
| A14 | mathematica | 85 | `expectApprox S_q recon from kernel` | S_q(Π_*) (independent route) | yes |
| A15 | mathematica | 89 | `expectApprox R_q^comp − 1/4` | R_q^comp definitional | partial (definitional) |

A9/A15 are knowingly definitional (the script itself says so at py:81-83 / wl:40-44) and are not relied on for branch discrimination; the falsifiable content lives in A3/A12 (g_- value) and A2/A11 (R_q^nat ≠ 1/4). No tautological or unanchored load-bearing rows.

## Findings

None.

## Independent-derivation check (Mathematica)

A `.wl` exists, so this is mandatory. Three corresponding sections:

1. Compensated branch root:
   - py:21 `g_minus = sp.N(rF - sp.sqrt(1 + rF**2)/2, 30)`
   - wl:45 `gMinus = N[rF - Sqrt[1 + rF^2]/2, 30]`
2. `S_q` kernel reconstruction (R1):
   - py:67-69 `S_kernel = (Pi_star * (kappa_q*sp.tanh(kappa_q) + Pi_star*(sp.exp(-Pi_star)/sp.cosh(kappa_q) - 1)) / ((1 - sp.exp(-Pi_star)) * (kappa_q**2 - Pi_star**2)))`
   - wl:82-84 `sQRecon = N[piStar*(kappaQ*Tanh[kappaQ] + piStar*(Exp[-piStar]/Cosh[kappaQ] - 1)) / ((1 - Exp[-piStar])*(kappaQ^2 - piStar^2)), 30]`
3. Natural loading ratio:
   - py:13 `Rq_nat = sp.N((1 - rF)**2 / (1 + rF**2), 30)`
   - wl:36 `rQNat = N[(1 - rF)^2/(1 + rF^2), 30]`

**Call: PARTIAL.** Structurally the `.wl` is a faithful symbol-for-symbol port of the `.py` algebra — same variable choreography, same intermediate steps, same expression order; the `.wl` comment at wl:42 even states it computes `g_-` as "a sanctioned mirror of the SymPy route." However, Stage 139 is a pure numerical-evaluation / gain-selection ledger step: it substitutes fixed upstream closed forms (Stage 121 `r_F1`, Stage 134 `Π_*`/`S_q`, Stage 138 `R_q`) into the gain expressions. There is no meaningful alternative algebraic route to a numerical substitution, so the legitimate second-engine value here is cross-CAS numerical agreement, which IS delivered independently: SymPy and Mathematica are different arithmetic engines and report distinct residual magnitudes (e.g. SymPy kernel residual region vs Mathematica `4.468…*^-16` at wl-output:30; `M_q^nat` residual `9.44e-16` at wl-output:24). The one genuinely route-independent check (R1, evaluating the Stage-134 kernel rather than reusing the outlet form `M_s = Π_*/(1−R_q S_q)`) is present in BOTH engines. For an evaluation stage this is acceptable and consistent with the established policy that a single algebraic route is permissible where independent re-derivation is genuinely impossible — so I do NOT raise a blocking `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce identical deliverable values to ~15 digits and all checks PASS:

| value | SymPy output | Mathematica output |
|---|---|---|
| r_F1 | 1.77799353547497811851563491229 (sympy:1) | 1.77799353547497811851563491229… (wl:5) |
| R_q^nat | 0.145454452260420126101421595368 (sympy:3) | 0.1454544522604201261014215953679… (wl:7) |
| M_s^nat,* | 1.66854252965623917495245119417 (sympy:4) | 1.6685425296562391506838… (wl:8) |
| g_-^F1 | 0.758035078944662826919680890414 (sympy:6) | 0.7580350789446628269196808904141… (wl:10) |
| R_q^comp | 0.25 (sympy:7) | 0.25 (wl:11) |
| M_s^comp,* | 1.80594111095635901074188994207 (sympy:8) | 1.805941110956358994454… (wl:12) |
| M_q^comp,* | -0.451485277739089752685472485519 (sympy:9) | -0.45148527773908974861… (wl:13) |
| shell shift | 0.082346466366924027407 (sympy:10) | 0.08234646636692403338800328712513354163 (wl:14) |
| mixed ratio | 1.8602841809701057886 (sympy:11) | 1.86028418097010579885369254932667945052 (wl:15) |

Differences appear only beyond the 15th significant digit (the SymPy literals are printed at 30 digits but anchored to ~15-digit Stage-134 imports; the agreement is at the level both engines claim). No engine disagreement.

## Verdict justification

`clean`. I read the card, the notes, and the relevant appendix section before the scripts. Attacks tried that failed: (1) the `168π²`-vs-`100π²` staleness — the radius uses `12/π²·(37/20)²−1 = (4107−100π²)/(100π²)`, i.e. `√(4107−100π²)/(10π)`, the canonical Family-1 form with the correct `100π²`, no stale `168π²`; (2) branch blindness — the `g_-^F1 = 0.758035…` check (A3/A12) genuinely discriminates the lower branch from the upper `g_+ ≈ 2.79`, and the `R_q^comp = 1/4` check is correctly self-labelled as branch-blind/definitional and not relied upon; (3) tautology in the outlet form — `Π_* = M_s + M_q S_q` is correctly demoted to a structural print (py:76-79) and is NOT counted as a literal verification; (4) tautology in the S_q check — the kernel reconstruction (A8/A14) recomputes `S_q(Π_*)` from the Stage-134 closed-form kernel without referencing `Sq_star`, and the kernel matches the Stage-134 notes form (notes134:52-60) exactly, so it is a real falsifiable cross-check (residual `~4.5e-16`); (5) value drift — every emitted deliverable matches the notes literals to within the `1e-12` tolerance (actual residuals `1e-15`–`1e-16`). The `.wl` is a structural port but for this numerical-evaluation stage that is acceptable (cross-CAS numerical confirmation; no independent algebraic route exists), so no transliteration finding. Outputs are fresh (both `.txt` mtime 2026-05-29 18:02, both scripts 16:51/16:52).

## Self-test notes

Checked: (1) variable-independence — no `diff`/`D` in either script (pure numerical substitution), so the zero-derivative trap does not apply. (2) symmetry/parity — no unbounded-domain integrals; the kernel is a closed-form evaluation at a single point. (3) trivial-case pre-check — I hand-verified each gain by substituting the rounded constants (`R_q^nat ≈ 0.1455`, `M_s^nat ≈ 1.6685`, `g_- ≈ 0.7580`, `R_q^comp = 0.25`, `M_s^comp ≈ 1.8059`, ratios `0.0823` / `1.8603`); all reproduce the asserted literals. (4)/(5) N/A — no missing-script finding and no prescribed fixes, so no paper round-trip needed.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit. Base record = script source + committed `.txt` outputs (both present and fresh; nothing executed).

| value | source (py / wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `r_F1 = 1.77799353547497811851…` | py:26 / wl:50; sympy-out:1, wl-out:5 | notes:22 `≈ 1.77799353547498` | MATCH |
| `S_q(Π_*) = 0.658075937605429` | py:27 / wl:51; sympy-out:2, wl-out:6 | notes:28 `≈ 0.658075937605429` | MATCH |
| `R_q^nat = 0.145454452260420…` | py:28 / wl:52; sympy-out:3, wl-out:7 | notes:64 `≈ 0.145454452260421` | MATCH |
| `M_s^nat,* = 1.66854252965623…` | py:29 / wl:53; sympy-out:4, wl-out:8 | notes:74 `≈ 1.66854252965624` | MATCH |
| `M_q^nat,* = -0.242696939724364…` | py:30 / wl:54; sympy-out:5, wl-out:9 | notes:81 `≈ -0.242696939724365` | MATCH |
| `g_-^F1 = 0.758035078944662826…` | py:31 / wl:55; sympy-out:6, wl-out:10 | notes:98 `≈ 0.758035078944663` | MATCH |
| `R_q^comp = 0.25` | py:32 / wl:56; sympy-out:7, wl-out:11 | notes:100 `R_q=1/4` | MATCH |
| `M_s^comp,* = 1.80594111095635…` | py:33 / wl:57; sympy-out:8, wl-out:12 | notes:110 `≈ 1.80594111095636` | MATCH |
| `M_q^comp,* = -0.451485277739089…` | py:34 / wl:58; sympy-out:9, wl-out:13 | notes:116 `≈ -0.451485277739090` | MATCH |
| `shell gain fractional shift = 0.0823464663669…` | py:35 / wl:59; sympy-out:10, wl-out:14 | notes:128 `≈ 0.0823464663669` | MATCH |
| `mixed gain magnitude ratio = 1.86028418097…` | py:36 / wl:60; sympy-out:11, wl-out:15 | notes:132 `≈ 1.86028418097` | MATCH |

INTERNAL (verification scaffolding, no prose expected; no finding):
`Sq_recon` (kernel reconstruction of S_q, drives A8/A14 — its target value `S_q` is already reconciled above), `outlet form residual (nat/comp, structural)` (= 0 by construction, py:76-79 / wl:87-88), `tol_literal`/`tol_algebraic` (tolerances), all `… residual = …` lines and `PASS:`/`all assertions passed` flags.

reconciliation: complete; 11 values checked, 0 misaligned.

(Note: the `.tex` card is intentionally terse and reports no numerics — per its own text it is "a derivation ledger entry," and the augmentation guard explicitly allows the card to omit values when the `.md` notes carry them correctly, which they do. Also noted: Stage-134 notes:86 quote `S_q(Π_*) ≈ …428` (last digit 8) vs Stage-139 notes/script `…429` (last digit 9) — a sub-ULP last-digit rounding of the same true value, upstream of this stage and within the 1e-12 tolerance; not a Stage-139 deliverable mismatch.)
