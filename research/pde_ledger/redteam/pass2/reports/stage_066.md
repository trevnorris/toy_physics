---
unit_id: 066
batch: III.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage066_wall_figure_of_merit.md]
  paper_appendix: present
---

# Audit unit 066 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_066.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage066_wall_figure_of_merit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows 110, 118)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage066_wall_figure_of_merit_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage066_wall_figure_of_merit_mathematica_audit.txt`

## What the paper claims

Stage 066 collapses the explicit thin-wall parent support/source threshold (from Stages 064–065) onto a single dimensionless wall control. `\stagefield{Output}` reads verbatim: "The wall control \eqref{eq:app-stage066-Wwall} and exact window \eqref{eq:app-stage066-Wwindow}." Deliverable 1 is the boxed figure of merit `W_wall = 4 pi a^2 L^2 J_1 V0^2 / (T_X ell)` (eq line 17). Deliverable 2 is the exact matched-branch window: `W_wall <= Pe_req/Delta_inf => fail` and `W_wall >= Pe_req/Delta_0 => succeed` (eq lines 23–25), with only the intermediate band needing the full operator solve. The notes add three supporting deliverables: the `W_wall = kappa G_eq^(tw)` link with `kappa = K_X L^2/T_X` (notes §1), strict monotonicity in six signed directions (notes §3), and the constant-compressibility form `W_H = 4 pi a^2 L^2 I_f V0^2 / (H_w T_X ell)` with the same window (notes §4). The appendix row 110 summarizes it as "Single wall control W_wall and exact window." This is a purely symbolic, dimensionless stage; the only numeric constant is the `4 pi` prefactor — there is no 100/168-type figure-of-merit number.

## What the script claims to verify

The SymPy script defines `W_wall`, `W_fail = Pe_req/Delta_inf`, `W_suff = Pe_req/Delta_0`, then checks that substituting the inverted Stage-065 wall-amplitude thresholds (`V0_fail_sq`, `V0_suff_sq`) into `W_wall` reproduces `W_fail`/`W_suff` exactly (the closure that makes the wall window exact). It then verifies the six monotonicity signs via signed simplification of the partials, and checks the `W_H` constant-compressibility reduction (`J1 -> I_f/H_w`) plus its two threshold round-trips. The Mathematica script mirrors all of this and additionally re-derives `W_wall` independently from the gain `G_eq^(tw) = 4 pi a^2 J_1 V0^2/(K_X ell)` via the `kappa = K_X L^2/T_X` substitution (`W_wall - kappa G_eq^(tw) = 0`). The bottom-line printed statements assert the two boxed paper deliverables.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| D1: `W_wall = 4 pi a^2 L^2 J_1 V0^2/(T_X ell)` | py:51 `W_wall = ...`; wl:40; wl:48 derives it from `kappa G_eq^(tw)` | match |
| D2: exact window `W_wall <= Pe_req/Delta_inf`, `>= Pe_req/Delta_0` | py:63–64 threshold round-trips; wl:53–54 | match |
| Notes §1: `W_wall = kappa G_eq^(tw)`, `kappa=K_X L^2/T_X` | wl:48 only (py has no kappa-route check) | match (mathematica only) |
| Notes §3: six signed monotonicity directions | py:87–92; wl:73–78 | match |
| Notes §4: `W_H = 4 pi a^2 L^2 I_f V0^2/(H_w T_X ell)` + same window | py:95,97–99; wl:82,85–92 | match |

Every paper-side deliverable has a faithful, non-tautological script-side check. `paper_alignment: aligned`. The `4 pi` prefactor matches across card, notes, and both scripts; no stale figure-of-merit number is present.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 63 | `simplify(W_wall.subs(V0**2,V0_fail_sq) - W_fail)==0` | D2 (fail threshold) | yes |
| A2 | sympy | 64 | `... - W_suff == 0` | D2 (succeed threshold) | yes |
| A3 | sympy | 87–92 | six `simplify(d... >/<0) is sp.true` | Notes §3 monotonicity | yes |
| A4 | sympy | 97 | `W_wall.subs(J1,If/Hw) - W_H == 0` | Notes §4 reduction | yes |
| A5 | sympy | 98–99 | `W_H` threshold round-trips == 0 | Notes §4 window | yes |
| A6 | mathematica | 48 | `wWall - (gEqTw/.kx->tx*kappa/len^2)*kappa == 0` | D1 + Notes §1 (kappa route) | yes |
| A7 | mathematica | 53–54 | `W_wall(V0_*)^2 - W_*` == 0 | D2 | yes |
| A8 | mathematica | 73–78 | six `expectTrue` sign checks | Notes §3 | yes |
| A9 | mathematica | 85–92 | `W_H` reduction + window round-trips | Notes §4 | yes |

No tautological or unanchored rows. The `V0_*_sq` substitutions are the inverted Stage-065 thresholds written as independent literals; the round-trip would break if the `4 pi` or any exponent were wrong, so the checks are substantive.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.txt:11,22-25,36`
- (also) `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage066_wall_figure_of_merit_mathematica_audit.txt:11`

**What's wrong:**
The committed SymPy `.txt` (mtime 2026-05-11) predates the current script (mtime 2026-06-03) and does not reflect it. The output banner reads "STAGE 49 — DIMENSIONLESS WALL FIGURE OF MERIT" (line 11) and "STAGE 49 AUDIT PASSED" (line 36), whereas the current script prints "STAGE 66" (py:45) and "STAGE 66 AUDIT PASSED" (py:101). More substantively, the output's MONOTONICITY block (lines 22–25) prints only four derivatives (`dW/d(V0^2)`, `dW/da`, `dW/dL`, `dW/dell`) and is missing the `dW/dJ1` and `dW/dT_X` lines that the current script prints (py:84–85) and asserts (py:91–92). The Mathematica output banner is likewise stale at line 11 ("STAGE 049"), though its body content already includes all six derivatives and matches the current `.wl`.

**Why this matters:**
The committed transcript no longer evidences the two monotonicity checks (J1, T_X) that the current script actually performs, and carries the wrong stage label. A reader trusting the transcript would not see the full set of asserted directions.

**Required change:**
Re-run both scripts and refresh both committed `.txt` outputs so the banner reads "STAGE 66"/"STAGE 066" and the SymPy MONOTONICITY block shows all six `dW/...` lines. No script-logic change required (the script already produces the correct content).

**Verification:**
After refresh, the SymPy `.txt` MONOTONICITY block contains `dW/dJ1 =` and `dW/dT_X =` lines, and both transcripts' banners read STAGE 66/066. The orchestrator's independent re-run regenerates these.

### F2 — stale numbering self-labels (paper_misalignment / notes_contradicts_script, label-only)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py:2` ("Stage 49 SymPy audit")
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py:14` ("Stage-48 thin-wall thresholds")
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py:59` ("# Stage-48 thresholds")

**What's wrong:**
The SymPy docstring titles the file "Moving-Throat PDE — Stage 49 SymPy audit" (line 2) — a stale +17 pre-renumber self-label; the canonical stage is 066 (filename, card, banner py:45 all say 66). Two inline comments reference "Stage-48" thresholds (lines 14, 59); the card states the inputs are Stage 065 (thin-wall branch) and Stage 059 (operator thresholds), and the notes attribute the wall thresholds to Stage 065. `48 + 17 = 65`, confirming these are stale pre-renumber labels for the Stage-065 source. These are unambiguous self/upstream labels (not genuine cross-stage references to a stage 48), so under the in-loop Reading-2 numbering policy they are corrected on a verdict:findings stage. The Mathematica `.wl` banner already uses the correct "STAGE 066" (wl:32); no stale labels were found in the `.wl` source.

**Why this matters:**
Stale stage labels in the docstring/comments mislead readers about which canonical stage the file belongs to and which upstream stage supplies the thresholds. Label-only; no math changes.

**Required change:**
In the SymPy script: line 2 `Stage 49 SymPy audit` -> `Stage 066 SymPy audit`; line 14 `the Stage-48 thin-wall thresholds` -> `the Stage-065 thin-wall thresholds`; line 59 comment `# Stage-48 thresholds` -> `# Stage-065 thresholds`. No assertion or numeric change; re-run to refresh the transcript.

**Verification:**
`grep -n "Stage 49\|Stage-48"` over the SymPy script returns nothing after the edit; banner/output unaffected except for the unrelated F1 refresh.

## Independent-derivation check (Mathematica)

The `.wl` is not a pure line-by-line transliteration. It uses its own symbol set (`len`, `kx`, `kappa`, `gEqTw`, `ifMom`) and, crucially, adds an independent derivation path the SymPy script lacks: at wl:43 it defines the thin-wall gain `gEqTw = 4 Pi a^2 j1 v0^2/(kx ell)` and at wl:48 checks `wWall - (gEqTw /. kx -> tx*kappa/len^2)*kappa == 0`, i.e. it re-derives `W_wall` from the gain via `kappa = K_X L^2/T_X` rather than asserting the closed form against itself. The remaining checks (threshold round-trips, monotonicity, `W_H` reduction) parallel the SymPy ones but operate on independently re-entered expressions under FullSimplify. Verdict: independent (the kappa/gain route is a genuinely distinct derivation of D1).

## Engine cross-check

Both engines produce identical closed forms. SymPy `.txt`: `W_wall = 4*pi*J1*L**2*V0**2*a**2/(T_X*ell)`; Mathematica `.txt`: `W_wall = (4*a^2*j1*len^2*Pi*v0^2)/(ell*tx)` — algebraically identical. `W_fail = Pe_req/Delta_inf`, `W_suff = Pe_req/Delta_0` agree. All threshold round-trips, the `J1 = I_f/H_w` reduction, and all six monotonicity derivatives match between transcripts (the SymPy `.txt` is merely missing the dJ1/dT_X *print* lines due to staleness, but the assertions are in the source and the Mathematica transcript shows them PASS). No sign or factor disagreement. `engines_agree: true`.

## Verdict justification

The scripts faithfully and non-tautologically verify both boxed paper deliverables (the `W_wall` form and the exact fail/succeed window) plus all three notes-level supporting deliverables (kappa link, monotonicity, constant-compressibility form). I attacked the threshold round-trips (could they be tautological by construction? — no, they cross-check two independently-written formulas and would break on a wrong `4 pi` or exponent), the monotonicity derivatives (variable-independence trap — all six variables genuinely appear in `W_wall`, so no zero-derivative trivial pass), the symbol domains (all positive/real, justified by the physical setup), and the `4 pi` prefactor against the prior stale-"168" warning (no numeric figure-of-merit number exists here; only `4 pi`, which matches everywhere). The only defects are non-math: the committed SymPy transcript is stale (wrong "STAGE 49" label, missing two monotonicity print lines) and the SymPy docstring/comments carry stale +17 pre-renumber labels ("Stage 49", "Stage-48"). Hence `verdict: findings`, both low-severity and label/refresh-only; the math is clean and aligned.

## Value Reconciliation (pass-2 augmentation)

All emitted deliverable values are symbolic (this is a dimensionless figure-of-merit stage); the only numeric constant is the `4 pi` prefactor.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `W_wall = 4 pi a^2 L^2 J_1 V0^2/(T_X ell)` | py:51 / wl:40 / sympy.txt:13, math.txt:13 | tex:14-18 (boxed eq:app-stage066-Wwall); md:53, md:87 | MATCH |
| `4 pi` prefactor of W_wall | py:51 / wl:40 | tex:17 `4\pi`; md:53 `4 pi` | MATCH |
| `W_fail = Pe_req/Delta_inf` | py:52 / wl:41 / sympy.txt:14, math.txt:14 | tex:23 (eq:app-stage066-Wwindow); md:69 (`W_fail`) | MATCH |
| `W_suff = Pe_req/Delta_0` | py:53 / wl:42 / sympy.txt:15, math.txt:15 | tex:25; md:71 (`W_suff`) | MATCH |
| `G_eq^(tw) = 4 pi a^2 V0^2 J_1/(K_X ell)` | wl:43 | md:33-34 | MATCH |
| `W_wall = kappa G_eq^(tw)`, `kappa=K_X L^2/T_X` | wl:48 / math.txt:16 | md:48-53 (§1) | MATCH |
| `W_H = 4 pi a^2 L^2 I_f V0^2/(H_w T_X ell)` | py:95 / wl:82 / sympy.txt:30, math.txt:48 | tex (window uses W_wall; W_H is notes-only); md:118 | MATCH (notes carrier) |
| `J1 = I_f/H_w` reduction | py:97 / wl:85 | md:111 | MATCH |
| six monotonicity sign directions | py:80-92 / wl:66-78 / math.txt:32-43 | md:91-96 (§3, six signed directions) | MATCH |

Internal scaffolding (no finding): `V0_fail_sq`, `V0_suff_sq` (inverted Stage-065 thresholds used only to drive the round-trip asserts), `Vp`/`v0Sq` (monotonicity substitution variable), the six `dW/...` partials (intermediate quantities for the sign asserts), all `PASS`/`FAIL`/residual flags, banner strings.

`reconciliation: complete; 9 deliverable values checked, 0 misaligned`
