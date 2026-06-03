---
unit_id: 228
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 228 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_228.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 68, 750-804 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

The stage card `\stagefield{Output}` states: "Numerator/denominator split, dynamic-window audit, and first branch-level comparison between the static $\Xi_1$ budget and resonant dynamic survival." The appendix (eq:app-part07-numerator-denominator-split, line 785) gives the exact split `\Xi_1=2(\pi_1-\delta_1)` with `\pi_1:=P_{01}/P`, `\delta_1:=\Delta_{01}/\Delta`, and the rigid-subcorridor theorem (lines 791-797): pure-transfer + `\pi_1=0` leaves a 1D numerator-rigid survivor; pure-transfer + `\delta_1=0` leaves a 1D denominator-rigid survivor; both rigidities kill the corridor. Both survivors pass the first wall-like dynamic-window audit, and the first active ceiling remains the transported static budget, not the dynamic window. The notes enumerate nine concrete deliverables (notes section 7): the split identity, the exact row formulas for `\pi_1`/`\delta_1`, the four rank/nullity counts, the reduced determinant, the positive-`\Xi_1` unit generators `v_num`/`v_den`, the undeformed wall-like pole census + `R_Q` values, the first-order log-slopes of `P_0` and the two `R_Q`, the resulting dynamic ceilings, and the comparison showing the static ceiling still binds first. The claim is `Mixed: ExactClosure + Numerical`.

## What the script claims to verify

The SymPy script builds the explicit finite-throat one-port compiler on the concrete sample branch (`kappa=2sqrt2/pi`, `lamB=1/2, lamU=3/10, lamW=2/5, lamR=1/4, OmU=1, OmW=7/5, varpi=2, M=1`), derives the compatibility wall stiffness `K_compat` in-script, then dresses each primitive with `exp(eps*x)` and differentiates at `eps=0` to get the first-order responses. It asserts: (1) `Xi_transfer - 2*(pi1-delta1) == 0` symbolically (line 191); (2) the exact coefficient rows of `pi1` and `delta1` against literal targets (lines 199-202); (3) the four rank/nullity counts via `Matrix.rank()`/`nullspace()` (lines 217-247); (4) the reduced determinant against a literal symbolic target (lines 249-254); (5) the unit generators against numeric targets and their `pi1`/`delta1`/`Xi1` projections (lines 279-301); (6) the numeric pole census and `R_Q` values via `np.roots` on the quartic and the residue-ratio formula (lines 366-373); (7) finite-difference log-slopes of `P_0` and the two `R_Q` (lines 406-421); (8) the dynamic and transported-static ceilings, plus the `dynamic > static` inequalities (lines 436-456). The assertions are substantive and non-tautological.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Split `Xi_1 = 2(pi_1 - delta_1)` | line 191 `assert simplify(Xi_transfer - 2*(pi1-delta1)) == 0` | match |
| Exact row formula for `pi_1` | lines 193, 199-200 vs `[3/19,16/19,3/19,32/19,0]` | match |
| Exact row formula for `delta_1` | lines 194, 201-202 vs `[0,0,50/(25-98pi^2),196pi^2/(98pi^2-25),196pi^2/(98pi^2-25)]` | match (script); **notes say 247pi^2 — see F2** |
| Rank/nullity counts (3/2,4/1,4/1,5/0) | lines 217-247 | match |
| Reduced determinant != 0 | lines 249-254 vs `196*(200+147pi^2)*(...)/(475*(...))` | match (script); **notes say 247(251+215pi^2)(...) — see F3** |
| Unit generators `v_num`, `v_den` | lines 276-282 numeric targets | match |
| Wall-like pole census + R_Q | lines 367-373 (`30.1999..`, `36.1711..`, `omega`s, `P_0`) | match |
| First-order log-slopes (P_0, R_Q+/-) | lines 406-421 | match |
| Dynamic + static ceilings, dynamic>static | lines 436-456 | match |
| Independent (Mathematica) re-derivation | none | missing (see F1) |

Set `paper_alignment: partial` — the script faithfully covers every paper-side deliverable and matches the appendix .tex exactly; the partial flag is for the two notes-side coefficient typos (F2/F3) and the absent second engine (F1).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 191 | `simplify(Xi_transfer - 2*(pi1-delta1)) == 0` | split identity | yes |
| A2 | sympy | 199-200 | `simplify(got - expected) == 0` (pi1 rows) | exact pi1 row | yes |
| A3 | sympy | 201-202 | `simplify(got - expected) == 0` (delta1 rows) | exact delta1 row | yes |
| A4 | sympy | 217-247 | `rank()==N`, `len(nullspace())==M` (4 blocks) | rank/nullity counts | yes |
| A5 | sympy | 249-254 | `simplify(det_reduced - expected_det) == 0` | reduced determinant != 0 | yes |
| A6 | sympy | 279-282 | `assert_close` v_num/v_den vs numeric targets | unit generators | yes |
| A7 | sympy | 292-301 | `assert_close` pi/delta/Xi projections incl `Xi=-2delta`,`Xi=2pi` | generator scalars | yes |
| A8 | sympy | 369-373 | `assert_close` wall poles, R_Q, P_0 | pole census | yes |
| A9 | sympy | 406-421 | `assert_close` p0_slope (=Xi), R_Q slopes; `abs(wall_slope)<5e-5` | log-slopes | yes |
| A10 | sympy | 436-456 | `assert_close` ceilings; `assert num_both_dyn>num_both_stat` etc. | dynamic/static ceilings | yes |
| A11 | mathematica | — | (absent) | — | — |

All sympy rows are anchored and non-tautological. Note A7's `pi_num`/`delta_den` "==0" checks are genuine: `v_num` comes from `M_num.nullspace()` (whose 4th row is the `pi1` coefficient row), and the assertion re-projects it through the independently-built `pi_coeff` vector, so it cross-checks nullspace-construction against coefficient-extraction rather than restating a definition.

## Findings

### F1 — missing_verification_script (subtype: missing_mathematica)

**Severity:** medium
**Files:**
- `(missing)` mathematica script for unit 228

**What's wrong:**
No `.wl` script exists for stage 228. The stage card line 11 records "Mathematica audit: none yet." The unit is not status-only (`is_status_only_candidate: False`) and the appendix marks the stage `ExactClosure` for the pure-transfer identities and rigidity splits — so it carries substantive, exact symbolic claims a second engine can independently verify. Every one of the nine deliverables is reachable with native Mathematica primitives that differ in kind from the SymPy route (notably the quartic pole census: SymPy uses `np.roots` on a numeric companion matrix, whereas Mathematica can use exact `Solve`/`Roots` or `NSolve` with a different algorithm). The dual-engine rule applies because verification is POSSIBLE here, not merely because it is desirable.

**Why this matters:**
Without a second engine, the `ExactClosure` claims (split identity, exact coefficient rows, rank/nullity, reduced determinant) rest on a single SymPy implementation. The two notes-side typos found below (F2/F3) show how a transcription slip can propagate; an independent re-derivation is the cross-check that would catch a SymPy-side error.

**Required change:**
Add an independent Mathematica audit (see directive F1 claim manifest M1-M9). It must re-derive each claim from the physical premises using native primitives and a different decomposition — NOT a line-by-line port of the `.py`.

**Verification:**
`redteam exec-mathematica 228` runs the new `.wl`, it exits 0, and all `expectZero`/`expect_close` checks pass; the verifier confirms it is not a transliteration.

### F2 — paper_misalignment (subtype: notes_contradicts_script)

**Severity:** medium
**Files:**
- notes `:151-152` vs script `:194` and output `:14`

**What's wrong:**
The notes state the `delta_1` row coefficient on `x_OmU` and `x_OmW` as `247π²/(98π²−25)`:
> notes line 151: `+\frac{247\pi^2}{98\pi^2-25}x_{\Omega_U},`
> notes line 152: `+\frac{247\pi^2}{98\pi^2-25}x_{\Omega_W}.`

The script (and its committed output) use `196π²`:
> script line 194: `196 * sp.pi**2 / (98 * sp.pi**2 - 25), 196 * sp.pi**2 / (98 * sp.pi**2 - 25)]`
> output line 14: `delta_1= ... + 196*pi**2*xOU/(-25 + 98*pi**2) + 196*pi**2*xOW/(-25 + 98*pi**2)`

I checked the math by hand. `Delta = OmU^2 OmW^2 − R^2`, `R = kappa*lamR`. The eps-derivative at 0 gives `Delta_01 = 2 OmU^2 OmW^2 (x_OU + x_OW) − 2 R^2 x_LR`. On the sample branch `OmU^2 OmW^2 = 49/25`, `R^2 = 1/(2π²)`, so `Delta = (98π²−25)/(50π²)`. The `x_OU` coefficient of `delta_1 = Delta_01/Delta` is `2*(49/25)/Delta = (98/25)*(50π²/(98π²−25)) = 196π²/(98π²−25)`. The script value `196π²` is correct; the notes value `247π²` is wrong. This is the same notes-side typo pattern flagged for sibling stages 222/223. The appendix .tex does NOT carry this coefficient (it states only the abstract `delta_1:=Delta_01/Delta`), so the public paper is unaffected; only the notes .md is wrong.

**Why this matters:**
The notes are the derivation source of record. A reader cross-checking the script against the notes would see a contradiction and could "correct" the script to the wrong value.

**Required change:**
`## Resolve before fix_loop` — user-routed. Do NOT auto-edit.

**Verification:**
After user resolution, notes lines 151-152 read `196π²` (matching the verified script), or the user confirms an alternative derivation.

### F3 — paper_misalignment (subtype: notes_contradicts_script)

**Severity:** medium
**Files:**
- notes `:196` vs script `:250-252` and output `:21`

**What's wrong:**
The notes give the reduced determinant as
> notes line 196: `\frac{247(251+215\pi^2)(80000+343225\pi^2+43218\pi^4)}{475(8670000+14894275\pi^2+2117682\pi^4)}`

The script (and committed output) give
> script lines 250-252: `196 * (200 + 147 * sp.pi**2) * (80000 + 343225 * sp.pi**2 + 43218 * sp.pi**4) / (475 * (8670000 + 14894275 * sp.pi**2 + 2117682 * sp.pi**4))`
> output line 21: `det[...] = 196*(200 + 147*pi**2)*(80000 + 343225*pi**2 + 43218*pi**4)/(475*(8670000 + 14894275*pi**2 + 2117682*pi**4))`

The two leading factors differ: notes `247(251+215π²)` vs script `196(200+147π²)`. Numerically these are not equal (`247*251=61997` vs `196*200=39200`; ratio ~1.81 at π²≈9.87), so it is a genuine mismatch, not a regrouping. The determinant is produced symbolically by SymPy and the script asserts equality with the literal target at line 254 (`assert simplify(det_reduced - expected_det) == 0`) — so the script's form is the machine-derived one and the notes' `247(...)` is the typo companion of F2 (same `247`-vs-`196` slip). The appendix .tex does not carry this determinant, so the public paper is unaffected.

**Why this matters:**
Same as F2: the notes are the source of record and would mislead a manual cross-check of the determinant's nonvanishing factor.

**Required change:**
`## Resolve before fix_loop` — user-routed. Do NOT auto-edit.

**Verification:**
After user resolution, notes line 196 reads `196(200+147π²)...` (matching the verified script), or the user confirms an alternative derivation.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration cannot be assessed. The directive's claim manifest (F1) requires native-primitive independence and explicitly forbids porting the `np.roots` choreography.

## Engine cross-check

Only one engine present; n/a.

## Verdict justification

The SymPy script holds up against the paper. I attacked: (1) the split identity — it is a genuine symbolic `simplify(...)==0`, not tautological, derived from independent `P01`/`Delta01`/`N01` differentiations; (2) the coefficient rows — I re-derived `delta_1`'s `x_OU` coefficient by hand and the script's `196π²/(98π²−25)` is correct; (3) the generator "==0" projections — they cross-check nullspace construction against coefficient extraction, not a restated definition; (4) the dynamic>static inequalities — only the finite cases are asserted (`num_nonempty` is correctly `inf` and omitted), all consistent with the carried Stage-241 budgets `0.367930328492646`/`0.737619063660757`, which I confirmed are identical upstream carry-forwards used in sibling stage 227 (lines 271-272), not values invented here. The `xi_val = (...)*0 + (...)` construct at line 268 has dead `*0` scaffolding but does not produce a false pass (orientation correctly uses the Xi sign, independently re-asserted at lines 294/298). The two notes-side typos (`247`-vs-`196` in F2/F3) are the only correctness discrepancies and they sit in the notes .md, not the script or the appendix .tex — the script and public paper are correct, so these route to the user as `notes_contradicts_script` rather than blocking the script. The verdict is `findings`: F1 (missing second engine, script-side, directive-actionable) plus F2/F3 (notes typos, user-resolution). No `stop_cold`: the math is internally consistent and self-correcting against the notes typos.

## Self-test notes

I checked variable-independence (the `sp.diff(..., eps).subs(eps,0)` derivatives all act on `eps`-dressed expressions that genuinely depend on `eps`, so no identically-zero derivative trap), the hand-derivation of `delta_1`'s coefficient (confirmed `196π²`, exposing the notes typo), the determinant numeric inequality (confirmed `247(251+215π²) != 196(200+147π²)`), upstream provenance of the budget literals (matched in stage 227), and the ceiling inequalities (only finite cases asserted, `num_nonempty=inf` correctly excluded). No self-test trap fires against the script; the directive's F1 claim manifest passes the path/primitive-independence self-test.
