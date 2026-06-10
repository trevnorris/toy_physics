---
unit_id: 228
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: ["moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 228 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_228.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows/sections at lines 68, 747-804, 837)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Numerator/denominator split, dynamic-window audit, and first branch-level comparison between the static $\Xi_1$ budget and resonant dynamic survival." The card carries the split `\Xi_1=2(\pi_1-\delta_1)` as an Input/structure and routes the detailed theorem-block formulas to the Part VII appendix. The appendix (eq. `app-part07-numerator-denominator-split`, theorem `thm:app-part07-pure-transfer-sieve`) states the deliverables explicitly: (1) the split `\Xi_1=2(\pi_1-\delta_1)` with `\pi_1:=P_{01}/P`, `\delta_1:=\Delta_{01}/\Delta`; (2) the rigid-subcorridor theorem (pure-transfer + `\pi_1=0` → 1D numerator-rigid survivor; pure-transfer + `\delta_1=0` → 1D denominator-rigid survivor; both → trivial); (3) both 1D survivors pass the first wall-like dynamic-window audit on the compatibility point; (4) the verdict that the first active ceiling is the transported static `|\varepsilon\Xi_1|` budget, not the dynamic window. The notes add the full numerics: the exact `\pi_1` and `\delta_1` row coefficients (incl. the corrected `\delta_1` coefficient `196\pi^2/(98\pi^2-25)` on `\Omega_U,\Omega_W` and the corrected reduced determinant `196(200+147\pi^2)(...)`), rank/nullity counts, the positive-`\Xi_1` unit generators `v_num`/`v_den`, the wall-pole census (`\omega_\pm`, `R_{Q,\pm}`, `P_0`), the first-order log-slopes, and the dynamic vs. transported-static ceilings.

## What the script claims to verify

Both scripts (M1-M9 / Sections 2-7) verify, on the explicit Stage-223 compatibility branch: the split identity `\Xi_1 - 2(\pi_1-\delta_1) = 0`; the exact `\pi_1` row `{3/19,16/19,3/19,32/19,0}` and `\delta_1` row `{0,0,50/(25-98\pi^2),196\pi^2/(98\pi^2-25),196\pi^2/(98\pi^2-25)}`; the rank/nullity ladder (3/2, 4/1, 4/1, 5/0) for pure-transfer and the three rigidity-augmented matrices; the exact reduced determinant `196(200+147\pi^2)(80000+343225\pi^2+43218\pi^4)/(475(8670000+14894275\pi^2+2117682\pi^4))` and its nonvanishing; the oriented positive-`\Xi_1` unit generators `v_num`/`v_den` and the identities `\Xi_1=-2\delta_1` (num) / `\Xi_1=2\pi_1` (den); the wall-pole census; the first-order dynamic log-slopes of `P_0`, `R_{Q,\pm}`, `\omega_\pm`; the dynamic ceilings and the transported-static ceilings, with the strict comparison that the dynamic ceiling exceeds the transported static ceiling on both rigid splits. The bottom-line claims are exercised by hard `assert ... == 0` (SymPy) and `expectZero`/`expectTrue`/`expectClose` with `Exit[1]` on failure (Mathematica).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `\Xi_1=2(\pi_1-\delta_1)` split | py:191 `assert simplify(Xi_transfer - 2*(pi1-delta1))==0`; wl:158 M1 | match |
| exact `\pi_1`, `\delta_1` row coefficients (notes §2) | py:199-202 vs `expected_pi`/`expected_delta` (194); wl:168-169 M2/M3 | match |
| rigid-subcorridor theorem (1D/1D/trivial) | py:217-247 rank/nullity asserts; wl:186-193 M4 | match |
| reduced determinant `196(200+147\pi^2)(...)` ≠ 0 | py:249-254; wl:199-206 M5 | match |
| positive-`\Xi_1` generators `v_num`/`v_den`, `\Xi_1=-2\delta_1`/`2\pi_1` | py:273-301; wl:214-234 M6 | match |
| wall-pole census `\omega_\pm,R_{Q,\pm},P_0` | py:366-373; wl:251-265 M7 | match |
| first-order dynamic slopes (`P_0`,`R_{Q,\pm}`) | py:403-421; wl:305-333 M8 | match |
| dynamic vs. transported-static ceilings + comparison | py:433-456; wl:355-375 M9 | match |

`paper_alignment: aligned`. Every appendix/notes deliverable maps to a substantive, non-tautological script check in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 191 | `simplify(Xi - 2(pi1-delta1))==0` | split identity | yes |
| A2 | sympy | 199-202 | per-coeff `simplify(got-expected)==0` | exact `\pi_1`,`\delta_1` rows | yes |
| A3 | sympy | 217-247 | `rank()==`/`len(nullspace())==` | rigidity theorem ranks | yes |
| A4 | sympy | 254 | `simplify(det_reduced-expected_det)==0` | reduced determinant | yes |
| A5 | sympy | 279-301 | `assert_close` generators+ identities | positive-`\Xi_1` survivors | yes |
| A6 | sympy | 369-373 | `assert_close` pole census | wall-pole census | yes |
| A7 | sympy | 406-421 | `assert_close` dynamic slopes | first-order dynamic slopes | yes |
| A8 | sympy | 436-456 | `assert_close`+`assert >` ceilings | dynamic vs static ceilings | yes |
| M1 | mathematica | 158 | `expectZero[xi1-2(pi1-delta1)]` | split identity | yes |
| M2/M3 | mathematica | 168-169 | `expectZero[Total[(row-expected)^2]]` | exact rows | yes |
| M4 | mathematica | 186-193 | `expectZero[rank-N]`/`nullity-N` | rigidity ranks | yes |
| M5 | mathematica | 205-206 | `expectZero[det-expected]`+`!=0` | reduced determinant | yes |
| M6 | mathematica | 219-234 | `expectVectorClose`/`expectZero`/`expectClose` | survivors+identities | yes |
| M7 | mathematica | 261-265 | `expectClose` pole census | wall-pole census | yes |
| M8 | mathematica | 321-333 | `expectClose`/`expectSmall` slopes | dynamic slopes | yes |
| M9 | mathematica | 365-375 | `expectClose`/`expectTrue` ceilings | ceilings + comparison | yes |

No row traces to "none." No orphaned scaffolding. The `\pi_1`/`\delta_1` rows and the reduced determinant are DERIVED in-script from the dressed-primitive linearization (py:136-186 / wl:97-151) and only THEN compared to the expected literals; the expected literals are the paper-side claim, not a self-referential definition — so A2/A4/M2/M3/M5 are non-tautological. (If the derivation produced a different coefficient, the comparison would fail.)

## Findings

None. Zero findings.

## Independent-derivation check (Mathematica)

The `.wl` is a GENUINELY INDEPENDENT re-derivation, not a transliteration. Three load-bearing route differences:

1. **Linearization.** SymPy uses `sp.diff(EXPR, eps).subs(eps,0)` (py:136-141). The `.wl` uses `Coefficient[Normal[Series[EXPR /. growthRules, {eps,0,1}]], eps, 1]` (wl:97-107) — a Taylor-series extraction, not symbolic differentiation.
2. **Pole census.** SymPy builds `sp.Poly`, lambdifies the quartic coefficients, and finds roots with `numpy.roots` on numeric coefficients in `y=\omega^2` (py:320-321, 350). The `.wl` uses `NSolve[poly==0, y, WorkingPrecision->60]` on the substituted polynomial (wl:244-249) — a different root finder at 60-digit precision.
3. **Dynamic slopes (the decisive divergence).** SymPy computes log-slopes by **symmetric finite difference** (`step=1e-8`, evaluate poles/`R_Q`/`P_0` at `±step`, take `(log(p)-log(m))/(2*step)`; py:383-393). The `.wl` computes them **analytically via the implicit-function theorem**: `dy = -(∂_ε F)/(∂_y F)` and `rqSlope = ∂_ε log(R_Q) + ∂_y log(R_Q)·dy`, all evaluated symbolically at `eps=0` (wl:282-303). These are fundamentally different mathematics, not the same algebra in two syntaxes.

Corroborating numerical evidence of independence: the two engines' `R_{Q,-}` numerator-rigid slope is `0.7135849466877175` (SymPy, finite-diff) vs `0.71358483642759108` (`.wl`, analytic) — agreement to ~7 digits with divergence in the 8th, exactly the finite-difference truncation-error signature a transliteration could not produce. Both scripts are self-contained (no `Get`/`Needs`/`Import`; SymPy imports only `math`/`numpy`/`sympy`). **INDEPENDENCE CALL: independent.**

## Engine cross-check

Both saved outputs end with the success banner and all PASS lines. The symbolic deliverables agree exactly:
- `\delta_1` row: SymPy out:6 `196*pi**2*xOU/(-25+98*pi**2)+196*pi**2*xOW/(-25+98*pi**2)`; Mathematica out:5 `(196*Pi^2*xOU)/(-25+98*Pi^2)+(196*Pi^2*xOW)/(-25+98*Pi^2)`. Identical. (The `.wl` "delta_1 row" print at out:9 shows the equivalent reduced form `2+50/(-25+98*Pi^2)` = `196\pi^2/(98\pi^2-25)`; M3 squared-difference test passes as 0.)
- Reduced determinant: SymPy out:13 and Mathematica out:31 both `196*(200+147*Pi^2)*(80000+343225*Pi^2+43218*Pi^4)/(475*(8670000+14894275*Pi^2+2117682*Pi^4))`. Identical.
- Numeric deliverables (`v_num`/`v_den`, pole census, slopes, ceilings) agree to ~7-8 digits, with the expected finite-diff-vs-analytic divergence on the dynamic slopes. No engine disagreement.

## Verdict justification

CLEAN. I attacked: (a) the corrected-value reconfirmation — both engines derive and assert `\delta_1` coefficient `196\pi^2/(98\pi^2-25)` and reduced determinant `196(200+147\pi^2)(...)`; no surviving `247`, `247\pi^2/(98\pi^2-25)`, `247(251+215\pi^2)`, `251+215\pi^2`, or `215\pi` in any stage-228 script, output, notes, or `.tex` (grep-confirmed empty); the notes read the corrected values at lines 152 and 196. (b) Tautology — the row coefficients and determinant are derived from the dressed-primitive linearization before being compared to the literals, so the asserts can fail. (c) Independence — confirmed three distinct mathematical routes plus a finite-diff truncation-signature in the numerics. (d) Symbol domains — both engines declare the mixed-drift `x*` symbols Real and the physical primitives positive, consistent with the setup; the positivity is used only in `FullSimplify`/`sqrt` normalization and is justified. (e) Paper alignment — every appendix/notes deliverable has a substantive matching check in both engines, and the published card/appendix carry only the high-level split/theorem/verdict (they never carried the explicit `\delta_1` coefficient or the determinant literal), so they are correctly unaffected by the notes-only value corrections. (f) Output freshness — SymPy `.py` mtime May 11 11:58, its `.txt` Jun 2 17:26 (output newer → not stale); `.wl` Jun 2 17:20, its `.txt` Jun 2 17:26 (fresh). I read paper card, notes, and appendix, and the scripts' verified claims match the paper's claims.

## Value Reconciliation (pass-2 augmentation)

Authoritative record = script source + committed `.txt` outputs (read-only; nothing executed). Both outputs are fresh.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| split identity `\Xi_1=2(\pi_1-\delta_1)` | py:205 / wl:154, py-out:4, wl-out:3 | tex card l.9; appendix l.785; notes l.130 | MATCH |
| `\pi_1` row `{3/19,16/19,3/19,32/19,0}` | py-out:5, wl-out:4 | notes l.140-146 | MATCH |
| `\delta_1` row, coeff `50/(25-98\pi^2)` on `\lambda_R` | py-out:6, wl-out:5 | notes l.150 | MATCH |
| `\delta_1` coeff `196\pi^2/(98\pi^2-25)` on `\Omega_U`,`\Omega_W` (CORRECTED) | py:194, wl:164, py-out:6, wl-out:5 | notes l.151-152 | MATCH |
| reduced determinant `196(200+147\pi^2)(...)/(...)` (CORRECTED) | py-out:13, wl-out:31 | notes l.196 | MATCH |
| rank/nullity ladder 3/2,4/1,4/1,5/0 | py-out:9-12, wl-out:15-30 | notes l.166-190; appendix l.792-796 | MATCH |
| `v_num` ≈ (-0.55551149,0.31814576,-0.65766801,-0.04533730,-0.39447126) | py-out:16, wl-out:36 | notes l.230-235 | MATCH |
| `v_den` ≈ (-0.26583993,0.18448137,0.94454459,0.04984499,-0.02543112) | py-out:20, wl-out:38 | notes l.254-260 | MATCH |
| `\delta_1(v_num)≈-0.86805617`, `\Xi_1(v_num)≈1.73611235` | py-out:18-19, wl-out:42-44 | notes l.241-243 | MATCH |
| `\pi_1(v_den)≈0.34646608`, `\Xi_1(v_den)≈0.69293215` | py-out:21-23, wl-out:48-50 | notes l.266-268 | MATCH |
| `\omega_-≈1.99753568`, `\omega_+≈4.94905432` | py-out:26-27, wl-out:57 | notes l.296,302 | MATCH |
| `R_{Q,-}≈30.19990756`, `R_{Q,+}≈36.17118648` | py-out:26-27, wl-out:58 | notes l.298,303 | MATCH |
| `P_0≈0.00206979232` | py-out:28, wl-out:59 | notes l.307 | MATCH |
| num dynamic slopes `R_{Q,+}≈-0.52346582`,`R_{Q,-}≈+0.71358484` | py-out:32, wl-out:71 | notes l.319,321 | MATCH |
| den dynamic slopes `R_{Q,+}≈-0.35245541`,`R_{Q,-}≈-0.23169484` | py-out:34, wl-out:74 | notes l.339,341 | MATCH |
| num dynamic ceilings `(0.96253269, ∞)` | py-out:37, wl-out:100 | notes l.370,375 | MATCH |
| den dynamic ceilings `(1.39592653, 1.42955095)` | py-out:38, wl-out:101 | notes l.383,388 | MATCH |
| num static ceilings `(0.21192772, 0.42486828)` | py-out:39, wl-out:102 | notes l.400,405 | MATCH |
| den static ceilings `(0.53097598, 1.06448959)` | py-out:40, wl-out:103 | notes l.412,417 | MATCH |
| `R_{Q,req}≈21.85456630` (10%-loss threshold) | py:432, wl:351 | notes l.361 | MATCH |
| `K_compat≈24.473754879` | py:166 (exact), wl-out:2 (exact) | notes l.79 | MATCH (decimal carrier in notes) |
| `\kappa=2\sqrt2/\pi` | py:144, wl:117 | notes l.70 | MATCH |

INTERNAL (scaffolding, no finding): transported static budgets `budget_both=0.367930328492646`, `budget_nonempty=0.737619063660757` (py:441-442 / wl:352-353 — carried-in Stage-224 transport used to compute the static ceilings, which ARE reported); finite-difference `step=1e-8` (py:383); `WorkingPrecision`/precision tags; `eps=0` linearization point; pass/fail flags and residual `diff` values; the four intermediate `D01/D21/D41/N01/P01/Delta01` linearized carriers; the omega-slope near-zero magnitudes (verification scaffolding, asserted `<5e-5`).

reconciliation: complete; 23 deliverable values checked, 0 misaligned.

## Self-test notes

Checked: (1) variable independence — the `.wl` `D[fBranch, eps]` and `D[Log[rqBranch], eps]` (wl:288-291) operate on `growthRules`-dressed expressions that genuinely depend on `eps`, so the slopes are not identically zero (confirmed nonzero in wl-out:71,74); (2) the reduced determinant and row-coefficient asserts are derivation-then-compare, not define-then-assert, so non-tautological; (3) trivial-case parity not applicable (no symmetric-domain integrals here — the dynamic slopes are pole-census derivatives, not integrals). Concluded no self-test trap is tripped; no directive written (zero findings).
