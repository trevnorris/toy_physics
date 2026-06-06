---
unit_id: 096
batch: IV.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage096_geometry_lane_check_verdict.md]
  paper_appendix: present
---

# Audit unit 096 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_096.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage096_geometry_lane_check_verdict.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows for MTDC-T8.1 read; this stage is `\input{stages/stage_096}` at line 1226, anchor MTDC-T8.1)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.txt`

## What the paper claims

Stage 096 is a "geometry-lane firewall" verdict / derivation-ledger step. The card's `Derivation ledger` states the computation "isolates the forced conservative carrier \(\widehat Y_Q^{\rm cons}=3/4+(1/4)/(1-\omega^2/\Omega_Q^2)\) or the obstruction variables \((\epsilon_2,\epsilon_4)\)." The card's `Checks` enumerate three deliverables: (1) "Check the static limit \(\epsilon_2=\epsilon_4=0\) returns \(c_{\rm pole}=1/4\)"; (2) "Check \(l=0\) and \(l=2\) orthogonality before applying the geometry firewall"; (3) any support/source success statement still carries the minimal-module hypothesis. The notes elaborate: on the isotropic branch the scalar/geometry \(l=0\) lane and the grouped real \(l=2\) bundle are block diagonal (proved upstream at Stage 77), so the contamination numbers `eps_2 = eps_4 = 0`; the obstruction formula then collapses to `c_pole = 1/4`, `c_geom = 3/4`, giving `Yhat_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)`, and the carried support/source identification reproduces `rho_alpha = 4/3`, `zeta_req = 1/3`. The card explicitly marks the result as conditional (status tag, not an unconditional theorem). Distinct deliverables: the orthogonality/decoupling premise (l=0 ⊥ l=2), eps_2=eps_4=0, c_pole=1/4, c_geom=3/4, Yhat_Q^cons closed form, rho_alpha=4/3, zeta_req=1/3.

## What the script claims to verify

The SymPy docstring states three goals: (1) the isotropic wall branch keeps the l=0 lane orthogonal to the grouped real l=2 bundle so the contamination numbers vanish; (2) the obstruction formula collapses to 3/4 + 1/4; (3) the carried contact-plus-pole identification reproduces rho_alpha=4/3 and zeta_req=1/3. SECTION I computes, for each of the five grouped real l=2 harmonics (Y20, Y21c, Y21s, Y22c, Y22s): the L²(S²) overlap with Y00 (claims =0), the spherical-Laplacian eigenvalue residual `(-Delta)Y2m - 6 Y2m` (claims =0), and the overlap of Y00 with `(-Delta)Y2m` (claims =0). SECTION II sets `eps_2 = eps_4 = 0` as integer literals, evaluates `c_pole = (1+eps_4)/(4(1+eps_2)^2)`, `c_geom = 1 - c_pole`, `rho_alpha = 1/c_geom`, `zeta_req = c_pole/c_geom`, `Yhat_cons = c_geom + c_pole/(1-omega^2/Omega_Q^2)`, and asserts each equals its target literal (1/4, 3/4, 4/3, 1/3, and the closed form). The Mathematica `.wl` mirrors this structure exactly.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| l=0 ⊥ l=2 orthogonality (Check 2) | SECTION I `<Y00\|Y2m> = 0` for all five 2m harmonics (both engines) | match |
| l=2 Laplace eigenvalue ℓ(ℓ+1)=6 | SECTION I `(-Delta)Y2m - 6 Y2m = 0` for all five (both engines) | match |
| firewall: geometry lane stays orthogonal under wall operator | SECTION I `<Y00\|(-Delta)Y2m> = 0` for all five (both engines) | match |
| eps_2 = eps_4 = 0 (carry-forward from Stage 77) | hardcoded `eps_2 = sp.Integer(0)`, `eps_4 = sp.Integer(0)` — asserted, not derived | partial (carry-forward, by design) |
| static limit returns c_pole = 1/4 (Check 1) | SECTION II `c_pole - 1/4 == 0` (both) | match (but tautological at eps=0 — see F1) |
| c_geom = 3/4 | SECTION II `c_geom - 3/4 == 0` (both) | match (tautological — see F1) |
| Yhat_Q^cons = 3/4 + (1/4)/(1-ω²/Ω_Q²) | SECTION II `Yhat_cons - Yhat_expected == 0` (both) | match (tautological — see F1) |
| rho_alpha = 4/3 | SECTION II `rho_alpha - 4/3 == 0` (both) | match (tautological — see F1) |
| zeta_req = 1/3 | SECTION II `zeta_req - 1/3 == 0` (both) | match (tautological — see F1) |

`paper_alignment: aligned` — every paper deliverable maps to a script-side check that targets the correct identity/value; no value or target mismatch. The weakness is verification *strength* in SECTION II (F1), not paper alignment.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 82 | `expect_zero(<Y00\|Y2m>)` ×5 | orthogonality (Check 2) | yes |
| A2 | sympy | 83 | `expect_zero((-Delta)Y2m - 6 Y2m)` ×5 | l=2 Laplace eigenvalue | yes |
| A3 | sympy | 84 | `expect_zero(<Y00\|(-Delta)Y2m>)` ×5 | firewall decoupling under wall op | yes |
| A4 | sympy | 109 | `expect_zero(c_pole - 1/4)` | static-limit c_pole (Check 1) | no (tautological at eps=0) |
| A5 | sympy | 110 | `expect_zero(c_geom - 3/4)` | c_geom | no (tautological) |
| A6 | sympy | 111–114 | `expect_zero(Yhat_cons - Yhat_expected)` | Yhat closed form | partial (ω-dependence real; coefficients tautological) |
| A7 | sympy | 115 | `expect_zero(rho_alpha - 4/3)` | rho_alpha | no (tautological) |
| A8 | sympy | 116 | `expect_zero(zeta_req - 1/3)` | zeta_req | no (tautological) |
| A9 | mathematica | 50 | `expectZero[<Y00\|Y2m>]` ×5 | orthogonality | yes |
| A10 | mathematica | 51 | `expectZero[(-Delta)Y2m - 6 Y2m]` ×5 | Laplace eigenvalue | yes |
| A11 | mathematica | 52 | `expectZero[<Y00\|(-Delta)Y2m>]` ×5 | firewall decoupling | yes |
| A12 | mathematica | 76–80 | `expectZero[c_pole-1/4 …]` ×5 | SECTION II values | no/partial (tautological at eps=0, mirrors A4–A8) |

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py:88-116`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.wl:57-80`

**What's wrong:**
SECTION II hardcodes the contamination numbers as integer literals (`eps_2 = sp.Integer(0)`, `eps_4 = sp.Integer(0)`; `.wl` `eps2 = 0; eps4 = 0`) and then evaluates the obstruction formula `c_pole = (1+eps_4)/(4*(1+eps_2)^2)` at exactly that point. With eps_2=eps_4=0 the formula degenerates to the literal `1/4` *by construction*, and everything downstream (`c_geom = 1 - c_pole = 3/4`, `rho_alpha = 1/c_geom = 4/3`, `zeta_req = c_pole/c_geom = 1/3`, and the ω-coefficients of `Yhat_cons`) is then fixed arithmetic. The five SECTION II `expect_zero` assertions (`c_pole - 1/4 == 0`, `c_geom - 3/4 == 0`, `rho_alpha - 4/3 == 0`, `zeta_req - 1/3 == 0`, and the constant part of `Yhat_cons - Yhat_expected`) therefore cannot fail no matter what the obstruction formula's structure is: they confirm `(1+0)/(4·(1+0)^2) = 1/4` etc. The load-bearing object the card names — the obstruction formula relating `(eps_2, eps_4)` to `c_pole` — is never exercised at a nonzero point, so the script never demonstrates that the formula actually *reduces correctly* to the static module from the general (contaminated) form. The only ω-dependent piece of A6/A12 (the `1/(1-ω²/Ω_Q²)` factor) is real and non-tautological; the numeric coefficients sitting in front of it are not.

This is a CHECKPOINT (is_checkpoint: True), which carries the higher bar: assertions must be substantive and non-tautological. SECTION I genuinely meets that bar (the orthogonality integrals and the Laplace eigenvalue residual are real computations that would fail for a wrong harmonic or wrong eigenvalue). SECTION II does not.

Note this is NOT `paper_misalignment`: the card's Check 1 ("Check the static limit \(\epsilon_2=\epsilon_4=0\) returns \(c_{\rm pole}=1/4\)") legitimately asks for the eps→0 evaluation, and the values all reconcile. The defect is that the script implements that check as a literal substitution rather than a symbolic limit of the general formula, leaving the formula's structure unverified at the checkpoint bar.

**Why this matters:**
At a checkpoint the eps→0 collapse is the actual physics deliverable (geometry lane is inert ⟹ contamination numbers vanish ⟹ the conservative module is forced). Verifying it only at the degenerate point where the formula reduces to a constant gives no protection against a transcription error in the obstruction formula itself (e.g., a wrong power on `(1+eps_2)`, a wrong factor of 4, or eps_2/eps_4 swapped) — all such errors are invisible at eps_2=eps_4=0. A substantive check must exercise the formula with eps_2, eps_4 as free symbols and confirm both (a) the static limit and (b) that the general structure is the documented `(1+eps_4)/(4(1+eps_2)^2)`.

**Required change:**
De-tautologize SECTION II in both engines by keeping `eps_2, eps_4` as free symbols and asserting the structural identity, then taking the static limit. Concretely (SymPy): introduce `eps_2, eps_4 = sp.symbols("eps_2 eps_4", real=True)`, define `c_pole_gen = (1 + eps_4) / (4*(1 + eps_2)**2)`, and add a non-tautological structural assertion that the general loading ratios match the documented composite form, e.g. assert `simplify(c_pole_gen*(4*(1+eps_2)**2) - (1+eps_4)) == 0` is NOT how to do it (that is tautological too) — instead assert the *physical* reduction: `expect_zero("c_pole|eps=0 - 1/4", c_pole_gen.subs({eps_2:0, eps_4:0}) - sp.Rational(1,4))` AND a can-fail sensitivity check that the formula is genuinely eps-dependent, e.g. `c_pole_gen.subs({eps_2:0, eps_4:1})` must equal `sp.Rational(1,2)` (a literal that differs from 1/4, proving eps_4 actually enters) and `c_pole_gen.subs({eps_2:1, eps_4:0})` must equal `sp.Rational(1,16)` (proving eps_2 enters with the documented power). Then derive c_geom/rho_alpha/zeta_req/Yhat from `c_pole_gen` at eps=0 as now. Mirror in `.wl`. This makes the obstruction-formula structure load-bearing (a wrong power or factor now fails a concrete literal) while keeping the eps→0 deliverable values exactly as the card states.

**Verification:**
After the fix, the SymPy output should show new lines confirming the eps-sensitivity literals (e.g. `c_pole|eps_4=1 = 1/2`, `c_pole|eps_2=1 = 1/16`) in addition to the existing `c_pole - 1/4 = 0` etc., and the `.wl` output the Mathematica equivalents. Both scripts must still exit 0 and all the card's stated values (1/4, 3/4, 4/3, 1/3, Yhat closed form) must continue to reconcile.

## Independent-derivation check (Mathematica)

The `.wl` is a close structural mirror of the `.py`: identical Y2m definitions (lines 39–44 vs py 58–63), an identically-shaped `lapS2` (line 28–30 vs py `lap_s2` 46–50), an identically-shaped `dOmegaIntegral` (line 32–34 vs py `domega` 39–43), the same five-harmonic Do-loop over the same three checks, the same `eps2=0; eps4=0` hardcoding, and the same five SECTION II assertions in the same order with the same target literals. This is borderline `mathematica_transliteration` — same variable choreography and same intermediate steps. I am NOT raising it as a separate finding because: (a) the underlying objects (spherical harmonics, the spherical Laplacian, the L²(S²) inner product) have essentially canonical definitions, so two independent authors would converge on nearly identical SECTION I code regardless of porting; (b) the engines use genuinely different symbolic backends (SymPy `integrate`/`simplify` vs Mathematica `Integrate`/`FullSimplify`), so the orthogonality and eigenvalue results constitute a real cross-engine corroboration of the same physics; and (c) the F1 de-tautologization, when applied symmetrically, will add independent eps-sensitivity content to both. The transliteration concern is noted but does not rise to a finding at this stage's scope.

## Engine cross-check

Both engines produce identical results. SECTION I: all 15 orthogonality/eigenvalue residuals = 0 in both (sympy output lines 5–34; mathematica output lines 9–38, all PASS). SECTION II: `c_pole=1/4`, `c_geom=3/4`, `rho_alpha=4/3`, `zeta_req=1/3` in both (sympy 37–46; mathematica 43–59). The Yhat closed form prints in slightly different but algebraically identical forms — SymPy `(Omega_Q**2 - 3*omega**2/4)/(Omega_Q**2 - omega**2)` (line 39) vs Mathematica `(3 + (1 - omega^2/omegaQ^2)^(-1))/4` (line 45) — and both confirm `Yhat_cons - Yhat_expected = 0`. Engines agree.

## Verdict justification

`verdict: findings` (one finding, medium). What holds up: SECTION I in both engines is substantive and non-tautological — the orthogonality integrals, the ℓ(ℓ+1)=6 Laplace eigenvalue, and the firewall decoupling under the wall operator are real computations that genuinely exercise the card's "l=0 ⊥ l=2 orthogonality" deliverable and the decoupling premise. Paper alignment is exact: every deliverable value (eps=0, c_pole=1/4, c_geom=3/4, rho_alpha=4/3, zeta_req=1/3, Yhat closed form) reconciles with the card and notes. What does not hold up at the checkpoint bar: SECTION II evaluates the load-bearing obstruction formula only at the degenerate eps_2=eps_4=0 point where it collapses to literals, so its five assertions are arithmetic-by-construction and would not catch a structural transcription error in the formula itself. Attacks tried that the script survived: I checked the spherical-harmonic normalizations and the eigenvalue constant (6 = ℓ(ℓ+1), correct), the parity/symmetry of the orthogonality integrands (Y00·Y2m and Y00·(-Δ)Y2m integrate to zero by genuine angular orthogonality, not by an assumption artifact), the symbol domains (theta/phi real, omega/Omega_Q positive — physically justified and not over-strong), and the value reconciliation (all match). The single residual weakness is the SECTION II tautology, addressed by F1. I read the paper card, the notes, and the appendix MTDC-T8.1 rows before reading the scripts; the script's verified claim matches the paper's claim in target and value.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| eps_2 = eps_4 = 0 | py:88-89 / wl:57-58; sympy out:35-36, math out:48-49 | notes:31, card:22 (`\epsilon_2=\epsilon_4=0`) | MATCH |
| c_pole = 1/4 | py:93/103 / wl:59/68; sympy out:37, math out:43 | notes:41, card:22 | MATCH |
| c_geom = 3/4 | py:94/104 / wl:60/69; sympy out:38, math out:44 | notes:43 | MATCH |
| Yhat_Q^cons = 3/4 + (1/4)/(1-ω²/Ω_Q²) | py:95/105 / wl:63/70; sympy out:39, math out:45 | notes:47, card:13 | MATCH |
| rho_alpha = 4/3 | py:100/106 / wl:61/71; sympy out:40, math out:46 | notes:51, card:27 | MATCH |
| zeta_req = 1/3 | py:101/107 / wl:62/72; sympy out:41, math out:47 | notes:53, card:27 | MATCH |
| l=2 Laplace eigenvalue 6 = ℓ(ℓ+1) | py:75/83 / wl:51; sympy out:6 etc., math out:11 etc. | derived in-script from ℓ=2 (docstring py:15); appendix MTDC-T8.1 cites the geometry-lane firewall | INTERNAL (derived scaffolding, not a prose deliverable) |

INTERNAL-only items (verification scaffolding, no finding): the per-harmonic orthogonality residuals `<Y00|Y2m> = 0` and `<Y00|(-Δ)Y2m> = 0` (these are pass/fail checks, not reported deliverable values); the intermediate `Yhat_expected` target form; pass/fail flags.

A note on provenance labeling (not a finding, deferred to the dedicated numbering pass per in-loop policy): the SymPy docstring (py:9) attributes the obstruction formula to "the Stage 092 obstruction formula" while the notes (notes:39) attribute the collapse to "the Stage-75 obstruction formula." This is a cross-reference to a *source* stage, not a self-label of 096 (all of 096's own self-labels correctly read "Stage 096"), so per the confirmed Reading-2 in-loop scope it is left for `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` rather than fixed here.

reconciliation: complete; 6 deliverable values checked, 0 misaligned

## Self-test notes

I checked: (1) Variable independence — the F1 directive's `c_pole_gen.subs(...)` substitutions are not derivatives, and the eps-sensitivity literals (eps_4=1 ⟹ 1/2; eps_2=1 ⟹ 1/16) were hand-evaluated from `(1+eps_4)/(4(1+eps_2)^2)` and differ from 1/4, so the new checks can fail (non-tautological). (2) Symmetry/parity — the orthogonality integrands Y00·Y2m carry net azimuthal/polar parity that makes the S² integral genuinely vanish by angular orthogonality, confirmed against the zero outputs; not an assumption artifact. (3) Trivial-case pre-check — at eps_2=eps_4=0 the general formula gives c_pole=1/4 (matches), and the proposed sensitivity literals give 1/2 and 1/16 (nonzero, distinct), confirming the directive's `assert`s land on the right values. (5) Paper round-trip — the F1 fix keeps all card/notes deliverable values unchanged (it only adds can-fail structural checks around the existing eps→0 evaluation), so it introduces no new paper_misalignment.
