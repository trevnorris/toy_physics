---
unit_id: 076
batch: III.4
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
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md]
  paper_appendix: present
---

# Audit unit 076 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_076.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 130 + `\input` at 270)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.txt`

## What the paper claims

Stage 076 expresses the remaining wall-depth datum `Theta_w` in normalized Family-1 density variables. The `\stagefield{Output}` is "The wall-depth lock \eqref{eq:app-stage076-theta}", whose boxed equation is verbatim `\boxed{\Theta_w=25\lambda_\mu^2\rho_w^2}` (tex:16–17). The notes derive this through a four-step chain: (1) the exact `n=5` GNLS thermodynamic identity `h(rho) = m c_s^2 / 4` (notes §1; also stated in `\stagefield{Inputs}`, tex:7); (2) the local enthalpy lock `mu_* = lambda_mu h_w` with `h_w = m c_{s,w}^2/4`, giving `Theta_w = 4 rho_w^2 mu_*^2/(hbar^2 c_{s,w}^2) = lambda_mu^2 m^2 rho_w^2 c_{s,w}^2/(4 hbar^2)` (notes §2); (3) the Stage-074 healing-lock `ell = hbar/(2 m c_{s,w})`, giving `Theta_w = lambda_mu^2 rho_w^2/(16 ell^2)` (notes §3); (4) the reference-branch convention `ell/a = 1/20` with `a=1` in normalized wall coordinates, so `1/ell^2 = 400` and `Theta_w = 25 lambda_mu^2 rho_w^2` (notes §4–5). The appendix row 130 echoes the single deliverable: `\(\Theta_w=25\lambda_\mu^2\rho_w^2\)`. The single deliverable is the boxed reference-branch identity; the intermediate forms (h-identity, enthalpy-lock form, healing-lock form) are the supporting closed forms.

## What the script claims to verify

Both scripts (docstring py:2–10; banner) verify the four-step chain as four `expect_zero`/`expectZero` residual checks plus one non-tautology guard. They (1) build `P = K rho^n`, derive `U = rho·∫(P/rho^2)dρ = K rho^n/(n-1)`, `h = dU/dρ`, `c_s^2 = (dP/dρ)/m`, specialize to `n=5`, and assert `h - m c_s^2/4 = 0`, with an explicit guard that the same residual is nonzero at `n=3`; (2) solve the enthalpy-lock `mu_* = lambda_mu m c_{s,w}^2/4` and assert `4 rho_w^2 mu_*^2/(hbar^2 c_{s,w}^2)` equals the independently-written closed form `lambda_mu^2 m^2 rho_w^2 c_{s,w}^2/(4 hbar^2)`; (3) substitute `c_{s,w} = hbar/(2 m ell)` and assert the result equals `lambda_mu^2 rho_w^2/(16 ell^2)`; (4) substitute `ell = a/20`, then `a=1`, and assert the normalized form equals `(1/(1/20))^2/16 · lambda_mu^2 rho_w^2 = 25 lambda_mu^2 rho_w^2`. The bottom line the script proves is exactly the boxed `Theta_w = 25 lambda_mu^2 rho_w^2`.

## Paper ↔ script cross-check

| paper-side deliverable | script-side check | status |
|---|---|---|
| `n=5` identity `h = m c_s^2/4` (tex:7, notes §1) | py50 / wl43 `expect_zero("n=5 enthalpy identity", h_n5 - m cs2_n5/4)` + n=3 non-tautology guard py52–55 / wl45–48 | match |
| enthalpy-lock form `Theta_w = lambda_mu^2 m^2 rho_w^2 c_sw^2/(4 hbar^2)` (notes §2) | py67 / wl58 `expect_zero("Theta_w under enthalpy lock", Theta_w - Theta_target)` | match |
| healing-lock form `Theta_w = lambda_mu^2 rho_w^2/(16 ell^2)` (notes §3) | py76 / wl66 `expect_zero("healing-lock reduction", ...)` | match |
| boxed reference form `Theta_w = 25 lambda_mu^2 rho_w^2` (tex:17, notes §4–5, appendix:130) | py91 / wl79 `expect_zero("normalized reference factor", Theta_ref_norm - ref_target)`, with `Theta_ref_norm` printed `= 25 lambda_mu^2 rho_w^2` (out L15) | match |

Every paper-side deliverable maps to a faithful script-side check. `paper_alignment: aligned`. The stage does not claim `alpha_r`, `Upsilon_w`, or `V0` results (those belong to Stage 075); the scripts correctly do not test them, so there is no `script_missing_paper_claim`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 50 | `expect_zero(h_n5 - m·cs2_n5/4)` | n=5 identity | yes |
| A1g | sympy | 52–55 | `assert residual_n3 != 0` (n=3 guard) | n=5 identity (non-tautology) | yes |
| A2 | sympy | 67 | `expect_zero(Theta_w - Theta_target)` | enthalpy-lock form | yes |
| A3 | sympy | 76 | `expect_zero(Theta_w_in_ell - Theta_heal_target)` | healing-lock form | yes |
| A4 | sympy | 91 | `expect_zero(Theta_ref_norm - ref_target)` | boxed `=25...` | partial (see F-note) |
| B1 | mathematica | 43 | `expectZero[hRhoN5 - mpsi·csSqN5/4]` | n=5 identity | yes |
| B1g | mathematica | 45–48 | `If[residualN3===0, fail]` | n=5 identity (non-tautology) | yes |
| B2 | mathematica | 58 | `expectZero[thetaW - thetaTarget]` | enthalpy-lock form | yes |
| B3 | mathematica | 66 | `expectZero[thetaWInEll - thetaHealReduced]` | healing-lock form | yes |
| B4 | mathematica | 79 | `expectZero[thetaRefNorm - refTarget]` | boxed `=25...` | partial (see F-note) |

F-note on A4/B4: `ref_target` is computed from the same `ref_factor = 1/20` that produced `Theta_ref_norm`, so this last assertion is the weakest link — it verifies that the substitution route (`subs ell=a/20`, `a=1`) agrees with the direct-formula route (`(1/ref_factor)^2/16`), both built from `1/20` and `1/16`. It is mildly self-consistent rather than fully independent. But it is NOT a finding: (a) it still exercises the `1/16` factor and the squaring, catching any transcription error in either route; (b) the substantive, non-trivial work is upstream in A3/B3 (the healing-lock reduction that produces `lambda_mu^2 rho_w^2/(16 ell^2)` from the enthalpy-lock form, a real algebraic cancellation), and the `ell = a/20`, `a=1` substitution producing `25` is genuinely printed (`out L14: 25 lambda_mu^2 rho_w^2/a^2`, `out L15: 25 lambda_mu^2 rho_w^2`); (c) `ref_factor = 1/20` is anchored to notes §4 (`ell/a = 1/20`). The chain as a whole non-tautologically delivers the boxed value.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.txt:3` (banner reads `STAGE 59 — EXACT n=5 WALL-DEPTH LOCK`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.txt:3` (banner reads `STAGE 059 — EXACT n=5 WALL-DEPTH LOCK`)

**What's wrong:**
Both committed transcripts predate their scripts and carry the pre-renumber banner. mtimes: sympy `.py` `2026-05-27 02:38`, sympy `.txt` `2026-05-27 02:15`; mathematica `.wl` `2026-05-27 02:38`, mathematica `.txt` `2026-05-27 02:19`. Both `.txt` are ~20 minutes older than their scripts. The current source emits `STAGE 076` (py:30 `banner("STAGE 076 — EXACT n=5 WALL-DEPTH LOCK")`, wl:26 `banner["STAGE 076 — EXACT n=5 WALL-DEPTH LOCK"]`), but the captured transcripts emit `STAGE 59` (sympy out L3) and `STAGE 059` (mathematica out L3). The numeric/symbolic content of the transcripts otherwise matches what the current scripts produce exactly (all closed forms, all `= 0` residuals, the `n=3` nonzero residual `3*K*rho**2/4`, all `PASS` lines), so this is a freshness/label drift, not a content disagreement.

**Why this matters:**
The committed transcript is the auditable record of the stage; a stale banner means the captured run does not correspond to the current source numbering. Low risk because all residual/closed-form content agrees, but the orchestrator's independent re-run should refresh both `.txt` so banners read `STAGE 076`.

**Required change:**
Re-run both scripts and overwrite the committed transcripts so the banner lines read `STAGE 076 — EXACT n=5 WALL-DEPTH LOCK`. No source edit is needed (the source banners already read `STAGE 076`); only the captured output needs refreshing.

**Verification:**
After re-run, `scripts/output/...stage076..._sympy_audit.txt` line 3 reads `STAGE 076 — EXACT n=5 WALL-DEPTH LOCK`; `mathematica/output/...stage076..._mathematica_audit.txt` line 3 reads `STAGE 076 — EXACT n=5 WALL-DEPTH LOCK`. All residual/closed-form/PASS lines unchanged.

## Independent-derivation check (Mathematica)

The `.wl` mirrors the `.py` choreography closely: same EOS `press = kConst*rho^nPoly` (wl:33 ≡ py:36), same energy-density convention `uPerMass = Integrate[press/rho^2]`, `uRho = rho*uPerMass` (wl:34–35 ≡ py:41–42), same `csSq = D[press,rho]/mpsi`, `hRho = D[uRho,rho]` (wl:36–37 ≡ py:43–44), same `n=3` non-tautology guard (wl:45–48 ≡ py:52–55), same enthalpy-lock `Solve` for `muStar` and same hand-written `thetaTarget` (wl:52–56 ≡ py:59–65), same healing `Solve` for `cSw`/`ell` and same `thetaHealReduced` (wl:60–65 ≡ py:69–75), same `refFactor = 1/20`, `ell->a*refFactor`, `a->1` route and same `refTarget` (wl:73–78 ≡ py:83–90). This is structurally a transliteration. However, stage 076's claim is a pure chain of closed-form thermodynamic identities and substitutions, where the only "derivation" freedom is the EOS convention (fixed: `P=K rho^5`) and which carried-forward lock to substitute (fixed upstream: `mu_*=lambda_mu h_w`, `ell=hbar/(2 m c_sw)`, `ell/a=1/20`). Both engines independently execute `FullSimplify`/`simplify`, `Integrate`/`integrate`, `Solve`/`solve`, and `D`/`diff` on each side rather than echoing a cached numeric result, so the algebra is genuinely re-run by each engine. Per the established batch policy for pure identity/substitution stages (cf. stage_072 report), the mirror is acceptable and not raised as a hard `mathematica_transliteration` finding — noted only (`independent_wl: transliteration`).

## Engine cross-check

Both engines agree on every emitted form:

| quantity | sympy (out) | mathematica (out) |
|---|---|---|
| `c_s^2 [n=5]` | `5*K*rho**4/m` | `(5*kConst*rho^4)/mpsi` |
| `h [n=5]` | `5*K*rho**4/4` | `(5*kConst*rho^4)/4` |
| n=5 identity residual | `0` | `0` (PASS) |
| n=3 residual | `3*K*rho**2/4` | `(3*kConst*rho^2)/4` |
| `Theta_w` (enthalpy lock) | `c_sw**2*lambda_mu**2*m**2*rho_w**2/(4*hbar**2)` | `(cSw^2*lambdaMu^2*mpsi^2*rhoW^2)/(4*hbar^2)` |
| `Theta_w` enthalpy residual | `0` | `0` (PASS) |
| `ell` from healing | `hbar/(2*c_sw*m)` | `hbar/(2*cSw*mpsi)` |
| healing-lock residual | `0` | `0` (PASS) |
| `Theta_w` (healing lock) | `lambda_mu**2*rho_w**2/(16*ell**2)` | `(lambdaMu^2*rhoW^2)/(16*ell^2)` |
| `Theta_w` (ref, general a) | `25*lambda_mu**2*rho_w**2/a**2` | `(25*lambdaMu^2*rhoW^2)/a^2` |
| `Theta_w` (normalized) | `25*lambda_mu**2*rho_w**2` | `25*lambdaMu^2*rhoW^2` |
| normalized-factor residual | `0` | `0` (PASS) |

Identical modulo symbol-name conventions (`K`↔`kConst`, `m`↔`mpsi`, `c_sw`↔`cSw`, `lambda_mu`↔`lambdaMu`, `rho_w`↔`rhoW`). `engines_agree: true`.

## Verdict justification

The scripts faithfully and non-tautologically verify exactly what the paper card and notes claim: the `n=5` enthalpy identity (with an explicit `n=3` guard proving the identity is special to `n=5`, not generic), the enthalpy-lock form, the healing-lock reduction, and the boxed reference-branch lock `Theta_w = 25 lambda_mu^2 rho_w^2`. Attacks tried and failed: (1) re-derived `U = K rho^n/(n-1)`, `h = n K rho^{n-1}/(n-1)`, `c_s^2 = n K rho^{n-1}/m` by hand and confirmed `h - m c_s^2/4 = 0` only at `n=5` and `=3K rho^2/4≠0` at `n=3` — the guard is real; (2) re-derived `Theta_w = 4 rho_w^2 mu_*^2/(hbar^2 c_sw^2)` under `mu_* = lambda_mu m c_sw^2/4` and confirmed the `c_sw^2` cancellation and `/4` factor give the printed `lambda_mu^2 m^2 rho_w^2 c_sw^2/(4 hbar^2)`; (3) substituted `c_sw = hbar/(2 m ell)` and confirmed the `m^2` and `hbar^2` cancellation yields `lambda_mu^2 rho_w^2/(16 ell^2)` — this is the one genuinely non-trivial algebraic step and it is non-tautological (target written independently); (4) confirmed `25` emerges as `(20)^2/16 = 400/16` from `ref_factor = 1/20` (notes §4), not hardcoded; (5) checked the weakest assertion (A4/B4) is self-consistent-only but still exercises `1/16` and squaring and is anchored to the upstream-non-trivial healing form; (6) checked symbol domains — all positive amplitudes/ratios, `n_poly != 1` correctly guards the `1/(n-1)` integration; (7) confirmed the dead `h_polytrope` (py:38, carries a spurious extra `m` factor) is never referenced or printed and cannot affect any assertion or output; (8) reconciled every emitted deliverable value against the notes/tex/appendix — all MATCH. The only defect is `stale_output`: both committed transcripts carry the pre-renumber `STAGE 59/059` banner and predate the current scripts; content otherwise agrees exactly. Verdict `findings` (one low-severity, non-blocking stale-output), `paper_alignment: aligned`, `needs_user_resolution: false`.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 8 deliverable values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `c_s^2 [n=5] = 5K rho^4/m` | py:48 / wl:41; sympy out L5, mma out L5 | notes §1 `c_s^2(rho)=(5K/m)rho^4` (md:31) | MATCH |
| `h [n=5] = 5K rho^4/4` | py:49 / wl:42; sympy out L6, mma out L6 | notes §1 `h(rho)=(5K/4)rho^4` (md:40) | MATCH |
| `h = m c_s^2/4` (n=5 identity) | py:50 / wl:43; out L7 (`=0`) | tex:7 `\(h(\rho)=mc_s^2/4\)`; notes §1 `h(rho)=m c_s(rho)^2/4` (md:44) | MATCH |
| `Theta_w (enthalpy lock) = lambda_mu^2 m^2 rho_w^2 c_sw^2/(4 hbar^2)` | py:66 / wl:57; sympy out L9, mma out L10 | notes §2 `Theta_w = lambda_mu^2 m^2 rho_w^2 c_(s,w)^2/(4 hbar^2)` (md:71) | MATCH |
| `ell = hbar/(2 m c_sw)` (healing lock) | py:71 / wl:63; sympy out L11, mma out L13 | notes §3 `ell = hbar/(2 m c_(s,w))` (md:80); tex:7 "healing lock" (Inputs) | MATCH |
| `Theta_w (healing lock) = lambda_mu^2 rho_w^2/(16 ell^2)` | py:77 / wl:67; sympy out L13, mma out L16 | notes §3 `Theta_w = lambda_mu^2 rho_w^2/(16 ell^2)` (md:89) | MATCH |
| `Theta_w (ref, general a) = 25 lambda_mu^2 rho_w^2/a^2` | py:86 / wl:76; sympy out L14, mma out L17 | notes §4 (`ell/a=1/20`, `Theta_w=25 lambda_mu^2 rho_w^2`, md:101–113) — general-a intermediate | MATCH |
| `Theta_w (normalized) = 25 lambda_mu^2 rho_w^2` | py:88 / wl:77; sympy out L15, mma out L18 | tex:16–17 boxed `\Theta_w=25\lambda_\mu^2\rho_w^2`; notes §4 (md:113) / §5 (md:125); appendix:130 | MATCH |

Load-bearing input check: `ref_factor = 1/20` (py:83 / wl:73), the source of the `25`, is anchored in notes §4 `ell/a = 1/20` (md:101) — so the constant is paper/notes-anchored, not an unsupported script literal. No MISMATCH, no MISSING-DELIVERABLE.

Internal scaffolding (accounted for, no finding): pass/fail flags and `= 0` residual lines for the four `expect_zero`/`expectZero` checks; the `n=3` non-tautology residual (`3K rho^2/4`, a deliberate guard value, not a stage deliverable); `mu_star_solved`, `csw_from_ell`, `ell_solved` intermediates; the unused `h_polytrope` symbol (py:38, dead code).

## Self-test notes

Traps checked: (1) **variable independence** — `sp.diff(P, rho)` and `sp.diff(U, rho)` both depend genuinely on `rho` (P, U are polynomials in rho), so neither derivative is identically zero; `sp.diff` is over the right variable. (2) **No unbounded integrals** in this stage (the only `integrate` is the indefinite `∫P/rho^2 dρ = K rho^{n-1}/(n-1)`, an antiderivative, parity-irrelevant). (3) **Trivial-case substitution** — at `n=5` the residual `h - m c_s^2/4 = 5K rho^4/4 - 5K rho^4/4 = 0` (asserted zero, confirmed), and at `n=3` it is `3K rho^2/4 ≠ 0` (the guard's `assert_nonzero` correctly gives a nonzero literal). (4) **Path spec** — n/a, no missing-script directive; only an output refresh is prescribed. (5) **Paper round-trip** — the single fix (refresh `.txt` banners to `STAGE 076`) introduces no new `paper_misalignment`; no source/constant changes, so all reconciled values are untouched.
