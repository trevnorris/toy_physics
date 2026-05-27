---
unit_id: 100
batch: IV.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 5
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage100_outgoing_normalization_factorization.md]
  paper_appendix: present
---

# Audit unit 100 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_100.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage100_outgoing_normalization_factorization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows around L60-86, L260-296 for the retarded factorization block; theorem at L65-83; main factorization eq at L73-74)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.txt`

## What the paper claims

The stage card's bottom-line Output is verbatim (stage_100.tex:15-17): "Full odd normalization factorizes as `mhat_0^2 chi_Q N_Q = 1`." The card lists three explicit Checks (stage_100.tex:21-25):
(i) the product `mhat_0^2 chi_Q N_Q` keeps source-map, conservative-normalization, and outgoing-normalization factors separate;
(ii) higher odd terms begin beyond the point-particle 2.5PN coefficient (i.e., next imaginary term is `O(omega^7)` per appendix L295);
(iii) the outgoing `l=2` DtN fingerprint matches the normalized `z = omega a/c_s` expansion `Yhat_2^out(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + O(z^6)` (appendix L317-322), which fixes `chi_Q = 1` (appendix L324-327).
The notes additionally enumerate: `Kbar_2 = Kbar_0/(4 Omega_Q^2)`, `Kbar_4 = Kbar_0/(4 Omega_Q^4)`, `Gammabar_5 = chi_Q * 9 Kbar_0/(32 Omega_Q^5)`, the three target invariants, the three ratios `Kbar_2/Kbar_2^target = N_Q`, `Kbar_4/Kbar_4^target = N_Q`, `Gammabar_5/Gammabar_5^target = chi_Q N_Q`, and the observable-closure condition `mhat_0^2 Gammabar_5 = 2G/(5c^5)` collapsing to `mhat_0^2 chi_Q N_Q = 1` (notes Sec. 3). The canonical sigma equivalence `sigma_Q^can = 9/(8 Omega_Q^5) = 4 a^5/(27 c_s^5)` (notes L22) is also asserted.

## What the script claims to verify

The SymPy script (and the structurally identical Mathematica mirror) (i) defines the retarded-branch `Y = 3/4 + 1/4 / (1 - omega^2/Omega^2 - I*chiQ*sigma_can*omega^5)` with `sigma_can = 9/(8 Omega^5)`, (ii) extracts the `omega^2`, `omega^4`, and `omega^5` coefficients to build `K2`, `K4`, `Gamma5`, (iii) defines the four targets `K0_t = 64 G Omega^5/(45 c^5)`, `K2_t = K0_t/(4 Omega^2)`, `K4_t = K0_t/(4 Omega^4)`, `Gamma5_t = 2 G/(5 c^5)`, and `NQ := K0/K0_t`, (iv) asserts the three ratio identities `K2/K2_t - NQ = 0`, `K4/K4_t - NQ = 0`, `Gamma5/Gamma5_t - chiQ*NQ = 0`, (v) asserts a fourth identity that is just `mhat0^2 *` the third (`mhat0^2*Gamma5/Gamma5_t - mhat0^2*chiQ*NQ = 0`), and (vi) solves the equation `mhat0^2*chiQ*NQ_sym = 1` for `NQ_sym` and confirms the SymPy result equals `1/(mhat0^2*chiQ)`. There is no assertion that the headline factorization `mhat_0^2 chi_Q N_Q = 1` actually follows from the observable normalization condition `mhat_0^2 Gamma_5 = Gamma_5^target`; no test that the next odd denominator term sits at `O(omega^7)`; no test of the DtN fingerprint against the `z = omega a/c_s` expansion; and no test of the equivalence `9/(8 Omega_Q^5) = 4 a^5/(27 c_s^5)`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Kbar_2 = Kbar_0/(4 Omega^2)` and `Kbar_2/target = NQ` | `K2/K2_t - NQ == 0` (sympy:29, math:53) | match |
| `Kbar_4 = Kbar_0/(4 Omega^4)` and `Kbar_4/target = NQ` | `K4/K4_t - NQ == 0` (sympy:30, math:54) | match |
| `Gammabar_5 = chi_Q * 9 Kbar_0/(32 Omega^5)` and `Gammabar_5/target = chi_Q NQ` | `Gamma5/Gamma5_t - chiQ*NQ == 0` (sympy:31, math:55) | match |
| Output: `mhat_0^2 chi_Q N_Q = 1` derived from observable closure `mhat_0^2 Gamma_5 = Gamma_5^target` (notes Sec. 3; appendix eq:app-part04-main-factorization L73-74) | A4 multiplies A3 by `mhat0^2` (no closure imposed); A5 inverts a freestanding equation (tautological) | mismatch |
| Check (ii): higher odd terms begin beyond 2.5PN (next imaginary at O(omega^7); appendix L295) | none | missing |
| Check (iii): DtN fingerprint vs `z = omega a/c_s` expansion (appendix L307-327); equivalence `9/(8 Omega^5) = 4 a^5/(27 c_s^5)` (notes L22); `chi_Q = 1` follows | none | missing |

Three deliverables match, one is mismatched (the central headline), two are missing. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 29 | `simplify(K2/K2_t - NQ) == 0` | even ratio (notes L57) | yes |
| A2 | sympy | 30 | `simplify(K4/K4_t - NQ) == 0` | even ratio (notes L59) | yes |
| A3 | sympy | 31 | `simplify(Gamma5/Gamma5_t - chiQ*NQ) == 0` | odd ratio (notes L77) | yes |
| A4 | sympy | 32 | `simplify(mhat0^2*Gamma5/Gamma5_t - mhat0^2*chiQ*NQ) == 0` | nominally headline closure | no — algebraically `mhat0^2 *` A3 |
| A5 | sympy | 33 | `simplify(solve(mhat0^2*chiQ*NQ=1, NQ)[0] - 1/(mhat0^2*chiQ)) == 0` | nominally headline closure | no — tests Solve, not physics |
| A6 | mathematica | 53 | `expectZero[K2/K2Target - nQ]` | even ratio | yes |
| A7 | mathematica | 54 | `expectZero[K4/K4Target - nQ]` | even ratio | yes |
| A8 | mathematica | 55 | `expectZero[Gamma5/Gamma5Target - chiQ*nQ]` | odd ratio | yes |
| A9 | mathematica | 56-59 | `expectZero[mHat0^2*Gamma5/Gamma5Target - mHat0^2*chiQ*nQ]` | nominally headline closure | no — `mHat0^2 *` A8 |
| A10 | mathematica | 61-62 | `expectZero[(nQSym /. Solve[mHat0^2*chiQ*nQSym==1, nQSym]) - 1/(mHat0^2*chiQ)]` | nominally headline closure | no — tests Solve |

A4/A9 and A5/A10 are the candidates for `tautological_check` / `insufficient_verification`; the latter pair are also `paper_misalignment` because the script presents them as the headline closure check but does not actually impose the observable condition.

## Findings

### F1 — paper_misalignment

**Severity:** high
**Subtype:** script_missing_paper_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_100.tex:15-17` (Output line)
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_100.tex:22-25` (Checks)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage100_outgoing_normalization_factorization.md:83-97` (Section 3, "Observable point-particle normalization condition")
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:73-76, 284-288` (main factorization equation)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:32-33` (A4/A5)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:56-62` (A9/A10)

**What's wrong:**
The paper card's Output reads (stage_100.tex:16): "Full odd normalization factorizes as `mhat_0^2 chi_Q N_Q = 1`." The notes derive this closure (Sec. 3) by *imposing* the observable condition `mhat_0^2 Gammabar_5 = 2 G/(5 c^5)` and substituting `Gammabar_5 = chi_Q * 9 Kbar_0/(32 Omega^5)` together with `Kbar_0 = N_Q * Kbar_0^target` to get `mhat_0^2 chi_Q N_Q = 1`. The appendix encodes the same statement at eq:app-part04-main-factorization (L73-74).

The script never imposes the observable-closure condition. The two assertions that nominally cover it are tautological:
- A4 (sympy:32): `simplify(mhat0**2*Gamma5/Gamma5_t - mhat0**2*chiQ*NQ)` is algebraically `mhat0**2 *` A3; if A3 is zero, A4 is automatically zero. It tests nothing new.
- A5 (sympy:33): `sp.solve(sp.Eq(mhat0**2*chiQ*NQ, 1), NQ)[0] - 1/(mhat0**2*chiQ)` solves the equation for `NQ` and confirms the SymPy result equals `1/(mhat0**2*chiQ)`. This tests SymPy's solver, not physics — `solve(x*y=1, y)` returns `1/x` by construction.
The headline `mhat_0^2 chi_Q N_Q = 1` is therefore unverified by the script.

In addition, two paper-side checks have no script-side counterpart:
- Check (ii) (stage_100.tex:23): "higher odd terms begin beyond the point-particle 2.5PN coefficient." The appendix (L295) makes this concrete: the next imaginary denominator term sits at `O(omega^7)`. The series in both engines truncates at `omega^5` (sympy series order `6`, Mathematica `Series[..., {omega, 0, 5}]`), so the script cannot even see the next odd term.
- Check (iii) (stage_100.tex:24): "outgoing l=2 DtN fingerprint against the normalized z = omega a/c_s expansion." The appendix gives the exact target (L317-322): `Yhat_2^out(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + O(z^6)`, which fixes `chi_Q = 1` (L324-327). The notes also state `sigma_Q^can = 9/(8 Omega_Q^5) = 4 a^5/(27 c_s^5)` (L22). The script uses `sigma_can = 9/(8 Omega^5)` but never (a) tests the equivalence with `4 a^5/(27 c_s^5)` under `Omega_Q = (3/2) c_s/a`, nor (b) compares `Yhat_Q^ret` against `Yhat_2^out(z)` to derive `chi_Q = 1`.

**Why this matters:**
The stage's Output is the central reduced-2.5PN factorization gate cited downstream in Stages 107-113 (appendix L86). The script as written cannot fail in a way that would reveal `mhat_0^2 chi_Q N_Q = 1` is wrong, because no actual constraint is imposed in the script. The checks the card explicitly enumerates as auditing this stage are absent.

**Required change:**
See `## Resolve before fix_loop` in the directive — direction depends on whether the script should be expanded to honor the paper's checklist, or the paper card's Checks list should be trimmed to match the script's current scope.

**Verification:**
Once user resolves direction, verifier confirms either (a) the new SymPy/Mathematica checks for `mhat_0^2 Gamma_5 = Gamma_5_target → mhat_0^2 chi_Q NQ = 1`, the higher-odd-term sentinel, and the DtN fingerprint comparison appear and exit 0; or (b) the paper card and notes/appendix are amended to remove the corresponding deliverables.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:33`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:61-62`

**What's wrong:**
The "NQ - 1/(mhat0^2*chiQ) after odd normalization" check is `sp.solve(sp.Eq(mhat0**2*chiQ*NQ, 1), NQ)[0] - 1/(mhat0**2*chiQ)`. Solving `mhat0^2 * chiQ * NQ = 1` for `NQ` yields `NQ = 1/(mhat0^2 * chiQ)` by construction. Subtracting `1/(mhat0^2 * chiQ)` is guaranteed zero independent of any physics; it certifies SymPy's `Solve` correctly inverts a linear equation. The Mathematica mirror at L61-62 has the identical structure.

**Why this matters:**
This is the closest the script comes to verifying the headline output `mhat_0^2 chi_Q N_Q = 1`, but it cannot fail. The "PASS" line in the transcript creates a false impression that the factorization closure has been verified.

**Required change:**
Resolution coupled to F1. If the user chooses to add a real closure check, A5/A10 are subsumed and should be replaced by a derivation of `mhat_0^2 chi_Q NQ = 1` from the observable condition `mhat_0^2 * Gamma_5 = Gamma_5_target` using the script's *derived* `Gamma_5` (from the series) and `NQ = K0/K0_target`. If the user chooses to trim the paper claim instead, A5/A10 should be deleted (they verify nothing).

**Verification:**
After fix, the residual the script prints comes from substituting derived `Gamma_5(K0, Omega, chiQ)` and `NQ(K0, G, Omega, c)` into `mhat0^2 * Gamma_5 / Gamma_5_target` and simplifying to `mhat0^2 * chiQ * NQ`, not from inverting an algebraic identity.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:32`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:56-59`

**What's wrong:**
The check `simplify(mhat0**2*Gamma5/Gamma5_t - mhat0**2*chiQ*NQ)` is algebraically `mhat0**2 * (Gamma5/Gamma5_t - chiQ*NQ)`. Because the prior assertion A3 (sympy:31) already establishes `Gamma5/Gamma5_t - chiQ*NQ == 0`, A4 reduces to `mhat0**2 * 0 == 0` automatically. It contributes no new test.

**Why this matters:**
A4/A9 looks superficially like an independent factorization test, but does not exercise anything beyond what A3/A8 already tests.

**Required change:**
Resolution coupled to F1. Either replace A4 with a substantive closure check (see F1/F2), or delete it as dead-weight scaffolding.

**Verification:**
After replacement, the residual must depend on quantities not present in A3 (e.g., the script must impose `mhat_0^2 * Gamma_5 = Gamma_5_target` and *derive* `mhat_0^2 chi_Q NQ - 1 = 0`).

### F4 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:33-62`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:10-33`

**What's wrong:**
The `.wl` script is a structural line-by-line port of the `.py` script with matched variable renames:

| sympy (.py) | mathematica (.wl) |
|---|---|
| L10 `sigma_can = sp.Rational(9, 8) / Omega**5` | L33 `sigmaCan = FullSimplify[(9/8)/omegaQ^5, ...]` |
| L11 `Y = sp.Rational(3,4) + sp.Rational(1,4)/(1 - omega^2/Omega^2 - sp.I*chiQ*sigma_can*omega^5)` | L34 `yRet = 3/4 + (1/4)/(1 - omega^2/omegaQ^2 - I*chiQ*sigmaCan*omega^5)` |
| L12 `Yser = sp.expand(sp.series(Y, omega, 0, 6).removeO())` | L35 `ySeries = Expand[Normal[Series[yRet, {omega, 0, 5}]]]` |
| L14-16 `K2 = K0*Yser.coeff(omega, 2)`, `K4 = K0*Yser.coeff(omega, 4)`, `Gamma5 = sp.im(Yser.coeff(omega,5))*K0` | L37-39 `k2 = k0*Coefficient[ySeries, omega, 2]`, `k4 = k0*Coefficient[ySeries, omega, 4]`, `gamma5 = (Coefficient[ySeries, omega, 5]/I)*k0` |
| L18-22 `K0_t, K2_t, K4_t, Gamma5_t, NQ` defs | L41-45 same defs in same order |
| L29-33: five `print(simplify(...))` calls in fixed order | L53-62: five `expectZero[...]` calls in identical fixed order |

There is no independent re-derivation in the Mathematica file — same Y, same series, same coefficient choreography, same five identities in the same order. The Mathematica banner header even reads `STAGE 083 — OUTGOING NORMALIZATION FACTORIZATION` (wl:26), consistent with copy-paste origin.

**Why this matters:**
Per the second-engine policy, the two engines must independently exercise the physics so that an engine-specific simplifier bug or normalization choice that lets a false identity slip through one engine has to also slip through the other. A transliteration provides no such cross-check; both will succeed or fail together.

**Required change:**
Re-derive the Mathematica side from the *physical premises*. Concretely: (a) start from the DtN side `Yhat_2^out(z)` (appendix eq:app-part04-Yout-dtn) or from the spherical-Hankel logarithmic derivative `Lambda_2^out(z) = z d/dz ln h_2^(1)(z)` and compare term-by-term against `Yhat_Q^ret(omega)` to fix `chi_Q = 1`; or (b) impose the observable closure `mHat0^2 * gamma5 = gamma5Target` symbolically and *solve* for the factorization condition, rather than transcribing the SymPy bookkeeping. This is design-level work and is not safe to mechanize. Flag for user direction.

**Verification:**
After rewrite, the `.wl` should not share the same step ordering as the `.py`. Final asserted identities may coincide (they must), but the intermediate algebra should differ visibly.

### F5 — symbol_assumption_error

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:7`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:30-31`

**What's wrong:**
Both engines declare `chiQ` (resp. `chiQ`) as `positive=True` / `chiQ > 0`. The paper card and notes treat `chi_Q` as a real renormalization factor whose canonical value is `1` (appendix L324-327) and whose deviation from `1` is the obstruction `Delta_Q = chi_Q - 1` (appendix eq:app-part04-DeltaQ-def L292-294). Nothing in the stage's algebraic ratios depends on the sign of `chi_Q`. Restricting `chi_Q > 0` quietly excludes the negative branch from `simplify` / `FullSimplify` even though the identities being tested do not require it.

**Why this matters:**
For the present assertions A1-A3 (the substantive ones), the positivity is harmless because the identities are linear in `chi_Q`. But it creates a latent trap if a future revision tests a quantity where positivity matters (e.g., a square-root branch of `chi_Q`).

**Required change:**
Replace `positive=True` with `real=True` for `chiQ` in the SymPy script line 7; in the Mathematica script L31, remove `chiQ > 0` from `$Assumptions`. Re-run; assertions remain zero (no algebraic dependence on the positivity of `chiQ`).

**Verification:**
After Codex applies, run sympy and mathematica audits; all five assertions still report `= 0`.

## Independent-derivation check (Mathematica)

Not independent. The `.wl` mirrors the `.py` step-for-step (see F4 table above). Both define `sigma_can = 9/(8 Omega^5)`, build the same retarded series, extract identical even/odd coefficients, and assert the same five identities in the same order. The Mathematica banner even cites the wrong stage number ("STAGE 083") on `wl:26`, consistent with a copy-paste origin from an earlier stage. The Mathematica side does not, for example, work from the spherical Hankel function `h_2^(1)(z)` and the DtN definition (appendix eq:app-part04-Lambda-out-dtn L307-313) to verify `chi_Q = 1` independently.

## Engine cross-check

Both engines report identical zero residuals on the five assertions (sympy.txt:14-18; math.txt:18-27). Headline-printed quantities also match modulo symbol names:
- `K2 = K0/(4*Omega^2)` ↔ `K2 = k0/(4*omegaQ^2)`
- `K4 = K0/(4*Omega^4)` ↔ `K4 = k0/(4*omegaQ^4)`
- `Gamma5 = 9*K0*chiQ/(32*Omega^5)` ↔ `Gamma5 = (9*chiQ*k0)/(32*omegaQ^5)`
- `NQ = 45*K0*c^5/(64*G*Omega^5)` ↔ `NQ = (45*cLight^5*k0)/(64*gNewton*omegaQ^5)`
Engines agree; but agreement is uninformative given F4 (the engines are running the same algebra).

## Verdict justification

`findings`. The three even/odd ratio identities (A1-A3 / A6-A8) hold up: they substantively test the form of the chosen retarded-branch coefficients against the targets. But the headline output of the stage card — `mhat_0^2 chi_Q N_Q = 1` — is not verified: A4/A9 are algebraic shadows of A3/A8, and A5/A10 invert an equation rather than imposing the observable closure. Two of the three Checks listed on the stage card (higher odd terms; DtN fingerprint against `z = omega a/c_s`) have no script-side counterpart. The Mathematica file is a transliteration. Symbol assumptions are slightly over-restrictive on `chi_Q`. The factorization stage cannot fail in this script in a way that would expose a wrong headline. Outputs are fresh.

Attacks tried that failed: (i) tried to break A1-A3 by reading the series coefficients off the chosen `Y` form — they reduce correctly to the canonical `1/(4 Omega^2)`, `1/(4 Omega^4)`, `9 chi_Q/(32 Omega^5)` and so the ratio identities are real; (ii) tried to break the engine agreement — both produce identical algebraic forms.

## Self-test notes

Walked through each Codex-applicable change before writing the directive: (1) F5 — replacing `positive=True` with `real=True` on `chiQ` (and dropping `chiQ > 0` from `$Assumptions`) leaves A1-A3 as linear-in-`chiQ` identities that simplify to zero regardless of sign; confirmed by inspection of the explicit forms in the saved output. (2) For F1/F2/F3, no mechanical Codex change is being prescribed — these are paper_misalignment-coupled and require user direction. (3) For F4, the rewrite is non-mechanical and flagged as Blocked. No `sp.diff`/`D[...]` differentiation patterns appear in this audit, so the variable-independence trap does not apply. No parity/symmetry integrals appear. Path specifications are unchanged (scripts/ and mathematica/ files already exist).
