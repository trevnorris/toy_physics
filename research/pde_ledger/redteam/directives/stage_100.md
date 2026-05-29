---
unit_id: 100
batch: IV.1
created_at: 2026-05-27T00:00:00Z
findings_count: 5
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 100

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_100.tex:16` quote: "Full odd normalization factorizes as `\widehat m_0^{\,2}\chi_Q N_Q=1`."
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_100.tex:22-24` quote: "Check the product `\widehat m_0^{\,2}\chi_QN_Q` keeps source, conservative, and outgoing factors separate. Check that higher odd terms begin beyond the point-particle 2.5PN coefficient. Check the outgoing l=2 DtN fingerprint against the normalized z=omega a/c_s expansion."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage100_outgoing_normalization_factorization.md:83-97` quote: "The actual point-particle observable branch includes the source-map factor `mhat_0`, so the full odd normalization condition is `mhat_0^2 Gammabar_5 = 2 G/(5 c^5)`. Substituting the renormalized one-pole branch gives `mhat_0^2 chi_Q N_Q = 1`."
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:73-74` quote: "`\widehat m_0^{\,2}\chi_Q N_Q=1.`"
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:295` quote: "Higher odd denominator terms beginning at `O(omega^7)` are invisible to the point-particle 2.5PN coefficient."
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:317-322` quote: "`\widehat Y_2^{\rm out}(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + O(z^6).`"
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage100_outgoing_normalization_factorization.md:22` quote: "`sigma_Q^can = 9/(8 Omega_Q^5) = 4 a^5/(27 c_s^5)`."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:32` quote: `print('mhat0^2*Gamma5/Gamma5_target - mhat0^2*chiQ*NQ =', sp.simplify(mhat0**2*Gamma5/Gamma5_t - mhat0**2*chiQ*NQ))`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:33` quote: `print('NQ - 1/(mhat0^2*chiQ) after odd normalization =', sp.simplify(sp.solve(sp.Eq(mhat0**2*chiQ*NQ, 1), NQ)[0] - 1/(mhat0**2*chiQ)))`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:56-62` quote: same five identities; A9/A10 mirror A4/A5.

## Resolve before fix_loop

The stage card's Output (`mhat_0^2 chi_Q N_Q = 1`) and three Checks have no substantive script-side verification:

1. The headline closure `mhat_0^2 chi_Q N_Q = 1`: the only scripted checks (A4 sympy:32, A5 sympy:33, A9 math:56-59, A10 math:61-62) are algebraically tautological — A4/A9 multiply A3/A8 by `mhat0^2`, and A5/A10 just confirm `solve(x*y=1, y)` returns `1/x`. Neither imposes the observable closure `mhat_0^2 * Gamma_5 = Gamma_5_target` on the *derived* `Gamma_5` and `NQ`.
2. Stage-card Check (ii) ("higher odd terms begin beyond the point-particle 2.5PN coefficient") is unverified — the series truncates at `omega^5` and never inspects the `omega^7` (or higher) odd content.
3. Stage-card Check (iii) ("outgoing l=2 DtN fingerprint against the normalized `z = omega a/c_s` expansion") is unverified — the script defines `sigma_can = 9/(8 Omega^5)` symbolically but never tests the equivalence with `4 a^5/(27 c_s^5)` (which fixes `Omega_Q = (3/2) c_s/a`), nor compares `Yhat_Q^ret(omega)` against `Yhat_2^out(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + O(z^6)` to conclude `chi_Q = 1`.

Possible directions (the user picks one or a combination):

- (a) **Strengthen the scripts to honor the card's Output and three Checks.** Replace A4/A5 (and A9/A10) with a derivation of `mhat_0^2 chi_Q NQ = 1` from the imposed observable condition `mhat_0^2 * Gamma_5 = Gamma_5_target`, using the script's series-derived `Gamma_5(K0, Omega, chiQ)` and `NQ = K0/K0_target`. Extend the series to `omega^7` (or higher) and assert the next imaginary coefficient sits at `omega^7` (or report which order is the first nonzero imaginary term beyond `omega^5`). Add a DtN-side derivation: define `Lambda_2^out(z) = z d/dz ln h_2^(1)(z)` (or hand-expand `h_2^(1)` to the needed order), compute `Yhat_2^out(z) = -3/Lambda_2^out(z)`, and assert `Yhat_2^out(z) - Yhat_Q^ret(omega)|_{omega = c_s z / a, chi_Q = 1} = O(z^6)`; also assert `9/(8 Omega_Q^5) == 4 a^5/(27 c_s^5)` under `Omega_Q = (3/2) c_s/a`.
- (b) **Trim the paper card's Checks list and Output to match the script's current scope.** Specifically, weaken stage_100.tex Output to "Ratios `K_2/K_2^target = K_4/K_4^target = N_Q`, `Gamma_5/Gamma_5^target = chi_Q N_Q`" and drop or relocate the Checks for higher-odd-term placement and the DtN fingerprint comparison (these would then have to be carried by another stage's audit). Notes Section 3 and appendix eq:app-part04-main-factorization would still need a citation to whichever stage *does* verify the closure.
- (c) **Split the burden:** (a)-style script extension for the closure (which is internal to this stage) but defer Checks (ii) and (iii) to a different audit unit. In that case the stage card's Checks list still needs editing so that this audit is graded against the deliverables it now claims.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. Findings F2 and F3 are coupled to F1 and will be resolved together; do not act on them in isolation.

## F2 — tautological_check (coupled to F1)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:33` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:61-62`

**Issue:**
A5/A10 invoke `solve` on the equation `mhat0^2 * chiQ * NQ = 1`, then assert the returned `NQ` equals `1/(mhat0^2 * chiQ)`. This is algebraically guaranteed by the inverse-linear structure of `solve`. The check cannot fail and certifies nothing about the physics; it is presented to the reader of the transcript as the headline closure check but is not.

**Required change:**
HOLD. Direction depends on F1 resolution. If direction (a): delete A5/A10 and replace with a substantive imposition of `mhat0^2 * Gamma_5 = Gamma_5_target` followed by simplification of `mhat0^2 * chi_Q * NQ - 1` (using the script's derived `Gamma_5` and `NQ`) to zero. If direction (b): delete A5/A10 outright with a comment that the closure is no longer claimed by this stage card.

## F3 — insufficient_verification (coupled to F1)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:32` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:56-59`

**Issue:**
A4/A9 (`mhat0^2 * Gamma5/Gamma5_t - mhat0^2 * chiQ * NQ`) is algebraically `mhat0^2 * (Gamma5/Gamma5_t - chiQ * NQ)`. Given A3/A8 already drives the parenthetical to zero, A4/A9 is trivially zero and contributes no new content.

**Required change:**
HOLD. Direction depends on F1 resolution. Under direction (a), A4/A9 is replaced by the closure derivation. Under direction (b), A4/A9 is deleted.

## F4 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:26-62`

**Issue:**
The `.wl` script is structurally a line-by-line port of the `.py` script (same `sigma_can = 9/(8 Omega^5)`, same `Y = 3/4 + (1/4)/(1 - omega^2/Omega^2 - I*chiQ*sigma_can*omega^5)`, same series order, same coefficient extraction, same five identities in the same order). The banner header even reads `STAGE 083 — OUTGOING NORMALIZATION FACTORIZATION` (wl:26), consistent with copy-paste from an earlier stage. The second-engine policy requires Mathematica to re-derive the result from physical premises (e.g., starting from the spherical Hankel function `h_2^(1)(z)` and the DtN definition) rather than echo SymPy's algebra.

**Required change:**
Re-derive the Mathematica side from a different physical starting point. Concretely: start from `Lambda_2^out(z) = z * D[Log[SphericalHankelH1[2, z]], z]` (use the `SphericalHankelH1` built-in or hand-expand `h_2^(1)(z)` to the needed `O(z^6)` order), compute `Yhat_2^out(z) = -3/Lambda_2^out(z)`, expand to `O(z^6)`, and verify the canonical form `1 + z^2/9 + 4 z^4/81 + i z^5/27`. Then identify `omega^2/Omega_Q^2 = z^2 * (c_s/(a Omega_Q))^2`, use `Omega_Q = (3/2) c_s/a` (from `sigma_Q^can`-equivalence), and confirm `chi_Q = 1` reproduces the same series as `Yhat_Q^ret`. The K2/K4/Gamma_5 ratio identities then follow from the canonical DtN side, not from re-extracting the SymPy series. Also fix the banner at wl:26 to read `STAGE 100 — ...`.

This is a design-level rewrite, not a mechanical edit. Codex should not attempt to apply it without the user's go-ahead because the choice of independent derivation path determines what intermediate quantities the new `.wl` will print and how the verifier judges "engines independently agree."

**Blocked: F4** — design-level rewrite. Please confirm the preferred independent derivation path before Codex acts.

## F5 — symbol_assumption_error

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:7` and `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:30-31`

**Issue:**
Both engines declare `chiQ` as positive. The paper treats `chi_Q` as a real renormalization factor whose canonical value is `1` and whose departures (positive or negative) are the obstruction `Delta_Q = chi_Q - 1`. The positivity is not required by any of the assertions A1-A3 (all linear in `chi_Q`) and creates a latent trap for future revisions.

**Required change:**

1. In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py:7`, change the line:
   ```
   K0, mhat0, chiQ = sp.symbols('K0 mhat0 chiQ', positive=True, real=True)
   ```
   to:
   ```
   K0, mhat0 = sp.symbols('K0 mhat0', positive=True, real=True)
   chiQ = sp.symbols('chiQ', real=True)
   ```
   `K0` and `mhat0` remain positive (they are amplitudes/source-map factors and the paper's setup keeps them positive); `chiQ` becomes simply real.

2. In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl:30-31`, change:
   ```
   $Assumptions =
     Element[{gNewton, cLight, omegaQ, k0, mHat0, chiQ, omega, nQSym}, Reals] &&
     gNewton > 0 && cLight > 0 && omegaQ > 0 && k0 > 0 && mHat0 > 0 && chiQ > 0 && nQSym > 0;
   ```
   to:
   ```
   $Assumptions =
     Element[{gNewton, cLight, omegaQ, k0, mHat0, chiQ, omega, nQSym}, Reals] &&
     gNewton > 0 && cLight > 0 && omegaQ > 0 && k0 > 0 && mHat0 > 0 && nQSym > 0;
   ```
   (Drop `chiQ > 0`; keep all other positivity constraints.)

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 100` and `redteam exec-mathematica 100` and confirms (i) the script files reflect the new symbol declarations and (ii) all five existing assertions still report `= 0` and `PASS` (the assertions are linear in `chi_Q` and do not depend on its sign).

---
## Applied: F1+F2+F3 (orchestrator-direct, post-user-resolution per batch-IV1-paper-alignment Cluster B direction (c))

- files_changed: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py, mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl
- summary: Both scripts now add docstring/comment annotations naming stage 097 (Check 3 DtN fingerprint pinning chi_Q = 1) and stage 102 (Check 2 higher-odd-term placement) as upstream Check anchors. The substantive closure for the headline Output mhat_0^2 chi_Q N_Q = 1 is derived by IMPOSING the observable condition mhat_0^2 * Gamma_5 = Gamma_5_target on the script-derived Gamma_5(K_0, chi_Q, Omega), then asserting closure_ratio - (mhat_0^2 chi_Q N_Q - 1) = 0. Tautological A4/A5 (and A9/A10) are replaced by this substantive derivation. F5 (chiQ symbol_assumption_error) also closed — chiQ is now declared real, not positive.
- deviation: F4 (mathematica_transliteration, design-level rewrite) remains blocked per the auditor's "BLOCKED" tag; the orchestrator did not extend Mathematica to a fully independent derivation route. The substantive closure derivation IS however a content difference from the prior tautological cross-check, making the two engines no longer pure transliterations of each other at the load-bearing assertion.

## Applied: F5
- files_changed: scripts/moving_throat_pde_stage100_outgoing_normalization_factorization_sympy_audit.py, mathematica/moving_throat_pde_stage100_outgoing_normalization_factorization_mathematica_audit.wl
- summary: chiQ now declared real (not positive) in both engines, per the directive's symbol-assumption fix; banners updated to STAGE 100.
- deviation: none
