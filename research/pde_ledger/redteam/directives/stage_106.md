---
unit_id: 106
batch: IV.2
created_at: 2026-05-27T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-29T16:49:38Z
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 106

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_106.tex:21-25` quote: "Check the product `\widehat m_0^{\,2}\chi_Q N_Q` keeps source, conservative, and outgoing factors separate. Check that higher odd terms begin beyond the point-particle 2.5PN coefficient. Check the outgoing l=2 DtN fingerprint against the normalized `z=\omega a/c_s` expansion."
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:308-313` quote: "`\Lambda_2^{\rm out}(z) = z\frac{d}{dz}\ln h_2^{(1)}(z) = -3+\frac{z^2}{3}+\frac{z^4}{9}+i\frac{z^5}{9}+O(z^6).`"
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex:317-322` quote: "`\widehat Y_2^{\rm out}(z) = \frac{-3}{\Lambda_2^{\rm out}(z)} = 1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6).`"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:25-71` — no z-expansion of `\Lambda_2^{\rm out}` or `\widehat Y_2^{\rm out}` anywhere; `\chi_Q=1` is imported as a substitution at line 36.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:27-75` — same gap.

## Resolve before fix_loop

The paper card's Checks block (stage_106.tex:21-25) names three verification items. The scripts cover only the first (and weakly): they display the factored product `m0hat^2·chi_Q·N_Q` and solve for `N_Q`, but they do not perform the z-expansion that derives `\chi_Q=1` (item iii) and they do not show that higher odd terms first appear at `O(\omega^7)` (item ii). Item (iii) in particular is the load-bearing derivation — without it, `\chi_Q=1` is an imported assumption, not a script-verified result, and the script's "closure" reduces to substitution algebra.

Possible directions (the user picks one):
- (a) **Script must own all three Checks.** Extend both scripts to (i) construct `\Lambda_2^{\rm out}(z) = z·d/dz·log(h_2^{(1)}(z))` from the spherical-Hankel definition, series-expand to `O(z^6)`, and assert the coefficients match `(-3, 0, 1/3, 0, 1/9, i/9)`; (ii) form `\widehat Y_2^{\rm out}(z) = -3/\Lambda_2^{\rm out}(z)`, expand to `O(z^6)`, and assert the coefficients match `(1, 0, 1/9, 0, 4/81, i/27)`; (iii) match the canonical retarded form `1/(1 - \omega^2/\Omega_Q^2 - i\chi_Q\sigma_Q^{\rm can}\omega^5)` against this expansion and verify `\chi_Q = 1` falls out; (iv) verify the next odd term is at `O(\omega^7)`. After this is added, Checks (ii) and (iii) become independently verified rather than imported.
- (b) **Card delegates Checks (ii) and (iii) to upstream Stage 88.** Per the notes ("From Stage 88, the canonical compact passive/outgoing grouped-`P2` DtN model gives `\chi_Q=1`"), the z-expansion may already be the responsibility of Stage 88's scripts. If so, update the stage_106.tex Checks block to read "carry forward from Stage 88" for items (ii) and (iii), and only list item (i) as a Stage-106-side check.
- (c) **Item (iii) belongs upstream but item (ii) belongs here.** Split: keep (ii) as a Stage-106-side script check (assert no `\omega^7` term in the expanded `\widehat Y_Q^{\rm ret}` at chi_Q=1, m0hat=1), delegate (iii) to Stage 88.

## RESOLVED — direction (b)-style, refined by Codex (consult 019e748e, 2026-05-29)

Codex CONCUR: chi_Q=1 is a LEGITIMATE carry-in (not a hidden assumption), so **no stage-106 Hankel z-expansion / new verification is needed.** Codex refined the upstream citation (R4 DISPUTE): **Stage 104 proves the DtN fingerprint but introduces no `chi_Q` symbol — chi_Q=1 is FIXED at Stage 105** (inheriting 104's fingerprint). Accurate ownership:
- item (ii) [higher odd at O(ω⁷)] → stage **102**
- item (iii) [outgoing l=2 DtN fingerprint] → stage **104**; **chi_Q=1 → stage 105**

Paper-card cross-reference (ii→102, iii→104+105) logged to PAPER_CLEANUP_TRACKER (manual paper pass).

**Codex: ONE script-doc accuracy edit rides this fix loop (the only F1 script action):** in BOTH the SymPy docstring and the `.wl` header comment, the line attributing chi_Q=1 to "Stage 104" must be corrected — Stage 104 provides the DtN fingerprint (item iii); **chi_Q=1 is fixed at Stage 105** from that fingerprint; Stage 102 covers item (ii). Keep chi_Q=1 as the carry-in. Then append `## Applied: F1`. Do NOT add any Hankel/z-expansion verification to 106 and do NOT edit paper/notes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl`
- summary: Corrected the script comments to cite Stage 104 for the DtN fingerprint and Stage 105 for fixing `chi_Q=1`.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:55`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:60`

**Issue:** The assertion `K4 - 4*K2**2/K0 == 0` is algebraically tautological because the script defines `K2 = K0/(4*(3*c_s/(2*a))**2)` and `K4 = K0/(4*(3*c_s/(2*a))**4)` immediately above. Substituting yields `4*K2^2/K0 = K0/(4Ω^4) = K4` identically, regardless of any numerical constant. The K0-K2-K4 canonical-target relation `K_0K_4 = 4K_2^2` (appendix line 114) is the intended physical claim, but the current check enforces this by construction rather than by independent numerical agreement among the four hardcoded `*_target` literals.

**Required change:**

In the SymPy script (`scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py`):

Replace lines 55-59 (the two `expect_zero` calls for the branch identities) with versions that test the four hardcoded `*_target` literals' mutual consistency, not the script-defined K2/K4. Specifically:

Before (lines 55-59):
```python
expect_zero("branch identity K4 - 4 K2^2 / K0", K4 - 4 * K2**2 / K0)
expect_zero(
    "branch identity Gamma5 - 9 K2^(5/2)/K0^(3/2)",
    Gamma5 - 9 * sp.sqrt(K2**5 / K0**3),
)
```

After:
```python
expect_zero(
    "target identity K0_target K4_target - 4 K2_target^2",
    K0_target * K4_target - 4 * K2_target**2,
)
expect_zero(
    "target identity Gamma5_target - 9 sqrt(K2_target^5 / K0_target^3)",
    Gamma5_target - 9 * sp.sqrt(K2_target**5 / K0_target**3),
)
```

These two assertions now test the four hardcoded literals' mutual algebraic consistency. Substituting (for example) `K2_target = 7*G*c_s**3/(5*a**3*c**5)` would cause both to fail.

In the Mathematica script (`mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl`):

Replace lines 60-64 with the Mathematica equivalent:

Before (lines 60-64):
```mathematica
expectZero["branch identity K4 - 4 K2^2 / K0", k4 - 4*k2^2/k0];
expectZero[
  "branch identity Gamma5 - 9 K2^(5/2)/K0^(3/2)",
  gamma5 - 9*Sqrt[k2^5/k0^3]
];
```

After:
```mathematica
expectZero["target identity K0_target K4_target - 4 K2_target^2",
  k0Target*k4Target - 4*k2Target^2];
expectZero["target identity Gamma5_target - 9 sqrt(K2_target^5 / K0_target^3)",
  gamma5Target - 9*Sqrt[k2Target^5/k0Target^3]];
```

After applying F3 (Mathematica re-authoring), this F2 fix may need to be re-applied inside the new derivation. Apply F3 first, then port the corrected target-identity assertions into the new structure.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 106` and `redteam exec-mathematica 106` and confirm: (1) both scripts still exit 0; (2) the new assertion lines print residuals = 0 (since the four target literals are mutually consistent); (3) replacing `K2_target = 6*G*c_s**3/(5*a**3*c**5)` with `K2_target = 7*G*c_s**3/(5*a**3*c**5)` causes the SymPy run to fail at the new assertion (manual spot-check; not part of automated verifier).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl`
- summary: Replaced tautological branch checks with target-literal consistency checks for the canonical even and odd coefficients.
- deviation: none

## F3 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl` (entire body, lines 27-75)

**Issue:** The Mathematica script reproduces the SymPy script's variable choreography line for line — same `constraint`, same `Solve` → `NQ` substitution, same `chi_Q → 1` substitution, same `K0 = NQ * K0_target`, same `K2 = K0/(4Ω^2)`, same final three `expectZero` calls in the same order, even the same prose in the RESULT block. This violates the second-engine policy: a SymPy bug would be reproduced verbatim.

**Required change:**

Re-author the Mathematica audit so that it derives the same bottom-line claim (`N_Q = 1` on the canonical point-particle branch, `\gamma_{\rm quad}^{\rm eff} = 2G/(5c^5)`) by a structurally independent path. Specifically, start from the appendix's retarded one-pole form `\widehat Y_Q^{\rm ret}(\omega) = 3/4 + (1/4)·1/(1 - \omega^2/\Omega_Q^2 - i\chi_Q\sigma_Q^{\rm can}\omega^5)` (appendix line 264-270) and proceed by symbolic series expansion in `\omega`:

1. Define `\Omega_Q = 3 c_s/(2 a)` and `\sigma_Q^{\rm can} = 9/(8\Omega_Q^5)` symbolically.
2. Define `Yret[omega_, chiQ_] := 3/4 + (1/4) · 1/(1 - omega^2/OmegaQ^2 - I*chiQ*sigmaQcan*omega^5)`.
3. Series-expand `Yret[omega, chiQ]` in omega to order 5; read off the `\omega^0`, `\omega^2`, `\omega^4`, `\omega^5` coefficients. The `\omega^5` coefficient is `i·chi_Q·sigmaQcan/4`.
4. The canonical odd-coefficient relation (appendix eq. app-part04-gamma5-chiN, line 279-281): `\overline\Gamma_5 = \chi_Q · 9\overline K_0/(32\Omega_Q^5)`. Substitute `\overline K_0 = N_Q · K0_target`, `\overline\Gamma_5 = (chiQ·NQ) · (2G/(5c^5))` (this should drop out from the series coefficient).
5. Independently impose `\widehat m_0^2 · \chi_Q · N_Q = 1` from the source-map relation `\widehat m_0 → 1` plus the canonical `chi_Q → 1` and derive `N_Q = 1`.
6. Assert: at `m0hat=1, chi_Q=1`, the `\omega^5` coefficient of `Yret` equals `i·sigmaQcan/4 = i·9/(32·\Omega_Q^5)`, and verify this matches `9·K0_target/(32·\Omega_Q^5)/K0_target = 9/(32·\Omega_Q^5)` from the K0_target normalization, i.e., the canonical Gamma5_target = 2G/(5c^5) is consistent.
7. Verify the next-order odd coefficient in `Yret` begins at `\omega^7` (this also covers paper Check (ii) — see F1).

Use entirely different intermediate variable names from the SymPy file: `Yret`, `seriesY`, `OmegaQ`, `sigmaQcan`, `omegaSym`, etc. Do NOT name any intermediate `nqGeneral`, `nqCanonical`, `k0`, `k2`, `k4`, `gamma5`.

**Claim manifest** (this is a re-derivation, not a missing-script, but listing required physical results explicitly):
- M1: `\chi_Q = 1` follows from matching the `\omega^5` coefficient of `\widehat Y_Q^{\rm ret}` against `\chi_Q · 9·\overline K_0/(32·\Omega_Q^5)` at the canonical `\overline\Gamma_5 = 2G/(5c^5)`.
- M2: `N_Q = 1` follows from `\widehat m_0^2 · \chi_Q · N_Q = 1` at `\widehat m_0 = 1, \chi_Q = 1`.
- M3: The next odd term in `\widehat Y_Q^{\rm ret}` begins at `O(\omega^7)`, so the point-particle 2.5PN coefficient is exhausted by the `\omega^5` term.
- M4: The canonical invariant coefficients `\overline K_0 = 54Gc_s^5/(5a^5c^5)`, `\overline K_2 = 6Gc_s^3/(5a^3c^5)`, `\overline K_4 = 8Gc_s/(15ac^5)`, `\overline\Gamma_5 = 2G/(5c^5)` satisfy the canonical-even identity `K_0 K_4 = 4 K_2^2` and `Gamma5 = 9·sqrt(K2^5/K0^3)`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 106` and confirm: (1) the script exits 0; (2) the script's variable names and assertion sequence are visibly different from the SymPy script (the verifier may diff structurally); (3) the bottom-line printed result still says `N_Q = 1` on the canonical point-particle branch and `\gamma_{\rm quad}^{\rm eff} = 2G/(5c^5)`.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl`
- summary: Re-authored the Mathematica audit around the retarded one-pole `Yret` omega-series witness and independent target-normalization checks.
- deviation: none

## F4 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py:61-65` (add new assertion immediately after the `canonical gamma_eff - target` check)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl:66-67` (add new assertion immediately after the `canonical gamma_eff - target` check; apply after F3 re-authoring)

**Issue:** The two non-tautological assertions (`NQ_canonical at m0hat=1 equals 1`, and `gamma_eff_canonical equals target`) both collapse to substitution identities and provide weak evidence. The off-canonical sensitivity `\Delta_Q := \chi_Q - 1` (appendix eq. app-part04-DeltaQ-def, line 293), which "measures the entire reduced failure at leading order" per appendix line 342, is never exercised. Adding a first-order sensitivity check exercises the algebra non-trivially.

**Required change:**

In the SymPy script, after line 65, append:

```python
Delta_Q = sp.symbols("Delta_Q", real=True)
gamma_eff_off = (m0hat**2 * Gamma5).subs([(m0hat, 1), (chi_Q, 1 + Delta_Q)])
gamma_eff_series = sp.series(gamma_eff_off, Delta_Q, 0, 2).removeO()
linear_coeff = sp.simplify(gamma_eff_series.coeff(Delta_Q, 1))
expect_zero(
    "Delta_Q first-order sensitivity coefficient",
    linear_coeff - (-2 * G / (5 * c**5)),
)
zeroth_coeff = sp.simplify(gamma_eff_series.coeff(Delta_Q, 0))
expect_zero(
    "Delta_Q zeroth-order coefficient equals Gamma5_target",
    zeroth_coeff - Gamma5_target,
)
```

In the Mathematica script (apply after F3 re-authoring; insert at the analogous location near the bottom, before the `RESULT` block), add the same content with Mathematica syntax:

```mathematica
gammaEffOff = (m0hat^2*gamma5) /. {m0hat -> 1, chiQ -> 1 + DeltaQ};
gammaEffSeries = Normal[Series[gammaEffOff, {DeltaQ, 0, 1}]];
linearCoeff = FullSimplify[Coefficient[gammaEffSeries, DeltaQ, 1], Assumptions -> $Assumptions];
expectZero["Delta_Q first-order sensitivity coefficient",
  linearCoeff - (-2*G/(5*c^5))];
zerothCoeff = FullSimplify[Coefficient[gammaEffSeries, DeltaQ, 0], Assumptions -> $Assumptions];
expectZero["Delta_Q zeroth-order coefficient equals Gamma5_target",
  zerothCoeff - gamma5Target];
```

Self-test (already in report's Self-test notes):
- N_Q on natural branch (m0hat=1) at chi_Q = 1 + Delta_Q expands as `1/(1+Delta_Q) = 1 - Delta_Q + Delta_Q^2 - ...`
- gamma_eff = m0hat^2 * N_Q * Gamma5_target = (1 - Delta_Q + ...) * (2G/(5c^5))
- Zeroth-order coefficient = 2G/(5c^5) = Gamma5_target ✓
- First-order coefficient = -2G/(5c^5) ✓
- A wrong N_Q (sign flip, factor error) would change the first-order coefficient, so the assertion is sensitive.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 106` and `redteam exec-mathematica 106` and confirm: (1) both scripts still exit 0; (2) two new assertion lines appear in each output transcript with residuals = 0; (3) the new printed lines reference `Delta_Q first-order sensitivity coefficient` and `Delta_Q zeroth-order coefficient equals Gamma5_target`.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage106_canonical_outgoing_reduced_closure_mathematica_audit.wl`
- summary: Added the first-order `Delta_Q` sensitivity and zeroth-order target checks to both script witnesses.
- deviation: none
