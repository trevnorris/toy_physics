---
unit_id: 177
batch: V.2
created_at: 2026-05-30T06:42:53Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-30T07:09:53Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 177

Apply the non-paper_misalignment finding below (only F1 carries edits). After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F2 and F3 are informational only (see report `/var/projects/toy_physics/research/pde_ledger/redteam/reports/stage_177.md`): the missing grouped-collapse check and the prefactor-slope identification are trivial-by-linearity / definitional, and any added assertion would be tautological. Do NOT add a `sum_r rho_r sigma_r == Xi_1` check and do NOT alter the prefactor block. No edits for F2/F3.

If F1's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes beyond what F1 names. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN both affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

## F1 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_mathematica_audit.wl:26,62`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.py:32`

**Issue:** The `.wl` is a line-by-line transliteration of the `.py` (identical variable choreography, perturbation construction, `sigma_r`, hardcoded lane pattern, projection weights, and byte-identical carry-forward prints; only the slope-extraction primitive differs). This defeats the second-engine independence policy: a setup error in the shared algebra (e.g. a wrong `lambda0` definition relative to the `(M,I,H)` decomposition) would pass in both engines because both encode the identical expression. Additionally both scripts mislabel the banner as "STAGE 160".

**Required change:**

(1) Add an independent, non-shared check to the Mathematica script that verifies the load-factor factorization the appendix actually asserts (eq:app-part05-load-factor-factorization, `Lambda_r^2/K = M_r^2 (1+I_r)^2/(1-H_r)^2`) — a fact neither script currently verifies and which is the bridge between the raw `lambda0` and the `(M,I,H)` decomposition. Insert immediately AFTER the `expectZero["weak-axisymmetric d ln H", ...]` line (currently wl:59), before the `banner["Portwise outgoing-defect amplitude"]` line (currently wl:61):

```
banner["Load-factor factorization (independent check)"];
expectZero["load-factor factorization lambda0^2/k = M^2 (1+I)^2/(1-H)^2",
  lambda0^2/k - mCal^2*(1 + iCal)^2/(1 - hCal)^2];
```

This is non-tautological: `lambda0 = (ou2*gw + r*gu)/(ou2*ow2 - r^2)`, `mCal = gw/(Sqrt[k]*ow2)`, `iCal = r*gu/(ou2*gw)`, `hCal = r^2/(ou2*ow2)` are already defined at wl:32-35; the identity holds only because the algebra works out (verified by hand in the report Self-test notes), and would fail if `lambda0` were defined inconsistently with `(M,I,H)`. Keep the existing checks; this is additive.

(2) Re-anchor the per-port amplitude line (currently wl:62) so the Mathematica engine derives the slope from the FACTORED form rather than echoing the SymPy raw-`lambda0` route. Replace the line

```
sigmaExact = FullSimplify[D[Log[(lambdaP^2/kP)/(lambda0^2/k)], eps] /. eps -> 0, Assumptions -> $Assumptions];
```

with a derivation that builds the perturbed load factor from the perturbed `(M,I,H)` factors:

```
lambdaSqOverKP = mCalP^2*(1 + iCalP)^2/(1 - hCalP)^2;
lambdaSqOverK0 = mCal^2*(1 + iCal)^2/(1 - hCal)^2;
sigmaExact = FullSimplify[D[Log[lambdaSqOverKP/lambdaSqOverK0], eps] /. eps -> 0, Assumptions -> $Assumptions];
```

where `mCalP, iCalP, hCalP` are the already-defined perturbed factors (wl:44-46). This makes the Mathematica per-port check derive `sigma_r` through the `(M,I,H)` factorization (a genuinely different intermediate route from the SymPy script's raw `lambdaP^2/kP`), while still comparing against the same `lam*sigmaR`. The result must still be `0`.

(3) Fix the banner string in BOTH scripts. In `.../mathematica/...wl:26` change `banner["STAGE 160 — WEAK-AXISYMMETRIC OUTGOING-SLIPPAGE COLLAPSE"];` to `banner["STAGE 177 — WEAK-AXISYMMETRIC OUTGOING-SLIPPAGE COLLAPSE"];`. In `.../scripts/...sympy_audit.py:32` change `banner("STAGE 160 — WEAK-AXISYMMETRIC OUTGOING-SLIPPAGE COLLAPSE")` to `banner("STAGE 177 — WEAK-AXISYMMETRIC OUTGOING-SLIPPAGE COLLAPSE")`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 177` and `redteam exec-sympy 177` and confirm: (a) the new `load-factor factorization ... = 0` / `PASS` line appears in the Mathematica transcript; (b) `Sigma_{A,r} = lambda_A sigma_r = 0` still passes via the factored-form derivation; (c) both transcripts print the "STAGE 177" banner; (d) both scripts exit 0.

## Applied: F1

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_mathematica_audit.wl`
  - `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.py`
- summary: Added the Mathematica load-factor factorization check, re-anchored the per-port slope through the factored `(M,I,H)` form, and corrected both Stage 177 banners.
- deviation: none

## F2 — insufficient_verification (informational, no edit)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:91-100`, `/var/projects/toy_physics/research/pde_ledger/mathematica/...mathematica_audit.wl:66-74`

**Issue:** The headline grouped collapse `Xi_1 = sum_r rho_r^(N) sigma_r` is not formed in-script (`Xi1` is an abstract free symbol). Per Self-test trap 3 this step is trivial-by-linearity and any added check would be tautological.

**Required change:** None. Do NOT add a collapse assertion.

## F3 — insufficient_verification (informational, no edit)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/...sympy_audit.py:111-122`, `/var/projects/toy_physics/research/pde_ledger/mathematica/...mathematica_audit.wl:85-95`

**Issue:** The prefactor-slope block asserts `P1/P0 - (n1-d1) == 0`, which is the definitional log-of-quotient identity and cannot fail. Strengthening it requires re-deriving Stage 241 (out of scope) or a tautological sum (F2).

**Required change:** None. Do NOT alter the prefactor block.
