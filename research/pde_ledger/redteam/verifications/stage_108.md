---
unit_id: 108
batch: IV.2
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: true
---

# Verification — unit 108

## Per-finding outcomes

### F1 — paper_misalignment (script_missing_paper_claim)

**Classification:** resolved

**What changed:**
Both engines now carry a full "Class D" block that exercises the general β-parameterized preservation submanifold.
- `scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:72-106` adds `Lambda_gen = S*Lambda_out(beta*z) + Sigma0 + Sigma2 z^2 + Sigma4 z^4 + I Sigma5 z^5`, re-solves `(Sigma2, Sigma4)` against the canonical even-moment fingerprint (the solutions now carry β), prints β-dependent forms, then asserts `general preservation submanifold = S(1 - beta^5)/9 - Sigma0/27` (L98-101), the residual `chi_gen(Sigma5 := chi_pres_gen) - 1 = 0` (L102), and the β=1 reduction `chi_pres_gen|_{β=1} - (-Sigma0/27) = 0` (L103-106).
- `mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:83-112` mirrors the same construction with independent variable names (`lambdaGen`, `solG`, `sigma5PresGen`), printing β-dependent `Sigma2(beta)`, `Sigma4(beta)`, `chi_gen(beta)`, and asserting the same three identities (L104-112).

**Assessment:**
Correct and non-tautological. The locus identity is `chi_pres_gen - (S(1 - β^5)/9 - Sigma0/27)` where `chi_pres_gen` is obtained by `solve(chi_gen == 1, Sigma5)`; the assertion would fail if e.g. the `z^5` coefficient in `Lambda_out` were perturbed or if the β-scaling were applied wrong. Exec logs confirm:
- SymPy prints `Sigma2(beta) = -S*beta**2/3 + S/3 - Sigma0/9`, `Sigma4(beta) = -S*beta**4/9 + S/9 - Sigma0/27`, `chi_gen(beta) = 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0)`, `Sigma5 general preservation locus = -S*beta**5/9 + S/9 - Sigma0/27`, all three Class D residuals = 0.
- Mathematica prints the algebraically equivalent factored forms (`(-sigma0 - 3*(-1 + beta^k)*sNorm)/k!`-style) and the same three residuals = 0 with explicit PASS lines.
The β=1 reduction assertion (`general locus reduces to beta=1`) is the requested sanity check that Class D collapses to Class C. The existing β=1 Class C block remains intact as an independent sanity reduction. No collateral edits beyond the Class D block.

### F2 — insufficient_verification (Mathematica precedence bug)

**Classification:** resolved

**What changed:**
`mathematica/.../stage108...wl:55` now reads `expectZero["chi_arg(beta=1) - 1", (chiArg /. beta -> 1) - 1];` — exactly the parenthesization the directive required. The `(chiArg /. beta -> 1)` is now evaluated before subtracting 1, so the assertion genuinely tests `chi_arg|_{β=1} - 1`.

**Assessment:**
Correct one-line edit. The exec log line 18 reads `chi_arg(beta=1) - 1 = 0` followed by `PASS: chi_arg(beta=1) - 1` (line 19). Since `chiArg = beta^5`, the residual is `1 - 1 = 0` rather than the previously-tautological `0^5 = 0`. Non-tautological: if `chiArg` were `2*beta^5` (a hypothetical perturbation), the residual would now be `1` and fail; under the old parsing it would still be `0`. The SymPy counterpart at L48 was already correct and is unchanged.

### F3 — stale_output (banner mismatch)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:25` now reads `banner('STAGE 108 — ROBUSTNESS CLASSES FOR chi_Q')`.
- `mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:26` now reads `banner["STAGE 108 — ROBUSTNESS CLASSES FOR chi_Q"];`.

**Assessment:**
Both banners now print `STAGE 108 — ROBUSTNESS CLASSES FOR chi_Q` in the transcripts (line 3 of each `.txt`). Em-dash preserved, surrounding whitespace unchanged.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `pure scale invariance = 0`
- `chi_arg(beta=1) - 1 = 0`
- `chi_add - 3(S+9Sigma5)/(3S-Sigma0) = 0`
- `general preservation submanifold = S(1 - beta^5)/9 - Sigma0/27 = 0`
- `general preservation locus check = 0`
- `general locus reduces to beta=1 (Class C) = 0`

**Mathematica:** exit=0. Notable lines:
- `PASS: pure scale invariance`
- `PASS: chi_arg(beta=1) - 1` (now genuine substitution at β=1)
- `PASS: chi_add - 3(sNorm + 9 sigma5)/(3 sNorm - sigma0)`
- `PASS: general preservation submanifold = S(1 - beta^5)/9 - sigma0/27`
- `PASS: general preservation locus check`
- `PASS: general locus reduces to beta=1 (Class C)`
- Final: `Stage 108 Mathematica audit passed.`

**Output freshness:** confirmed. SymPy output mtime 1779916684 > script mtime 1779916168; Mathematica output mtime 1779917061 > script mtime 1779916185. Both transcripts post-date the script edits.

## Material-change assessment

`material_change`: true.

The scripts now verify a strictly stronger claim (general β-parameterized preservation submanifold) than before. The original `chi_Q = 3(S + 9 Sigma_5)/(3S - Sigma_0)` Class-C result is unchanged; the β=1 reduction `Sigma_5 = -Sigma_0/27` is unchanged. The new general locus `Sigma_5 = S(1 - β^5)/9 - Sigma_0/27` is the unifier the notes claimed. No downstream unit's *result* changes; this is a coverage extension, not a correction.

Downstream units that might reference the chi_Q robustness classification should be flagged as `upstream_stale: true` per orchestrator policy, but the user should not expect any of them to actually break: the new check strengthens the existing alignment, it does not invalidate any prior numerical or symbolic identity.

## Side observations (non-blocking)

- Mathematica's `Sigma4(beta=1) = (3*sigma0^2 - 9*sigma0*sNorm)/(-81*sigma0 + 243*sNorm)` prints in an un-cancelled form (transcript L21); the subsequent `Sigma4(beta=1) + sigma0/27 = 0` PASS line confirms it simplifies to `-sigma0/27`, so this is a cosmetic display artifact only — not a defect. No action needed.
- Both engines' Class C and Class D blocks now have parallel pipelines; if future maintenance simplifies, the β=1 block could be removed as redundant with `D|_{β=1}`. Out of scope for this verification.

## Verdict justification

All three findings are resolved. F1's Class D block is a substantive, non-tautological extension that verifies the notes' boxed general preservation submanifold in both engines, with an explicit β=1 reduction sanity check tying it back to Class C. F2's parenthesization fix makes the Mathematica β=1 check discriminating rather than vacuously true. F3's banner strings now match the unit number. Both exec logs exit 0 with all PASS lines printed; output mtimes confirm the transcripts were regenerated post-fix. Engines remain in agreement on every shared assertion.
