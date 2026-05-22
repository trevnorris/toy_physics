---
unit_id: 033
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 033

## Per-finding outcomes

### F1 — tautological_check (Stage 16.6 gate identity)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py:102-130`: replaced the block that defined `gate_num = simplify(expand((alpha_crit_mic - alpha0_mic) * gate_den))` and asserted `alpha_crit_mic - alpha0_mic - gate_num/gate_den == 0`. The new code computes
  - `gate_diff = sp.cancel(sp.together(alpha_crit_mic - alpha0_mic))`,
  - `gate_num_actual, gate_den_actual = sp.fraction(gate_diff)` (independent of the claim),
  - `den_ratio = sp.simplify(gate_den_actual / gate_den_claim)`,
  - guard `assert den_ratio.is_number` so a residual symbolic factor would fail the check,
  - then reconstructs `gate_num_target = sp.simplify(gate_num_actual / den_ratio)` and asserts the final residual zero.
- `mathematica/...stage033...mathematica_audit.wl:111-131`: mirror edit using `Cancel[Together[...]]`, `Numerator`/`Denominator`, `FullSimplify[gateDenActual/gateDenClaim]`, and a `NumericQ[denRatio]` guard before the final `expectZero`.

**Assessment:**
The new check is non-tautological: `gate_den_actual` is derived from `together(alpha_crit_mic - alpha0_mic)` without ever referencing `gate_den_claim`, so a wrong claimed denominator (e.g. one missing a `varpi^2` factor or with a wrong polynomial in K0) would leave free symbols in `den_ratio`, triggering the `is_number` / `NumericQ` failure. The actual ratio produced is `den_ratio = 9*pi**2` (SymPy) / `-9*Pi^2` (Mathematica, sign artifact of `Numerator/Denominator` factoring), which are `is_number=True`/`NumericQ=True` parameter-free constants. The Pi^2 factor traces to `Delta0 = OmegaU^2*OmegaW^2 - 88*g_R^2/(9*pi^2)`: when `Together` is applied, the `1/(9*pi^2)` factor of Delta0 is cleared from the denominator into the rational scaling. Codex's deviation from "rational" (in the original directive) to "parameter-free numeric constant" is correctly documented in the `## Applied: F1` block and is the right call — the structural property the test needs is "no residual symbols," not "rational coefficient." The final identity `alpha_crit(mic) - alpha_0(mic) - gate_num_target/gate_den_claim = 0` then verifies that, after the constant ratio correction, the claimed denominator is indeed (up to sign and the Pi^2 factor) the canonical denominator of the combined difference.

The directive's stipulated structural property — "the ratio test can fail if the claim is wrong" — holds. No collateral edits outside the cited block.

### F2 — tautological_check (Mathematica k0Onset hardcoded)

**Classification:** resolved

**What changed:**
- `mathematica/...stage033...mathematica_audit.wl:95-97`: replaced
  ```
  k0Onset = FullSimplify[gU^2/OmegaU^2 + kappa0Sq*chi^2/(NQ*delta0^2), ...];
  ```
  with
  ```
  k0OnsetSolutions = Solve[n0Mic == NQ, K0];
  If[Length[k0OnsetSolutions] == 0, fail[...]];
  k0Onset = FullSimplify[K0 /. First[k0OnsetSolutions], Assumptions -> $Assumptions];
  ```

**Assessment:**
Edit is exactly as the directive specified. The saved Mathematica output (line 27) shows `K0_onset` printed as the genuine `Solve` result:
`(7744*gR^4*gU^2*NQ + 72*OmegaU^2*(9*(gR*gU + gW*OmegaU^2)^2 - 22*gR^2*gU^2*NQ*OmegaW^2)*Pi^2 + 81*gU^2*NQ*OmegaU^4*OmegaW^4*Pi^4)/(NQ*(88*gR^2*OmegaU - 9*OmegaU^3*OmegaW^2*Pi^2)^2)`
— a non-trivial polynomial form rather than the hand-stated closed form. The assertion at lines 108-111 then PASSES (`K0_onset - [gU^2/OmegaU^2 + kappa0^2 Chi^2/(NQ Delta0^2)] = 0`), which is now a genuine consistency check between solved value and closed-form target. The back-substitution check at line 107 (`N_-(0) at K0_onset - NQ = 0`) also still passes. No collateral edits.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/...stage033...mathematica_audit.wl:133-160`: inserted the numeric cross-check block exactly as specified in the directive. Two rational substitution rules (`numericRule1`, `numericRule2`) cover Stage 16.1 (`dN - dNFormula`) and Stage 16.6 (`alphaCritMic - alpha0Mic - gateNumTarget/gateDenClaim`) at 30-digit floating-point precision. `Do` loop iterates both rules; each iteration uses `If[Abs[...] > 10^-20, fail[...], pass[...]]`.

**Assessment:**
Insertion location and content match the directive exactly. The saved Mathematica output shows four `PASS` lines for the numeric cross-checks (two assertions × two rules):
- `monotonicity numeric residual = 0``78.83717778732962` → PASS
- `gate-identity numeric residual = 0``78.70733917439959` → PASS
- `monotonicity numeric residual = 0``78.99802761920512` → PASS (rule 2)
- `gate-identity numeric residual = 0``78.57895656247477` → PASS (rule 2)

Each residual is reported as a `0``<precision>` value — a zero with reduced trailing precision, but still numerically `0`, so `Abs[res] > 10^-20` evaluates False and PASS is taken. The `N::meprec` warnings (Internal precision limit `$MaxExtraPrecision = 50` reached) are expected because `dN`/`dNFormula` contain nested `Sqrt[7744*alpha0^2 + ...]` radicals that require extra precision when evaluated at rational points — but the final residual still resolves to zero. This is genuinely structurally distinct from the analytic `FullSimplify` path used earlier in the script, so the cross-check exercises an independent derivation chain (numeric substitution + arbitrary-precision arithmetic vs. symbolic simplification). A latent algebra error in `dNFormula`, `gateNumTarget`, or `gateDenClaim` that happened to cancel under `FullSimplify` would not cancel at the rational test points.

## Exec log assessment

**SymPy:** dedicated `stage_033_sympy.log` was not captured by the orchestrator under `redteam/exec_logs/` (only the diff patch is present). The canonical saved transcript at `scripts/output/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.txt` (mtime 2026-05-21 17:31:40, ~2 minutes newer than the edited script at 17:29:39) is the authoritative log. Notable lines:

- L33: `alpha_crit - closed finite-throat form = 0`
- L58: `K0_onset - [gU^2/OmegaU^2 + kappa0^2 Chi^2/(NQ Delta0^2)] = 0`
- L132: `denominator ratio (must be parameter-free) = 9*pi**2`
- L133: `alpha_crit(mic) - alpha_0(mic) - gate_num_target/gate_den_claim = 0`
- L182: `All Stage 16 checks passed.`

Since the script terminates only by completing `print("All Stage 16 checks passed.")` and no `AssertionError` was raised, sympy_exit is recorded as 0.

**Mathematica:** dedicated `stage_033_mathematica.log` was not captured under `redteam/exec_logs/`. The canonical saved transcript at `mathematica/output/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.txt` (mtime 2026-05-21 17:31:50, ~2 minutes newer than the edited script at 17:29:39) is the authoritative log. Notable lines:

- L29: `PASS: N_-(0) at K0_onset - NQ`
- L31: `PASS: K0_onset - [gU^2/OmegaU^2 + kappa0^2 Chi^2/(NQ Delta0^2)]` (now non-tautological)
- L34: `denominator ratio (must be parameter-free) = -9*Pi^2`
- L35: `PASS: gate denominator matches claim up to parameter-free constant`
- L38: `PASS: alpha_crit(mic) - alpha_0(mic) - gate_num_target/gate_den_claim`
- L52, L66, L85, L86: four PASS lines from the F3 numeric cross-checks
- L89: `Stage 033 Mathematica audit passed.` and the script terminates via `Exit[0]`.

The `N::meprec` warnings noted in the output do not cause failure — Mathematica issues them as informational notes when extra precision is added during evaluation of the radical-bearing expressions; the residuals still evaluate to zero (precision-tagged `0``78.x`). Mathematica_exit = 0 (script reaches `Exit[0]`).

**Output freshness:** confirmed. Edited scripts have mtime 2026-05-21 17:29:39; saved `.txt` outputs have mtime 2026-05-21 17:31:40 (sympy) and 17:31:50 (mathematica). Outputs are 2 minutes newer than the corresponding edited scripts.

## Material-change assessment

`material_change`: false.

The fixes change *how* the gate identity and `k0Onset` are verified, not the verified results themselves. The Stage 16.6 gate denominator is still `8 varpi^2 OmegaU^2 Delta0 (11 A_mic + 9 DeltaK)` (up to the universal `9*pi^2` rescaling absorbed by `Delta0`'s `9*pi^2` factor); the Stage 16.5 `K0_onset` value is still `g_U^2/Omega_U^2 + 648*pi^2*(Omega_U^2*g_W + g_R*g_U)^2 / (NQ*(9*pi^2*Omega_U^2*Omega_W^2 - 88*g_R^2)^2)` (sympy output L57). The numeric cross-check (F3) is purely additional verification. No symbolic results that downstream units (034+) consume have shifted.

## Side observations (non-blocking)

1. The displayed sympy `gate numerator =` block at the end of the output is now `gate_num_target = sp.simplify(gate_num_actual / den_ratio)`, which still equals `(alpha_crit_mic - alpha0_mic) * gate_den_claim` up to simplification. This is fine for display purposes only and does not re-introduce the tautology — the *assertion* is the final `expect_zero(...)` which now hinges on the independently-derived `gate_num_actual / den_ratio` reconciling with `gate_den_claim`.
2. In the Mathematica saved transcript, the F3 numeric residuals are shown as precision-tagged zeros like `0``78.83717778732962`. These are bona fide zeros; the trailing precision number just reflects how many digits of accuracy Mathematica believes it has computed. The `Abs[res] > 10^-20` guard correctly evaluates False on such inputs.
3. The original auditor's directive for F1 asked for `is_number and is_rational`; the actual ratio is `9*pi**2` (transcendental), so Codex correctly relaxed to `is_number` only and documented the deviation in the `## Applied: F1` block. This is the substantively correct relaxation — the test purpose is "no free symbols remain," which `is_number` enforces.
4. No `stage_033_sympy.log` / `stage_033_mathematica.log` artifacts under `redteam/exec_logs/`; only `stage_033_diff.patch` is present. The orchestrator should consider also storing the canonical `.txt` outputs there for future stages, but this is an infrastructure note, not a verification blocker — the freshness check via mtime on the output files confirms the post-fix scripts produced the saved transcripts.

## Verdict justification

All three findings are correctly addressed by non-tautological constructions. F1 replaces a `(x*d)/d == x` identity with an independent denominator extraction via `together`/`fraction` followed by an `is_number`/`NumericQ` guard on the ratio; F2 replaces a hardcoded Mathematica `k0Onset` with a real `Solve[n0Mic == NQ, K0]` inversion mirroring the SymPy path; F3 inserts a structurally distinct numeric cross-check at two rational test points covering both the Stage 16.1 monotonicity identity and the Stage 16.6 gate identity. The saved post-fix transcripts (both fresher than the edited scripts) show every assertion passing, both engines reaching their respective success footers, and the gate-denominator structural claim (`8 varpi^2 Omega_U^2 Delta0 (11 A_mic + 9 DeltaK)` up to the `9*pi^2` Delta0 normalization) genuinely verified. Verdict: `verified`, `material_change: false`.
