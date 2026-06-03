---
unit_id: 220
batch: VII.1
created_at: 2026-06-02T13:10:00-06:00
findings_count: 2
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
applied_at: 2026-06-02T14:03:36-06:00
findings_applied: 2
findings_blocked: 0
---

# Codex directive — unit 220

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. Do NOT touch paper.tex, notes/, or any prose document.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

## F1 — missing_verification_script (missing_mathematica)

**Target:** `mathematica/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_mathematica_audit.wl` (new file)

**Issue:** Unit 220 is `is_status_only_candidate: false` and `is_checkpoint: false`, with only a SymPy engine. Every deliverable — the `3x3` symbolic determinant identity, the six entries of the symbolic inverse, the susceptibility quadratic form, the collinear-source and primitive-source family factorizations, the symbolic `Pi`-derivative identity, the linear-outgoing correction, and the real/imaginary split of the phase-lag term — is fully and independently verifiable in Mathematica. The dual-engine contract therefore requires an independent `.wl`. (The stage card already notes "Mathematica audit: none yet" at `paper/stages/stage_220.tex:11`; that is the gap to close, not a justification for single-engine.)

**Required change:** Create the file at the Target path. It must independently verify the claim manifest below using native Mathematica primitives and a DIFFERENT decomposition from the `.py`. Each claim must be guarded so the script exits nonzero on failure — use an `expectZero[expr_]` / `expectTrue[cond_]` helper that prints a per-claim `PASS`/`FAIL` line and calls `Exit[1]` on failure (mirror the project's other `_mathematica_audit.wl` files). The script must print a clear per-claim PASS line and exit 0 only when all of M1–M9 hold.

**Symbol setup (use the same physical definitions as the paper, NOT the .py's `together` choreography):**
- Real structural constants `K, M, C, varpi, OmegaU, OmegaW, R, GU, GW`; real `omega`; unrestricted `Pi`; real nonnegative `Gamma`; real sources `Jq, JU, JW`, `sq, sU, sW, S`; positive `x, kappa`; real `betaQ, betaU, betaW`.
- `KB = K - M omega^2 - C^2/(varpi^2 - omega^2)`, `Aw = OmegaU^2 - omega^2`, `Ww = OmegaW^2 - omega^2 - Pi`.
- `Kdyn = {{KB, -GU, -GW}, {-GU, Aw, -R}, {-GW, -R, Ww}}`.
- `DeltaPi = Aw Ww - R^2`, `QPi = GU^2 Ww + 2 GU GW R + GW^2 Aw`, `DPi = KB - QPi/DeltaPi`.

**Anti-transliteration guard (mandatory):** Derive each result independently. In particular:
- M1/M3: use `Det[Kdyn]` and `Inverse[Kdyn]` directly (Mathematica's cofactor/adjugate route), NOT a hand-typed transcription of the SymPy closed forms; then compare against the paper's stated forms with `FullSimplify[... == 0]`. The independent direction is "engine inverts the matrix, we check it reproduces the paper's chi entries," the reverse of the `.py` (which builds chi by hand and checks `Kinv - chi`). Both are legitimate but must not be the same keystrokes.
- M6 (spatial families): extract the three coefficients with `Coefficient`/`Collect` on the monomial basis `{x^-6, E^(-2 kappa x)/x^4, E^(-4 kappa x)/x^2}` (or `SeriesCoefficient` after the substitution), and ALSO assert no OTHER power of `x` survives — do NOT reproduce the SymPy `together` then `simplify(diff)==0` regrouping verbatim. Verifying the *absence* of any 5th family is the substantive independent content.
- M9 (phase-lag): use `ComplexExpand[Re[...]]` / `ComplexExpand[Im[...]]` to split `-I/2 Gamma TJ0^2`, NOT `as_real_imag`. A line-by-line port of the SymPy variable choreography (same `together`/`inv` sequence rewritten in WL syntax) is REJECTED as transliteration.

**Claim manifest** (the new `.wl` must independently verify each):

- **M1 — determinant identity.** `FullSimplify[Det[Kdyn] - DeltaPi*DPi] == 0`. (Compute `Det[Kdyn]` natively; `DPi` carries the `1/DeltaPi`, so clear denominators or `Together` before `FullSimplify`.)
- **M2 — static reduction.** Substituting `omega -> 0, Pi -> 0`: `KB -> K - C^2/varpi^2` (call it `Kstar`), `Aw -> OmegaU^2`, `Ww -> OmegaW^2`, `DeltaPi -> OmegaU^2 OmegaW^2 - R^2`, `QPi -> GU^2 OmegaW^2 + 2 GU GW R + GW^2 OmegaU^2`, and `DPi -> Kstar - (that Q)/(that Delta)`. Assert each difference `FullSimplify`s to 0. (This reproduces the carried Stage-219 one-port bundle.)
- **M3 — inverse entries.** Let `Kinv = Inverse[Kdyn]`. Assert `FullSimplify` of each of the six paper entries minus the corresponding `Kinv` entry is 0: `chi_qq = 1/DPi` (Kinv[[1,1]]); `chi_qU = (GU Ww + R GW)/(DeltaPi DPi)` (Kinv[[1,2]]); `chi_qW = (Aw GW + R GU)/(DeltaPi DPi)` (Kinv[[1,3]]); `chi_UU = (KB Ww - GW^2)/(DeltaPi DPi)` (Kinv[[2,2]]); `chi_UW = (KB R + GU GW)/(DeltaPi DPi)` (Kinv[[2,3]]); `chi_WW = (KB Aw - GU^2)/(DeltaPi DPi)` (Kinv[[3,3]]).
- **M4 — susceptibility law.** With `Jvec = {Jq, JU, JW}`, `Vmix = -1/2 (Jvec . Kinv . Jvec)`. Assert `FullSimplify[Vmix - VmixExpected] == 0` where `VmixExpected = -1/2 (chi_qq Jq^2 + 2 chi_qU Jq JU + 2 chi_qW Jq JW + chi_UU JU^2 + 2 chi_UW JU JW + chi_WW JW^2)`.
- **M5 — collinear-source factorization.** With `Jcol = S {sq, sU, sW}`, assert `FullSimplify[(-1/2 (Jcol . Kinv . Jcol)) + 1/2 chiS S^2] == 0`, where `chiS = Ns/(DeltaPi DPi)` and `Ns = DeltaPi sq^2 + 2 (GU Ww + R GW) sq sU + 2 (Aw GW + R GU) sq sW + (KB Ww - GW^2) sU^2 + 2 (KB R + GU GW) sU sW + (KB Aw - GU^2) sW^2`.
- **M6 — primitive product-family theorem.** With `Jprim = {betaQ x^-3, betaU E^(-2 kappa x)/x, betaW E^(-2 kappa x)/x}`, form `Vprim = -1/2 (Jprim . Kinv . Jprim)`. Define `C6 = chi_qq betaQ^2`, `C4 = chi_qU betaQ betaU + chi_qW betaQ betaW`, `C2 = chi_UU betaU^2 + 2 chi_UW betaU betaW + chi_WW betaW^2`. (i) Assert `FullSimplify[Vprim - (-1/2 (C6 x^-6 + 2 C4 E^(-2 kappa x)/x^4 + C2 E^(-4 kappa x)/x^2))] == 0`. (ii) INDEPENDENT CONTENT: after substituting `y = E^(-2 kappa x)`, assert that `Vprim` is, in `x` and `y`, supported ONLY on `{x^-6 y^0, x^-4 y^1, x^-2 y^2}` — i.e. the coefficient of every other `x^a y^b` monomial in the collected form is 0. Implement by `Collect`/`CoefficientRules` and `expectTrue` that the support set equals exactly those three monomials. This is the "no new spatial family" claim and must NOT be a re-statement of (i).
- **M7 — outgoing-port derivative identity.** Define `TJ = chi_qW Jq + chi_UW JU + chi_WW JW` (the `e_W^T Kinv J` transfer factor). Assert `FullSimplify[D[VmixExpected, Pi] + 1/2 TJ^2] == 0`. (Verify `VmixExpected` genuinely depends on `Pi` through `Ww` before trusting the derivative — see self-test.)
- **M8 — linear outgoing correction.** Let `TJ0 = TJ /. Pi -> 0`. Assert `FullSimplify[(D[VmixExpected, Pi] /. Pi -> 0) I Gamma - (-1/2 I Gamma TJ0^2)] == 0`, i.e. `deltaV^(1) = -I/2 Gamma TJ0^2`.
- **M9 — phase-lag no-go + absorbed power.** Let `dV1 = -1/2 I Gamma TJ0^2`. Using `ComplexExpand` (all symbols real, `Pi` already removed): assert `ComplexExpand[Re[dV1]] === 0` (general, symbolic). Define `Pabs = -omega ComplexExpand[Im[dV1]]`. Assert `FullSimplify[Pabs - omega Gamma/2 TJ0^2] == 0` (the perfect-square / always-dissipative form, giving `Pabs >= 0` for `omega, Gamma >= 0`). Optionally confirm on the same numeric off-pole sample used by the `.py` (`K->11, M->2, C->1, varpi->5, OmegaU->3, OmegaW->4, R->2, GU->1, GW->2, omega->1/2, Jq->1, JU->2, JW->1, Gamma->1/10`) that `DeltaPi != 0`, `DPi != 0`, `Re[dV1] == 0`, `Im[dV1] != 0`, `Pabs > 0`.

**Verification command:** After Codex applies, the verifier runs `redteam exec-mathematica 220` and confirms the new `.wl` exists at the Target path, contains M1–M9 as guarded checks, and exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_mathematica_audit.wl`
- summary: Created an independent Mathematica audit covering M1 through M9 with guarded PASS/FAIL checks and hard failure on any residual.
- deviation: none

## F2 — insufficient_verification

**Target:** `scripts/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.py:195`

**Issue:** The notes (deliverable 9, lines 441–443) and appendix eq. app-part07-phase-lag-no-real (lines 456–457) state the absorbed-power result as the general inequality `P_abs^(1) = (omega Gamma / 2) T_J^2 >= 0`. The script verifies the `Re = 0` half symbolically at line 194 (good), but the non-negativity half is exercised ONLY at the single numeric sample `P_abs_sample > 0` (line 239). The symbolic `P_abs` (line 195) is `omega*Gamma*(num)^2/(2*(den)^2)`, so its non-negativity for all admissible `omega, Gamma >= 0` follows from it being `omega*Gamma/2` times a perfect square — but the script never asserts that structural fact. One sample does not exercise the general "off-pole, any drive" claim.

**Required change:**
Add ONE symbolic assertion immediately after the `P_abs = sp.factor(...)` line (line 195), pinning `P_abs` to its perfect-square form:

```python
    # General (not sample-only) non-negativity: P_abs = (omega*Gamma/2) * T_J0^2,
    # i.e. omega*Gamma times a perfect square, hence >= 0 for omega, Gamma >= 0.
    assert sp.simplify(P_abs - omega * Gamma / 2 * T_J0**2) == 0
```

Do NOT remove or alter the existing numeric sample checks at lines 235–239; this assertion is added alongside them. Do not introduce new symbols (`omega`, `Gamma`, `T_J0` already exist).

Self-test confirmation (already done by auditor): at the script's sample point this gives `(1/2)(1/10)/2 * T_J0^2 = 0.0003378`, exactly matching the printed `P_abs(sample)`, so the assertion is true and the residual `simplify`s to 0.

**Verification command:** The verifier runs `redteam exec-sympy 220` and confirms a new symbolic assertion `P_abs == omega*Gamma/2 * T_J0^2` appears near line 195 and the script exits 0 with all checks passing.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.py`
- summary: Added the symbolic perfect-square assertion for `P_abs = omega*Gamma/2*T_J0**2` immediately after the absorbed-power expression.
- deviation: none

## Authorized notes renumber (USER-AUTHORIZED 2026-06-02)

The user authorized renumbering stale VII.1 notes-prose stage labels to canonical. The audit logged a notes-side renumbering drift here (notes reference pre-renumber stage numbers; the card/appendix/scripts use canonical). Notes-only cleanup in THIS fix loop (Codex applies notes/*.md; Claude reviews):
- In `notes/stages/moving_throat_pde_stage220_..._sympy_audit.md`, renumber every stale stage-number reference to match the canonical numbering used in this stage's SymPy script comments + the paper card (self-reference → Stage 220; cited upstream stages → the numbers the .py comments use). Math/content unchanged.
- Do NOT touch scripts, paper.tex, or appendix. Acceptance: notes stage labels match the .py comments + card. Append `## Applied: notes-renumber`.

## Applied: notes-renumber

- files_changed:
  - `notes/stages/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.md`
- summary: Renumbered stale notes-side VII.1 labels so Stage 237 self-references now read Stage 220 and Stage 253 upstream static-bundle references now read Stage 219.
- deviation: none
