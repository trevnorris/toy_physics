---
unit_id: 168
batch: V.1
created_at: 2026-06-08T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-08T16:33:04-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 168

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected script (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. Getting the script to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl:71`

**Issue:** The `.wl` is a verbatim line-for-line port of the `.py` and adds no independent route to the load-bearing reduction result `δ_⊥ = -ε_⊥`. In particular the boxed slippage combination `epsPerp` is hand-typed identically in both engines (py:77, wl:71), so a transcription slip in those weights would be copied into both and pass in both. The achievable independence is limited (pure substitution-and-cancellation of imported carry-forward formulas — no Solve/Series/matrix route exists), so this is a LOW-severity policy-fidelity flag, not a correctness flag; the engines already agree exactly (all residuals 0). The fix is to make the Mathematica path to `ε_⊥` *derive* the slippage weights instead of transcribing them.

**Required change:**
Replace the hand-typed `epsPerp` literal at wl:71 with a quantity DERIVED by reading the `{epsT, epsv, epsL}` weights out of `-deltaPerpSlip` (which is already built upstream at wl:64-70 by substituting the slippages into `deltaPerp`), then add a one-line assertion that the derived weights equal the boxed Family-1 forms. Concretely:

1. At wl:71, replace
   ```
   epsPerp = g*epsT + (g + bCoeff)*epsv + cCoeff*epsL;
   ```
   with a derived definition that extracts the coefficients of the three slippages from the negated substituted normal coordinate:
   ```
   wT = Coefficient[Expand[-deltaPerpSlip], epsT];
   wv = Coefficient[Expand[-deltaPerpSlip], epsv];
   wL = Coefficient[Expand[-deltaPerpSlip], epsL];
   epsPerp = wT*epsT + wv*epsv + wL*epsL;
   ```
2. Immediately after, add an `expectZero` that confirms the DERIVED weights match the boxed Family-1 forms from notes §3 (so a wrong upstream `deltaPerp`/transport law would now be caught by an independent coefficient read, not by re-typing the same target):
   ```
   expectZero["epsPerp weights match boxed form (g, g+1/(2s), 2g+3/(4s))",
     (wT - g)^2 + (wv - (g + bCoeff))^2 + (wL - (2*g + 3/(4*s)))^2];
   ```
   (Sum-of-squares so a single `expectZero` covers all three; FullSimplify under the existing `$Assumptions` reduces it to 0.)
3. Leave the existing `expectZero["delta_perp + eps_perp", deltaPerpSlip + epsPerp]` at wl:73 unchanged — it now compares the substituted coordinate against an INDEPENDENTLY DERIVED `epsPerp`, so it is no longer transcription-vs-transcription.
4. In the "Numeric Family-1 coefficients" block, de-transcribe the hand-typed decimal radius so it is computed, not pasted from the SymPy literal. Replace `rNum = SetPrecision[1.77799353547498, 30];` (wl:100) with the canonical Family-1 closed form:
   ```
   rExact = Sqrt[4107 - 100*Pi^2]/(10*Pi);   (* canonical Family-1 radius = 1.77799353547498... *)
   rNum = N[rExact, 30];
   ```
   Leave `gNum = SetPrecision[0.758035078944663, 30];` as a decimal (no simpler closed form available in this card's inputs; do not fabricate one). This mirrors the de-transcription applied at stage 169.
5. Leave the rest of the script (mouth-bias, outlet, carry-forward blocks) unchanged.

**Self-test (already done by auditor):** From the script's own printed `delta_perp with slippages = -2*epsL*g - epsT*g - epsv*g - (3*epsL)/(4*s) - epsv/(2*s)`, negating gives `epsT`→`g`, `epsv`→`g + 1/(2s)` (since `1/(2s)=bCoeff`), `epsL`→`2g + 3/(4s)`. So `wT=g`, `wv=g+bCoeff`, `wL=2g+3/(4s)=cCoeff` — the new `expectZero` residual reduces to 0, and the existing `delta_perp + eps_perp` check still gives 0. No new constant is introduced (the weights are the same boxed forms the notes state), so no new paper_misalignment.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 168` and confirms (a) wl:71 no longer holds the literal `epsPerp = g*epsT + (g + bCoeff)*epsv + cCoeff*epsL;` verbatim, (b) the new "epsPerp weights match boxed form" `expectZero` PASSES, (c) all seven original `expectZero` checks still PASS, (d) the bare decimal `1.77799353547498` no longer appears as the radius source (replaced by `Sqrt[4107 - 100*Pi^2]/(10*Pi)`) and the numeric coefficient prints are unchanged, and (e) the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl`
- summary: Derived `epsPerp` from coefficients of `-deltaPerpSlip`, added the boxed-weight check, and computed the Family-1 radius from its closed form.
- deviation: none
