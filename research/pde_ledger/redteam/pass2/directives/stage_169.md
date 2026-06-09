---
unit_id: 169
batch: V.1
created_at: 2026-06-08T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-08T22:33:11Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 169

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected script (`math -script <path>`) and iterate until it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl:96-119`

**Issue:** The `.wl` is a line-by-line port of the `.py` (identical block choreography, identical hand-typed weight expressions, and — most importantly — the identical hand-typed decimal literals `rNum = 1.77799353547498` and `gNum = 0.758035078944663` that drive the load-bearing Stage-168 transport coefficients). Because both engines source the three Family-1 coefficients from the same transcribed decimals, the second engine provides no independent check on exactly the quantities most likely to carry a transcription error. The fix makes the Mathematica numeric evaluation independent by deriving the Family-1 radius `r_*` from its canonical closed form instead of pasting the SymPy decimal. The rest of the script (the integral-based `sphereAvg`, the structural `iXY = (dX . gMat . dY)/5`, the symbolic factorization check) is already structural and may stay.

**Required change:**
In the "Stage 168 weighted transport" block, replace the hand-typed decimal radius with the canonical Family-1 closed form so the coefficients are computed, not transcribed. The Family-1 radius is `r_* = Sqrt[4107 - 100*Pi^2]/(10*Pi)` (numerically `1.77799353547498`).

Before (current `wl` L101-103):
```
rNum = SetPrecision[1.77799353547498, 30];
gNum = SetPrecision[0.758035078944663, 30];
Print["Numeric Xi_perp combination = ", fmt[N[xiPerp /. {g -> gNum, r -> rNum}, 20]]];
```
After:
```
rExact = Sqrt[4107 - 100*Pi^2]/(10*Pi);   (* canonical Family-1 radius *)
rNum = N[rExact, 30];
gNum = SetPrecision[0.758035078944663, 30];
Print["Numeric Xi_perp combination = ", fmt[N[xiPerp /. {g -> gNum, r -> rNum}, 20]]];
```

Leave the three paper-comparison target literals (`0.758035078944663`, `1.00314310113848`, `1.88373219118005`) in the `checks` table at L110-112 exactly as-is — they are the paper-side values the script compares AGAINST, and must remain hand-typed so the check is a genuine paper-comparison.

Note: `gNum` (= `g_*`) has no simpler closed form available in this card's inputs, so it legitimately remains a decimal; this directive only de-transcribes the radius `r_*`. If you can source `g_*` symbolically from an upstream Stage-16x definition without inventing a value, you may, but do not fabricate a closed form.

**Self-test (already done by auditor):** `Sqrt[4107 - 100*Pi^2]/(10*Pi) = 1.77799353547498…` (hand-checked: `100*Pi^2 ≈ 986.9604`, `4107 - 986.9604 = 3120.0396`, `Sqrt ≈ 55.857`, `/31.4159 ≈ 1.77799`). With `g_*=0.758035` and `s=Sqrt[1+r_*²]≈2.03992`: coeff_T `= g_* = 0.758035`, coeff_v `= g_*+1/(2s) = 1.003143`, coeff_L `= 2g_*+3/(4s) = 1.883731` — all still pass the `1e-12` checks. The substitution changes nothing numerically; it only removes the duplicated transcription.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 169` and confirms: (a) the bare decimal `1.77799353547498` no longer appears as the radius source (replaced by `Sqrt[4107 - 100*Pi^2]/(10*Pi)`), (b) the three coefficient checks still print `PASS`, and (c) the script exits 0. The SymPy script and its output are unchanged.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl`
- summary: Replaced the transcribed Family-1 radius decimal with the canonical closed form before numeric evaluation.
- deviation: none
