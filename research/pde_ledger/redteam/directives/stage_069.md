---
unit_id: 069
batch: III.3
created_at: 2026-05-22T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-22T20:04:22-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 069

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py:74-90`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py:99-110`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl:57-73`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl:81-94`

**Issue:** Every `expect_zero` / `expectZero` assertion in these ranges reduces to an algebraic identity that follows directly from the definitions on the preceding lines (sympy:54-66, wl:35-50). In particular, `Cres2 := 1/Pres`, `Wfail_res := Pres * Wfail_match`, and `delta_fail := Wfail_res - Wfail_match` make assertions like `Pres - 1/Cres2 == 0`, `Wfail_res - Pe_req/(Cres2*Deltainf) == 0`, and `delta_fail - Pe_req*(1-Cres2)/(Cres2*Deltainf) == 0` tautological. A checkpoint stage requires substantive consolidation checks, not algebraic self-identities.

**Required change:**

In the SymPy script, before line 74 (after the print block at line 71), add a substantive derivation block:

```python
# Stage-49 matched-window derivation as a generating function.
# Encodes the theorem "W_match(Delta_eff) = Pe_req / Delta_eff" with
# Delta_eff ranging over [Delta_0, Delta_inf] on the matched branch.
Delta_eff = sp.symbols("Delta_eff", positive=True, real=True)
W_match_generator = Pe_req / Delta_eff
expect_zero(
    "matched fail edge from W_match(Delta_inf)",
    W_match_generator.subs(Delta_eff, Deltainf) - Wfail_match,
)
expect_zero(
    "matched success edge from W_match(Delta_0)",
    W_match_generator.subs(Delta_eff, Delta0) - Wsuff_match,
)
expect_positive(
    "W_match decreasing in Delta_eff",
    -sp.diff(W_match_generator, Delta_eff) * Delta_eff**2 / Pe_req,
)
```

The third assertion's expression simplifies to `1`, which is provably positive — it confirms that `dW_match/dDelta_eff = -Pe_req/Delta_eff^2 < 0`, anchoring the monotonicity of the matched-window edges.

Then, before line 81 (after the existing line 80 assertion), add the resonance penalty derivation:

```python
# Stage-51 / Stage-68 resonance penalty as a band-edge ratio.
# Encodes the theorem "P_res = W_fail_res / W_fail_match" rather than
# defining P_res and C_res^2 to satisfy it by construction.
Pres_from_ratio = sp.simplify(Wfail_res / Wfail_match)
expect_zero(
    "P_res from band-edge ratio matches (1 + Pres_gap)",
    Pres_from_ratio - (1 + Pres_gap),
)
Pres_from_success_ratio = sp.simplify(Wsuff_res / Wsuff_match)
expect_zero(
    "P_res from success-band ratio agrees with failure-band ratio",
    Pres_from_ratio - Pres_from_success_ratio,
)
```

Both new assertions exercise the *consistency* of the resonance penalty across the two band edges, which is a non-trivial consolidation claim.

In the Mathematica script, mirror these changes with **independent algebraic structure** (see F3). Specifically, insert after line 55 (the threshold-printing block):

```mathematica
(* Stage-49 matched-window derivation via parameterized effective gap. *)
Clear[DeltaEff];
$Assumptions = $Assumptions && DeltaEff > 0;
WMatchGen[de_] := PeReq/de;
expectZero[
  "matched fail edge from W_match(Delta_inf)",
  WMatchGen[Deltainf] - WfailMatch
];
expectZero[
  "matched success edge from W_match(Delta_0)",
  WMatchGen[Delta0] - WsuffMatch
];
expectPositive[
  "W_match decreasing in Delta_eff",
  -D[WMatchGen[DeltaEff], DeltaEff] * DeltaEff^2 / PeReq
];
```

And after line 63 (between the existing `expectPositive` blocks):

```mathematica
(* Stage-51 / Stage-68 resonance penalty via band-edge ratio extraction. *)
PresFromRatio = FullSimplify[WfailRes/WfailMatch, Assumptions -> $Assumptions];
expectZero[
  "P_res from band-edge ratio matches (1 + Pres_gap)",
  PresFromRatio - (1 + PresGap)
];
PresFromSuccessRatio = FullSimplify[WsuffRes/WsuffMatch, Assumptions -> $Assumptions];
expectZero[
  "P_res from success-band ratio agrees with failure-band ratio",
  PresFromRatio - PresFromSuccessRatio
];
```

Do not delete the existing assertions in lines 74-90 and 99-110 (sympy) / 57-94 (wl); leave them as redundant checks. The point is to *add* substantive consolidation assertions so the suite is not exclusively tautological.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 069` and `redteam exec-mathematica 069` and confirm the new lines `matched fail edge from W_match(Delta_inf) = 0`, `matched success edge from W_match(Delta_0) = 0`, `W_match decreasing in Delta_eff > 0`, `P_res from band-edge ratio matches (1 + Pres_gap) = 0`, and `P_res from success-band ratio agrees with failure-band ratio = 0` appear in both `.txt` outputs AND both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl`
- summary: Added matched-window generator checks and resonance band-edge ratio checks to both audit scripts.
- deviation: Mathematica checks were placed after the F3-ordered resonance checks so the first six Mathematica assertions remain in the required F3 order.

## F2 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py:18-24` (docstring provenance block)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py:54-66`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl:35-50`

**Issue:** The docstring claims that `[Pe_req/Delta_inf, Pe_req/Delta_0]` is the "Stage 066 theorem surface carried forward verbatim" and `P_res = 1/C_res^2` is the "Stage 068 resonance penalty," but the script introduces all of these (`Pe_req`, `Delta_0`, `Delta_gap`, `Pres_gap`) as free symbols on lines 54-55 with no link to the upstream stages. The provenance comment is unsupported by the script body.

**Required change:**

Step 1: In the SymPy script, after line 23 (before the closing `"""` on line 24), append a precise reference to the source-script lines:

```
- Source for the matched-window form `Pe_req/Delta_inf`:
  `scripts/moving_throat_pde_stage066_*.py` produces this form from the
  matched-branch theorem at its file's `expect_zero("Stage 049 matched
  window", ...)` block. This script imports the form symbolically as
  `Pe_req/Deltainf` and `Pe_req/Delta_0`; it does not re-derive it.
- Source for the resonance penalty `P_res = 1 + Pres_gap`:
  `scripts/moving_throat_pde_stage068_*.py` produces this form from the
  resonance amplification factor. This script imports the form symbolically
  as `Pres = 1 + Pres_gap` and asserts the band-edge ratio matches; it does
  not re-derive the amplification factor.
- The assertions in this script are conditional on Stages 066 and 068
  being correct. The new ratio-based assertions added per F1 verify the
  consolidation step but do not re-verify the upstream derivations.
```

Step 2: In the SymPy script, after line 66 (the `Wsuff_res` definition), add a comment block stating the upstream-anchor assumption explicitly:

```python
# --- Upstream-anchor assumption (see docstring) ---
# Wfail_match, Wsuff_match, Wfail_res, Wsuff_res are constructed here from
# free positive symbols. Their interpretation as the Stage 066 / Stage 068
# threshold values is asserted by the F1 consolidation block below
# (matched fail edge from W_match(Delta_inf) = 0, P_res from band-edge
# ratio matches (1 + Pres_gap) = 0), but the upstream derivations
# themselves are not re-checked here.
```

Step 3: In the Mathematica script, mirror Step 2 — after line 50 (the `WsuffRes` definition), add the analogous Mathematica comment block:

```mathematica
(* --- Upstream-anchor assumption ---
   WfailMatch, WsuffMatch, WfailRes, WsuffRes are constructed here from
   free positive symbols. Their interpretation as Stage 066 / Stage 068
   threshold values is asserted by the F1 consolidation block below;
   the upstream derivations themselves are not re-checked here. *)
```

The Mathematica script's header has no docstring; the comment block above the definitions is sufficient.

**Verification command:**
After Codex applies, the verifier will inspect the docstring of the SymPy script and the comment blocks in both scripts, and confirm both new comment blocks are present. Both scripts must still exit 0 after the new F1 assertions land.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl`
- summary: Added upstream provenance notes to the SymPy docstring and upstream-anchor assumption comments to both scripts.
- deviation: none

## F3 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl:35-50` (definition block)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl:57-73` (first assertion block)

**Issue:** The Mathematica script transliterates the SymPy script line-by-line, sharing identical variable names, the same definition order, and the same assertion sequence. Because both scripts make the same definitions (`Cres2 := 1/Pres`, `Wfail_res := Pres * Wfail_match`), the Mathematica engine cannot catch errors built into those definitions.

**Required change:**

Step 1: Reverse the resonance-penalty parameterization in the Mathematica script. Replace the `Cres2 = FullSimplify[1/Pres, ...]` definition (currently line 43) with the following derivation block. Edit lines 35-50 of the Mathematica script as follows:

**Before** (lines 35-43):
```mathematica
Clear[PeReq, Delta0, DeltaGap, PresGap, Deltainf, Pres, Cres2, uFail, uSucc, vFail, vSucc];
$Assumptions =
  Element[{PeReq, Delta0, DeltaGap, PresGap, vFail, vSucc}, Reals] &&
  PeReq > 0 && Delta0 > 0 && DeltaGap > 0 && PresGap > 0 &&
  vFail > 0 && vSucc > 0;

Deltainf = FullSimplify[Delta0 + DeltaGap, Assumptions -> $Assumptions];
Pres = FullSimplify[1 + PresGap, Assumptions -> $Assumptions];
Cres2 = FullSimplify[1/Pres, Assumptions -> $Assumptions];
```

**After:**
```mathematica
Clear[PeReq, Delta0, DeltaGap, PresGap, Cres2Prim, Deltainf, Pres, Cres2, uFail, uSucc, vFail, vSucc];
$Assumptions =
  Element[{PeReq, Delta0, DeltaGap, PresGap, Cres2Prim, vFail, vSucc}, Reals] &&
  PeReq > 0 && Delta0 > 0 && DeltaGap > 0 && PresGap > 0 &&
  Cres2Prim > 0 && Cres2Prim < 1 &&
  vFail > 0 && vSucc > 0;

Deltainf = FullSimplify[Delta0 + DeltaGap, Assumptions -> $Assumptions];
(* Parameterize via Cres2 as primitive; derive Pres = 1/Cres2 and verify it
   equals (1 + PresGap) by Solve.  This routes the resonance-penalty
   identity through a different algebraic operation than the SymPy script. *)
Cres2 = Cres2Prim;
Pres = FullSimplify[1/Cres2, Assumptions -> $Assumptions];
PresFromSolve = Pres /. First[Solve[Pres == 1 + presGapFree, presGapFree]];
(* PresFromSolve evaluates to 1/Cres2; verify it equals 1 + PresGap when
   PresGap = 1/Cres2 - 1.  This anchors the (Pres, PresGap, Cres2) triple
   via Solve rather than direct algebraic substitution. *)
expectZero[
  "Pres-PresGap consistency via Solve",
  (Pres /. presGapFree -> (1/Cres2 - 1)) - (1 + (1/Cres2 - 1))
];
```

Step 2: Rename `WfailRes` to `WfailResViaCres2` and define it via division rather than multiplication, so the Mathematica resonance threshold goes through `1/Cres2` directly rather than `Pres`:

**Before** (lines 49-50):
```mathematica
WfailRes = FullSimplify[Pres WfailMatch, Assumptions -> $Assumptions];
WsuffRes = FullSimplify[Pres WsuffMatch, Assumptions -> $Assumptions];
```

**After:**
```mathematica
(* Resonance thresholds via direct 1/Cres2 scaling rather than Pres
   multiplication.  This is algebraically equivalent but routes the
   threshold construction through a different intermediate quantity. *)
WfailResViaCres2 = FullSimplify[WfailMatch / Cres2, Assumptions -> $Assumptions];
WsuffResViaCres2 = FullSimplify[WsuffMatch / Cres2, Assumptions -> $Assumptions];
(* Compatibility aliases so downstream assertions still reference the
   `Res` names but the values were derived via Cres2 rather than Pres. *)
WfailRes = WfailResViaCres2;
WsuffRes = WsuffResViaCres2;
```

Step 3: Reorder the first assertion block (lines 57-73) to interleave the resonance and matched checks differently from the SymPy order. After the change, the first six `expectZero` / `expectPositive` calls in the wl file should be (in order):

```
1.  expectPositive["1 - C_res^2", 1 - Cres2]                      (was wl:72)
2.  expectPositive["P_res - 1", Pres - 1]                          (was wl:73)
3.  expectZero["P_res - 1/C_res^2", Pres - 1/Cres2]                (was wl:57)
4.  expectPositive["Delta_inf - Delta_0", Deltainf - Delta0]       (was wl:62)
5.  expectZero["resonance fail threshold - Pe_req/(C_res^2 Delta_inf)",
              WfailRes - PeReq/(Cres2 Deltainf)]                   (was wl:64)
6.  expectZero["resonance success threshold - Pe_req/(C_res^2 Delta_0)",
              WsuffRes - PeReq/(Cres2 Delta0)]                     (was wl:68)
```

Then continue with the remaining assertions (matched window width, matched success - matched fail) in their existing positions but reordered to appear after the resonance block. Concretely, edit lines 57-73 of the original to produce the above order, preserving all existing `expectZero` / `expectPositive` calls (just permuted).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 069` and confirm:
- The script's first `expectZero` line in the output is now `1 - C_res^2 > 0 -> ...` (not `P_res - 1/C_res^2 = 0` as before).
- A new line `Pres-PresGap consistency via Solve = 0` appears after the threshold definitions.
- The variable `WfailResViaCres2` (or `Cres2Prim`) appears in the script body.
- The script exits 0.
- Final assertions still produce zero residuals matching the SymPy output's final residuals.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl`
- summary: Reparameterized the Mathematica resonance path through primitive `Cres2`, built resonance thresholds via `1/Cres2`, and reordered the first assertion block.
- deviation: The `Pres-PresGap consistency via Solve` residual is defined in the derivation block but asserted after the first six reordered checks to preserve the required initial assertion order.
