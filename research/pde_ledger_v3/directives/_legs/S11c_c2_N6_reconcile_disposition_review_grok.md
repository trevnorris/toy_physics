I'll start from the source-of-truth files (tally, engines, question/consult docs) and form my own view before reading the disposition.Tally and framing docs are in; next I’ll read both engines, the comparator, and the consult/collapse reviews before looking at the disposition.The EL counterexample matches PATH_B on paper; I’ll now check the consults, collapse reviews, and the cited engine lines for reason (b).The consults live under `directives/_legs/`; I’ll read those plus the v3 reviews and `S11c_decisions.md` before opening the disposition.Independent picture of the comparison is in place. I’ll read the disposition last and score it against that.`48a0b4e7` is cited for the Wolfram per-engine result but a later commit marks it stale — I’ll check which engine actually produced the `R_cov` vanishing.# Independent review — S11c-c2 N6 reconcile DISPOSITION

Reading order: tally → engines → question/PATH_B/RESOLVED → consults/v3 reviews → **then** the artifact.

## Independent picture (before the artifact)

**What the comparison established**

- Within each engine, `R_cov` and the carrier bridge `C_E−C_M` give no nonzero at the retained rectangle (PIT-qualified). Cross-engine, every *matched* key of those vanishing objects is 0: that is `(0)−(0)`, a dual-engine confirmation of the vanishing statement, not operand agreement.
- The operands themselves are the open object: `CARRIER_EULERIAN`/`MATERIAL` 40 nonzero, `SOURCE_ACTUAL`/`BASELINE`/`PREDICTED` 76, `FROZEN_PHI` 18. `SOURCE_BASELINE` is a shipped-settings duplicate of `SOURCE_ACTUAL`. Schema-unmatched families (`R_N6`, channels, RC sources, guards) were never compared.
- The collapse dictionary’s load-bearing map 4 is a correct *ungraded density* table and the wrong object to apply to already-EL’d, already-graded source leaves. Path B (stop, carry a debt) is the structurally indicated disposition, not a cost call.

**What it did not establish:** operand AGREE on `C_E`/`C_M`/sources/Φ; that the leftover shape is “just thickness”; that Φ is physically correct; that `V` transforms; that omitted blocks do not leak.

---

## 1. Numerical fidelity

Tally `_measurements/S11c_c2_N6_comparator_run_tally.txt` vs artifact §1, §2, §8:

| Family | Tally | Artifact |
|---|---|---|
| `N6COV_R_COV` / `_BASELINE` / `_CONTROL_DELTA` | 160/0 each (ll. 8–10) | §1 l. 33–34 |
| `N6COV_SOURCE_CONTROL_DELTA` | 160/0 (l. 14) | §1 l. 34 |
| `N6RC_CARRIER_BRIDGE_RESIDUAL` | 320/0 (l. 17) | §1 l. 35 |
| `N6RC_ADVECTION_ABSENCE` | 6/0 (l. 16) | §1 l. 35 |
| `N6RC_FROZEN_RELATIONS` | 8/0 (l. 24) | §1 l. 35 |
| support total | 400/0, 1680 UNDEC (l. 53) | §1 l. 35–36 |
| `N6RC_CARRIER_EULERIAN`/`MATERIAL` | 280/40 each (ll. 19–20) | §2 ll. 49–50 |
| `N6COV_SOURCE_ACTUAL`/`BASELINE`/`PREDICTED` | 16/76 each (ll. 12, 13, 15) | §2 ll. 52–53 |
| `N6COV_FROZEN_PHI` | 42/18 (l. 6) | §2 l. 54 |
| `N6COV_ACTUAL_CONTROL_PARAMETERS` | 0/4 (l. 5) | §8 l. 136 |
| `N6COV_PHI_DOMAIN_CENSUS` | 0/4 + 44 BOOL (l. 7) | §8 ll. 136–137 |

No mismatch.

---

## 2. The `(0)−(0)` claim

Within-engine zeros (already 0 before the comparator subtracts):

- SymPy `R_cov = actual − predicted` (`covariance_sympy.py:194–204`); WL `covDelta = mapCombine[sMMap, pMap, sub]` (`.wl:878`).
- SymPy carrier bridge `residual(ce, cm)` (`reconcile_sympy.py:210,254`); WL `bridge = mapCombine[eMap, mMap, sub]` (`.wl:874,876`).
- `SOURCE_CONTROL_DELTA`: shipped `kappa_a=1, kappa_j=0` ⇒ `ACTUAL≡BASELINE` (`covariance_sympy.py:34–35,140–145`; `.wl:18–20,845–848`).

Question doc §2 (ll. 70–75) is the right reading: each object already 0 inside each engine ⇒ cross-engine 0 is trivial and says nothing about `C_E`/`C_M`/ACTUAL/PREDICTED/Φ.

Artifact §1 ll. 37–40 names `R_cov` and `C_E−C_M` as the `(0)−(0)` objects and refuses operand AGREE. §2 keeps the 40/76/18 as UNADJUDICATED and does not call any surfaced-nonzero family agreed.

The 280 carrier zeros and 16 source zeros in §2 are genuine *operand* zeros on those keys (the carriers/sources are not within-engine zeros). The artifact reports the split and withholds family-wide AGREE because of the 40/76 leftover — that is the conservative family-level claim PATH_B asked for, not a lump into `(0)−(0)`.

`FROZEN_RELATIONS` 8/0 and `ADVECTION_ABSENCE` 6/0 sit in the §1 vanishing list (question §1.i did the same; PATH_B’s EARNED list omitted them). The `(0)−(0)` *explanation* still names only `R_cov` and `C_E−C_M`. Those 8/6 are metadata joins, not a second covariance channel, and do not move the N6 claim. Not a must-fix.

---

## 3. Structural obstruction — EL calc (reason 3a)

Standard 1-D EL in θ, no naked θ: `EL_θ L = −d_x(∂L/∂θ')`.

`L = a θ' e' + b e W' θ'`, with `T(a)=k R`, `T(b)=−k R/W`, `R=W_0/W` (`W=W_bg`).

**After EL, then T** (coefficients constant through the derivative):

`EL_θ L = −[a e'' + b(e' W' + e W'')]`

`T(EL_θ L) = −k (W_0/W) e'' + k (W_0/W²) e' W' + k (W_0/W²) e W''`

**T first, then EL** (`T(a)`, `T(b)` are x-dependent):

`∂(TL)/∂θ' = k (W_0/W) e' − k (W_0/W²) e W'`

Differentiating:

`EL_θ(TL) = 2 k W_0 (W'/W²) e' − (k W_0/W) e'' + (k W_0/W²) e W'' − 2 k W_0 e (W')²/W³`

**Commutator:**

\[
EL_θ(TL) − T(EL_θ L) = \frac{k W_0}{W^2} W' e' − \frac{2 k W_0}{W^3} e (W')^2
\]

which is astra’s identity (`strategy_consult_astra.md:17–21`; PATH_B ll. 18–21; artifact §3 ll. 79–81).

Retained rectangle, from the engines: undifferentiated `WBg → W_0(1+η w1)`; a spatial `WBg` jet → `σ_W · (profile jet)` (`.wl:127–132`; `diagnostic:245–248`). So `W' = O(σ_W^1)` with no extra η on the jet. The first term is `η^0 σ_W^1` and survives. The second is `O(σ_W^2)` and is dropped. γ_14→0 only isolates κ; it is not a retained-order truncation.

μ really is EL of the energy density in both engines: WL `el` (`.wl:144`) then `muE = el[energy["DENSITY"]]` (`.wl:843`) and `MU -> el[pulled]` (`.wl:247`); SymPy `variation` (`diagnostic:345–349`). Map 4 sends constant WL coefficients to `R_W = W_0/W_bg` (position-dependent). Applying that table to already-differentiated sources misses the product-rule term above. No grading-stage rewrite restores a derivative that was never taken.

Reason (b): both engines expand profiles *before* grade extraction — `.wl:127–141` (`profileRules` + `finish`) and `diagnostic:245–248` (`grades()` `xreplace(profiles)`), not `reconcile_sympy.py:245–248`. Graded leaves carry `W_0`, not live `W_bg`. `R_W · W_bg = W_0` is a product of two η-series; coefficientwise reproduction is the convolution `(TF)_1 = T_0 F_1 + T_1 F_0`, which the no-grade-mixing tripwire forbids. That is Grok’s v3 MUST (`collapse_v3_review_grok.md:34–52`) and the G4 note in `collapse_directive_gate.md:121–128`.

The “upstream of EL; a replay does not reconcile the already-emitted `.out`” line is not stronger than the two reasons. Astra (`strategy_consult_astra.md:48–50`): a canonicalized replay validates the replay, not the original streams. Grok: a sound collapse needs a new emit path in both constructors.

---

## 4. Over-claim / under-claim

Honored:

- `(0)−(0)` is not upgraded to operand AGREE (§1 ll. 37–40).
- The 19-row ungraded table is not upgraded to “the sources agree” (§4 ll. 102–105).
- UNADJUDICATED is not upgraded to “just thickness”; leftover SHAPE is recorded as uninspected (§5 ll. 107–114).
- B is not justified as “c2 already has everything”; S11c-d mixing/leakage is cited (`S11c_decisions.md:83–91`, N5: off-diagonal coupling `O(η)`, leakage `O(η²)`) and B is an explicit debt (§6 ll. 116–121).
- Not called “weak N6”; headline is dual-engine covariance + components not in a common thickness coordinate + operand debt (§1 ll. 42–44).

---

## 5. Completeness + carry-open

§7 carries: the operand debt; uninspected SHAPE; the three RESOLVED caveats (Φ physical-correctness, `V_E≡V_M` as builder agreement, extracted-block leakage); 2 S11c-b signs / 6 §3d / c1 ENERGY, with an explicit do-not-pre-adjudicate. §8 surfaces the 8 census leaves as a fact-lookup. Nothing the disposition needs to say is missing. (`I_{M→E}` is a step-record item, not this record.)

---

## 6. Scope

Reading B is restated as standing (`d21c8ff5`), not re-opened. c1 is left standing; c1 ENERGY stays UNDECIDED. No carried-open item is pre-adjudicated. Collapse-directive rebuild / v4 is forbidden (§ provenance ll. 157–158).

`48a0b4e7` in §1 l. 31 is the *first* WL clearance; the compared `.out` is `ae73b884` from certified engine `e11f2f82` (Grounding l. 20 is the right hash). The physics claim (WL reproduces `R_cov`, own primes, import-free) still holds of the compared stream — the 160/0 matched `R_cov` zeros would not be 0 if that `.out` had a live residual. Provenance nit, not a claim change.

---

**SOUND** — nothing outstanding changes what the disposition may claim.
