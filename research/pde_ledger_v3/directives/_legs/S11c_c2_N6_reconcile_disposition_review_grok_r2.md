I'll start with the source-of-truth files (tally, engines, question, consults) before reading the disposition, so the verdict is grounded in those rather than the artifact's framing.Tally and framing docs are in. Next I’ll read the engine sources, consults, and collapse reviews at the cited lines before looking at the disposition.The sympy grade-extraction citation looks off. I’ll locate the actual `grades()` / `xreplace(profiles)` sites and read the consults plus the v3 Grok review before the disposition.I’ll verify the EL commutator in sympy, then inspect how `FROZEN_RELATIONS` keys actually match before opening the disposition.Independent picture of the comparison is in place. I’ll read the disposition last and score it against the tally, engines, and EL calc.# Independent review — S11c-c2 N6 reconcile DISPOSITION

Reading order: tally → engines/comparator → question / PATH_B / RESOLVED → consults / v3 reviews → **then** `_measurements/S11c_c2_N6_reconcile_disposition.md`.

## Independent picture (before the artifact)

**Established.** Within each engine, `R_cov` and `C_E−C_M` give no nonzero at the retained rectangle (PIT-qualified). Cross-engine, every *matched* key of those already-zero objects is 0: that is `(0)−(0)`, a dual-engine confirmation of the vanishing statement, not operand agreement. `SOURCE_CONTROL_DELTA` is the same kind of triviality (`ACTUAL≡BASELINE` at shipped `κ_a=1, κ_j=0`). Structural support agrees on 400 matched keys, and only when both sides witnessed nonzero support.

**Not established.** Operand AGREE on `C_E`/`C_M` (40 leftover), `SOURCE_*` (76), Φ (18). Schema-unmatched families (`R_N6`, channels, RC sources, guards) were never compared. Map 4 is a correct *ungraded density* table and the wrong object to apply to already-EL’d, already-graded source leaves. Path B (stop, carry a debt) is the structurally indicated disposition.

---

## 1. Numerical fidelity

Committed tally `_measurements/S11c_c2_N6_comparator_run_tally.txt` vs artifact §1, §2, §8:

| Family | Tally | Artifact |
|---|---|---|
| `N6COV_R_COV` / `_BASELINE` / `_CONTROL_DELTA` | 160/0 each (ll. 8–10) | §1 ll. 37–38 |
| `N6COV_SOURCE_CONTROL_DELTA` | 160/0 (l. 14) | §1 l. 38 |
| `N6RC_CARRIER_BRIDGE_RESIDUAL` | 320/0 (l. 17) | §1 l. 39 |
| `N6RC_ADVECTION_ABSENCE` | 6/0 (l. 16) | §1 l. 40 |
| `N6RC_FROZEN_RELATIONS` | 8/0 (l. 24) | §1 ll. 40–43 (split 4+4; tally ZERO = residual-0, which includes both) |
| support total | 400/0, 1680 UNDEC (l. 53) | §1 l. 45 |
| `N6RC_CARRIER_EULERIAN`/`MATERIAL` | 280/40 each (ll. 19–20) | §2 ll. 62–63 |
| `N6COV_SOURCE_ACTUAL`/`BASELINE`/`PREDICTED` | 16/76 each (ll. 12, 13, 15) | §2 ll. 65–66 |
| `N6COV_FROZEN_PHI` | 42/18 (l. 6) | §2 l. 67 |
| `N6RC_SOURCE_{EULERIAN,MATERIAL,BRIDGE_RESIDUAL}` | 0/0/640 each (ll. 29–31) | §2 ll. 70–71 |
| `N6COV_ACTUAL_CONTROL_PARAMETERS` | 0/4 (l. 5) | §8 l. 159 |
| `N6COV_PHI_DOMAIN_CENSUS` | 0/4 + 44 BOOL (l. 7) | §8 ll. 159–160 |

No mismatch.

---

## 2. The `(0)−(0)` claim

Within-engine zeros, before the comparator subtracts:

- SymPy `R_cov = actual − predicted` (`covariance_sympy.py:194–204`); WL `covDelta = mapCombine[sMMap, pMap, sub]` (`.wl:878,883`).
- SymPy carrier bridge `residual(ce, cm)` (`reconcile_sympy.py:210,254`); WL `bridge = mapCombine[eMap, mMap, sub]` (`.wl:874,876`).
- `SOURCE_CONTROL_DELTA`: shipped `ACTUAL_A_RHO=1`, `ACTUAL_JUNK=0` (`covariance_sympy.py:34–35,140–145`); WL `actualAdvectionCoefficient=1`, `actualJunkCoefficient=0` (`.wl:18–20,845–848`).

Question doc §2 (ll. 70–75) is the right reading: each of those objects is already 0 inside each engine, so the cross-engine 0 is trivial and says nothing about `C_E`/`C_M` / ACTUAL / PREDICTED / Φ.

The artifact does **not** lump every matched residual-zero into `(0)−(0)`. §1 splits:

- **(a)** the within-engine-zero families above, plus `ADVECTION_ABSENCE` 6/0 and four zero-valued `FROZEN_RELATIONS` leaves. The 6 are the joined `density_gradient` components on the two `RHO4_CONSTANT` cases (`reconcile_sympy.py:241–244`; `.wl:971–974`); PY `a_rho` vs WL `ADVECTION` do not join by `field_name`. Those gradients are the zero object. Correctly `(0)−(0)`.
- **(b)** the other four `FROZEN_RELATIONS` matched leaves, plus support 400/0.

The 8 `FROZEN_RELATIONS` residual-zeros are exactly 4 cases × `{a_rho, h_alpha}` (the only top-level `field_name` joins; `extract_meta` at `comparator.py:687–715`). Constructors:

```52:52:research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py
    h = b.dot(b.u, b.grad_W) / b.W_bg if alpha == 'LAB_HELD' else sp.S.Zero
```

```841:842:research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl
  a = displacement.Table[td[density4, i]/density4, {i, 3}];
  h = If[anchor === "LAB_HELD", displacement.Table[td[WBg, i]/WBg, {i, 3}], 0];
```

So `h_α ≠ 0` on the two `LAB_HELD` cases and `a_ρ ≠ 0` on the two `RHOBR_CONSTANT` cases (`density4 = rhoBr/WBg` at `.wl:839`). Those four residual-zeros are genuine nonzero-operand agreement, and the artifact says so. Support 400/0 is also not `(0)−(0)`: `support_residual` returns UNDECIDED if either side is `NO_NONZERO_FOUND` (`comparator.py:537–540`), so a ZERO support residual is two `NONZERO_WITNESSED` observations agreeing.

§2 keeps the 40/76/18 as UNADJUDICATED and does not call any surfaced-nonzero family agreed. The 280/16/42 residual-zeros on those families are reported as the split, not as family-wide AGREE.

---

## 3. Structural obstruction — EL calc (reason 3a)

1-D EL in θ, no naked θ: `EL_θ L = −d_x(∂L/∂θ')`.

`L = a θ' e' + b e W' θ'`, with `T(a)=k R`, `T(b)=−k R/W`, `R=W_0/W` (`W=W_bg`).

**After EL, then T** (coefficients constant through the derivative):

`EL_θ L = −[a e'' + b(e' W' + e W'')]`

`T(EL_θ L) = k W_0 (−W e'' + e W'' + W' e') / W²`

**T first, then EL** (`T(a)`, `T(b)` are x-dependent):

`∂(TL)/∂θ' = k (W_0/W) e' − k (W_0/W²) e W'`

`EL_θ(TL) = k W_0 [(e W'' + 2 W' e') W − W² e'' − 2 e (W')²] / W³`

**Commutator** (sympy, expanded):

\[
EL_θ(TL) − T(EL_θ L) = \frac{k W_0}{W^2} W' e' − \frac{2 k W_0}{W^3} e (W')^2
\]

This is astra’s identity (`strategy_consult_astra.md:17–21`; PATH_B ll. 18–21; artifact §3 ll. 95–97).

Retained rectangle: undifferentiated `WBg → W_0(1+η w1)`; a spatial `WBg` jet → `σ_W · (profile jet)` (`.wl:127–132`; `profile_definitions` at `brane_operator_sympy_audit.py:900–905`; `grades()` at `diagnostic.py:245–248`). So `W' = O(σ_W^1)` with no extra η on the jet. The first term is `η^0 σ_W^1` and **survives**. The second is `O(σ_W^2)` and is dropped by `retainTerm` / `GRADES`. `γ_14→0` only isolates κ; it is not a retained-order truncation.

μ really is EL of the energy density in both engines: WL `el` (`.wl:144`) then `muE = el[energy["DENSITY"]]` (`.wl:843`) and `MU -> el[pulled]` (`.wl:247`). Map 4 sends constant WL coefficients to `R_W = W_0/W_bg` (position-dependent). Applying that table to already-differentiated sources misses the product-rule term above. No grading-stage rewrite restores a derivative that was never taken.

**Reason (b).** Both engines expand profiles *before* grade extraction — `.wl:127–141` (`profileRules` + `finish`) and `diagnostic.py:245–248` (`grades()` `xreplace(self.inputs.profiles)`), **not** `reconcile_sympy.py:245–248` (that span is object emission). The disposition cites the diagnostic site correctly (§3 l. 103). Graded leaves carry `W_0`, not live `W_bg`. `R_W · W_bg = W_0` is a product of two η-series; coefficientwise reproduction is the convolution `(TF)_1 = T_0 F_1 + T_1 F_0`, which the no-grade-mixing tripwire forbids. That is Grok’s v3 MUST (`collapse_v3_review_grok.md:34–52`) and the G4 note in `collapse_directive_gate.md:121–128`.

The “upstream of EL; a replay does not reconcile the already-emitted `.out`” line is **not** stronger than the two reasons. The artifact explicitly says the sound instrument is not unique (§3 ll. 108–116): cleaner path = thickness/basis map before EL; acknowledged alternative = a derivative-aware chain-rule bridge (`strategy_consult_astra.md:48`); neither reaches the frozen streams without new construction; a canonicalized replay validates the replay unless the old-stream relation is separately established (`strategy_consult_astra.md:50`). That matches the consults. The two reasons are scoped to channel-(b) source operands, not over-applied to the carrier 40.

---

## 4. Over-claim / under-claim

Honored:

- `(0)−(0)` is not upgraded to operand AGREE (§1 ll. 47–53).
- The 19-row ungraded table is not upgraded to “the sources agree” (§4 ll. 125–128). Opus v3 independently re-derived the table as a constructor identity (`collapse_v3_review_opus.md:10–22`); the artifact keeps that distinct from operand collapse.
- UNADJUDICATED is not upgraded to “just thickness”; leftover SHAPE is recorded as uninspected (§5 ll. 132–137).
- B is not justified as “c2 already has everything”; S11c-d mixing/leakage is cited (`S11c_decisions.md:52` table + `:83–91` N5: off-diagonal coupling `O(η)`, leakage `O(η²)`) and B is an explicit debt (§6 ll. 141–144; `RESOLVED.md:36,55`).
- Not called “weak N6”; headline is dual-engine covariance + components not in a common thickness coordinate + operand debt (§1 ll. 55–57).

The four `FROZEN_RELATIONS` nonzero-premise equalities are called a “small positive cross-engine signal” and immediately bounded as premise/support, not carrier/source/Φ AGREE. That is not an upgrade.

---

## 5. Completeness + carry-open

§7 carries: the operand debt; uninspected SHAPE; the three `RESOLVED.md:36–45` caveats (Φ physical-correctness; `V_E≡V_M` as builder agreement, prediction uses `V_E` not `Φ(V_E)`; extracted-block leakage); 2 S11c-b signs / 6 §3d / c1 ENERGY, with an explicit do-not-pre-adjudicate. §8 surfaces the 8 census leaves as a fact-lookup (6-row crosswalk, WL-only `THICKNESS`/`MATERIAL_NORMAL`, `n6cov_J_mu`↔`junkMu`). Nothing the disposition needs to say is missing. (`I_{M→E}` is a step-record item, PATH_B execute-2, not this record.)

---

## 6. Scope

Reading B is restated as standing (`d21c8ff5`), not re-opened. The compared WL stream is identified as `ae73b884` from certified engine `e11f2f82` (file sha `e5cea55b`), superseding first-clearance `48a0b4e7` — that matches `STATUS.md` and the harness `production_sha256`. c1 is left standing; c1 ENERGY stays UNDECIDED. No carried-open item is pre-adjudicated. Collapse-directive rebuild / v4 is forbidden (provenance ll. 180–181).

---

**SOUND** — nothing outstanding changes what the disposition may claim.
