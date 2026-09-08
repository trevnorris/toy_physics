I'll execute the Wolfram FORM ablation myself: copy the harness to `/tmp`, read the artifacts, then run each mode one kernel at a time.Harness, directive, and astra transcript are in hand. Next I’ll copy the files to `/tmp`, pin the engine construction sites, and start the sequential kernel runs.Copies are in `/tmp`. I'll make the three knives no-ops, confirm the engine sites, then start BASELINE.Engine sites are unique and the no-op harness is in place. Starting BASELINE (one kernel, `timeout 600`).BASELINE finished in ~2 min with the same operator digest as astra. Running NO-OP `K_A` next.NO-OP `K_A` left the operator digest unchanged. Running NO-OP `K_T`.NO-OP `K_T` also matched baseline. Running NO-OP `K_W`.All three no-op knives matched baseline. Running live `RESCALE_K_A`.`RESCALE_K_A` doubled only the affinity term. Running `DEAD_PATH`.All six kernel runs finished. Writing a comparison script to import the WXF records and print slot-wise diffs.Comparison finished. Extracting the slot-wise diffs and remaining static checks.# WL carrier-ablation harness — FORM execution review (Grok)

Artifact: `research/pde_ledger_v3/mathematica/S11c_b_carrier_ablation_harness.wl` (SHA-256 `9cdce88952726e8e297b102720c37f60190f6ce20c622d39b764b152873c6470`, matches astra’s committed digest). Engine wrapped, not re-derived: `S11c_b_brane_operator_mathematica_audit.wl`. All kernel runs used `timeout 600`, one seat at a time, copies under `/tmp/gwl` only. Comparison instrument: `/tmp/gwl/compare_carriers.wl`.

---

## Static (read, no execution)

**Carrier is the live emitted operator.** Worker loads definitions only (prefix above the `(* Main variable-coefficient objects. *)` marker), then:

```183:186:research/pde_ledger_v3/mathematica/S11c_b_carrier_ablation_harness.wl
  operator = Global`evaluatedModel["EULERIAN", "MATERIAL_ADVECTED",
    "RHO4_CONSTANT"]["OPERATOR"];
  If[!AssociationQ[operator], stop["operator evaluation"]];
  carrier = extractCarrier[operator];
```

`extractCarrier` (L121–132) does `D[row, atom] /. (atoms → 0)` on `scalarRows` of that operator (`U_MOMENTUM_ROWS`, `MASS_EVOLUTION_ROW`, `THICKNESS_ROW`). No reimplementation of `faceSources` / `pressureField`, no hand-typed carrier. `modelRecord` lives after the marker (`engine:1572`) and is not loaded.

**Knives = directive’s one site × one FORM, exact-string patch.** `patches` (L66–80) + `StringReplace` of `spec[[2]] → spec[[3]]` (L165–166):

| Knife | Function | Before (engine) | After |
|---|---|---|---|
| K_A | `faceSources` | `flux = lambdaAResponse affinity + lambdaVResponse normalVelocity;` (`engine:1080`) | drop affinity addend |
| K_T | `faceSources` | `virtualWork = virtualWork tractionPressure virtualNormalDisplacement;` (`engine:1083`) | `virtualWork = 0;` |
| K_W | `pressureField[-1]` | `pressureField[-1] := pressureLower[…];` (`engine:1015`) | `:= pressureUpper[…]` |

`validateSites` (L91–112) requires `StringCount[source, spec[[2]]] === 1` **and** the same unique count inside the named-function window. Engine counts (mechanical): `fluxBefore`, `workBefore`, `pressureBefore`, `faceSources[…]`, `projectedFaceFlux[faceAssociation_Association]`, `evaluatedModel[…]`, `activateSpatialDivergences[expression_Association]`, marker — each **1**. Tower flux at `engine:2862` is a different string (`velocity` not `normalVelocity`) and sits after the marker, so it is not loaded.

**PRINT-not-PASS.** No `PASS`/`FAIL`/verdict/`bite` in executable code (only the L14–17 prohibition comment). `emitTriple` prints `{baseline, corrupted, diff}` then guards **shape** and `$Aborted|$Failed|Indeterminate` (L144–151). The only `SameQ` stop is the canonical-copy **drift** guard (L212–213), not a knife-bite assertion.

---

## Execution ablation

NO-OP copy `/tmp/gwl/harness_noop.wl` replaced the third element of `patches["K_A"|"K_T"|"K_W"]` with `fluxBefore` / `workBefore` / `pressureBefore` (so `StringReplace` is identity). `validateSites` still sees the original before-strings.

Worker CONSTRUCTION_PATCH on the NO-OP runs (baseline_source **equals** corrupted_source):

```
CONSTRUCTION_PATCH: <|"MODE" -> "K_A", ..., "baseline_source" -> "  flux = lambdaAResponse affinity + lambdaVResponse normalVelocity;", "corrupted_source" -> "  flux = lambdaAResponse affinity + lambdaVResponse normalVelocity;"|>
CONSTRUCTION_PATCH: <|"MODE" -> "K_T", ..., "baseline_source" -> "  virtualWork = virtualWork tractionPressure virtualNormalDisplacement;", "corrupted_source" -> "  virtualWork = virtualWork tractionPressure virtualNormalDisplacement;"|>
CONSTRUCTION_PATCH: <|"MODE" -> "K_W", ..., "baseline_source" -> "pressureField[-1] := pressureLower[xOne, xTwo, xThree, time];", "corrupted_source" -> "pressureField[-1] := pressureLower[xOne, xTwo, xThree, time];"|>
```

BASELINE `OPERATOR_SHA256` = `8cba37be5fe30437e6bf09c69a6e51767a635d809d4f803995cebd1808f4334a` (matches astra’s transcript). NO-OP K_A/K_T/K_W share that digest. Live `RESCALE_K_A` = `4519e9b9…` (matches astra). Live `DEAD_PATH` = `6f101c16…` (operator **did** change).

Comparison script `/tmp/gwl/compare_carriers.wl`, command `timeout 600 wolframscript -file /tmp/gwl/compare_carriers.wl`, exit 0. Literal slot-wise `Expand[mode − baseline]` (in-plane `d_w_*` slots printed `NA`):

### NO-OP K_A / K_T / K_W — diffs identically 0

```
K_A_MINUS_BASELINE: <|"U_MOMENTUM_ROWS[1]" -> <|"delta_p_plus" -> 0, ..., "delta_p_minus" -> 0, ...|>, ... all five rows 0 ...|>
K_A_DIFF_IDENTICALLY_ZERO: True
K_T_MINUS_BASELINE: <|"U_MOMENTUM_ROWS[1]" -> <|"delta_p_plus" -> 0, ..., "delta_p_minus" -> 0, ...|>, ... all five rows 0 ...|>
K_T_DIFF_IDENTICALLY_ZERO: True
K_W_MINUS_BASELINE: <|"U_MOMENTUM_ROWS[1]" -> <|"delta_p_plus" -> 0, ..., "delta_p_minus" -> 0, ...|>, ... all five rows 0 ...|>
K_W_DIFF_IDENTICALLY_ZERO: True
```

A no-op patch cannot fabricate a bite. **Pass.**

### RESCALE_K_A — coefficient of the same live term, not a structural change

```
RESCALE_K_A_MINUS_BASELINE: <|"U_MOMENTUM_ROWS[1|2|3]" -> 0, "THICKNESS_ROW" -> 0,
  "MASS_EVOLUTION_ROW" -> <|"delta_p_plus" -> -(lambdaAZero/(rhoM*(1 - I*frequency*tauA))),
                            "delta_p_minus" -> -(lambdaAZero/(rhoM*(1 - I*frequency*tauA)))|>|>
RESCALE_K_A_DIFF_MINUS_AFFINITY_CARRIER: <|"..." -> 0  (all rows/slots)|>
RESCALE_K_A_DIFF_PLUS_ASTRA_K_A_BITE: <|"..." -> 0  (all rows/slots)|>
```

The extra term is exactly the baseline `MASS_EVOLUTION_ROW` affinity carrier (`−λ_A0 / (ρ_M (1 − i ω τ_A))` on both faces). Astra’s live K_A bite (committed transcript L123) is that same term with coefficient **−1** (`+λ_A0 / …`); rescale is **+1 ×** the affinity carrier. No vanished/new coupling on momentum or thickness. **Pass.** (`RESCALE_K_T` / `RESCALE_K_W` not run; not required.)

### DEAD_PATH — pressure-free deletion does not move the carrier

Operator digest **changed** (`8cba37be…` → `6f101c16…`; kinetic `THICKNESS_ROW` addend actually removed), but:

```
DEAD_PATH_MINUS_BASELINE: <|"..." -> 0  (all rows/slots)|>
DEAD_PATH_DIFF_IDENTICALLY_ZERO: True
```

Extractor responds only to pressure-bearing structure. **Pass.**

### EXTRACTOR_ORDER — harness order is the physical one

Default extractor (L127): `D[row, atom] /. zeroRules` (∂ then → 0). `zeroFirst` (L126): `D[row /. zeroRules, atom]` (→ 0 then ∂). Literal `CARRIER − ZERO_FIRST_CARRIER`:

```
EXTRACTOR_ORDER_CARRIER_MINUS_ZERO_FIRST: <|"U_MOMENTUM_ROWS[1]" -> <|"delta_p_plus" -> -1/2*(sigmaW*w1JetOne) + (lambdaXZero*sigmaW*w1JetOne)/(2*rhoM*(1 - I*frequency*tauX)), "delta_p_minus" -> (same)|>, … momentum 2/3 analogous …,
  "MASS_EVOLUTION_ROW" -> <|"delta_p_plus" -> -(lambdaAZero/(rhoM*(1 - I*frequency*tauA))), "delta_p_minus" -> (same)|>,
  "THICKNESS_ROW" -> <|"delta_p_plus" -> WZero/2 - (lambdaXZero*WZero)/(2*rhoM*(1 - I*frequency*tauX)), "delta_p_minus" -> (same)|>|>
ZERO_FIRST_CARRIER: all pressure slots 0
EXTRACTOR_ORDER_IDENTICALLY_ZERO: False
```

`∂/∂p |_{p→0}` is differentiate first, then evaluate at 0. Zeroing first kills the linear response (all 0). The harness’s chosen order is the physically correct one; the self-test prints the disagreement rather than silently swapping it. **Pass.**

---

**SOUND**
