# S11c-c2 blind Wolfram N6 engine — BUILD REVIEW: both legs CLEAR (2026-09-07)

Artifact: `mathematica/S11c_c2_N6_mathematica_audit.wl` (astra-built per CLEAR-TO-BUILD directive `33c45297`; 1004 lines,
40 tags). Codex-written → **two build legs: a fresh Claude agent + Grok**, identical prompt
`directives/_legs/S11c_c2_N6_wl_build_review_prompt.md`, **serialized** (both ablate Mathematica; 2-seat licence).
Ablation artifacts under `/tmp/.../scratchpad/n6wl/` (agent) and `/tmp/n6_wl_review_grok/` (Grok) — ephemeral.

## Both legs: BUILD CLEAR (convergent, independent methods)
Each formed its own view of `C_E−C_M` and `R_cov` from §5c + route-2 + the sibling specs BEFORE reading the engine
(Grok saved an independent derivation `derive_n6_identities.py`), then line-traced every load-bearing object and ran the
five mandatory FORM ablations (each a /tmp copy, single case, one kernel, `timeout 600`; agent 7 runs 451–513 s, Grok
5 runs 2–4 min skip-PIT). Working tree untouched by both.

**Load-bearing objects — all reached by computation (no hand-typed payload), same line-traces from both legs:**
`μ_E = el[energy.DENSITY]` (825); `μ_M` = phiMap-on-energy + Jacobian×wave-deg-2 + EL (233–236, 828); `μ_E∘Φ = μ_E /.
predictedMap` (835) — three DISTINCT constructions (leaf counts 669 / 4179 / 3344). `C_E = carrier[eulerianSlabFace]`
(848), `C_M = carrier[materialFaceFold]` with covector inside (840, 264–294) — two independent routes. `es`/`ms` =
`sourceBind` stripping Z/resolvent at `pFace→0` (355–370). Residuals = un-cancelled `arithmeticPlus[a,neg[b]]` circuits
evaluated only at PIT sampling. No `assert`, no `VERDICT`, no residual-zero exit (all `Quit[]` operational).

**Baseline (computed findings, ⛔ NOT asserted):** `SPLIT_CHECK` = 0/288 (affine split exact, samplewise, not a
literal-0 node); `R_N6` = 18/288 and `SOURCE_BRIDGE_RESIDUAL` = 18/72 (genuine nonzero); `CARRIER_BRIDGE_RESIDUAL` = 0
and `R_COV` = 0 (carrier + source-naturality hold at these cases — a computed measurement, disposition deferred).

**The five ablations — all bite one-sided (both legs):**
1. **Carrier FORM knife** (contaminate the material covector `invT[[1,4]]`, hold `ms`): `CARRIER_BRIDGE_RESIDUAL`
   0→8 nonzero, `C_M`/`CARRIER_CHANNEL`/`MATERIAL_OPERAND`/`R_N6` move; **`C_E`, Eulerian operand, `es`, `ms`,
   `SOURCE_BRIDGE_RESIDUAL` byte-identical.** (Baseline `C_E` hash = `C_M` hash; knife splits them.)
2. **Φ-coefficient knife** (`a_ρ→2·a_ρ` on the ACTUAL path, prediction held): `R_COV` 0→18, `SOURCE_ACTUAL` moves;
   `SOURCE_PREDICTED`, `FROZEN_PHI` truth table, all carriers fixed. Bites where the `a_ρ+h_α` truth table is blind.
3. **θ-junk knife** (`κ_j·J_μ·e_W` into `μ_M` only, `MATERIAL_ADVECTED.RHO4` where `a_ρ=h_α=0`): `R_COV` 0→6,
   prediction fixed; `ADVECTION_ABSENCE` computed (not A−A).
4. **Non-circularity** (force `ms_pred` from the material pullback): `R_COV` collapses to structural zero while
   `SOURCE_ACTUAL` stays byte-identical across the Φ-knife and the circular version ⇒ the shipped `ms_pred` (μ_E+Φ) is
   genuinely independent; breaking the prediction route does not move the material route.
5. **PIT** (`k_in=k_out`): the shipped engine rejects `kOut===kIn` (line 763, `EXCLUDED_LOCI`); an off-diagonal
   `(k_out−k_in)` jet identity would be hidden if kept. **δ is WL-derived** — `D`,`E` per-object from the interned
   circuit `degreePair`/`nodeDegree` (WL primes `{1000000009, 998244353, 1004535809}` ≠ the SymPy set), real branch
   charts selected before modular reduction, `min(1,D/(N−E))` with an `Exists[goodPrime]` bad-prime condition.

**Operational/discipline (both):** import-free (no Get/Import/abs-path/.py/.out); byte-identical isolated-vs-in-repo
(blindness); flush observed; light RSS (~240–420 MB); Φ domain census rank-2 `UNCOVERED={}` (all μ_E jet atoms covered);
slot-linearity/closure guards computed (`G_lin`,`G_cross`,denominator); able-to-fail dimensions (extra-`W_0` →
`Inactive[Equal][2,1]`, flagged inconsistent, not coerced). **Three caveats correctly OPEN:** Φ is the *declared* map
(not derived from motion); prediction uses `V_E` not `Φ(V_E)` (`V_E≡V_M` = builder agreement); block-leakage unchecked.

## Adjudication (G4)
Both legs demonstrated the controls BITE (reproducible ablation scripts + stdout), independence by one-sided corruption,
and no hand-typed payloads — positive evidence, ⛔ not weak "found nothing." Convergent across two independent harnesses
+ Grok's independent derivation. No physics-wrong defect survived either filter. ⇒ **WL N6 engine per-engine BUILD
CLEAR.** ⚠ Scope: this clears the ENGINE (correct construction, biting controls); the computed results' cross-engine
disposition (carrier/`R_cov` matching, the T7 join) is the comparator/reconcile stage — ⛔ NOT claimed here, ⛔ neither
engine subtracts its own discrepancy and calls it covariance.

## ⚠ BLINDNESS — PRECISE POSITION (corrected 2026-09-07, compact-prep Codex-verify)
"Blindness byte-identical" above = **import-freedom + isolated-vs-in-repo** (astra imported nothing; both runs
byte-identical) — TRUE and unaffected. ⚠ BUT the compact-prep verify caught that the build directive was **not fully
outcome-value-free**: two residual soft-leaks of the expected-zero outcomes (`ΔC=0`≡`C_E=C_M`; `R_cov=0`) survived my
round-3 lint (scrubbed post-hoc; see the directive-review adjudication's POST-HOC CORRECTION), so **"outcome-blind
authorship" was weaker than claimed** — astra saw that the carrier + source-naturality residuals are expected to vanish.
⭐ **This did NOT compromise the engine, and the build legs are the direct proof:** the risk of leaking the expected
zeros is a *designed-to-agree / fix-until-matches* engine (forced zeros). A forced/hardcoded zero **cannot move under
corruption** — yet both legs showed the carrier knife moves `C_E−C_M` off zero and the Φ knife moves `R_cov` off zero
(one-sided, others fixed), the non-circularity ablation proves the routes independent, and both legs independently
DERIVED the physics and line-traced every object to a computation. ⇒ the residuals are **genuinely computed from the
physics**, not tuned to the leaked zeros. **No rebuild is required** — the ablations verify exactly what the (breached)
value-free-directive control was meant to protect; the leak was a directive/records defect (now fixed), not a
substantive blindness failure. (User to confirm this disposition; the alternative for maximal rigor is a fresh-builder
rebuild + repeat legs.)

## Next
Commit reviewed baseline → generate + commit the reviewed `.out` (git-annex/GIN) → c2 T7 comparator + reconcile.
