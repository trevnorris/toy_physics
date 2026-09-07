# S11c-b carrier ablation-harness — ORCHESTRATOR adjudication key (⛔ NOT handed to the blind builder)

**Role.** This is the orchestrator-side companion to the blind builder directive
`directives/S11c_b_carrier_ablation_harness_directive.md`. It holds the case rationale, the per-knife
**impact cones** (must-MOVE / must-NOT-move), the **expected bites**, and the **self-test expected outcomes** —
i.e. the *interpretation* the harness must NOT contain (it prints; the orchestrator adjudicates). ⛔ The blind
builder (`gpt-6-astra`) gets ONLY the directive. This key goes to **me** and to the **harness-review legs**
(fresh Claude + Grok), who need the cones to verify the harness. It is the durable record of the 2026-09-07
review-round-2 legs (codex-sol + Grok, both **NOT-SOUND** on the round-1 draft; raw logs
`scratchpad/s11cb_reauthor/{sol,grok}.log`, and sol's grounded stdout is quoted below).

## What the round-1 legs confirmed CORRECT (do not re-litigate)
The **central re-author is right**: the carrier is **face-law-sourced**; the §3a energy basis and the
virtual-constraint fold are **pressure-independent**, so the old K1 (energy spurion) / K2 (constraint fold)
knives correctly stay dropped. Grounded (sol `/tmp/s11cb_review_sympy.py` literal stdout, `LAB_HELD × EULERIAN ×
RHO4_CONSTANT`):
```
ENERGY_PRESSURE_DEPS      (False, False, False, False)
CONSTRAINT_PRESSURE_DEPS  (False, False, False, False)
UNFOLDED_MASS_PRESSURE_DEPS (False, False, False, False)   # bulk mass balance pressure-free before the closure fold
CLOSURE_PRESSURE_DEPS     ((True,False,True,False),(False,True,False,True))   # θ enters via closure  → K_A
VIRTUAL_WORK_PRESSURE_DEPS (True, True, True, True)         # U/e_W enter via virtual work → K_T
```
WL check confirmed `ENERGY_PRESSURE_FREE=True`, `CONSTRAINT_PRESSURE_FREE=True`.

## Case pin — MATERIAL_ADVECTED × EULERIAN × RHO4_CONSTANT (a deliberate change from the gate record's LAB_HELD)
The prior gate record pinned `LAB_HELD`. **Change it to `MATERIAL_ADVECTED`** (route stays EULERIAN, so C_E is
still the certified object; density stays RHO4_CONSTANT). Reason, grounded: in `LAB_HELD × EULERIAN` the **U
carrier is identically zero in BOTH engines** — the LAB_HELD virtual normal displacement has no virtual-U
dependence (`WL:1040-1052`), so `faceGeneralizedRows` puts no pressure into U. sol's stdout:
```
RETAINED_ROW U0/U1/U2 CARRIER_NZ (False,False,False,False)   # U pressure-free in LAB_HELD
RETAINED_ROW E_W      CARRIER_NZ (True, True, True, True)
RETAINED_ROW THETA    CARRIER_NZ (True, True, True, True)
```
So under LAB_HELD, K_T's U channel has nothing to bite. In `MATERIAL_ADVECTED` the virtual map carries
`virtualUVector` (`WL:952`, `:1046-1048`), so U becomes a **live** carrier in both engines (sol: "material-case
probes also confirmed nonzero U carriers"). MATERIAL_ADVECTED therefore exercises **all three** channels (θ, U,
e_W); it strictly dominates LAB_HELD for coverage, and the construction code is shared across branches, so
certifying it certifies the code path that also serves LAB_HELD. ⚠ Harness-review legs: challenge this if c2 in
fact binds the LAB_HELD carrier specifically — I did not find that it does.

## Per-knife cones (CARRIER level — ∂(row)/∂(pressure atom)|_{→0}; NOT the full row)
The harness prints the **complete** carrier (every row × every slot) for baseline + each knife + diff. I read the
printed triples against these:

- **K_A — Λ_A response → θ.** must-MOVE: the **θ** carrier. must-NOT-move: the **U** and **e_W** carriers.
  ⛔ Do NOT list "Λ_V channel" or "the bulk base" as must-NOT — at the carrier level they are pressure-free and
  trivially zero (differentiating a pressure-free object is 0 regardless of damage — sol F13), so they carry no
  information. The Λ_A-only FORM (structurally delete the expanded `Lambda_A_0`-bearing addends, keeping Λ_V) is
  what makes this cone honest.
- **K_T — traction/virtual-work → U/e_W.** must-MOVE: the **U** carrier AND the **e_W** carrier (both live under
  MATERIAL_ADVECTED). must-NOT-move: the **θ** carrier. The FORM is the **whole traction channel** (SymPy: drop
  all four `face_u`/`face_e` additions; WL: `virtualWork = 0`) — this closes the Λ_X coverage hole (dropping only
  bare-`p` leaves `lambdaXResponse affinity`, which sol confirmed is nonzero:
  `EW_BARE_PRESSURE_SUBCHANNEL_NONZERO={True,True}`, `EW_LAMBDAX_PRESSURE_SUBCHANNEL_NONZERO={True,True}`).
- **K_W — face COLLAPSE (not exchange).** The collapse applies the slot-map `C_plus ← C_plus + C_minus`,
  `C_minus ← 0` to each pressure-bearing row's carrier columns. Adjudicate the printed collapse diff (nonzero)
  against this **column-support** transform — ⛔ NOT a full carrier-matrix rank change (round-2 legs computed the
  base and collapsed symbolic column ranks are **both 1**; the rank change is of the **minus-column block** only,
  1→0). Grounded by BOTH legs (independent):
  - **minus-face carrier columns → 0** (the identification kills the separate minus face).
  - **plus VALUE column moves** — the ±-value coefficients are **symmetric** (`p_minus − p_plus ≡ 0`), so collapse
    **doubles** the plus value column.
  - **plus ∂_w column → 0** — the ±-`d_w` coefficients are **antisymmetric** (`d_w_minus + d_w_plus ≡ 0`), so
    collapse **cancels** the plus `d_w` (⛔ it does NOT "absorb both"). sol: `COLLAPSED_DW_PLUS_ALL_ZERO=True`;
    grok: `EW_COLLAPSE_CARRIER_NZ (True,False,False,False)`.
  - **exchange is the weak alternative** the collapse replaces — dead in WL (equal ± carriers,
    `MATERIAL_SWAP_COLUMN_DIFF_NZ={False,False}`), and in SymPy it moves only the `d_w` columns
    (`(False,True,False,True)`); the collapse bites on all slots.
  ⛔ Delete the vacuous "pressure-free rows' carriers" must-NOT (they are identically 0 under any damage — no
  information). WL has no ∂_w pressure slot, so its `d_w` columns are N/A, not zero-by-cancellation.

## Disjointness (holds — keep)
K_A and K_T patch **disjoint** downstream consumers and **neither** patches the shared `affinity` (WL `:1079`) /
shared source: SymPy closure modifies θ at `:2834`, face additions modify only U/e_W at `:2998`; WL flux
(`:1080`) and traction (`:1081`) are separate fields. A knife that moves the *other* channel's carrier is a
non-disjointness finding.

## Self-test expected outcomes (the harness RUNS + PRINTS these; expectations live HERE, not in the directive)
- **Extractor-order self-test** — apply `→0` (P=0) **before** `∂/∂p` instead of after. Expected: **every** carrier
  (baseline and every knife's corrupted) becomes 0, so every diff is 0. If any knife still shows a nonzero diff
  under this order, the "carrier" is reading a stored value, not differentiating the live object.
- **Dead-path self-test** (one per engine, pinned pressure-free site: SymPy e_W kinetic addend `:2970`; WL kinetic
  addend `:1347`). Expected: the carrier diff is **identically 0** (proves the harness reports no spurious motion
  from a pressure-free structural change).
- **Live-rescale contrast** (×2 at each knife's own site: K_A `Lambda_A_0` addends; K_T `face_u`/`face_e` /
  `virtualWork`; K_W minus pressure slot). Expected: a **nonzero raw diff** that is NOT a FORM bite (arithmetic,
  not physics). ⚠ This only discriminates because K_A's FORM is now a structural addend **deletion**, not a
  `Lambda_A_0→0` coefficient zeroing — the ×2 rescale and the deletion are then genuinely different operations.

## Round-2 findings folded (both legs, 2026-09-07) — traceability
1. WL entrypoint OOM: `S11CB_PRIMARIES_ONLY` (`:2206`) tests AFTER `extractCouplingData` (`:2199`) +
   `kernelOriginsFromOrigins` (`:2203`) already ran (verified by me). ⇒ definitions-only load + direct
   `evaluatedModel`; never the emit `Do`. **2.** K_T must be one whole-channel FORM on both engines. **3.** K_A
   SymPy retargeted to Λ_A-only (structural deletion of the expanded `Lambda_A_0` addends — round-3 item 2
   sharpened this from a coefficient zeroing) so it no longer kills Λ_V. **4.** K_W exchange → collapse. **5.**
   pinned branch LAB_HELD → MATERIAL_ADVECTED (U live). **6.** SymPy `build_operator` returns `casify(...)` nested
   Tuple (`:3261`/`:603-613`, verified) ⇒ `named_tuple_row`, not `["…"]`. **7.** SymPy θ obs citation
   `:2367/:2372/:2377` are energy templates (overwritten at `:3042`) ⇒ delete. **8.** WL obs `:1257/:1259` are
   `frozenEvaluatedModel`, `:1189-1190` `rawModel` ⇒ pin `evaluatedModel :1344-1351`, corrupted=False; ∂ the
   APPLIED `pressureUpper[…]`, hit the EXPRESSION slot; `truncateBackground` is not P→0. **9.** SymPy drift guard:
   compare retained single-case object from imported `build_operator` vs an unablated temp copy (no single-case
   tag exists). **10.** self-tests pinned (above). **11.** must-MOVE/must-NOT + zero/nonzero interpretation moved
   here, out of the builder packet; harness prints the complete carrier with NO selection. **12.** removed the
   builder-facing G2 pointer. **13.** dropped the vacuous "bulk base" must-NOT.
Non-blocking: N6 atom **order** is `(δp+, ∂_wδp+, δp−, ∂_wδp−)` (`N6:321-322`) — pin that order in the manifest.

## Round-3 findings folded (both legs, 2026-09-07; the whole design PASSED — only precision items remained)
Round-2 legs PASSED: WL entrypoint/OOM, MATERIAL_ADVECTED pin (U live in MAT / zero in LAB, both computed), K_T
whole-channel, SymPy access & citations, drift guard, extractor-order self-test, complete-carrier, disjointness,
slot order. Round-3 folds:
1. **Directive blindness** — removed the cone-map from the target-object section, the "changes the carrier"
   disclosure, the K_W "exchange-is-dead / rank-changing" rationale, and neutralized the K_A/K_T headings that
   named the target row; reframed the coefficient-exception clause to implementation-fidelity ("if a site can't be
   patched as specified, stop") so it no longer discloses an expected bite.
2. **K_A SymPy FORM** — `Lambda_A_0→0` reads as a coefficient zeroing (the c→0 limit of a rescale) and would make
   the FORM-vs-rescale self-test indistinguishable; changed to a **structural addend deletion** (`sp.Add.make_args`
   filter dropping the Λ_A-bearing terms at `:2815-2828`), leaving Λ_V addends. sol computed the addends are
   cleanly separable (per face: 4 Λ_A, 4 Λ_V, 0 shared). ⛔ Not the `:408` registration bind.
3. **K_W cone** — corrected above (`C_plus ← C_plus+C_minus`, `C_minus ← 0`; minus columns → 0; plus-value doubles;
   plus-`d_w` **cancels**; ⛔ no full-matrix rank change).
4. **Self-tests** — pinned the dead-path to one pressure-free site per engine (`:2970` / `:1347`) and the
   live-rescale to ×2 at each knife's own site.
5. **K_W one function** — pinned to `substrate_substitutions` (`:1995`, single consumer `filtered_substrate`
   `:2032-2048`), not two.
