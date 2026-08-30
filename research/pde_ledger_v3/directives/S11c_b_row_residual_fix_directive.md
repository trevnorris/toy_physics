# Fix directive (v3) — S11c-b row-residual instrument: coupling in-scope weak verdict + deterministic witness

Two build legs (fresh Claude agent + Grok) found the instrument sound except for the **coupling** verdict; two
decision-list legs (Codex + Grok) then rejected two earlier drafts of this fix and **converged on the exact
construction below**. Implement that construction. Leave the truncation function, assembly accounting,
ε-metadata, strong slab rows, mass row, and admissibility componentwise comparison untouched (ablation-verified
sound). Records: `directives/_measurements/S11c_b_row_residual_instrument_build_review.md`.

## Background — why the current coupling path is wrong (do not repeat it)
The current build classifies the **post-`_bridge_d`** residual with `apply_bridge_d=False`
(`scripts/S11c_b_row_residual.py:304, 406-409`). `_bridge_d` substitutes the profile relations (introducing
`η,σ_W` and relating the jets), and it does **not** commute with `formal_divergence`/`formal_dx` (layer fixture
`scripts/S11c_b_adjudicated_comparison.py:1304-1308`). Running the Euler operator on the bridged residual treats
related jets as independent and manufactures false bulk (e.g. `div(W_bg u_T)` is a pure divergence pre-bridge but
post-bridge/`apply_bridge_d=False` yields a spurious Euler `−W_0 η ∂w1 + σ_W ∂w1`). Therefore: **the Euler /
homotopy operators must run on the PRE-bridge residual (independent jets); the requested truncation must run on
the BRIDGED result; and the in-scope weak representative must be the bridged-and-truncated homotopy REMAINDER,
never a bridged Euler signature** (the signature is a field-keyed obstruction, is early-exited/incomplete
`…adjudicated_comparison.py:664`, and is absent on the `RESIDUAL_BULK` short-circuit `:735-739`).

## The construction to implement (per coupling case: both cross-sector blocks + both relabelled adjoints)

All of `A._arithmetic_residual` (L593), `A.classify_total_divergence` (L726), `A._bridge_d` (L230),
`A.formal_divergence` (L291), `A._homotopy_vector` (L686), `A._normalise_exact` (L617),
`A.PRODUCTION_FIELD_REGISTRY` (L205) are in the committed layer `A = S11c_b_adjudicated_comparison`; reuse them,
do **not** modify the layer. Let `T = requested_truncation` (the instrument's existing L171 function).

1. `R_pre := A._arithmetic_residual(left_pre_bridge_d, right_pre_bridge_d)` — the **pre-bridge** cross-engine
   A−B residual (the same operands the layer's weak route uses at L985-992: the aligned coupling operands
   **before** `_bridge_d`).
2. **Full-order diagnostic** (emit, label clearly as full-order, NOT the in-scope verdict):
   `full := A.classify_total_divergence(R_pre, A.PRODUCTION_FIELD_REGISTRY, apply_bridge_d=True)`. Emit
   `FULL_PREBRIDGE_ROUTE := full.route` and `EULER_SIGNATURE := full.euler_signature` (deterministic-serialized,
   see §Determinism). This is the layer's **reviewed** convention for the full-order classification.
3. **`ROW_RESIDUAL` (the A−B measurement, §382):** `ROW_RESIDUAL := T(A._bridge_d(R_pre))` — the requested-
   truncated, bridged, in-scope cross-engine **A−B density**. It stays a genuine `operand_A − operand_B` object
   (emit the truncated bridged `left`/`right` operands too, so A, B, A−B are all present).
4. **`IN_SCOPE_WEAK_REMAINDER` (the modulo-total-in-plane-divergence representative):** compute the homotopy
   remainder **in the instrument**, on the **pre-bridge** residual, so it is available **even when the layer
   short-circuits** (`RESIDUAL_BULK`) or its own homotopy is not invoked:
   ```
   try:
       V   = A._homotopy_vector(R_pre, A.PRODUCTION_FIELD_REGISTRY)
       rem = A._normalise_exact(R_pre - A.formal_divergence(V, A.PRODUCTION_FIELD_REGISTRY))
       IN_SCOPE_WEAK_REMAINDER := T(A._bridge_d(rem))
   except Exception:
       emit NO_CLEAN_QUOTIENT flag; IN_SCOPE_WEAK_REMAINDER unavailable (fall back to ROW_RESIDUAL only)
   ```
   `rem = R_pre − div(V)` is the part of the residual that is **not** a total in-plane divergence; bridged and
   truncated it is the in-scope genuine-bulk-modulo-divergence representative. ⛔ Do **not** substitute the Euler
   signature for it. ⛔ Do **not** truncate the pre-bridge operand. ⛔ Do **not** run Euler/homotopy on the
   bridged residual.

Apply steps 1-4 to both cross-sector blocks and both relabelled adjoint operands (coverage unchanged).

⛔ The strong slab rows, the mass row, and the admissibility comparison keep their exact residuals and must
**not** be routed through `classify_total_divergence`/homotopy (§1d: their first-jet content is physics, not a
representational divergence). ⛔ Do **not** modify the committed layer.

## Determinism
Canonicalize the emitted `EULER_SIGNATURE`/witness **presentation only** (the layer's `sp.cancel`-based
`_normalise_exact` is hash-order-sensitive; the route is stable). Serialize deterministically (sort terms /
pinned `signsimp`) on an **emission copy**, **after** route/residual/remainder are computed. ⛔ Never apply the
sign/serialization rewrite to `ROW_RESIDUAL`, `IN_SCOPE_WEAK_REMAINDER`, or any potential `V` (flipping a
remainder's sign changes the emitted residual). Emit a determinism demonstration (a stable per-witness hash or a
runnable two-invocation byte-identity check).

## Method, guards, obligations (rule 2, rule 5)
- Compute → emit → then guard. ⛔ Assert **no** route, **no** residual value, **no** in-scope verdict, **no**
  which-engine judgement. This directive states **no** expected route and **no** expected residual (rule 5): for
  each coupling case, `FULL_PREBRIDGE_ROUTE`, `ROW_RESIDUAL`, and `IN_SCOPE_WEAK_REMAINDER` are the measurement,
  adjudicated off-instrument.
- Keep every non-coupling emission **byte-identical** to the prior build (a reviewer will diff the strong-row and
  admissibility output and require it unchanged).
- Deliverable: edited `research/pde_ledger_v3/scripts/S11c_b_row_residual.py`, runnable as before, exit 0.
- Report to `directives/_measurements/S11c_b_row_residual_fix_directive.md`, each claim carrying its command:
  the file:line of (a) the pre-bridge `classify_total_divergence(..., apply_bridge_d=True)` full-order call, (b)
  the `ROW_RESIDUAL := T(bridge(R_pre))` A−B density, (c) the instrument-side homotopy remainder
  `T(bridge(R_pre − div V))` with the throw fallback, (d) the witness determinism canonicalization; a two-run
  byte-identical-witness demonstration; and a literal diff summary showing the strong-row and admissibility
  output is byte-identical to the pre-fix build.
