# S11c-b carrier ablation-harness knife-list — G2 decision gate: NOT-SOUND → re-author owed (2026-09-07)

**Artifact:** `directives/S11c_b_carrier_ablation_harness_directive.md` (orchestrator-written knife-list; UNCOMMITTED).
**Legs (G1 orchestrator-written → Codex + Grok; G2 one-pass decision gate):** codex-sol (gpt-5.6-sol xhigh) + Grok
(grok-4.6 high), parallel, identical prompt. Raw logs: `scratchpad/s11cb_gate/{sol,grok}.log` (session-scratch —
the exact line-level site corrections live there). **Both verdicts: NOT-SOUND / remove-and-retarget.**

## The load-bearing finding (both legs, independent; orchestrator-VERIFIED G4)
My knife-list targeted the **wrong object.** The pressure-slot carrier c2 binds — `∂(slab rows)/∂(δp±,∂_wδp±)|_{P=0}`
— is **face-law-sourced**, NOT energy-basis or constraint-sourced. Verified in N6's own extractor
`scripts/S11c_c2_N6_diagnostic_sympy.py`: `pressure_coefficients` (`:409-411`) = `∂(rows)/∂p .xreplace(p→0)`, and
`rows` (`:311`) = `{U, E_W, THETA: mass−correction}` built from `face_generalized_force_rows` — the traction /
virtual-work / #90-closure face laws. `N6_BASE_SOURCE_DEPENDENCY` (`:325`) measures `energy.has(p)` (energy is
P-independent). ⇒ **K1 (energy-spurion freeze) and K2 (constraint fold) move P-independent objects that are NOT in
the carrier** — freezing them leaves the extracted carrier byte-identical, and a builder would green-stamp the
wrong object (a 26=26-class energy freeze sold as a bite on `C_E`).

## Re-author spec (convergent across both legs — for the post-compact re-author; ⛔ do NOT build the current K1–K4)
- **DROP K1 (energy spurion) + K2 (constraint) from THIS carrier harness.** Energy-basis genuineness is already
  covered by the #89/#89a form ablation; pin-B by its own form ablation. (Retain K1/K2, corrected, only in a
  separate energy-basis/whole-slab harness if ever wanted — not here.)
- **KEEP K3 (θ / Λ_A closure fold, #90), corrected sites:** PY the fold is `closure_residuals` +
  `mass_balance_source − closure_residual_sum` in `build_operator` (`S11c_b_brane_operator_sympy_audit.py:2815-2834`)
  → writes `THETA_BALANCE` (`:3042-3046`). WL: `flux = lambdaAResponse·affinity + lambdaVResponse·normalVelocity`
  in `faceSources` (`mathematica/S11c_b_brane_operator_mathematica_audit.wl:1080`) → `projectedFaceFlux` (`:1094-1096`)
  → `MASS_EVOLUTION_ROW` (`:1349-1350`). FORM = change the flux STRUCTURE (drop the affinity channel / swap flux↔a
  non-response operand) — ⛔ patch the response at `:1080`, NOT the shared `affinity` at `:1079` (shared with Λ_X).
  ⛔ NOT sites `:389` (import-key), `:2214` (provenance), `:408-409` (knob binds = coefficient theater), `:2135`
  (traction/Λ_X = wrong channel), `:115-116` (dimension metadata). CONE: only **Λ_A** is pressure-dependent (Λ_V is
  absent from the pressure derivative); MUST-MOVE = extracted θ-row carrier; MUST-NOT = extracted U/e_W carrier, Λ_X,
  bulk ACCUMULATION+ADVECTIVE, energy-EL.
- **RETARGET K4 (slot index wiring) to native pressure atoms + face-swap only:** ⛔ the cited PY `delta_p_{face}`
  (`:484-493`) are symbol REGISTRATION (`bind_additional_inherited`) — patching them rewrites nothing (theater); the
  slots live inside the imported S11c-a expressions. WL has NO discrete δp slots (continuum `pressureField[±1]`
  `:1014-1015`). ⛔ The ∂_w-slot↔value-slot swap is dimensionally invalid (pressure vs pressure/length) AND would
  require re-deriving PY's linearization (forbidden). FORM = a FACE swap only: PY `delta_p_plus↔delta_p_minus` (with
  the matching `d_w_*` pair) substituted in the CONSUMED expressions (after `face_generalized_force_rows` /
  `closure`, `:2149-2220`/`:2815-2828`); WL `pressureUpper↔pressureLower` (`:1014-1015`) on the `evaluatedModel` path.
  Pin ONE, name it, ⛔ no builder "or".
- **ADD K_trac (traction / virtual-work — the U/e_W channel, 4 of 5 carrier rows; a coverage HOLE):** PY
  `face_generalized_force_rows` (`:2135-2175`) + the add into `U_BODY_BALANCE`/`E_W_BALANCE` (`:2998-3040`); WL
  `tractionPressure`/`VIRTUAL_WORK` (`:1081-1083`) + `faceGeneralizedRows` (`:1139-1146`) into `U_MOMENTUM_ROWS`/
  `THICKNESS_ROW` (`:1345-1348`). FORM = drop the bare-`p` traction term, or skip adding the face rows to U/e_W.
  MUST-MOVE = extracted U/e_W carrier; MUST-NOT = θ carrier, Λ_X-only-vs-bare-p (per leg detail), bulk. (K3+K_trac =
  the two disjoint CONTENT channels; K4 = the index structure.)

## Architecture fixes (both legs)
- **Observation object = `∂(emitted SLAB_OPERATOR tags)/∂(native pressure atoms)`** — `PY_S11CB_SLAB_OPERATOR`
  (`:4142`, rows `U_BODY_BALANCE`/`THETA_BALANCE`/`E_W_BALANCE`) and `WL_S11CB_SLAB_OPERATOR` (`:2281-2282`,
  `U_MOMENTUM_ROWS`/`MASS_EVOLUTION_ROW`/`THICKNESS_ROW`). ⛔ Do NOT reimplement N6's `face_factory`/`carrier[rows]`;
  ⛔ do NOT patch emit. Pin the exact source tag, row order, pressure-slot order, `∂/∂p`, `P=0` substitution.
- **Operator-only entrypoint (avoid the ≥64 GB OOM):** call `build_operator`/`evaluatedModel` for ONE case; skip
  `COUPLING_KERNEL`/tower/heavy controls. Still wrap-live-engine (call the production function), not a re-derivation.
- **⚠ IMPLEMENTABILITY BLOCKER:** both engines enumerate all 4 cases (PY `:4117-4120`, WL `:2190-2197`);
  `S11CB_PRIMARIES_ONLY` skips CONTROLS, not cases. "Canonical engine unchanged" + "single-case" are incompatible
  until a **production-supported case selector** is added. Pin the case **EULERIAN × LAB_HELD × RHO4_CONSTANT** (C_E
  is Eulerian; WL's `corrupted` flag is Eulerian-only) — ⛔ not an "e.g.".
- **Drift guard:** direct canonical run vs an unmodified temp copy under the SAME selector, scoped to the
  operator-only tags — ⛔ not "the committed tags" (stale + forces full emit).
- **Self-ablation SPECIFIED, not just named:** a dead-path mutation AND a carrier-extractor mutation must drive the
  printed must-MOVE `diff` identically 0 (a live coefficient rescale yields a nonzero diff, so nonzero ≠ FORM);
  print both transcripts, interpretation to review.
- **Delete the value-leaks + discretion:** ⛔ remove `26=26`, `15+15`, `40→10`, "10 uniform-sector terms" (expected
  counts of a different object — the compact-verify flagged the same); pin ONE function + ONE FORM per knife (no
  "or"/"e.g."); fix "…are G2-cleared" → "Once G2-cleared…".
- Path prefix: engine/deliverable paths are repo-root-relative (`research/pde_ledger_v3/…`).

## Disposition
The G2 gate did its job (caught a wrong-object knife-list before any build). ⛔ Do NOT build the current K1–K4.
NEXT = re-author the knife-list per this spec (physics-bearing targets → review-until-clear), then build. The
current directive stays UNCOMMITTED.
