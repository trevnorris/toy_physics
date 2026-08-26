# S11c-a WL fix — the free-premise background bulk current (§3c/§1b) — directive (rev 2, folded from 2 legs)

## Why this exists
T7 cross-engine reconciliation found the two blind engines diverge in the T-f projection: the Wolfram (WL)
engine carries the rest-frame background bulk current as independent free-premise symbols
(`currentWBackground`, `currentXBackground{i}`), while the SymPy engine constructs it as `ρ_4D⁰·v_bulk⁰` and so
obtains zero.

The defect: WL's treatment violates §3c ("*none may be introduced as a free premise*"), the current law §1b
(`j = ρ_4D v_bulk`), and — internally — WL's OWN background velocity, which WL already sets to zero
(`bulkVelocityZero`, lines 115-119, substitutes every component including the normal one to `0`). The engine
computes its background bulk current as unconstrained free symbols rather than from `ρ_4D⁰·v_bulk⁰`.

⛔ Which trace terms survive, and whether any downstream object vanishes, is computed from the premises and is
NOT stated here — §3c: "Which trace terms then survive is computed from these premises, not stated here." Do
NOT seek or engineer any particular downstream survival or cancellation; the only task is to correct the
current's provenance and let every downstream object recompute.

This directive re-scopes WL to construct the background bulk current faithfully to §1b/§2d/§3c. §3c already
forbids the free-premise treatment; no spec change is needed — derive from the existing spec.

⚠ SUPPLIED / UNFALSIFIABLE WITHIN THIS BUILD: the §3c/§1b/§2d readings and the background-state premises are
verified physics, supplied as premises. ⛔ Do NOT tune the output to match any recorded value or the sibling
engine, and ⛔ do NOT infer a "correct answer" to iterate toward — compute and print; the cross-engine diff and
the ablation review happen on the orchestrator's side, where a residual is a finding, not a build failure
(rule 5).

## The PROPERTY that must hold (state what must be TRUE; reach everything by computation)
P1 — NO FREE-PREMISE BACKGROUND CURRENT. The rest-frame background bulk current used anywhere in the engine
   (the `Zero`/background part of `currentW` and each `currentX{i}`, and everything derived from them) must be
   the CONSTRUCTED object `j⁰ = ρ_4D⁰ · v_bulk⁰` per §1b (`j = ρ_4D v_bulk`), formed from the engine's supplied
   background density and its supplied rest-frame background bulk velocity — the same `rhoBulkZero` (line 434)
   and `bulkVelocityZero` (lines 115-119) the engine already carries. It must be REACHED BY COMPUTATION from
   that construction, never introduced as the independent free symbols `currentWBackground` /
   `currentXBackground{i}` (defined at lines 438-446). Where the supplied background velocity is zero, the
   background current is reached AS zero by that computation, not asserted. (§3c: "Every background face value
   or normal derivative … is obtained by differentiating a member of the supplied background state `𝔅⁰`; none
   may be introduced as a free premise.")
   ⭐ This is the identical property already enforced on the SymPy engine's background jets; it is now applied
   to WL's background current.

## Scope — the background-current construction and every consumer
Change site: the `currentWZero` / `currentXZero` definitions (lines 438-446). Because these feed
`traceFieldInventory` (lines 612-618) and the T-f projection assembly (`projectionTermsSource`, lines ~798-822),
correcting the construction at the definition propagates to every consumer.
Re-emit (recompute automatically when the engine re-runs, but they MUST be present in the output): every object
that consumes the background current — the projection family (`PROJECTION_SHAPE_DERIV`,
`PROJECTION_DYNAMIC_OPERAND`, `PROJECTION_STATIC_OPERAND`, `PROJECTION_RESIDUAL`, `PROJECTION_TERM_ORIGINS`),
the shifted-face-trace / face-shift objects that read the current entries of `traceFieldInventory`, and the
representation/form-control/control-independence/uniform-limit packages that consume any of the above. Every
branch, both faces, both DOFs, both density representatives, both routes, and the controls must be recomputed
for the objects in scope.
⛔ Do NOT touch: the perturbation current (`currentWWave` / `currentXWave` = `currentWPerturbation` /
`currentXPerturbation{i}`), the nonzero background density (`rhoBulkZero` = `rhoBulkBackground`, which is
correctly nonzero — §3c: the density background does not vanish), the background velocity construction
(`bulkVelocityZero`), or the pure face-geometry objects. Only the background CURRENT construction changes.

## Script obligations — non-negotiable (rule 2 / §6)
1. The script PRINTS computed objects; it never states a conclusion. Every emit payload is a CAS object.
2. For every COMPARISON it emits operand A, operand B, and `A−B` BEFORE any guard. A physics disagreement
   emits and CONTINUES with exit status 0; a nonzero exit is reserved for operational failure only (§6).
   ⛔ Do not add a guard that asserts a physics value — before OR after the emit — and exits nonzero on a
   physics residual.
3. Interpretation belongs to the step record, not the script.
⛔ The only place physical symbols may be combined by hand is in constructing the supplied law
`j⁰ = ρ_4D⁰·v_bulk⁰` from the supplied `𝔅⁰` members; every downstream object must be REACHED BY COMPUTATION,
never typed as a result. A control re-enters the chain at the field construction (an ansatz/map/level-set/law),
never at a result.

## Acceptance — VALUE-FREE and CONSTRUCTION-VERIFIED (rule 5); correctness is verified DOWNSTREAM, not here
- The engine runs to completion and re-emits every in-scope object AND its derived packages, each in its
  EXISTING emit schema (a comparison object emits A/B/residual; a non-comparison object — e.g. `FACE_SHIFT`'s
  exact trace + shape derivative, or the projection operand and its already-separate residual — emits its
  computed object as it already does, over all branches/faces/DOFs/representatives/routes/controls).
  ⛔ Do NOT invent residual payloads for non-comparison objects or alter the established tag schema.
- P1 holds by CONSTRUCTION, verified at SOURCE — NOT by the emitted value. ⚠ Because the engine's background
  velocity is already zero, `j⁰ = ρ_4D⁰·v_bulk⁰` evaluates to the SAME `0` in the output as a bare `:= 0`
  would, so an output diff CANNOT distinguish a faithful construction from an assertion. P1 is therefore
  verified two ways, both of which a `:= 0`, a free-symbol left in place, or a post-hoc `/. freeSymbol -> 0`
  patch FAILS:
    (a) SOURCE-DIFF REVIEW — the `currentWZero`/`currentXZero` definitions must be BUILT from the engine's
        supplied background density and background velocity (`rhoBulkZero`, `bulkVelocityZero`) — the objects
        §1b's `j = ρ_4D v_bulk` names. The review reads the diff and confirms the construction.
    (b) BACKGROUND-VELOCITY PROBE (a FORM ablation the review runs on a copy): replace the background bulk
        velocity by a formal NONZERO probe; a faithful `j⁰ = ρ⁰·v⁰` construction makes the background current
        TRACK the probe (nonzero), while a `:= 0` assertion or a post-hoc substitution leaves it zero — so the
        probe distinguishes them. ⛔ The existing `W_bg`-jet form controls do NOT exercise the
        velocity→current dependency and cannot verify P1; do not rely on them for this.
  ⛔ Do NOT self-certify by grepping for the old symbol names.
- ⛔ Agreement with any recorded value or the sibling engine is NOT a build acceptance item. The orchestrator
  diffs cross-engine over all cases and the review legs ablate the deliverable; do not target any value.

## Target
`research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` (this file at HEAD; last
touched a7459cb8). The relevant spec is already in `S11c_a_SHARED_PHYSICS.md` §1/§1b/§2d/§3c — derive from it.
