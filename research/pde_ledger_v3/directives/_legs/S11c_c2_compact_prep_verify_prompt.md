# Verify — the S11c-c2 compact-prep state (docs + the load-bearing export-membership decision)

The orchestrator is about to compact context. Before that, verify the state docs are accurate and — most importantly
— that one **load-bearing decision** about the c2 export is correct. ⛔ Do not rubber-stamp; if something is wrong,
say so plainly. This is a document + light-computation verification; ⛔ do not modify the working tree.

## Background (what just happened)
- STEP 1 (c1 record corrections) is DONE + committed `efae0337` (review-until-clear, both legs CLEAR; c1 STANDS).
- The c2 SymPy build directive `directives/S11c_c2_sympy_build_directive.md` was GATED (2 decision legs) + committed
  `c1de32b0`, then an astra build ran it (19 min, EXIT 0). The build's physics looks faithful (its builder report
  `_measurements/S11c_c2_sympy_builder_report.md`), but it (a) OVERSTEPPED — ran its own invalid self-reviews +
  ~40 files (quarantined out of the repo) because the directive §7 embedded the orchestrator's review-launch; and
  (b) OVER-EXPORTED — `scripts/S11c_c2_exports.py` is 60 MB from 3 fully-expanded closed objects. The directive was
  then fixed (`1ae6c336`): §7 bounds the builder to build→verify→report; §5 splits EMIT (→`.out`) from EXPORT.

## 1. THE LOAD-BEARING DECISION — the c2 export membership (verify with the real files)
The orchestrator decided the c2 export delta should contain **ONLY the "coupled nonlocal self-energy operator"** that
S11c-d binds — concretely `S11CC2_CLOSED_COUPLING_KERNEL` (the re-extracted closed off-diagonal coupling) — and that:
- `S11CC2_SELF_ENERGY_INCREMENT` is a **comparator representation** (spec §3c) the T7 comparator reads from stdout,
  ⛔ NOT a downstream binding ⇒ EMIT-only (→ `.out`), NOT in the ledger export.
- `S11CC2_CLOSED_SLAB_OPERATOR` (the full closed operator) is likely NOT a downstream binding either (default: drop
  it from the export; S11c-d binds the coupling operator, not the full slab operator).
- term_origins, parity blocks, the six §3d re-adjudication objects, and every §5 control operand/residual are
  EMIT-only.

**Verify:** Is this correct? Ground it in `directives/S11c_decisions.md` (esp. the sub-step table row for S11c-c —
its "exports" column — and the S11c-d row: what does profile-conditioned spectrum/scattering actually consume?), the
c2 spec `directives/S11c_c2_SHARED_PHYSICS.md` (§0 scope, §2 the three named objects, §3b/§3c, §4, §7), and the
bind-closure design `directives/export_ledger_bind_closure_design.md` (D1 membership = what a LATER step binds).
Specifically settle:
- (a) Does S11c-d bind the **closed coupling kernel**, the **closed slab operator**, or **both**? Is defaulting to
  export only the coupling kernel right, or does d need the full closed slab operator to solve its spectrum?
- (b) Is the **self-energy increment** genuinely a comparator-only artifact with no downstream binding (so it belongs
  in `.out`, not the export)? Or does something downstream bind it?
- (c) Given D1 says "export a row iff a LATER step binds it" but S11c-d is not built yet — is the right membership the
  step's **named deliverable** (the coupled self-energy operator), with the rest EMIT-only? Or should c2 export more
  defensively?
Report the correct export membership set (which S11CC2_* objects belong in the delta), with your evidence.

## 2. The directive fixes (are §5 and §7 correct?)
Read the two fixed sections of `directives/S11c_c2_sympy_build_directive.md` (§5 EMIT≠EXPORT + held/factored
representation; §7 builder-bounds-and-no-orchestrator-process). Do they correctly encode the fix, and is the
"held/factored, not sp.expand'd; the comparator canonicalizes A−B" representation instruction sound (i.e. will a
held/factored export still be consumable by S11c-d and the comparator)?

## 3. The state docs (accuracy)
Check the new top clause of `STATUS.md` (2026-09-05, "STEP 1 DONE; c2 build directive gated; astra build RAN —
OVERSTEPPED + OVER-EXPORTED") against reality: `git log --oneline -8`, the committed directive, and the salvaged
artifacts. Any inaccuracy or overstatement? In particular: is "physics faithful per its builder report" appropriately
hedged (the builder report is the builder's own account, not yet independently verified)?

## 4. The plan (sound?)
The plan is: salvage the artifacts, run REAL independent legs (fresh Claude agent + Grok) on the PHYSICS, and
separately REPAIR the export (membership + held/factored). Is running the physics legs on the current script valid
given the export needs a repair — i.e. are the physics (the emitted .out objects + the audit computations) separable
from the export membership/representation, or is the held/factored representation change entangled with the physics
computation enough that a clean rebuild would be better? Give a recommendation.

## Output
For each of 1–4: your finding + evidence (commands/file:line for the real files). End with: **CLEAR TO COMPACT** or
the exact list to fix first, and your one-line recommendation on §4 (physics-legs-then-export-repair vs rebuild).
