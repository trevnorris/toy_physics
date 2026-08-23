# Independent physics review — the S11c-a SymPy engine

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`
(the engine; ~1900 lines). Its second product is
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_a_exports.py` (the frozen ledger it writes).

## What to check
This engine takes the **closed** S11c-a spec and computes the **first-order shape derivative of every
S11b interface law** once the uniform brane is un-frozen: `W₀ → W_bg(x)` (a non-uniform slab thickness)
and `μ_R → μ_R,bg(x)`. It emits, for **both anchoring branches** (`LAB_HELD` / `MATERIAL_ADVECTED`),
**both faces** `s∈{+1,−1}`, **both independent face DOFs** `{δW, ζ_c}`, and **both density representatives**,
the objects **T-0 … T-i**: background density map, outward face normal, conormal, true-area measure and its
shape derivative, face velocity, relative flux, kinematic balance, traction + virtual work, shifted-face
trace, dynamic-window projection (+ static operand, residual, term-origin inventory), virtual material
constraint, physical sourced mass balance, and the face-closure shape derivative. It also runs four control
families (§5a representation-invariance two-route + one-sided independence corruption; §5b per-direction
source-level form ablation; §5c uniform-limit regression; §6 homogeneity).

Your job: **derive these objects yourself from the spec and decide whether the engine computed them
correctly** — and, above all, whether any emitted payload is a **hand-typed CAS object with no data
dependence on the derivation** (a genuine SymPy expression, authored, that would not move if the derivation
above it were deleted). That is the one defect this review exists to catch, and only a **form ablation**
catches it. Everything else — a wrong sign, a wrong normal orientation, a dropped face-shift term, a
mis-bound bulk density, a broken two-route control — is in scope insofar as it makes the **physics** wrong.

## What you are handed
- The engine script (the artifact above) and its export `S11c_a_exports.py`.
- **The sole physics authority: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`** (HEAD-closed).
  §1 = inherited S11b setup (supplied, unfalsifiable here); §2 = the non-uniform background ansatz, the two
  anchoring branches, the two density representatives, the reserved admissibility computation; §3 = the
  supplied face maps and laws you take the shape derivative *of* (§3a maps/normal-orientation law/measure,
  §3b displacement/velocity/flux/traction/balance/closure, §3c shifted trace + dynamic window); §4 = the
  objects to compute (T-0…T-i); §5 = the routes and controls; §6 = method, dimensions, the three script
  clauses; §7 = names/F9/tag grammar; §8 = supplied vs computed.
- `scripts/S11b_exports.py` (the imported LEDGER) is in the tree if you want to check inherited binds.
- ⛔ You are **not** handed the build directive or any orchestrator commentary. The spec is the authority;
  derive from it, not from the engine's own comments.

## Required method — this is a SCRIPT
Derive independently. **Write your own derivation script BEFORE opening the artifact in earnest, and save
BOTH the script AND its literal stdout to named absolute paths; report those paths. ⛔ Without a runnable
script and its stdout, your derivation claims are discarded** — a prose re-derivation is the exact defect
(a typed conclusion with no computation) this rebuild exists to remove, relocated into the review. You do
not need to reproduce all 13 tasks × branches × faces × DOFs; pick the load-bearing ones and derive them
honestly (at minimum: the outward face normal and its orientation law §3a; the true-area measure shape
derivative; the relative flux `J_s^α` shape derivative §3b; the traction `t_s^α` and virtual work; the
closure shape derivative T-i; and the dynamic-window projection residual T-f).

Probe for: a value verified using the very predicate/definition that produced it; a conclusion emitted as
an unconditional literal; a check whose expected value lives inside the artifact; and a residual that is
zero **by construction** (operand B derived from operand A) rather than by two independent routes.

⛔⛔ **A FORM ABLATION IS MANDATORY — it is not optional, and it is the only thing that has ever caught the
worst defect.** Change the **structure** of a load-bearing object in a `/tmp` COPY of the engine — e.g.
flip a sign **and** an off-diagonal in a normal/traction/Jacobian construction, or **collapse two
independent symbols into one** (e.g. force `W_bg` and `μ_R,bg` to the same profile, or set the two
independent DOFs `δW` and `ζ_c` equal, or collapse the upper/lower face parameter `s`) — then re-run and
report the **literal** diff in the affected `PY_S11CA_*` tag lines. ⭐ A **coefficient** rescale tests only
arithmetic; only a **FORM** change tests physics, because a rescale never leaves the family. If the diff on
a load-bearing tag is **byte-identical** under a structural change, that tag is hand-typed — report it with
the tag name and the two outputs.

⭐ **Ask of every load-bearing claim: WHICH LINE COMPUTED THIS?** Give the line number, or report the object
as uncomputed. ⚠ **Report any `assert` (or raising guard) that precedes the `emit` of the value it guards** —
an assert-before-emit hides a hand-typed payload, because a perturbation strong enough to flip it kills the
process and the leg sees only PASS-or-crash. The spec's three script clauses (§6) require emit-then-compare,
never assert-before-emit.

**One-sided corruption of the engine's own two-route controls (§5a).** The engine claims a direct-Eulerian
route and a material-coordinate route agree (residual emitted). ⛔ A zero residual proves nothing if route B
is derived from route A. In your `/tmp` copy, **corrupt exactly one route at its source** (the spec §5a
names the exact source mutations) and confirm the residual moves and **only that route's operand moves**. If
breaking route A also moves route B, the two routes were never independent and the residual is decoration —
report it.

## Physics filter
Report a finding only if it catches a way the **physics** could be wrong: a hand-typed/undependent payload,
a wrong normal orientation or sign, a dropped or spurious term in a shape derivative, the physical `u`
substituted where the virtual `δ_vu` is required (T-g), a face-shift term evaluated-then-discarded (§3c), a
mis-bound inherited object (e.g. the bulk density `ρ_m`) that corrupts a flux/traction/closure, or a control
that cannot fail. ⛔ Do **not** report "the script would be wrong on a different input", and ⛔ do not report
style. "Nothing survives the filter" is weak evidence — state exactly what you derived, what you ablated,
and the literal diffs.

## Ablation sandbox and operational limits
⛔ Copy the engine to `/tmp` and ablate the **COPY**. ⛔ Never modify the working tree.
This is a **Python/SymPy** engine (no Mathematica licence seat), but a full run is ~**220 s**. Keep your own
derivations small and targeted; when you re-run an ablated full copy, allow it time rather than killing it
early, and save each ablation script + its literal stdout to named absolute paths. As a triage cross-check you may classify each emitted row as computed-vs-declared (a mostly-CONSTANT
engine has not been built) — ⛔ triage only, not a verdict; the FORM ablation above is the real control.

## Output
Findings most-serious first. For each: the spec section you derived from, the exact `PY_S11CA_*` tag and
engine line, the concrete way the physics is wrong, and the **literal** ablation diff or derivation-script
stdout that establishes it. Then answer: (1) is any load-bearing payload hand-typed / byte-identical under a
form ablation? (2) are the §5a two-route controls genuinely independent under one-sided corruption? (3) is
the engine safe to commit as one of the two S11c-a engines?
