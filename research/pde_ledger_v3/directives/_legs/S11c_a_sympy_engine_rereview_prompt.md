# Independent physics review — the S11c-a SymPy engine

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`
(the engine; ~2100 lines). Its second product is
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

Derive these objects yourself from the spec and decide whether the engine computed them correctly. Above
all: (a) is any emitted payload a **hand-typed CAS object with no data dependence on the derivation** (a
genuine SymPy expression, authored, that would not move if the derivation above it were deleted)? — only a
**form ablation** catches this; and (b) does the **§5a representation-invariance package actually compute
each object by two genuinely independent routes**, so its residual is a real test and not zero-by-construction?

## What you are handed
- The engine script (the artifact above) and its export `S11c_a_exports.py`.
- **The sole physics authority: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`** (HEAD-closed).
  §1 = inherited S11b setup (supplied, unfalsifiable here); §2 = the non-uniform background ansatz, the two
  anchoring branches, the two density representatives, the reserved admissibility computation; §3 = the
  supplied face maps and laws you take the shape derivative *of* (§3a maps/normal-orientation law/measure,
  §3b displacement/velocity/flux/traction/balance/closure, §3c shifted trace + dynamic window); §4 = the
  objects to compute (T-0…T-i); **§5 = the routes and controls — §5a specifies route 2 as the
  material-coordinate face-flattening derivation using `w′ = [w−ζ_c]/[W_bg+δW]`, mapped back to Eulerian,
  and a one-sided corruption that mutates only the direct route**; §6 = method, dimensions, the three script
  clauses; §7 = names/F9/tag grammar; §8 = supplied vs computed.
- `scripts/S11b_exports.py` (the imported LEDGER) is in the tree if you want to check inherited binds.
- ⛔ You are **not** handed any build directive or orchestrator commentary. The spec is the authority;
  derive from it, not from the engine's own comments (a comment claiming "an independent construction" is
  not evidence the construction is independent — you must test that).

## Required method — this is a SCRIPT
Derive independently. **Write your own SymPy derivation script BEFORE relying on any engine claim, and save
BOTH the script AND its literal stdout to named absolute paths under `/tmp/s11ca_rr/`; report those paths.
⛔ Without a runnable script and its stdout, your derivation claims are discarded** — a prose re-derivation
is the exact defect (a typed conclusion with no computation) this rebuild exists to remove, relocated into
the review. Pick the load-bearing objects and derive them honestly (at minimum: the outward face normal and
its orientation law §3a; the true-area measure shape derivative; the relative flux `J_s^α` shape derivative
§3b; the traction `t_s^α` and virtual work T-d; the closure T-i; the dynamic-window projection residual
T-f; the T-g virtual constraint; and the T-c′ kinematic-balance operands).

Probe for: a value verified using the very predicate/definition that produced it; a conclusion emitted as an
unconditional literal; a check whose expected value lives inside the artifact; and a residual that is zero
**by construction** (operand B derived from, or algebraically identical to, operand A) rather than by two
independent routes.

⛔⛔ **A FORM ABLATION IS MANDATORY — it is the only thing that has ever caught the worst defect.** Change
the **structure** of a load-bearing object in a `/tmp` COPY of the engine — flip a sign **and** an
off-diagonal in a normal/traction/Jacobian construction, or **collapse two independent symbols into one**
(force `W_bg` and `μ_R,bg` to one profile; set the two DOFs `δW` and `ζ_c` equal; collapse the face sign
`s`) — then re-run and report the **literal** diff in the affected `PY_S11CA_*` tag lines. A **coefficient**
rescale tests only arithmetic; only a **FORM** change tests physics. If a load-bearing tag is
**byte-identical** under a structural change, that payload is hand-typed — report it with the tag name and
both outputs. A full engine run is ~270 s — allow it time.

⭐ **Ask of every load-bearing claim: WHICH LINE COMPUTED THIS?** Give the line number, or report the object
as uncomputed. ⚠ **Report any `assert`/raising guard that precedes the `emit` of the value it guards** — it
hides a hand-typed payload, since a perturbation strong enough to flip it kills the process and the leg sees
only PASS-or-crash. The spec's §6 script clauses require emit-then-compare, never assert-before-emit.

⛔⛔ **MANDATORY: prove or refute the §5a two-route independence by ONE-SIDED CORRUPTION.** The engine emits a
representation-invariance residual (`PY_S11CA_REP_INVARIANCE_*`) claiming each of T-c, T-d, T-i, T-g (and
T-h) is computed by two independent routes (Eulerian direct vs material-coordinate `w′` face-flattening). ⛔
A zero residual is worthless if route 2 is derived from, or is algebraically the same expression as, route 1.
In a `/tmp` copy, **corrupt exactly one route at its source** (per §5a: reverse only the upper-face `x¹`
first jet of `W_bg` in the direct-route source; or mutate only the material-route source) and confirm that
**only that route's operand moves and the other route's operand stays fixed**. If corrupting route 1 also
moves route 2 (or the two operands are byte-identical expressions before any corruption), the two routes are
not independent and the residual is decoration — report it per object, with the literal operands. Do this
for each of T-c, T-d, T-i, T-g, and T-h separately; they need not all behave the same.

Also confirm the **§5a one-sided independence control the engine itself emits**
(`PY_S11CA_CONTROL_INDEPENDENCE_*`) mutates only the direct route and leaves the material route uncorrupted
— not a shared source that both routes read.

## Physics filter
Report a finding only if it catches a way the **physics** could be wrong: a hand-typed/undependent payload,
a wrong normal orientation or sign, a dropped or spurious term in a shape derivative, the physical `u`
substituted where the virtual `δ_vu` is required (T-g), a face-shift term evaluated-then-discarded (§3c), a
mis-bound inherited object (e.g. the bulk density `ρ_m`) corrupting a flux/traction/closure, a two-route
control whose routes are not independent, or a control that cannot fail. ⛔ Do not report "wrong on a
different input" or style. "Nothing survives" is weak evidence — state exactly what you derived, what you
ablated, and the literal diffs.

## Ablation sandbox and operational limits
⛔ Copy the engine to `/tmp` and ablate the **COPY**. ⛔ Never modify the working tree.
This is a **Python/SymPy** engine (no Mathematica licence seat), but a full run is ~**270 s**. Keep your own
derivations small and targeted; when you re-run an ablated full copy, allow it time rather than killing it
early, and save each ablation script + its literal stdout to named absolute paths under `/tmp/s11ca_rr/`.

## Output
Findings most-serious first. For each: the spec section you derived from, the exact `PY_S11CA_*` tag and
engine line, the concrete way the physics is wrong, and the **literal** ablation diff or derivation-script
stdout that establishes it. Then answer explicitly: (1) is any load-bearing payload hand-typed / byte-identical
under a form ablation? (2) for EACH of T-c, T-d, T-i, T-g, T-h — are the §5a two routes genuinely independent
under one-sided corruption (only the corrupted route moves)? (3) is T-c′'s kinematic-balance emitted as two
genuinely-computed operands + residual, or a bare/hand-typed zero? (4) are the primary T-0…T-i payloads
correct and form-ablation clean? (5) is the engine safe to commit as one of the two S11c-a engines? Report
the absolute paths of every derivation and ablation script + stdout you saved.
