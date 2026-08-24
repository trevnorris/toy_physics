# Independent physics review — the S11c-a blind Wolfram engine

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl`
(~1570 lines). It is a **blind** engine: it imports nothing and re-derives every object from the spec's
equations, streaming `WL_S11CA_<QUANTITY>` tags to stdout in Mathematica `InputForm`. A full run is heavy
(~20 minutes, peaks ~9-10 GB RSS).

## What to check
This engine takes the **closed** S11c-a spec and computes the **first-order shape derivative of every S11b
interface law** once the uniform brane is un-frozen: `W_0 → W_bg(x)`, `mu_R → mu_R,bg(x)`. It emits, for both
anchoring branches (`LAB_HELD` / `MATERIAL_ADVECTED`), both faces, both face DOFs `{delta_W, zeta_c}`, and
both density representatives, the objects **T-0 … T-i**: background density map, outward face normal,
conormal, true-area measure and its shape derivative, face velocity, relative flux, kinematic balance,
traction + virtual work, shifted-face trace, dynamic-window projection (+ static operand, residual,
term-origin inventory), virtual material constraint, physical sourced mass balance, and the face-closure
shape derivative — plus four control families (§5a representation-invariance two-route + one-sided
independence; §5b per-direction form ablation; §5c uniform-limit; §6 homogeneity).

Derive these yourself from the spec and decide whether the engine computed them correctly. Above all:
(a) is any emitted payload a **hand-typed CAS object with no data dependence on the derivation** (a genuine
Mathematica expression, authored, that would not move if the derivation above it were deleted)? — only a
**form ablation** catches this; and (b) does the **§5a representation-invariance package actually compute
each object by two genuinely independent routes** (a direct Eulerian route and a material-coordinate `w'`
face-flattening route), so its residual is a real test and not zero-by-construction?

## What you are handed
- The engine `.wl` (the artifact above).
- **The sole physics authority: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`** (HEAD-closed).
  §1 = inherited setup (supplied, unfalsifiable here); §2 = the non-uniform background ansatz, the two
  anchoring branches, the density representatives, the reserved admissibility premise; §3 = the supplied
  face maps, orientation law, measures, laws you take the shape derivative *of* (§3a maps/normal/measure,
  §3b displacement/velocity/flux/traction/balance/closure, §3c shifted trace + dynamic window); §4 = the
  objects to compute (T-0…T-i); **§5 = routes and controls — §5a specifies route 2 as the
  material-coordinate face-flattening derivation using `w' = [w−zeta_c]/[W_bg+delta_W]`, mapped back to
  Eulerian, and a one-sided corruption that mutates only the direct route**; §6 = method, dimensions, the
  three script clauses; §7 = names/tag grammar; §8 = supplied vs computed.
- ⛔ You are **not** handed any build directive, the sibling SymPy engine, or orchestrator commentary. Derive
  from the spec, not from the engine's own comments (a comment claiming "an independent construction" is not
  evidence — you must test it).

## Required method — this is a SCRIPT
Derive independently. **Write your own derivation script (SymPy or Wolfram, your choice) BEFORE relying on any
engine claim, and save BOTH the script AND its literal stdout to named absolute paths under `/tmp/s11ca_wl/`;
report those paths. ⛔ Without a runnable script and its stdout, your derivation claims are discarded** — a
prose re-derivation is the exact defect this rebuild exists to remove, relocated into the review. Pick the
load-bearing objects and derive them honestly (at minimum: the outward normal and its orientation law §3a;
the true-area measure shape derivative; the relative flux shape derivative §3b; the traction and virtual
work T-d; the closure T-i; the dynamic-window projection residual T-f; the T-g virtual constraint; and the
T-c′ kinematic-balance operands).

Probe for: a value verified using the very predicate/definition that produced it; a conclusion emitted as an
unconditional literal; and a residual that is zero **by construction** (operand B derived from, or
algebraically identical to, operand A) rather than by two independent routes.

⛔⛔ **A FORM ABLATION IS MANDATORY — it is the only thing that has ever caught the worst defect.** Change the
**structure** of a load-bearing object in a `/tmp` COPY of the `.wl` — flip a sign **and** an off-diagonal in
a normal/traction/Jacobian construction, or **collapse two independent symbols into one** (force `W_bg` and
`mu_R,bg` to one profile; set the two DOFs `delta_W` and `zeta_c` equal; collapse the face sign) — then
re-run and report the **literal** diff in the affected `WL_S11CA_*` tag lines. A **coefficient** rescale
tests only arithmetic; only a **FORM** change tests physics. If a load-bearing tag is **byte-identical**
under a structural change, that payload is hand-typed — report it with the tag name and both outputs.

⭐ **Ask of every load-bearing claim: WHICH LINE COMPUTED THIS?** Give the line number, or report the object
as uncomputed. ⚠ **Report any `Abort`/`Quit`/`Throw`/`Check` (or a guard that halts) that precedes the emit
of the value it guards** — it hides a hand-typed payload, since a perturbation strong enough to flip it kills
the process and you see only PASS-or-crash.

⛔⛔ **MANDATORY: prove or refute the §5a two-route independence by ONE-SIDED CORRUPTION.** The engine emits a
representation-invariance residual (`WL_S11CA_REP_INVARIANCE_*`) claiming each of T-c, T-d, T-i, T-g (and
T-h) is computed by two independent routes. ⛔ A zero residual is worthless if route 2 is derived from, or is
algebraically the same expression as, route 1. In a `/tmp` copy, **corrupt exactly one route at its source**
(per §5a: reverse only the upper-face `x^1` first jet of `W_bg` in the direct-route source; or mutate only
the material-route source) and confirm that **only that route's operand moves and the other stays fixed**. Do
this for each of T-c, T-d, T-i, T-g, T-h separately; they need not all behave the same. If corrupting route 1
also moves route 2, or the two operands are byte-identical before any corruption, the routes are not
independent and the residual is decoration — report it per object with the literal operands. Also confirm the
engine's own `WL_S11CA_CONTROL_INDEPENDENCE_*` mutates only the direct route, not a shared source.

Also check **T-c′ `KINEMATIC_BALANCE`**: two genuinely-computed operands + residual, or a bare/hand-typed
zero? And regression-check the primary T-0…T-i payloads are correct and form-ablation clean.

## Physics filter
Report a finding only if it catches a way the **physics** could be wrong: a hand-typed/undependent payload, a
wrong normal orientation or sign, a dropped/spurious shape-derivative term, the physical `u` used where the
virtual `delta_v u` is required (T-g), a face-shift term evaluated-then-discarded (§3c), a mis-derived
inherited object (e.g. `rho_m` in a flux/traction/closure), a two-route control whose routes are not
independent, or a control that cannot fail. ⛔ Not "wrong on a different input", not style. "Nothing survives"
is weak evidence — state exactly what you derived, what you ablated, and the literal diffs.

## ⛔ Mathematica operational constraints — Wolfram has a TWO-SEAT licence
⛔ Copy the `.wl` to `/tmp` and ablate the **COPY**. ⛔ Never modify the working tree.
⛔ **Wrap EVERY kernel run in `timeout 1800`** (a full engine run is ~20 min at ~9-10 GB; 1800 s leaves margin).
A 1800 s hit is a FAILED ablation — report it and move on; ⛔ never raise the timeout.
⛔ **Run at most ONE kernel at a time** (the licence has two seats and the other is reserved). ⛔ Never launch
two `wolframscript` runs concurrently.
⭐ To keep an ablation fast you MAY comment out heavy tasks you are not testing (the `CONTROL_FORM` §5b task
dominates runtime) — but **report exactly what you commented out**, and never comment out the object you are
ablating or its inputs.
⭐ Save every ablation script AND its literal stdout to named absolute paths under `/tmp/s11ca_wl/`, and
report those paths. ⛔ If a run dies with a healthy-looking log, check `free -h` and `ps -eo pid,rss,comm
--sort=-rss | head` for an orphaned kernel before anything else.

## Output
Findings most-serious first. For each: the spec section you derived from, the exact `WL_S11CA_*` tag and
engine line, the concrete way the physics is wrong, and the **literal** ablation diff or derivation-script
stdout that establishes it. Then answer explicitly: (1) any load-bearing payload hand-typed / byte-identical
under a form ablation? (2) for EACH of T-c, T-d, T-i, T-g, T-h — are the two §5a routes genuinely independent
under one-sided corruption? (3) is T-c′ two genuinely-computed operands + residual, or a bare zero? (4) are
the primary T-0…T-i payloads correct and form-ablation clean? (5) is the engine safe to commit as one of the
two S11c-a engines? Report the absolute paths of every derivation and ablation script + stdout you saved.
