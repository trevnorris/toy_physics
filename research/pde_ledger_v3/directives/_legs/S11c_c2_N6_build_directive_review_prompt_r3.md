# Independent review (round 3) — S11c-c2 N6 build directive + §5c

## Artifacts (read LAST)
- Directive: `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_sympy_build_directive.md`
- §5c (revised): `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md` (§5c ≈ 303-378)

Revised after two prior review rounds. Fresh, independent review of the CURRENT versions. Working dir
`/var/projects/toy_physics`; repo-root-relative paths unless absolute. Physics-bearing; review by reading (companion
script does not exist — defer executable script-control tests to the build).

## Settled (⛔ do not re-litigate unless you find NEW code evidence overturning it)
**Map-then-extract is correct** — `extract` (`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:325-342`) is a lossy
Eulerian projection; the N4 map mixes sectors; they do not commute; parent does map-then-extract. Two prior rounds
independently confirmed this from the code.

## The DECISION this round turns on
The route-2 CAS pipeline has been **deliberately delegated to the builder (astra)** — the directive now hands the
**object + verified refs + a definition-of-done** (§3) instead of a hand-authored recipe (two prior recipes were
type-wrong). Your job: decide whether that object + DoD is **correct and complete enough that a competent builder
will produce the RIGHT route-2**, with the build legs able to gate it — or whether a gap remains that would let a
wrong/vacuous route-2 pass. Reason from the code:
- native material δp-slot carrier from **S11c-a MATERIAL face sources** (traction / `closure_shape_deriv` / `V_s`,
  `scripts/S11c_a_interface_geometry_sympy_audit.py:703-845`) folded into the S11c-b face-fold δp slots; that S11c-a
  builder already maps its covector to Eulerian (`:742-755`);
- the `N4` field map defined by `material_pullback` (`scripts/S11c_b_brane_operator_sympy_audit.py:1916-1981`) — a
  bulk-density quadratic map (2nd scale-derivative annihilates linear rows; never touches δp slots);
- the trial-only `representation_pullback` (`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:1136-1171`).

Check specifically:
1. **Is the DoD sufficient + self-consistent?** The 6 DoD clauses (native not re-pullback; scalar-level map + Eulerian
   variation after + single extract; no annihilation; reverse u-row channel present; no double-map; carrier-first).
   Is any load-bearing constraint missing (e.g. how the imported same-`α` c1 response binds while staying opaque to
   the map; which face velocity closes route 2; how `SLAB_M` δp slots line up with route 1's for a meaningful diff)?
2. **Could a builder satisfy every DoD clause and still compute the WRONG object?** If yes, name the loophole.
3. **Is the open design question** (whether a bulk field map is needed at all, vs native material faces + material
   `V_s` already yielding the material increment) correctly left to the build, or does it hide a spec ambiguity that
   must be resolved BEFORE the build?

## Also verify the fixes landed + nothing regressed
- **F2 advection:** the probe omits `u·∇ρ⁰/ρ` inside the map; structural absence is **RHO4_CONSTANT** (`ρ₄=ρ_br/W₀`,
  `brane:1885-1890`), live probe is **RHOBR_CONSTANT**. Confirm §5c + directive now say this (⛔ not the reverse).
- **F3 tilt:** built through an injectable Eulerian carrier factory, PIT-checked vs the imported carrier (⛔ not the
  c1-DtN-jet override, `selfenergy_fold:384-387`).
- **Audit edit:** replace `:1087-1107` with raw `ANCHORING_L_MINUS_M = increment[LAB] − increment[MATERIAL]`
  (§5c orientation, no `representation_pullback`); keep density loop `:1085-1086`.
- **Tractability / leak / dimensions:** carrier-first both routes, no `build_case` end-to-end, PIT with FN bound;
  residual a finding with no target, never a builder exit; builder barred from `_measurements/`; `[L,T,M]` +
  able-to-fail + `(ε,η,σ_W)` on every emit.

## Physics filter
Report a finding only if it catches a way the built N6 control could be **wrong, vacuous, intractable, or
answer-leaking** — not a stylistic preference.

## Output
Findings each with file:line quote, the source-of-truth it violates, why it matters, minimal fix. If nothing
outstanding changes what will be computed or claimed, say **CLEAR TO BUILD**. Evidence-first, brief.
