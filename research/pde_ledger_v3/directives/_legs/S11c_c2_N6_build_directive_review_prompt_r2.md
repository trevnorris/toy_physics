# Independent review (round 2) — S11c-c2 N6 build directive + §5c route-2

## Artifacts (read them LAST — see method)
- Directive: `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_sympy_build_directive.md`
- Governing spec §5c (revised): `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md`
  (§5c ≈ 303-372)

These were revised after a round-1 review. This is a fresh, independent review of the CURRENT versions. Working dir
`/var/projects/toy_physics`; paths repo-root-relative unless absolute. Both are physics-bearing; review by **reading**
(the companion script does not exist yet — defer executable script-control tests to the build, but flag any control
that cannot be built as specified).

## What you are checking
Whether the directive + §5c **faithfully and tractably** commission the corrected N6 control **without leaking an
answer or a recipe**, so the build can proceed with no further change to what will be computed or claimed.

## Sources of truth (read FIRST, form your OWN account, THEN open the directive)
1. §5c itself (the governing object) + §3a (close), §3c (weak increment + slot-linearity + **`extract`**), §1d
   (routing), §5d (density), §6 (dimensions/multigrade). ⛔ You are NOT handed, and there is NOT, any expected value
   for the residual — §5c fixes no target.
2. The actual code the object depends on — read and reason from it, do not take the prose on trust:
   - the §3c weak extraction `scripts/S11c_c2_selfenergy_fold_sympy_audit.py:325-342` (`def extract`);
   - the N4 material map `scripts/S11c_b_brane_operator_sympy_audit.py:1916-1981` (`material_pullback`) and how
     `build_operator(route="MATERIAL")` uses it `:2762-2783`;
   - the trial-field-only `representation_pullback` `scripts/S11c_c2_selfenergy_fold_sympy_audit.py:1136-1171`;
   - the S11c-a material face-flattening `scripts/S11c_a_interface_geometry_sympy_audit.py:703-815`;
   - the audit's cross-anchoring loop `scripts/S11c_c2_selfenergy_fold_sympy_audit.py:1087-1107` and the density loop
     `:1085-1086`;
   - the carrier scaffold `scripts/S11c_c2_FG_diagnostic_sympy.py:99-134`.

## The ONE load-bearing physics question — vet it independently from the code, not from the prose
§5c route 2 is written **map-then-extract**: `I_{M→E} = extract( T_{M→E}[ close(SLAB_M) − SLAB_M ] )`. Decide, from
the actual `extract` (`:325-342`) and `material_pullback` (`:1942-1981`), whether that ordering is correct:
- Is `extract` a **lossless** read of the full operator, or a **lossy projection** (does it `restrict` to
  Helmholtz/transverse-longitudinal components, take div/curl, pair with fixed Eulerian test functions)?
- Does `T_{M→E}` **mix sectors** (e.g. `θ → θ + u·∇ρ⁰/ρ`)?
- Therefore: do `extract` and `T_{M→E}` **commute**? If they do not, is map-then-extract (map the carrier, then the
  single Eulerian extract) the correct order — and is extract-then-map (`T[extract(...)]`) wrong because it drops the
  advective off-diagonal? Does the parent (`build_operator("MATERIAL")`) do map-then-extract or extract-then-map?
- Is it right that there is **no pre-map `SLAB_M` callable** (that `build_operator(route="MATERIAL")` already folds
  the map in), so `SLAB_M` must be built natively unmapped and `T_{M→E}` applied once as `material_pullback`?

Report your independent conclusion with the code lines you relied on. If you conclude the ordering is wrong in either
direction, that is a spec finding.

## Also substantiate
1. **Native second route / no double-map / no tautology:** route 2 built natively (S11c-a face-flattening + material
   sources), `T` applied exactly once, closed with the SAME imported same-`α` c1 response; routes structurally
   independent so `R_N6` is not zero by construction; the one-sided controls (tilt on Eulerian; omit `u·∇ρ⁰` inside
   `T` only) each move one route; `A−A` → emit `RHOBR_CONSTANT` structural absence.
2. **Tractability:** carrier-first BOTH routes, ⛔ no `build_case` end-to-end; coefficients truncated/evaluated BEFORE
   closure expansion; residual + slot-linearity guard + controls all tested by exact finite-field PIT (several primes,
   shared generic samples, denominators cleared, singular/branch cells handled, FN bound). Is anything (building
   `SLAB_M`, applying `T`) still forcing whole-object symbolic materialization?
3. **Audit correction:** an EDIT (not companion-only) replacing `:1087-1107` with raw `ANCHORING_L_MINUS_M =
   increment[LAB] − increment[MATERIAL]` (§5c:366 orientation, no `representation_pullback`), density loop
   `:1085-1086` kept.
4. **Leakage / iterate-to-target:** residual framed as a computed finding, no supplied target, never a builder
   exit/assert/pass; builder barred from the `_measurements/` records; tag names / supplied text do not encode the
   sign/shape of the answer; anti-examples placeholder-only.
5. **Dimensions:** `[L,T,M]` restored + able-to-fail dimensional consistency + `(ε,η,σ_W)` multigrade on every
   emitted object (§6 / `N12`).
6. **Missing premise:** anything the builder needs that is prose where it should be an equation or a named source.

## Physics filter
Report a finding only if it catches a way the built N6 control could be **wrong, vacuous, intractable, or
answer-leaking** — not a stylistic preference.

## Output
Findings each with: the quote (with file:line), the source-of-truth quote it violates, why it matters, minimal fix.
If nothing outstanding changes what will be computed or claimed, say **CLEAR TO BUILD**. Evidence-first, brief.
