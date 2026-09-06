# Final directive review (round 4) — S11c-c2 N6: integration of the cleared route-2 spec

## Artifacts (read LAST)
- Directive: `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_sympy_build_directive.md`
- §5c: `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md` (§5c ≈ 303-393)

Working dir `/var/projects/toy_physics`; repo-root-relative paths unless absolute. Physics-bearing directive review by
reading (companion script does not exist — defer executable script-control tests to the build).

## What already CLEARED (⛔ do not re-litigate; verify only that it was carried faithfully)
The **route-2 construction spec** `research/pde_ledger_v3/_measurements/S11c_c2_N6_route2_spec_astra.md`
(astra-authored) was reviewed by a fresh Claude agent + Grok and is **SPEC CLEAR**. Its load-bearing calls — verified
against the code, incl. the grade fact on the real imported operator — are settled:
- route 2 has **NO `T` on the increment**; it is the parent `task_rep_invariance` pattern applied to the c2 face-slot
  increment (native material face sources + material μ_θ, differenced directly);
- **μ_θ binding = material** (route 2), Eulerian (route 1); geometry-only would make the advection probe vacuous;
- the **reverse-u channel is grade-suppressed** out of the retained rectangle — a mandatory-survival requirement is
  wrong; report reverse-block sensitivity + permit computed absence.

## What THIS round checks — the INTEGRATION (did the fold land faithfully + consistently?)
1. **§5c ↔ astra spec consistency.** Does §5c (formula `I_{M→E}=extract(close(SLAB_M)−SLAB_M)`, the "no T" text, the
   μ_θ-binding bullet, the reverse-u grade-suppression bullet, the advection probe = material-μ_θ tag t) faithfully
   match the cleared astra spec? Any residual **stale `T_{M→E}`-on-increment / map-then-extract / "reverse channel
   must survive" language** anywhere in §5c or the directive that now contradicts the cleared construction?
2. **Directive ↔ §5c ↔ astra spec consistency.** The directive defers route 2 to the astra spec and carries a 6-clause
   DoD. Does the DoD match the cleared spec (esp. clause 4 = reverse-u grade-suppressed/permit-absence, NOT mandatory)?
   Do the §1 definitions (SLAB_M, μ_θ binding, close opaque), §2 tractability, §4 controls (tilt factory; advection
   tag t at both sites, RHO4 absent / RHOBR live), §5 audit edit (`ANCHORING_L_MINUS_M` = LAB−MATERIAL raw), §7 output
   + dimensions, all agree with §5c and the astra spec?
3. **Leak / iterate-to-target.** The directive now points the builder at the astra spec (exempted from the
   `_measurements/` bar). Confirm the astra spec **withholds** the R_N6 target (no expected value, no "must vanish"),
   so exempting it does not leak; residual-zero is never a builder exit; other `_measurements/` records stay barred.
4. **Buildability gap correctly flagged.** The "no ready-made `SLAB_M` callable — a thin material face-fold adapter
   must be built" gap (astra §8) is carried, so the builder is not sent to a nonexistent function.

## Physics filter
Report a finding only if it catches a way the built N6 control could be **wrong, vacuous, intractable, or
answer-leaking**, OR a directive/§5c internal contradiction that would misdirect the builder — not a stylistic nit.

## Output
Findings each with file:line quote, what it contradicts, why it matters, minimal fix. If nothing outstanding changes
what will be computed or claimed, say **CLEAR TO BUILD**. Evidence-first, brief.
