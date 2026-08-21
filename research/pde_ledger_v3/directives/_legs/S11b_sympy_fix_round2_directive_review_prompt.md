# Directive review — S11b SymPy engine fix-round-2

You are one of two independent legs reviewing an **orchestrator-written repair directive** BEFORE any builder
runs it. Document review: does the directive specify the right repair, completely, without leaking a value,
without over-specifying, without breeding new defects, and without disturbing correct physics?

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_sympy_engine_fix_round2_directive.md`
and its rule-2 twin
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11b_sympy_engine_fix_round2_directive.md`.

## Context
The engine `research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py` (baseline `6d57b27e`)
had a round-1 emission-fidelity repair that fixed the impedance (now computed) and the big §6 checks (now
bite). This round-2 directive makes genuine a **tail of four checks** that still cannot fail under a one-sided
corruption of what they police. The physics is already correct and must not change.

## What you are handed (verify against these, ⛔ do not take the directive on faith)
- The engine (baseline `6d57b27e`) and its product `S11b_exports.py`.
- The physics authority `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md` — esp. §6 causality
  diagnostics + kernel-orientation (~381–406), the Onsager conversion + cross-check (~733–742), B0a parity
  (~597–603), B2d two-port (~704–711), B8 controls.
- `CLAUDE.md` (rules 2, 3, 5, 11, 15) + `.claude/skills/build/SKILL.md` (three script clauses; anti-tautology
  corollaries; the "hardening a check adds no-physics lines" measurement).

## The four fixes (verify each is real and correctly specified, against the engine + spec)
1. `parity_even_jw` (~669): even branch `integral_coefficient*(Omp*je + (-Omp)*je) ≡ 0`, typed pair →
   rebuild from the live `source_finite` under the parity/interval assumptions, with a **narrow** claim
   (odd-integrand corruption legitimately leaves even=0 on a symmetric interval).
2. `control_no_reciprocal_traction` (~1919–1942, `export=False`, a B8 control): both routes drop `Λ_X` →
   genuinely test the dropped-channel trap via independent routes (cut vs full B2a after `Λ_X⁰=0`; the
   full-minus-cut effect vs the independent placeholder constructor ~901–938; power vs the live B2d identity
   ~1332–1343).
3. `onsager_reciprocity` cross-check (~1370–1385): `force_res ≡ −flux_res` ⇒ cross-check is `2·flux_res` →
   derive the two conversions separately from the live mixed law, normalise each up to a scalar/sign, compare
   (spec-mandated, ~733–742).
4. `kernel_orientation_identities` (~942–957): both operands use the module-global `Λ_I` → compare a
   separately-constructed retarded reference (`Λ_I⁰, τ_I, ω`) against `K_I` from the live equations
   (spec-decisive, ~381–406); `causality_check` (~1123–1127) + growth/decay (~1518–1522) carry the same
   payload and repair with it.

## What to check — report a finding only if it changes the repair or its outcome
1. **Accuracy.** Is any fix mis-described, any line reference wrong, or any listed defect actually already
   genuine? For each, is a genuine independent second route actually **available** in the engine (verify), or
   would the honest disposition be to report an obstruction (§13) rather than force a check the spec does not
   ask for?
2. **Completeness.** Is any same-class exported tautology (residual inert under a one-sided corruption of its
   policed object) present in the engine but MISSING from these four? Name it with a line. Is the
   `causality_check`/growth-decay propagation correctly handled (they must not be treated as separate)?
3. **Leak (rule 5).** Does the directive state or let a builder infer any expected physics VALUE, sign, ratio
   or count?
4. **Over-specification (rule 3).** The fixes cite the spec's own construction (separate reference, two
   conversions, the placeholder constructor). Is that naming the object + property, or a line-by-line rewrite
   that removes the builder's independent construction?
5. **Breeding (rule 15 / step-run).** Could a builder satisfy this by ADDING no-physics lines or a new
   tautology (e.g. a moving ornament beside an inert residual)? Is the one-sided-corruption acceptance
   decisive, and does the "no shared constructor" rule close the loophole?
6. **Scope.** Is "LEAVE homogeneity" correct — verify `trace_dimension` (~1598–1624) recurses into nested
   Adds so `homogeneity_{thickness,mass,affinity}` can already fail? Are the value/wiring preservation guards
   right? Could any fix inadvertently change a physics VALUE or the export wiring?
7. **Spec fidelity.** Do the four fixes match what §6 / B2d / B0a / B8 actually mandate, or would any
   manufacture a check the spec does not ask for (esp. the parity "narrow claim" and the onsager two-route
   comparison)?

## Method
Read the sources first, form your own view of each defect and its available fix, THEN judge the directive.
Quote the exact engine/spec line for any claim of fact. ⛔ A prose assertion with no citation is discarded.

## Output
Findings most-severe first (directive location · engine/spec citation · the fix), then per-axis (1–7), then a
one-line verdict: safe to build from as-is, or must something be fixed first?
