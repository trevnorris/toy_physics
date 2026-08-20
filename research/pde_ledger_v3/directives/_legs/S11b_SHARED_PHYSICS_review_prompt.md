# Independent review — the unified S11b SHARED PHYSICS spec

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md`
Its rule-2 substantiation (the command + literal output behind each provenance claim):
`/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11b_SHARED_PHYSICS.md`

## What this is
The **shared physics specification** for step **S11b — the interface coupling law**, written to be inserted
byte-identically into two independent engine build directives (a SymPy engine and a **blind** Wolfram
engine that imports nothing and re-derives). It is **ONE unified step** that subsumes two historical
execution stages — **S11b-A** (the bulk's response to moving faces) and **S11b-B** (the homogeneous
interface assembly) — under a set of already-settled decisions. It supplies setup, governing equations,
supplied constitutive laws/routes, premises, and a task list; it must supply **no result**. Because both
engines read it, **an error here defeats dual-engine agreement by construction** — both engines compute the
same wrong thing and agree. Your job: find where the spec is **wrong in a way that would corrupt the
build**.

## Sources of truth — read these FIRST, form your OWN view, THEN read the spec
- **The physics** (the two historical stages the spec re-expresses):
  - `research/pde_ledger_v3/directives/S11bA_SHARED_PHYSICS.md`, `research/pde_ledger_v3/steps/S11bA_interface_response.md`
  - `research/pde_ledger_v3/directives/S11bB_SHARED_PHYSICS.md`, `research/pde_ledger_v3/steps/S11bB_interface_assembly.md`
- **The settled decisions the spec must honour** (structure / naming / provenance / scope):
  `research/pde_ledger_v3/directives/S11b_unified_decisions.md` (the "G-series", G1–G14) and its twin
  `research/pde_ledger_v3/directives/_measurements/S11b_unified_decisions.md`.
- **The ownership boundary S11 → S11b:** `research/pde_ledger_v3/steps/S11_stray_longitudinal.md` (near the
  end — what S11 defers, and to whom).
- **The chain rules the spec inherits** (⛔ for context on what the spec should POINT AT rather than
  restate): `research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` (F1–F9);
  `research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md` (the blind-Wolfram control);
  `research/pde_ledger_v3/directives/S11_C17_C18_spec_repair_decisions_v2.md` (the frozen `T7` comparator
  contract). A structural reference for the SHARED_PHYSICS pattern: `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`.

## Method (DOCUMENT review)
Read the physics sources first and form an independent view of what S11b-A and S11b-B actually established,
how their objects relate, and which of their statements are **results** (values the reduction is expected
to reproduce). **THEN** read the spec. For every finding, quote BOTH sides (source and spec) with
`file:line`. ⛔ A prose claim with no `file:line` citation is discarded. ⛔ **Do NOT trust the spec's own
citations or the twin's** — open each cited line yourself and confirm it says what is claimed. (Citation
errors have been found and fixed during authoring of the predecessor artifacts; assume more may remain.)

## The specific questions to settle — independently (read/derive; do not take the spec's word)

1. **The withhold discipline (the highest-value hunt).** This spec's predecessor drafts leaked answers. Find
   **every** place the spec **states a result the engines are meant to compute**, or **presupposes the FORM
   of an answer**. In particular, from your reading of the sources, identify the objects the sources treat
   as **acceptance values or derived results** — e.g. the permeable-face impedance and its regime forms; the
   per-face inertial-loading magnitude; the coefficient mapping fixed by the affinity; the reduced permeable
   response the μ_s=0 slice must reproduce; the breathing-mode stability quantity and its boundary; the
   closed energy-invariant count; the passivity region; the discarded convective correction's order — and
   check that the spec **withholds each** (names the object as an obligation-to-compute, diffs on the
   orchestrator's side) rather than printing it. ⭐ Flag any leak, and separately flag any task that fixes a
   power law, sign, functional form, or region the engine is supposed to discover. ⚠ One historical spec
   defect was ordering an engine to *"report the order in `v₀/c_s0`"* — a presupposed form both engines then
   obeyed and agreed on; check the spec does not repeat that class of defect anywhere.

2. **Provenance (G2–G6).** For each object, is its provenance right — **imported** (`c_s0`, `μ_R`),
   **originating here** (`ρ_m`, the background drain), **derived** (the bulk impedance and its regimes, the
   face response, the EOM, the dispersion), **supplied** (the affinity closure and its normalization, the
   `Λ_X` reciprocal-traction channel, the §1b branch prescription, the balance-law route), or a
   **reduction-slice** (the pressure-only, equal-time μ_s=0 case)? Specifically: is the canonical closure
   genuinely the **affinity** law with the pressure-only closure recovered as a `μ_s=0` slice — or has the
   spec inverted that, or aliased the two coefficients? Is the historical impermeable/permeable bulk
   response correctly **derived here** rather than smuggled in as "established input"? Is the background
   drain named **apart** from any in-plane background velocity, and is its role confined to a scope limit?

3. **Supplied physics faithfully carried?** Compare the spec's supplied blocks — the branch prescription,
   the affinity/closure, the balance-law virtual-displacement rule and its forbidden routes, the causality
   diagnostics, the two mandatory convention cross-checks, the three energy-accounting discriminators, the
   two-port power identity — against the sources. Is any **load-bearing** physics **dropped, weakened, or
   altered** in the re-packaging? Is any check turned **tautological** (an operand differenced against its
   own definition), or any "residual" that would be zero by construction for any input?

4. **Scope (G14).** Are the questions the spec places **in scope** genuinely answerable on a **uniform**
   background — specifically the assembled longitudinal mode's **fate** and its behaviour **at/near the
   sound-cone grazing threshold**? Is anything **non-uniform** (the variable-coefficient slab spectrum,
   unconditional confinement, non-uniform leakage rates) wrongly pulled **into** this step, or anything
   **uniform** wrongly **deferred out**? Does any tag or sentence imply a **uniform** result settles a
   **non-uniform** question? Is the DC/harmonic/sideband **radiation** audit correctly assigned to the
   nonlinear programme and **not** to S11b-C?

5. **Point at the obligation, don't restate it.** Find any place the spec **restates** an inherited
   obligation (F1–F9, the blind-Wolfram control, the frozen `T7` comparator contract) in **weaker or
   different words** instead of pointing at it — a known failure mode where the re-wording comes out weaker.
   (⚠ The three script clauses and the no-VERDICT rule are **required verbatim** in every script directive —
   those being present in full is correct, not a restatement defect.)

6. **Self-containment and the cut denylist.** The blind engine reads only this shared block. Is the spec
   **self-contained** — does any task depend on a file or object the spec does not itself supply? Separately:
   the historical stages carried a "WHAT NEITHER ENGINE MAY READ" **denylist**, which is cut (a denylist
   means the architecture is wrong; blindness is by **absence**). Confirm the spec carries **no** such
   denylist, and that its blindness claim rests only on the blind engine importing nothing.

7. **Object naming and one-tag-per-object.** Any object named such that an object-comparison across the
   export chain would prove a **false equal** (the reused-key hazard)? Any bundled or split named object that
   would make the two engines' tag sets uncomparable?

## Physics filter
Report a finding only if it catches a way a **decision or a physics statement in the spec is wrong** — a
leaked value or presupposed form, a mis-assigned provenance, dropped/weakened supplied physics, a wrong
scope boundary, a tautological check, a restated-not-pointed obligation, a self-containment gap, or a
name/tag hazard. ⛔ Do NOT report style, phrasing, or "you could add more."

## Hygiene
Read-only. ⛔ Do not modify the working tree. If you run any check (e.g. a SymPy script to test a
symbol-equality or a dimensional claim), save the script AND its literal stdout to a named absolute path and
report that path. ⛔ A prose "I checked and it's fine" is discarded.

## Output
A list of findings, each: the section/task of the spec, the defect, both-sided `file:line` quotes, and why
it would corrupt the build. If nothing survives the filter for a question, say so plainly. ⛔ Do not report
what any physical quantity "comes out to be." End with the single most important thing you would change
before this spec is folded.
