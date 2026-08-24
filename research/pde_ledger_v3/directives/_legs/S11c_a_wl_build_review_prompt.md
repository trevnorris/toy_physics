# Review — the S11c-a blind Wolfram (WL) engine BUILD DIRECTIVE (wiring + blindness, not physics)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_wl_build_directive.md`

This is an **orchestrator-written build directive** for the S11c-a **blind Wolfram engine**. ⚠ It is **not**
the physics spec — the physics lives in the **closed** `directives/S11c_a_SHARED_PHYSICS.md` (HEAD
`2926c71c`), already reviewed and locked. This directive adds only the **blindness mandate, the WL tag
grammar for the downstream comparator, run discipline, and one guard**. Your job is to find a defect that
would (a) corrupt the WL engine, (b) break its **blindness** (the design's one blindness control — the
engine must import nothing and re-derive from the spec alone), (c) break the tag join the cross-engine
comparator depends on, or (d) leak an expected value or a sibling-engine result. ⛔ Not to re-review the
spec's physics (closed). It gets two legs before the builder launches (rule 7).

## Read first
- `directives/S11c_a_SHARED_PHYSICS.md` — esp. §4 (T-0…T-i objects), §5 (routes/controls — **§5a route 2 =
  material-coordinate `w'` face-flattening, mapped back to Eulerian; one-sided mutates only the direct
  route**), §6 (three script clauses, no VERDICT), §7 (the `<ENGINE>_S11CA_<QUANTITY>` grammar; "the Wolfram
  engine re-derives the supplied §§1–3 inputs without an import"), §8 (supplied vs computed, builder report).
- `directives/S11b_wl_build_directive.md` — the S11b blind-WL build directive this adapts.
- `directives/S9_export_chain_rebuild_directive.md:16-18` — the blindness control anchor.
- `directives/S11_C17_C18_spec_repair_decisions_v2.md:53-60` — the frozen T7 comparator contract (rejects a
  native boolean as a residual operand).
- The rule-2 twin `directives/_measurements/S11c_a_wl_build_directive.md` (spot-check the anchors).

## What to check
1. **Faithful pointing, ⛔ not restatement/paraphrase.** The directive must POINT at the spec's §§ and the
   inherited anchors, not restate/narrow them. ⚠ The measured defect class is a paraphrase that silently
   drops or weakens a clause. Confirm §5a, §6's script clauses, §7's tag grammar, and the blindness mandate
   are pointed at whole, not re-rendered in a way that changes them.
2. **Blindness is STRUCTURAL, not a prohibition.** The engine must import nothing and re-derive §§1–3;
   blindness is enforced by *absence* (importing nothing), ⛔ never by a do-not-read sentence (rule 12). Does
   the directive correctly mandate no import / no LEDGER / no import tag, and declare the inherited constants
   as the engine's own symbols per §1/§2a? Is there any smuggled denylist?
3. **The §5a guard — correct AND blind-safe.** The directive requires route 2 to be the spec's genuine
   material-coordinate `w'` construction, not a re-expression of route 1. (a) Is this a faithful enforcement
   of spec §5a (⛔ not new physics, ⛔ not a recipe that over-specifies the derivation)? (b) ⭐ CRITICAL:
   does the guard leak any **sibling-engine (SymPy) construction, result, or defect**? A blind engine must
   not learn anything about the other engine's actual computation — flag any sentence that reveals it.
4. **Tag grammar enables the comparator.** Does the directive emit `WL_S11CA_<QUANTITY>` for exactly the
   spec's `S11CA_<QUANTITY>` base names (so the comparator can join the two engines on the shared stem)? Are
   per-case branch/face/DOF/density values kept as a keyed map in one tag's payload, not separate tags? Is
   the `_LOCAL_` convention correct? Are booleans required to be emitted as CAS objects (T7 rejects a native
   boolean)?
5. **No leaked value / no VERDICT / no self-review machinery (rule 5, rule 2, rule 12).** Any expected value,
   acceptance criterion referencing a value, PASS/FAIL/VERDICT, completeness registry the spec's task
   structure doesn't carry, or a check that would audit its own input?
6. **Mathematica idioms + run discipline + acceptance.** ConditionalExpression stripped / poles handled
   symbolically; one kernel at a time (2-seat licence); kill criteria (600 s no-output / 6 GB) by PID;
   demonstration runs under scratch, ⛔ never under `mathematica/out/`; `--sandbox danger-full-access`. Are
   the three acceptance tests executable and value-free (blindness copy-to-scratch byte-identical + no stray
   write; flush observed under redirect; clean `mathematica/out/`)? Is the comparator correctly SEPARATE
   (not built here)?

## Physics filter
Report a finding only if it catches a way the **WL engine, its blindness, or the comparator join** would be
wrong (a paraphrase that narrows the spec, a blindness break, a leaked sibling result, a wrong/unjoinable
tag, a leaked value, added physics, a smuggled denylist). ⛔ Not spec physics (closed), not style. "Nothing
survives" is weak evidence — say what you checked.

## Output
Findings most-serious first (source quote + directive line + the concrete failure). Then answer: (1) is the
blindness mandate structural and complete? (2) is the §5a guard a faithful, blind-safe enforcement of the
spec (no sibling leak)? (3) does the tag grammar enable the cross-engine join? (4) is the directive safe to
hand to the builder?
