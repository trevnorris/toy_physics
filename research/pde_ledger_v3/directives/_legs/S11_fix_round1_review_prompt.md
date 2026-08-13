# Independent review — a DECISION LIST, before any builder touches the engine

Read-only. ⛔ Change no file. This is `CLAUDE.md` rule 7's trigger: no builder launches until its decision
list has had two legs. The list is the one artifact the builder trusts, and it is otherwise checked zero
times.

## Artifact

`/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11_engine_fix_round1_brief.md`

Its evidence: `research/pde_ledger_v3/directives/_measurements/S11_engine_fix_round1_brief.md`
(regenerate with `docs/_measurements/gen_s11_fix_round1.sh`).

## Read the authorities FIRST, then the list

`CLAUDE.md` (rules 2, 3, 5, 7, 10, 11, 14). `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`.
`research/pde_ledger_v3/directives/S11_sympy_build_directive.md` (closed).
The engine: `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py`.

Reading the list first anchors you to its framing, which is the thing under test.

## Context

Codex built the engine from the closed directive. It compiles and emits genuine `srepr` objects. Four runs
have failed to publish `S11_exports.py`; `MAIN` D=5 raises `MemoryError` inside `sympy`'s row reduction,
is swallowed by a bare handler, and never enters `completed_pairs`, so the publish guard can never pass.
An isolated one-cell reproduction is at `/home/trevnorris/.s11_build/repro_d5.py` (~4 min).

## What to check

1. **IS ITEM 1'S DIAGNOSIS ACTUALLY RIGHT?** The list blames the `iszerofunc`/`simplify=False` pair at
   `:870` and `:920`. Test that. Is the structural zero test the *cause*, or a *symptom* of something
   upstream — a matrix that is already swollen when it arrives, roots carrying unreduced radicals, a
   coefficient domain choice? **Run the reproduction yourself** and say what you measured. If the cause is
   upstream, the list sends a builder to fix the wrong line.
2. **DOES ITEM 1 NAME AN OBJECT OR SPECIFY A RECIPE?** Rule 3. A "fix" that computes a *cheaper but
   different* object would look like success and corrupt the physics. Say whether the property as stated is
   implementable without inventing the routine, and whether it is tight enough to exclude a wrong object.
3. **ITEM 3.** `run_pairs_payload`/`skipped_pairs_payload` are built after the loop from `completed_pairs`
   and are consumed by `merged_export`. Can "a completed `MAIN` survives an interruption" and "no row
   claims a run record it does not have" both hold? If not, is the list right to say report it under §10
   and leave the placement alone?
4. **THE BARE HANDLERS** at `:756 :775 :955 :1317 :1450 :1668`. Which mask a real failure? Is asking the
   builder to *report* rather than *remove* them the right call?
5. **THE D2/D3/D4 CLAUSE.** The list says a changed value there after the fix is a finding, not a
   regression. Is that right, and is it stated strongly enough that a builder will not tune to reproduce
   prior output?
6. **WHAT IS MISSING OR WRONG**, including the premise. If this is not fixable in place — if the honest
   outcome is a §10 unavailable construction, or a spec question rather than an engine question — say so.

## Method

- **Every claim carries the command that produced it and its literal output.** A claim with no computation
  behind it is discarded (rule 2). This binds reviewers exactly as it binds the engine.
- ⛔ Do not propose a new registry, checklist or mechanism. Three of the four items should make the engine
  simpler or louder, not larger.
- ⛔ Do not edit any file. ⛔ Do not run the full engine — use the one-cell reproduction; the full run takes
  hours and has OOM-killed this machine.
- Report a finding only if it bears on whether the physics could come out wrong, or a wrong result survive.

## Deliverable

A verdict per numbered item, every additional defect with file:line and the command that found it, and a
plain statement of anything the list gets wrong.
