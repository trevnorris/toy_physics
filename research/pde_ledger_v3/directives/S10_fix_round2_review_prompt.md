# Independent physics review — S9/S10 export chain, fix round 2

## Artifact

Uncommitted working-tree change in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
git diff -- research/pde_ledger_v3/scripts/
git status --short research/pde_ledger_v3/
```

Both engines changed (`S9_light_requires_shear_sympy_audit.py`,
`S10_brane_mode_spectrum_sympy_audit.py`), possibly `extract_knob_inventory.py`, both **generated** modules
`scripts/S9_exports.py` and `scripts/S10_exports.py`, both regenerated transcripts under `scripts/out/`,
and a new directive under `directives/`.

Codex, one build pass. **You are one of two independent legs with this identical prompt. Earlier rounds
exist; do not look for them.**

## What this is

Each step's SymPy engine imports the previous step's flat `LEDGER`, binds its objects rather than
re-declaring them, adds its own derivations across several component counts, overwrites what it
re-derives, and writes the merged module the next step imports. This round repairs five properties of that
export. The physics was independently reproduced by an earlier leg and is **not** what this round touches.

## What to check — by mutation, not by reading

1. **Can a consumer bind every symbol it needs?** The stated aim is that a downstream step never
   re-declares a symbol that appears in an exported value. Take the position of that consumer: import the
   generated module and try to substitute into the objects it carries. ⚠ A same-named symbol with
   different assumptions substitutes **silently and changes nothing** — check that substitution actually
   bites, and do not accept a name match as an object match.
2. **Is the distinction between a bindable object and an authored field name structural?** It is claimed
   it is not a hand-maintained list. Test that claim: add a new emitted object and see whether it is
   classified without anyone editing a list.
3. **Can an imported entry's provenance change without the run saying so?** Construct a collision where a
   downstream record lands on an upstream key with an equal value and a different class or origin, and
   report whether the run rejects it.
4. **Does the symbol-identity scan actually reject?** It is claimed the written module is checked so that
   each symbol name maps to one object. Introduce a same-named symbol with different assumptions and show
   whether the run fails. ⛔ A passing run on unmodified input demonstrates nothing.
5. **Is every exported value traversable?** Walk every stored value with the ordinary structural
   accessors. Report anything that raises, and anything that stores a non-CAS container inside a CAS
   object.
6. **Is the same object exported twice under two keys?** A downstream step differencing two copies of one
   object gets zero by construction. Report any pair that would do that.
7. **Is the generated module a deterministic function of the run?** Run under several hash seeds and
   compare bytes.
8. **Did the derivation move?** Report what changed in each transcript beyond the export block, and say
   whether any computed physics value moved.
9. **Did this round introduce anything new?** It is a repair on files that have been repaired before. A
   round that breeds a fresh defect is the specific failure this question exists to find, and it has
   happened on this project.

## Required method

⛔⛔ **A FORM ABLATION IS MANDATORY** — change the *structure* of a load-bearing object, not a coefficient.
A coefficient rescale tests arithmetic; only a form change tests physics. Run a control at the unchanged
structure first and show your harness reproduces the artifact.

⛔⛔ **A test that passes on a weaker fix is not a test.** For each of the five repairs, construct the
weakest change it should reject and show whether it does.

⚠ **Watch for the recurring defect:** a residual whose two operands descend from one source is zero by
construction. **Four** such checks have been built and deleted on this project; one of them passed while
the change it policed was reverted entirely.

Report with line numbers: anything typed where it could be read; a check whose expected value lives inside
the artifact it checks; an `assert` preceding the value it guards; a redundant emitted object.

### Write your own script first

⛔ Write your own verification script **before** opening the artifact; save it and its literal stdout to
named absolute paths and report them. **Without them your claims will be discarded.** A prose
re-derivation is the same defect class this rebuild exists to remove, relocated into the review.

## Do not read

Anything under `directives/` beginning `S10_fix_round2_decisions`, `S10_fix_round1_`,
`S10_export_chain_`, `S10_comparator_decisions`, or `S9_export_` / `S9_structural_`.
Anything under `/tmp/claude-1000/` — that is other people's evidence. ⭐ Produce your own.

## Operational constraints

- Ablate **copies** under `/tmp`. ⛔ Never modify the working tree.
- ⚠ **The S10 engine is slow.** Wrap every run in `timeout 1800`; a hit is a failed ablation — report it
  and move on. ⛔ Never raise it further. Prefer a few decisive ablations over many shallow ones.
- ⛔ Do not start `wolframscript`.
- Save every ablation script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings

Carried limits, out of scope unless this round altered one: the field-dimension premise is unfalsifiable
within this engine now that the registry comparison is gone; the third overwritten record is the
difference of the other two and adds no falsification power; S9's assumption set hand-types its
predicates; `field_dimension` is an alias of `length_dimension`; the relational round-trip repair is inert
because no allowed stratum occurs on the export path; cross-engine (py↔wl) naming and the comparator are a
separate pass and the Wolfram engines are untouched.
