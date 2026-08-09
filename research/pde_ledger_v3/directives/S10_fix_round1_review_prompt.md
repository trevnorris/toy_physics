# Independent physics review — S10 export chain, fix round

## Artifact

Uncommitted working-tree change in `/var/projects/toy_physics`, branch `ledger-v3-rebuild`:

```
git diff -- research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py
git diff -- research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out
```

plus untracked `research/pde_ledger_v3/scripts/S10_exports.py` (**generated**) and
`research/pde_ledger_v3/directives/S10_export_chain_directive.md`.

Codex, two passes. You are one of two independent legs with this identical prompt. **Earlier rounds exist;
do not look for them.**

## What this is

S10 imports S9's generated flat `LEDGER`, binds its objects rather than re-declaring them, adds its own
derivations across several component counts, overwrites what it re-derives, and writes the merged
`S10_exports.py` for the next step. A previous review found the physics correct and six defects in the
export chain; this round repairs four and records three.

## What to check — by mutation, not by reading

1. **Is the generated module a deterministic function of the run?** Nothing on the export path may depend
   on set or hash iteration order. Run under several hash seeds and compare bytes. Then look for any
   *remaining* unordered iteration the seeds happened not to expose.
2. **Is every symbol one object?** A same-named symbol with different assumptions in two places produces a
   difference that **prints as though it vanished**. Scan the whole merged ledger yourself — ⛔ do not
   trust a scan the artifact performs on itself. Check substitution actually bites: substituting a
   reported solve variable into the object solved must change it.
3. **Does a failed run leave a generated module behind?** Make a guard fail and look on disk. Also check a
   stale module from a *previous* successful run cannot survive a later failure.
4. **Does reconstruction preserve what was stored?** An unevaluated numeric relational must survive the
   round trip with its operands and its type. Demonstrate on a construction that actually produces an
   allowed stratum. ⚠ Check the fix is general, not special-cased to the one relational found before.
5. **Did the derivation move?** It is claimed the stdout is byte-identical to the pre-fix artifact apart
   from the export block. Verify against the working tree and say what moved.
6. **Did this round introduce anything new?** It is a repair of a repair on the same file. A round that
   breeds a fresh defect is the specific failure this question exists to find, and it has happened on this
   project before.
7. **Are the three recorded items honestly stated in the directive?** One deleted check the run no longer
   has, one redundant exported record, and a movement inventory naming two moved values. ⚠ A record that
   overstates what remains policed is a finding.

## Required method

⛔⛔ **A FORM ABLATION IS MANDATORY** — change the *structure* of a load-bearing object, not a coefficient.
Run a control at the unchanged structure first and show your harness reproduces the artifact.

⛔⛔ **A test that passes on a weaker fix is not a test.** For each repair, construct the weakest change it
should reject and show whether it does.

⚠ **Watch for the recurring defect:** a residual whose two operands descend from one source is zero by
construction. Four such checks have been built and deleted on this project.

Report with line numbers: anything typed where it could be read; a check whose expected value lives inside
the artifact it checks; an `assert` preceding the value it guards; a redundant emitted object.

### Write your own script first

⛔ Write your own verification script **before** opening the artifact; save it and its literal stdout to
named absolute paths and report them. **Without them your claims will be discarded.**

## Do not read

Anything under `directives/` beginning `S10_fix_round1_decisions`, `S10_export_chain_decisions`,
`S10_export_chain_review_prompt`, `S10_decisions_review_prompt`, or `S9_export_` / `S9_structural_`.
Anything under `/tmp/s10_fix_round1_evidence/`, `/tmp/s10_export_chain_evidence/`, `/tmp/review_leg_s10/`
or `/tmp/s10_export_chain_review_leg/` — that is other people's evidence. ⭐ Produce your own.

## Operational constraints

- Ablate **copies** under `/tmp`. ⛔ Never modify the working tree.
- ⚠ **This engine is slow.** Wrap every run in `timeout 1800`; a hit is a failed ablation — report it and
  move on. ⛔ Never raise it further. Prefer a few decisive ablations over many shallow ones.
- ⛔ Do not start `wolframscript`.
- Save every ablation script and its literal stdout to named absolute paths and report them.

## Physics filter

Report a finding only if it catches a way the physics, or the ledger's record of it, could be wrong.

## Not findings

Carried limits, out of scope unless this round altered one: the field-dimension premise is unfalsifiable
within this engine now that the registry comparison is gone; the third overwritten record is the difference
of the other two and adds no falsification power; S9's assumption set hand-types its predicates;
`field_dimension` is an alias of `length_dimension`; cross-engine (py↔wl) naming and the comparator are a
separate pass and the Wolfram engine is untouched.
