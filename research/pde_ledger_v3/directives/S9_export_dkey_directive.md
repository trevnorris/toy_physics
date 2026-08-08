# S9 export re-key — computed production dimension in the ledger key

**User decision, 2026-08-08.** Read `CLAUDE.md` and
`research/pde_ledger_v3/S9_REWRITE_PLAN.md` first. The rewrite plan is CLOSED. This directive applies only
the export-key repair below. If an exported object does not fall cleanly on one side of the stated
production boundary, name it and stop; do not choose a key for it.

Repository root: `/var/projects/toy_physics`. Ledger root: `research/pde_ledger_v3`.

## Decision

A ledger key records the spatial component count at which its object was computed. A fixed-component
product has key `<name>_d<component-count>`; a product solved symbolically in general `D` has the
unsuffixed key `<name>`. Exported knobs and supplied dimension premises that are shared inputs rather than
products of a selected-component computation retain their unsuffixed keys.

The component count is read from `spatial_coordinates`, the same structure used by the fixed-component
action and ansatz. Every fixed-product collection reads that structure when it records its production
tag. A fixed specialization of a symbolic solve uses that same readout both for the substitution and for
the production tag. No component-count literal is permitted between the computation and the key.

The formatter accepts the production tag it is given. It first removes a legacy component-looking tail
from the internal producer name, then derives the complete ledger suffix from the tag alone. Thus an
internal name cannot supply, preserve, or override the component count in the ledger key. A general-`D`
tag yields no suffix; a fixed tag yields exactly the suffix formed from that tag.

## The dimension-solve boundary

Classification is a fact about construction, never expression inspection. Do not inspect `free_symbols`,
coordinates contained in a payload, expression shape, or payload text to classify an output.

Inside `dimension_block`, the cut is between the solve's inputs and outputs:

- The action terms, their field multiorders and dimensions, and the linear system are solve inputs. They
  are constructed from the fixed-component action and carry the component count read from
  `spatial_coordinates`.
- The solution and coefficient-dimension products are solve outputs. They are symbolic in general `D`
  and carry the general-`D` tag.
- Fixed specializations evaluate those symbolic outputs at the count read from `spatial_coordinates` and
  carry that same fixed tag.

The supplied energy-density, field, wavevector-norm, and squared-velocity dimension premises are not
products of the action or the solve. They are constructed together in one supplied-reference group and
receive its general/shared production tag together. They use the same metadata path and key formatter as
the other exported objects; none bypasses classification or receives a tag by name after construction.

This boundary is expressed by constructing the input and output collections separately and joining them;
it is not a post-hoc list of names. The ansatz/spectrum, directional, zero-wavevector, and standard-name
product paths are fixed-component paths. Production metadata is preserved through standard-name
substitution.

If any non-internal `MAIN` output reaches the writer without production metadata, the run must stop and
name that output. An object that cannot be assigned by the construction boundary above is a finding about
the export boundary and is not to be resolved in this round.

## Emitted partition

The writer emits the production-side D-partition from the recorded tags and the key-side partition read
from the generated key suffixes, after the existing export tallies. It emits those two operands and stops:
there is no partition residual, assertion, or replacement guard. The operands are a visible computed
object containing the export keys for a reader and an independent review leg to check. They are not an
in-engine audit of the classification, because the recorded classification is their own input.

The exact-`srepr` reconstruction check continues to look up every live record by its generated key and to
emit both counts and its residual before guarding. The entry population and every record payload field
remain unchanged; only keys move.

## Files and invariants

Change only:

- `scripts/S9_light_requires_shear_sympy_audit.py`, for production metadata, key formatting, partition
  emission, and removal of the partition residual and guard;
- `scripts/S9_exports.py`, only by running the engine;
- `scripts/out/S9_light_requires_shear_sympy_audit.out`, only as literal stdout from that run;
- this directive, as the additional authored deliverable.

Do not change the derivation, action, ansatz, assumptions, computed values, or Wolfram engine. Do not add
a registry, YAML, runner, test framework, or repository file. Do not address the open S9 naming and
dimension items reserved outside this round.

Verification must include a control rebuild at the unchanged spatial structure, followed by a rebuild
with only that structure's component count changed. The control must reproduce the generated export; the
changed structure must move the fixed keys to the suffix computed from the changed structure while the
general-`D` keys remain unsuffixed. For every supplied reference, also mutate the shared group's tag and
emit the resulting key. Build and run the specified weaker implementations: a partial supplied-reference
group, deletion of the named override without the group repair, and partition operands derived from a
single source. Save every ablation script and its literal stdout at a named absolute path outside the
repository.

Re-run S9 and capture literal stdout. Compare it with the committed
`scripts/out/S9_light_requires_shear_sympy_audit.out`, and report the complete diff.
