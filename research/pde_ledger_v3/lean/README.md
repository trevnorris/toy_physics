# v3 PDE ledger Lean formalizations

One pinned Lean environment serves the step-specific source directories:

- [s9/](s9/README.md): the original D=3 action, integrated variation, and mode census.
- [s10/](s10/README.md): the arbitrary-dimensional baseline, all five controls,
  expression-tree dimensions, complete basis constructions and Levi-Civita comparisons.

See [CHECKPOINT.md](CHECKPOINT.md) for the original S9/S10 checkpoint and
[s10/CAS_CHECKPOINT.md](s10/CAS_CHECKPOINT.md) for the subsequent CAS bridge
checkpoint, validation evidence and resume plan.
The subsequent [S10 CAS bridge](s10/CAS_BRIDGE_RESULT.md) connects 916 scalar
expressions from the two anisotropic D3 transcripts to checked Lean objects,
including complete minor families. Its [locus extension](s10/MINOR_LOCUS_RESULT.md)
also certifies 12 exceptional predicates and four targeted points. The
[rerun extension](s10/EXCEPTIONAL_RERUN_RESULT.md) connects the exceptional
matrices and complete bases, with an explicit physical scale convention. The
[count extension](s10/COUNT_RESULT.md) certifies all 112 generic and exceptional
rank, nullity and count records against those matrices and bases, with
explicit chart assumptions on the generic counts.
The [root-list extension](s10/ROOT_RESULT.md) certifies solution and distinct-root
lists, algebraic multiplicities, and syntactic filter counts, including the
parallel double root. The [coincidence extension](s10/COINCIDENCE_RESULT.md)
connects the primary emitted equations, full loci, allowed regions, decisions
and witnesses on the positive-coefficient domain. The
[metadata extension](s10/METADATA_RESULT.md) checks aggregate/Q8 copies, root
signs, solve operands and statuses, and the retained/skipped stratum records.

`lakefile.toml`, `lake-manifest.json`, `lean-toolchain`, `setup.sh`, and the ignored
`.lake/` dependency/build cache live here. Each Lean library declares its source
directory explicitly. This avoids duplicating the dependency installation while
keeping each step's proofs and reports under its own directory.

## Build

Run from `research/pde_ledger_v3/lean/`:

```sh
LAKE_CACHE_DIR=.lake/cache lake build
lake build S9Pilot
lake build S10Pilot
lake build S10Controls
lake build S10Anisotropic
LAKE_CACHE_DIR=.lake/cache lake build S10Audit
```

The first command checks all libraries. The others select one target. Their
modules treat warnings as errors. To print the root theorem audits directly:

```sh
lake env lean s9/S9Pilot.lean
lake env lean s10/S10Pilot.lean
lake env lean s10/S10Controls.lean
lake env lean s10/S10Anisotropic.lean
lake env lean s10/S10Audit.lean
```

## Pinned environment

- Lean 4.33.0, commit `d8b18978322de05a8f3dba51ef03cf5461676c17`.
- Lake 5.0.0-src+d8b1897, bundled with Lean.
- Physlib commit `8b2b23701c07409f882c49e0bd290a214e44450c`.
- Mathlib v4.33.0, commit `db584cd6d46c92f209a44c0f1c829460d327499d`.

`bash setup.sh` installs the project-selected toolchain, fetches cached library
artifacts, and builds the project. It uses the locked dependencies when the
manifest exists and does not change elan's global default. Setup needs network
access and writes to elan's installation and library cache locations. Git, curl,
and zstd are needed; VS Code with the Lean extension is optional for editing.

The proof reports state their assumptions and exclusions. These formalizations
use Mathlib directly and Physlib's dimensional algebra and Levi-Civita symbol;
its wave theorem is also used in the original installation check. The CAS
bridge report identifies exactly which emitted expressions are certified.
The CAS engines, exports, and ledger prose remain separate artifacts; their
full certification is not implied by a Lean build.
