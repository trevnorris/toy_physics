# v3 PDE ledger Lean formalizations

**Required before Lean work:** read and follow
[FORMALIZATION_POLICY.md](FORMALIZATION_POLICY.md), as directed by
[AGENTS.md](AGENTS.md). The deliverable is proof, compact object identification,
an exhaustive coverage contract and mutation controls, with fidelity review.
The earlier plan to extend the full S10 CAS bridge is superseded.

One pinned Lean environment serves the step-specific source directories:

- [s9/](s9/README.md): the original D=3 action, integrated variation, and mode census.
- [s10/](s10/README.md): the arbitrary-dimensional baseline, all five controls,
  expression-tree dimensions, complete basis constructions and Levi-Civita comparisons.
- [s11/](s11/README.md): the completed bounded homogeneous, invariant, bulk and
  analytic contracts, plus the separately identified work still in progress.

For a fresh checkout, start with [INSTALL.md](INSTALL.md). It covers elan/Lean,
the pinned Python environment, dependency caches, and the portable proof/control
runner. The original author verification scripts retain historical workspace
guards; the portable runner creates new evidence without changing those records.

See CHECKPOINT.md (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/CHECKPOINT.md`) for the original S9/S10 checkpoint and
s10/CAS_CHECKPOINT.md (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/CAS_CHECKPOINT.md`) for the subsequent CAS bridge
checkpoint, validation evidence and resume plan.
The later [S10 contract completion](s10/FIDELITY_REVIEW.md) records the compact
coverage contract, both independent CLEAR fidelity reviews, mutation evidence
and the stopping point for S10 Lean work. Separate CAS production and paper
obligations remain in [COVERAGE.md](s10/COVERAGE.md).
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
bash setup.sh
export PATH="${ELAN_HOME:-$HOME/.elan}/bin:$PATH"
.venv/bin/python verify.py --doctor
.venv/bin/python verify.py all
```

`verify.py all` freshly builds the local proofs, checks the recorded controls
and runs compact native checks for all **completed** contracts, including S11.
It uses one worker and separate objects/logs under `_scratch/lean_portable/`.
Use `--list`, `--plan`, or a contract name such as `analytic-error` to select a
smaller run. The completed D5 density classification is included as `d5`;
its fresh portable replay passed, as recorded in
INSTALL_D5_VALIDATION.json (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/INSTALL_D5_VALIDATION.json`). The new D5 bulk
contract is complete with its own recorded suite and two CLEAR fidelity
reviews; it is registered as `d5-bulk`. Registration checks are recorded in
`INSTALL_D5_BULK_REGISTRATION.json`; no new full D5 bulk portable replay is
claimed. Observe the host resource safeguards in INSTALL.md. See
[s11/D5_BULK_FIDELITY_REVIEW.md](s11/D5_BULK_FIDELITY_REVIEW.md).

For ordinary shared-cache builds after installation:

```sh
LAKE_CACHE_DIR=.lake/cache lake build
lake build S9Pilot
lake build S10Pilot
lake build S10Controls
lake build S10Anisotropic
LAKE_CACHE_DIR=.lake/cache lake build S10Audit
lake build S11D4Bulk
```

Bare `lake build` checks the default **S9/S10** targets; it does not build S11.
The other commands select one target. They treat warnings as errors but do not
run the mutation/native suites. To print the root theorem audits directly:

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

`bash setup.sh` installs elan if needed, the project-selected toolchain and a
pinned Python virtual environment, fetches cached library artifacts and builds
the required Physlib imports. It requires the committed dependency manifest
and does not change elan's global default. Run `verify.py` afterwards for proofs
and controls. Setup needs network access; see [INSTALL.md](INSTALL.md) for system
packages, exact commands, resource limits and troubleshooting.

The proof reports state their assumptions and exclusions. These formalizations
use Mathlib directly and Physlib's dimensional algebra and Levi-Civita symbol;
its wave theorem is also used in the original installation check. The CAS
bridge report identifies exactly which emitted expressions are certified.
The CAS engines, exports, and ledger prose remain separate artifacts; their
full certification is not implied by a Lean build.
