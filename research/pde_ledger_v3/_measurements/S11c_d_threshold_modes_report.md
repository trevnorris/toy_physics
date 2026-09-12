# S11c-d generalized threshold modes and local sheet connections

The finite-threshold construction extends the existing engine and follows the
[plan](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_modes_plan.md).
It retains the three-parent import fold, direct-lookup manifest, Fourier reduction,
previous spectrum/resolvent records and independent bulk-exception family.

[NormalTaylorChains.construct](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3471)
computes full block-Taylor kernels and their ranks in an exact algebraic field.
Kernel-nullity increments determine chain lengths; quotienting the leading
spaces of longer chains selects independent chains. Both left and right germs
are computed. A Taylor cap without stabilization is explicitly unresolved.
[at_point](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3534)
derives the Taylor matrices along the computed radical relation.
[plane_wave](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3690)
constructs polynomial-times-plane-wave modes and independently applies the
matrix differential series to their coefficients. Block systems use the stated
unit frame; restored row/column coefficient units and physical chain units
accompany their multigrade metadata.

[ThresholdModeAudit.construct](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3569)
uses the enumerated frequency loci, with separate zero-frequency records and
finite-denominator/nonzero-radical guards. The
[connection](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3719)
joins each radical lift to the actual reduced Fourier seed and matrix, evaluates
independent determinant-factor valuations, and derives the local normal-square
frequency unfolding by elimination. The opposite radical sheet remains a
separate record; its nonzero source-root difference is not a failed residual.

[local_chart](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3784)
evaluates both normal lifts above and below threshold at two radii. Full left
and right SVD spaces, gauge ranks, equation residuals and secants are retained.
Both starting normal signs traverse upper/lower semicircles and a loop, with
full mode checks at every one of seventeen nodes. Existing joint-root transport
checks normal and bulk-radical continuation separately. The bulk segments are
affine between mode nodes, as explicitly recorded. SVD arithmetic has 53-bit
mantissa and relative rank tolerance 1e-10; selected approach matrices are also
evaluated at 40/60 digits. Secant differences near roundoff do not establish
monotonic convergence as the radius shrinks.

[exception_disk](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3940)
checks all fifteen already-derived end-exception polynomials using exact
rational Taylor dominance on a complex-frequency disk. Target multiplicities
come from division by the target minimal polynomial. The emitted exact sign,
bound digest and reconstruction residual accompany displayed bounds. Paths
lie within half the accepted disk radius. This excludes additional enumerated
exceptional frequencies locally; it is not a global parameter/sheet atlas.

The full four-case producer completed in **10,051.411 seconds**, with
**1,733,148 KiB** peak child RSS, exit zero, empty stderr and all **27 source/input
pins unchanged**. The final engine SHA is
`5612475906982732e71b899190efe6bf833108041e9d694cf453b820efda7c6b`.
The three focused LAB_HELD/RHO4_CONSTANT runs use the same source:

| End | Seconds | Peak RSS (KiB) | Published bytes |
| --- | ---: | ---: | ---: |
| [Reference](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_threshold_modes_reference.out) | 131.815 | 173,212 | 3,706,891 |
| [Left](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_threshold_modes_left.out) | 133.750 | 174,272 | 3,668,304 |
| [Right](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_threshold_modes_right.out) | 185.115 | 175,132 | 3,662,063 |

The [full census](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_full_summary.json)
contains **24 finite threshold points**, **48 left/right chain spaces** and
**96 individual chains**. Every chain space has kernel counts and full-basis
ranks `[2,4,4]`, leading-space rank 2 and chain lengths `[2,2]`. All 7,320 scalar
block-kernel, chain, polynomial-mode, radical-Taylor and determinant-multiplicity
residuals are exactly zero. All 24 local connections and all 24 fifteen-condition
exception disks are defined. Twelve joins match the reduced real-axis seed;
twelve retain the opposite sheet. Twelve zero-frequency intersections remain
individually unresolved.

The 192 approach points and 2,448 path nodes all have rank 3 and full nullity 2
under the stated SVD tolerance. All 144 normal paths and 336 bulk paths are
defined. Across 48 normal loops, the maximum absolute end-plus-seed residual
is `2.082e-17` at `[-1,0,0]`. Selected mode-equation residual maxima are
`2.949e-15` at `[0,0,0]` and `1.536e-15` at `[-1,-2,1]`; a secant residual maximum
is `3.839e-13` at `[2,0,0]`. The
[dimension-separated inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_full_modes_inventory.json)
retains every residual category and approach radius. These are finite-precision
observations on the recorded slice, not a convergence theorem.

The [published main transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
is **83,919,481 bytes** with **229,636 unique tags**. The threshold family adds
52,416 computed objects, including 14,508 fingerprints. Its inventory found no
duplicate tags, missing metadata, nonfinite objects or unmatched objects.
Solved dimensional residuals and Fourier carrier/action reconstruction residuals
remain zero. The
[preservation inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_full_payload_preservation.json)
accounts for all 124,804 previous tags: 124,408 payloads are byte-identical;
396 changes are classified as run provenance, mapping/metadata ordering or
injective Dummy renaming. There are no unclassified changes, missing tags or
changed existing native values. The prior 432 modes, 432 residues, 1,728 contours
and 504 joint-sheet paths remain, including their unresolved domains. The
[codec check](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_full_codec_inventory.json)
restores every payload and all 229,634 source-index assignments. Its expanded
314,180,252-byte transcript is retained only in scratch space.

The main and three focused outputs were published atomically. The prior
189,492,142-byte annex payload remains unchanged. The
[run provenance](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_modes_runs.json)
records producer commands, source/input hashes, verification commands, artifact
hashes, resource measurements and development attempts. Those attempts include
replacing expensive whole-expression algebraic simplification with exact
coefficient-tree arithmetic, fixing empty sparse rows, and successive additions
of path-node coverage, unresolved-cap handling and exact disk certification.
The retained solver/export contract remains byte-identical.


The synthetic [chain controls](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_threshold_chain_controls.json)
produce a simple chain of length 1, a Jordan chain of length 2, and mixed lengths
3 and 2 with determinant valuation 5. Their exact block/chain residuals are
zero. Truncating the mixed example below the required Taylor order yields an
explicit unresolved status. These are dimensionless algorithm controls, not
completion of the physics controls in section 5.

Zero-frequency, singular-denominator, remaining mixed/defective domains,
global sheet coverage and physical current channel labels remain unresolved.
The earlier End/Bulk stages retain their original scope records; the new
threshold family provides the finite-point generalized-chain follow-up.
Constant-end normal-momentum poles remain distinct from the profile-dependent
frequency poles required by section 3b. No profile-frequency pole search or
bound-state existence claim is made.

All ten broad TODO categories remain, and `scripts/S11c_d_exports.py` remains
absent. The next construction is the reduced nonlocal S11b current and mode/flux
normalization, retaining explicit mode and sheet domain boundaries. Section 1
premises and the c2 operand/sign and shear-normalization debts remain supplied.
The two named generic triage helpers (`reduction/derived_or_declared.py` and
`reduction/engine_output_checks.py`) were unavailable after searches including
hidden/ignored repository files, Codex skills and Claude configuration; they
were not run. Dedicated dimension, coverage, preservation and codec inventories
are recorded instead. No new premise, upstream repair or change of work
breakdown was required. No S10/Lean changes, review, comparator, Wolfram,
downstream run, export, commit or push was performed in this checkpoint.
