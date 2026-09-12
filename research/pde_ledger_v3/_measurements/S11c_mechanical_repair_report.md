# S11c mechanical-load repair — execution record

The user authorized the [repair plan](S11c_mechanical_sign_repair_plan.md) on
2026-09-12. The mechanical-load repair and dependent regeneration are complete.
The pinned baseline is in
[S11c_mechanical_repair_baseline.json](S11c_mechanical_repair_baseline.json).

S11c-b now derives the external-work row multiplier from the supplied action's
stiffness coefficient and the actual stored mechanical row. It applies that
multiplier to the face contribution in both local and expanded rows and in weak
restriction provenance. The separately exported physical generalized force
retains its meaning. No authority was changed.

The native LAB_HELD/RHO4_CONSTANT rebuild computes multiplier `−1`. Its four
action-load comparisons, four non-face mechanical comparisons, mass comparison,
four physical-force comparisons, and chemical-derivative comparison are all
exactly zero against the pinned baseline. The run took 991.02 seconds with
1,643,780 KiB peak RSS and empty stderr.

S11c-c2's power check now uses b's separately computed, constrained per-term
energy variations and a functional time variation of the supplied kinetic
energy. This construction is independent of assembled slab/face routing; it is
not an independent derivation of b's energy basis or material constraint. The
energy-basis import supplies the stiffness orientation anchor. Pressure closure
and the physical traction reference remain native calculations.

The focused c2 run uses the newly computed native b case and the pinned existing
c1 response. It writes no export and makes no claim about the other three cases.
Its scalar-carrier-normalized power residual is exactly zero, as are all four
kinetic-normalization residuals. Separate controls reverse native traction and
the incoming face-work load before closure. Each computed control residual has
three nonzero PIT samples. These fingerprints demonstrate detection at those
samples; they do not prove nonvanishing on every parameter or spectral stratum.
The run took 489.18 seconds with 1,662,352 KiB peak RSS and empty stderr.

The two focused transcripts are under `scripts/out/`, with 15 and 11 unique
keys and sizes 118,120 and 308,291 bytes. The
[focused inventory](S11c_mechanical_repair_focused_checks.json) records source
and payload hashes, literal short residuals, control samples, and complete order
and restored-dimension metadata. An unnecessary post-run expansion of the
nonzero controls was interrupted; the recorded inventories use the already
computed carrier fingerprints and canonical nominal residual. Neither native
calculation was interrupted.

The full four-case b primary run completed in 13,850.05 seconds with
1,942,636 KiB peak RSS and empty stderr. Its 2,441-row export changes exactly
`slab_operator`, `slab_operator_term_origins`, `coupling_kernel`, and
`coupling_kernel_term_origins`. No c1 direct input changed. Source/export pins
and the native execution inventory have no discrepancies. Its ten previously
deferred heavy control families remain deferred. The four-case check completed in 514.07 seconds with 1,722,504 KiB peak RSS
and empty stderr. All 188 physical and 12 source-pin comparison scalars are
exactly zero, including eight separate uniform S11b face-load comparisons.
Whole mass, chemical, kinetic, and constrained energy provenance hashes are
unchanged; the zero-source control is zero in every component. The 202-key
check transcript has 596 metadata leaves with no gaps, unknown dimensions,
nonfinite objects, or duplicate keys. The native b and check transcripts are
published under `scripts/out/`. Three earlier helper attempts stopped on face
labels or legacy-coordinate metadata; their temporary logs are retained. The
completed check binds only the computed S11b face-work symbols to live inputs.
The native producer was run once and was unchanged during these helper fixes.

The c1 refresh completed in 497.28 seconds with 1,652,512 KiB peak RSS and
empty stderr. All 44 exported value serializations are unchanged; source and
export pins are consistent and the native transcript is published.

The full c2 run completed in 3,967.96 seconds with 1,911,220 KiB peak RSS
and empty stderr. Its 70-row export changes the closed slab and closed coupling
roots; the self-energy increment is unchanged in the value census. All 431
publication checks completed (387 `True`, 44 literal zero), and the 35 direct
lookups equal the import list, including the new energy-basis input. The
four canonical power residuals and all 16 kinetic-normalization entries are
exactly zero. Both controls have three finite, nonzero PIT samples in every
case (24 samples total); these are sampled detection results, not a global
nonvanishing theorem. The final 24-key inventory took 88.94 seconds with
619,608 KiB peak RSS and empty stderr. Its metadata has no gaps or unknown
units. The 6,456 infinity occurrences are integration/limit endpoints, with
none elsewhere and none in the numeric fingerprints; this is not a convergence
proof. The earlier inventory was refined to distinguish domain endpoints from
invalid nonfinite values without changing its computed output. Both native c2
and check transcripts are published under `scripts/out/`.

The required scoped dependency triage completed in 1,420.25 seconds with
1,648,232 KiB peak RSS. Its baseline emitted four tags over 2,046 symbols;
four of six collapse trials ran, and two reproduced the documented
`StrictGreaterThan`/`StrictLessThan` subtraction error. The three physical roots
were classified derived; the invariant dimension-binding tag leaves the literal
verdict `DERIVED_OR_DECLARED: FAIL`. No premises sidecar was added, and the
c2 export remained unchanged. This is the same scoped triage boundary recorded
by the prior build, not a full clearance.

The fresh d reference completed in 297.26 seconds with 1,714,256 KiB peak
RSS. All 39 reconstruction scalar digests match literal zero, and dimensions
are resolved. The original current check then completed in 180.97 seconds
with 122,636 KiB peak RSS and its post-emission zero-residual guard enabled.
All five mechanical-load coefficients now agree; the five mass coefficients,
remaining energy/current residuals, and independent stiffness anchor also
agree. All 30 recorded residual scalars are zero, with clean metadata. The
reference and current transcripts are published; the previous annex payload
is intact. The existing d formulas were unchanged.

The final four-case d rebuild completed in 9,989.25 seconds with
1,734,840 KiB peak RSS, stable source pins, exit zero, and empty stderr. All nine
serial recheck stages completed. The [full inventory](S11c_mechanical_repair_d_full_checks.json)
and [run records](S11c_mechanical_repair_d_rechecks.json) preserve their inputs,
residuals, domains, and source hashes. No historical root count was used as a
target.

The 24 native spectrum packets each contain nine isolated radical roots of a
degree-11 polynomial and both normal-momentum lifts: 432 candidates in total,
with 336 one-dimensional and 96 two-dimensional nullspaces. The computed full
basis and normal-derivative pairing ranks satisfy the emitted regularity checks
at these points. The separate rectangular-jet family contains 528 records
(336 with nullity one, 192 with nullity two). Its coverage inventory has no gaps.
Joint-sheet records retain 480 transported paths and 24 branch-locus paths with
an unresolved continuation status. These are point/path results, not global
parameter-variety or sheet coverage.

The constant-end resolvent census has 432 Laurent residues and 1,728 contour
records. The known near-pole double-precision inverse-jump cancellation remains
visible: `2.3508087365723638` at `[3,2,-1]`; the corresponding 60-digit identity
residual maximum is `4.0192237586793586e-52`. The threshold family computes 24
right and 24 left chain spaces, 96 individual chains, 144 normal paths and
336 bulk paths. All emitted plane-equation residuals are zero. Twelve threshold
points match the reduced source sheet and twelve lie on the opposite sheet;
the source-root comparison is not silently treated as zero on that other sheet.
Earlier locus-census unresolved labels remain as stage-local records preceding
these dedicated chain calculations. Zero-frequency, singular-denominator,
remaining mixed/defective domains, and global coverage remain open. These
normal-momentum data do not supply the profile-frequency bound poles of §3b.

Fourier reconstruction retains zero literal residuals and zero numeric
fingerprint projections in every case, with complete carrier/slot metadata.
The full output contains 229,636 unique tags and no inventoried metadata gaps,
unresolved dimensions, nonfinite objects, or duplicate tags. Lossless expansion
and re-encoding preserve every payload and all 229,634 indexed source-line
assignments. The published output is 83,848,737 bytes, SHA-256
`312935010ef56ce0eabe7b1685e906dfb785c40e02978c7b93fdbb4af11bae52`.
The [publication record](S11c_mechanical_repair_d_publications.json) verifies
that the previous annex payload remains unchanged.

The b/c1/c2 exports and their production transcripts, the focused diagnostics,
the repaired current transcript, and the complete four-case d transcript are
now regenerated and published under their expected paths. The original
five-coefficient mismatch is resolved. The remaining S11c-d program resumes
with the reduced nonlocal current and mode/flux normalization on one case,
then the complete two-ended matching/scattering construction. All ten broad
engine TODOs remain explicit; no incomplete d export was written.

The retained d solver contract, S10/Lean work, and authority files are unchanged.
Fourteen changed/new repair Python files compile, and the code/report diff
check is clean. No commit, review leg, comparator, Wolfram run, or downstream
stage occurred. The `.out` files are ready for DataLad/git-annex at the next
requested checkpoint; code and reports remain for Git.
