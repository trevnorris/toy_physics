# S11c-d wider-box three-momentum refinement

The complete fixed-box refinement is validated: four workers and the supervisor
exited zero with empty stderr after 158,602.48 seconds (44h03m). The run evaluated
1,048,960,640 new momentum nodes across twelve grids; all sixteen full-action
records retain every 40 single, 30 pair and 10 triple row. Single/pair values
remain explicitly held at the independently checked Gauss rules.

Raw integral changes from outer144 to216, innermost24 to32 and middle24 to32
reach 3.67575e-10, 3.91669e-18 and 1.55497e-14 respectively. Corresponding maximum
native-term changes are 3.67575e-15, 3.91619e-23 and 1.55497e-19; complete-action
changes are 3.06712e-15, zero in saved arithmetic and 1.69407e-21. The refined
cutoff3-to4 complete-action change remains 1.13467310948e-7. The zero rounded
complete-action difference does not erase the nonzero raw integral/term changes.
These are finite samples in the declared unit frame, not uniform error bounds.

Saved-operand acceptance checks 87 current/frozen sources, all 64,048 worker
artifacts including 64,016 partial sums, complete original row/source/profile/
field/ordered-limit identities, all actual measure controls and finite masses,
24 exact read-only caches and 5,040 zero held-term scalars. Native emission replay
passed 19,762 tags, 9,879 write keys and 134,254 metadata paths. Checks/stdout and
all pre/post packet hashes match; no numerical integration was repeated for
acceptance. Peak worker RSS was 133,848 KiB; coordinator peak was 316,972 KiB.

The canonical 5,135,198-byte transcript is `scripts/out/S11c_d_wide_three_momentum.out`,
SHA256 `61c531c071907d199d68a7665b55e72d5ba087ce637aba4bec811c37d76a39f2`.
Publication and annex identity are recorded in the execution/checkpoint files.
The earlier 14.6–16.4-hour timing estimate came from short prefixes; the full
measured runtime supersedes it for future cost planning. All original operands,
frozen inputs, worker groups/records/partials and logs remain in repository scratch.

Next: independent adaptive outer integration of all ten triple rows on both
cutoffs3/4 and both fields, holding accepted inner32/32 and source/profile256/512,
position bounds48/14 and regulator0.2. Keep the complete single/pair contributions
explicit. A tested conditioned native recursion already exists; validate its
wide-triple application before production. Do not repeat the completed grids.
No uniform/independent-grade/global exceptional coverage, infinite-tail or Abel
limit, scattering or physical pole result is established by this refinement.

The user-approved nonlinearPoleV2 addendum governs new work. The current run's
v10 authority, engine and imported exports remain byte-identical. New stages
must pin the addendum with an explicit dependency disposition; old results are
not relabeled. The retained solver/export suffix remains unchanged.
