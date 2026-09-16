# S11c-d SymPy builder checkpoint

RIGHT current/adjoint normalization is validated and published for the supplied
LAB_HELD / RHO4_CONSTANT case: all 18 root/lift candidates, 22 basis directions,
18 invertible field maps, and two current-normalized two-dimensional subspaces.
Validation covers 5,990 original tags and 5,850 numerical residual scalars.
The nine nonzero raw balance/reconstruction families are preserved. Independent
contractions of the saved eta-squared pairing remainder account for them to
2.15e-13; 275 decomposition/projection/limit residual scalars are zero. No new
upstream physical repair was needed. This is retained-order accounting, not an
exact finite-contrast balance or a higher-order continuum prediction.

The thickness-coordinate repair and native b/c1/c2/d regeneration are committed
through 8b2e3cf2. All eight endpoint source/frequency/pairing prerequisites remain
committed through 432db7e7; RIGHT, LEFT and REFERENCE each have 802 zero retained
pairing residual scalars. RIGHT is committed at 0da0f746. LEFT normalization
is now validated and published: 18 candidates, 22 basis directions and 18
adjoint field maps, with no raw norm above the diagnostic threshold. Its
independently computed discarded balance matrices vanish; the largest
residual-minus-remainder norm is 8.21e-13.

The [Mathematica audit](S11c_wolfram_repair_audit_report.md) found a duplicate
pressure shift in native c2 N6. The three mechanical issues are absent on the
audited b domain (200 zero residuals). The [repair plan](S11c_wolfram_pressure_trace_repair_plan.md)
is complete on its validated domain: 218 zero focused residuals and 9,968,256
zero native numerator evaluations across four cases under the adopted covariance
criterion. Raw R_N6 remains a representation diagnostic. The repaired main
transcript is committed at `5acfdf30` and its annex payload hash is verified.
REFERENCE normalization is validated: 18 candidates, all 22 basis directions,
18 invertible maps, 5,586 tags and 5,850 numerical residual scalars. Its discarded
balance matrices vanish; all 275 exact remainder checks are zero and the largest
residual-minus-remainder norm is 8.21e-13. Pre/post-emission packet joins agree.
RIGHT/LEFT/REFERENCE normalization is complete on the supplied case.
LEFT is committed at b271c71b; the accepted SymPy chain remains unchanged.
See [the native checkpoint](S11c_wolfram_pressure_trace_native_checkpoint.json).
See the [normalization report](S11c_d_end_normalization_report.md) and
[plan](S11c_d_end_normalization_plan.md) for complete evidence and boundaries.

All working runs and source snapshots live in the repository's _scratch/s11c/.
Published .out files use DataLad/git-annex with post-save full-hash verification;
ordinary sources, reports and inventories use Git. The user controls continuation.
Owned jobs use user-authorized local completion/error wake-ups, with no recurring
checks or model polling.

The [matching-channel construction](S11c_d_matching_channels_report.md) now
is validated: all 36 candidates and 44 basis directions are retained, with two
incoming and two outgoing directions per end. All cross-mode currents were
evaluated; reconstruction residual norms are below 6.34e-16. The 500-tag packet
and 2,763 metadata paths passed full replay. The packet is committed and annex-verified at 3a5d3d25. Reduced operator
assembly now has validated full five-slot source rows and all five native
probe-action columns. Exact assembly is validated: four local derivative
matrices, 80 distinct nonlocal integrals and 281 zero reconstruction/extraction
residuals; 1,226 metadata paths passed replay. The user approved the proposed values
for all 30 inherited free gradient-energy coefficients. The separate development
input preserves every original parameter and profile. Numerical action evaluation
continues with the accepted end-dependency proof and unchanged symbolic operators.

The complete two-ended variable-profile S-matrix, profile-frequency bound
poles/residues/overlap, survival, flux bookkeeping, weak coefficients, section 5
controls and final own-row export remain program work. Per-root and full-subspace
checks retain their stated sheet and exceptional-domain limits. Supplied physical
premises and c2 operand debt remain inherited premises. No S11c_d_exports.py exists.

The numerical action check completed: all 300 components agree within
2.49e-16, retaining 1,920 nonlocal contributions and passing mutation checks.
All 106 tags and 3,611 metadata paths replayed, and source/packet hashes agree.
Grid changes reach 0.04670; quadrature refinement, physical tails and Abel
weak limits remain work before boundary matching. Publication is annex-verified
at 009667f7. Quadrature-domain recovery passed with empty stderr: all source
and operand joins agree. Width-based split rules resolve the elementary Abel
mass to 5.33e-15; profile transforms agree to 1.25e-12 at the refined order.
The domain transcript is annex-verified at e514b704. Bounded source Fourier
factorization is now validated: all 80 original operators, 35 distinct source
integrals, 324 zero normalized residuals and 264 zero certificate proofs.
All 58 original expanded residuals remain beside exact saved-pair checks;
all denominator/branch restrictions are retained. The source packet is
byte-identical and its constructor unchanged. Full replay covers 1,670 tags,
4,752 metadata paths and 833 fresh write-keys. The transcript is published and annex-verified at bed088be. Source Fourier quadrature is now validated:
35 distinct source integrals on both approved Gaussian fields, 140 interval
records and 9,240 frequency evaluations. All 320 native occurrence joins and
104 binding proofs pass, with 42 additional zero replay proofs and all coefficient
and measure mutations detected. Scaled source/adaptive residuals reach only
1.73e-14 / 4.11e-13. Gauss changes are 1.43e-6 / 4.11e-13 and the largest finite
source-interval change is 4.40e-9. The 8.69 MB transcript has 10,400 tags and
151,318 metadata paths; source and numerical packet hashes are unchanged.
Original live raw forms remain distinct from their canonical pickle forms,
with twenty exact metadata-support transitions checked and forty altered-unit/
order controls rejected. See the source quadrature report and accepted checkpoint.
Next is complete-action momentum quadrature with all three Abel transfer pairs
and nested profile factors. Full-action convergence, physical tails and Abel
limits remain work before boundary matching; finite source tests do not replace
these requirements.

Finite momentum integration is validated on the declared finite grids: all
80 native integral rows, 70 bound sources, six nested profiles and three Abel
pairs are evaluated. All 300 native action and 1,920 term comparisons agree to
1.11e-16 scaled / 3.60e-17 absolute. Refinement changes fall roughly fortyfold,
to 3.64e-5 / 1.05e-5 for the two fields; further resolution is needed before
physical tails and Abel limits. The 2.22 MB transcript replays all 5,398 tags,
2,697 keys and 22,822 metadata paths, with unchanged operand packets. See the
momentum action report and accepted checkpoint. No scattering/pole claim follows.

The targeted one-momentum refinement is validated: final outer/source changes
are below 6.64e-15 / 3.81e-15; adaptive outer quadrature agrees within 1.90e-15
in complete actions. The accepted baseline and all held native terms replay
exactly. The 3.84 MB transcript has 13,230 tags and 64,385 metadata paths, with
unchanged packet hashes. Next is two-momentum refinement with the accepted
one- and three-momentum terms held explicitly. Full physical limits remain work.

Two-momentum refinement is validated: final outer/panel/source/profile raw
integral changes are 8.33e-14 / 1.50e-12 / 2.62e-16 / 8.55e-15. Adaptive outer
quadrature agrees within 2.45e-14 in integrals and 6.21e-17 in complete actions.
All 30 rows, 22 records and 91,047 metadata paths pass, with unchanged held terms
and packet hashes. The 5.65 MB transcript is published and annex-verified.
Next: the remaining ten three-momentum rows, before physical tails and Abel limits.

Initial three-momentum refinement is validated in 3h22m: all ten rows and both
fields, with final outer/innermost/middle raw changes 2.68e-12 / 8.22e-15 /
5.61e-14. All ten records, 10,184 partials and 34,691 metadata paths pass, with
unchanged held terms, source hashes and packet identities. The 2.24 MB transcript
is published and annex-verified. Next are source/profile order refinements on
the fixed final momentum grid; physical tails and Abel limits remain work.

Three-momentum source/profile refinement is validated in 2h18m with four
single-thread workers. The preserved 11,616,256-node prefix was resumed without
changing native summation. Raw source/profile order changes are at most
6.14e-19 / 3.64e-17. All ten rows, six records and 26,705 metadata paths pass;
source, worker, partial and pre/post-packet hashes agree. The 1.52 MB transcript
is published and annex-verified. Next is independent adaptive outer quadrature,
then physical domain/tail and regulator checks; finite agreement is not a limit.

Independent adaptive outer production is validated in 2h12m. All ten rows
and both fields agree with Gauss-144 within 1.74e-16 in raw integrals; complete
actions round equal. All 420 conditional points, 6,764 partials, 67 sources and
20,678 metadata paths pass with unchanged packets and held single/pair terms.
The 1.18 MB transcript is published and annex-verified. Next are finite-domain
and tail checks. Inner/source/profile rules were held, so independent outer
agreement does not establish uniform/grade coverage, tails or Abel limits.

Finite position-domain production is validated in 2h26m. All 80 rows and both
fields were reevaluated at source/profile bounds 48/10 and 48/14. Source/profile
domain changes reach 8.47e-14 / 4.40e-15 in complete actions. All 12 layouts,
8,484 partials, 71 sources and 40,672 metadata paths pass with unchanged packets.
The 2.31 MB transcript is published and annex-verified. Next is momentum-domain
expansion with freshly computed source-frequency/profile-transfer ranges and
resolution checks. Finite changes do not establish physical tail or Abel limits.

Momentum-domain preparation is validated in 7m37s: 152 source/profile records,
all 80-row coefficient probes, 74 sources and 263,184 metadata paths pass.
Source/profile adaptive residuals reach 7.17e-13 / 1.15e-12. At the new wider
transfer range, profile 256-to-384 changes reach 0.00404, falling to 1.80e-12
for 384-to-512. Next compute source256/profile512 complete actions with a
matching-rule cutoff-2 baseline before cutoffs 3/4. The 7.30 MB preparation
transcript is published and annex-verified; wider actions and physical limits remain.

Momentum-domain action preflight 649a98d2 passed 152 selected transforms,
six exact prefixes and six exact coarse full actions, with 37,614 metadata
paths. The coarse smoke remains instrument evidence only.

The separately discovered S11 Q9 coefficient-action orientation defect has
been dependency-traced. Its 72-row family is carried in accumulated exports,
but c1/c2/d import manifests exclude it; the actual d binder retains identical
73 inputs when all Q9 rows are removed. The current 78-source production hashes
match and its numerical operands remain unchanged. The run completed against
those immutable inputs. Refresh carried rows/provenance after repair acceptance,
with explicit consumed-root/closure joins. D3-D5 Q9 validation and parity-odd
extra-action checks remain owned by that repair. See the Q9 dependency report.

The complete momentum-domain production is now accepted: six clean workers,
240,460,912 new nodes in 9 h 17 m, 18 saved layouts, 14,672 partials and
66,798 metadata paths. The matching profile-order baseline changes the action
by 6.94e-18; cutoff 2→3 and 3→4 changes reach 2.25e-7 and 1.14e-7.
Wider-box paired-momentum refinement is next before tail interpretation.
The 3.41 MB canonical transcript and checkpoint retain every source/packet join.

Wider-box one-/two-momentum production is accepted in 10m28s: four clean
workers, 4,674,896 new nodes, 28 records, 274 partials and 236,526 metadata
paths. Final single and paired outer/inner changes reach 2.87e-15 and
4.84e-16 / 3.40e-16 in raw integrals. The refined cutoff 3-to-4 action change
remains 1.134674e-7. All held three-momentum values remain explicit. The
10.64 MB transcript is published and annex-verified. Next is independent
outer quadrature on boxes 3/4, followed by needed three-momentum checks;
physical tails, regulator limits and scattering remain open.

Independent wider-box outer production is accepted in4m20s: all40 single
and30 paired rows agree with refined Gauss within4.83e-15/7.71e-16; all eight
adaptive solves finish. The box3-to4 action difference remains1.134674e-7.
Four clean workers,84 sources,3066 conditional points and101909 metadata paths
pass. The6.75MB transcript is published and annex-verified. Next refine the
held wider-box three-momentum terms at fixed source/profile/domain/regulator
settings. Finite agreement does not establish physical tails or Abel limits.

Wider-box triple instrument preflight is accepted in 114.51 seconds: four exact
saved prefixes, four finest-rule cost prefixes, twelve isolated settings and
twelve exact coarse serial/worker plus independent cell comparisons. Four clean
workers, 87 sources, 5040 held-term scalars and 75886 metadata paths pass.
The 5.07 MB smoke stays instrument evidence in durable scratch. Production launched after acceptance 8ac0caaf at
216/24/24→216/32/24→216/32/32 with all single/pair values held. Four single-thread
2 GiB workers, 87 current/frozen sources and the silent watcher are verified;
no production result is accepted yet.
Prefix extrapolation suggests roughly 16 hours, with concurrent-load uncertainty.

## Retained user-approved solver/export contract


1. Preserve `EdgeReduction`, the positional three-parent fold, and the exact
   direct-lookup manifest. All numerical assembly must consume the computed
   reduced rows, including the full nonlocal terms and full coupling vertex.
2. Accept explicit, independent dimensionless profile functions w(xi), m(xi),
   their derivatives, asymptotic limits and tail information. A selected smooth
   step with an independently adjustable localized modulus bump is a numerical
   instance; it does not replace the interface class. Store the profile formula
   and digest in every case record. The current preflight input is recorded in `S11c_d_channel_preflight_input.json`.
3. Require a complete parameter map in a declared L/T/M unit frame, real
   continuum frequency and tangential momentum, small contrast, and
   sigma_W = eta_bg W_0/L_W for evaluations on the physical homotopy. Retain
   independent eta/sigma grades in the symbolic calculation. Test actual
   reference/end channel availability before attempting flux normalization.
   Do not manufacture an incident channel by assigning a sector label.
4. Compute both-end modes, left/right normalization, the S11b-derived current,
   and the variable-profile matching problem. Re-expand the continuum response
   to the retained rectangle; do not present a finite-contrast numerical
   solution as a higher-order continuum prediction. Retain evanescent matching
   modes, channel degeneracies and domain failures explicitly.
5. Numerical pole searches have an explicit profile, parameter map, sheet,
   bounded search region and isolating contours. Evaluate the retained operator
   without the continuum re-expansion. Record boundary/quadrature resolution,
   domain size, precision, root residuals, contour-count evidence, and changes
   under refinement. A bounded search does not establish a global pole set.
   An unsuccessful or inconclusive search is unresolved, not an empty pole set.
   Compute residues/projectors and sheet/decay/width/closure tests only for
   actually resolved candidates; emit spectral overlap, not capture probability.
6. Separate transparent symbolic expressions from evaluated numerical records.
   Symbolic operator/continuum/weak-coefficient exports remain differentiable
   SymPy expressions, compacted with algebraic equivalence checks. Numerical
   mode and pole datasets retain their input bindings, domain, convergence
   evidence, dimensions and truncated-model status. They are not stand-ins for
   a generic symbolic profile-dependent root function. This is the export
   distinction motivating the user-approved contract; the downstream consumer
   will need to bind the appropriate representation explicitly.
7. Fingerprints summarize already constructed/evaluated objects. They do not
   evaluate nonlocal integrals or replace a spectral solve. Algebraic PIT and
   physical numerical evaluation are separate records. All completed roots use
   fresh lowerCamel write-keys and the existing bind-closure/minimal-delta guards.
8. Finish the one-case path, then implement controls and bookkeeping, then run
   all four cases once and write the complete export. No review legs,
   comparator, Wolfram engine, downstream stage, or commit belongs to this lane.
