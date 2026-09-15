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
factorization remains pending. Phase extraction was repaired at 91f56ec6;
the retry exposed ten native operands with middle momentum outside source
position. The finite-domain adapter now locates source position by identity
and preserves the other limits. All 80 layout checks and 32 focused residuals
pass; the native engine prefix is unchanged. Per-row hashed saves now retain
completed downstream work. Repair 50c51500 is committed; the full
factorization/metadata retry was launched with a local completion/error watcher.
See the domain and factorization reports for evidence and scope limits.

All 80 bounded Fourier rows are now saved, with 35 distinct source
integrals. Full metadata/source replay completed, but 58 expanded
reconstruction residuals require exact normalization. A saved-pair recovery
is prepared; representative pairs pass and a coefficient mutation is
detected. The complete source packet, original residuals and all denominator
restrictions remain preserved. No factorization publication is accepted yet.

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
