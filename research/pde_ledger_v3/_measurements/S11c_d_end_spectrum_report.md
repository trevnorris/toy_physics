# S11c-d physical-field end-spectrum coverage

2026-09-10–11. The focused checks, fresh four-case regeneration, inventories
and canonical output publication completed. This extends checkpoint `c5af3181`
and remains unreviewed. The complete S11c-d engine and export remain unfinished.
The [plan](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_plan.md)
records the authorized scope and stop condition.

## Computed construction

[EndSpectrumCoverage](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2747)
in the existing engine consumes
both computed constant-end representations: the five-field closed operator
and the full canonical sector pencil. The actual field ansatz constructs the
coordinate lift and its dual test map. The matrix pullback residual compares
these separately built representations before spectral elimination.

At a bound input, the constructor clears the physical pencil's rational
row denominators, computes its determinant, and reduces its numerator on the
computed radical relation. It emits the eliminated-momentum remainder instead
of assuming evenness. Polynomial gcds identify intersections with denominator,
normal-threshold and radical-branch loci. The quotient determinant and coordinate
factor remain separate operands; their pullback residual is computed.

The [isolate computation](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:2776) square-free-factorizes the resulting
polynomial and computes roots at 50 and 80 decimal digits. Exact rational
Taylor bounds around rational centers establish one-root disks when the linear
term strictly dominates the remainder. Exact disk-separation comparisons and
the polynomial degree account for all finite polynomial roots with their
multiplicities. The certificate concerns the algebraic radical curve at the
bound input; membership in the Fourier sheet chart is a separate record.

The original rational five-field matrix is then evaluated at both normal-root
lifts. It supplies left/right nullspaces, literal field-equation residuals and
oblique nullspace projectors. Classifiers use field lifts derived from the
existing sector ansatz, choosing among the three axial charts. Mixed spaces
and singular charts retain explicit domain statuses. These are nullspace
projectors, not bound-pole Riesz projectors or flux normalization.

Heavy matrices use fingerprints of already constructed objects and SHA digests.
The input unit frame, parameter binding, source grade support, local grade
origin, physical mode/root dimensions and restored equation dimensions accompany
the records. Polynomial arithmetic and SVD diagnostics are coefficients in the
declared numerical unit frame. Finite-contrast end records evaluate the retained
operator; they are not re-expanded continuum predictions.

The normal full run adds one native algebraic-PIT packet at the reference and
both ends of each case. With an explicit channel input it also adds the three
physical-input packets. `--channel-input-scope spectrum` selects those input
spectra without duplicating the optional legacy local-jet packets. The default
`modes-and-jets` setting preserves that earlier input path.

## Focused evidence

All five final checks use the pinned symbol caches and transcript from the
inverse-Fourier run. The source/parent/input hashes stay fixed. Each exits 0
with empty stderr and empty dimensional constraints. The instrument and exact
commands are retained in the run record.

| Check, LAB_HELD/RHO4_CONSTANT | Wall seconds | Peak RSS, KiB |
|---|---:|---:|
| Physical reference | 19.2944 | 147512 |
| Physical left end | 21.4160 | 149504 |
| Physical right end | 25.7554 | 152376 |
| First tangential momentum set to zero | 20.6175 | 150804 |
| Algebraic PIT sample 0, right end | 25.7182 | 152396 |

Each physical polynomial has degree 11, nine distinct radical roots and
18 distinct `(k,q)` candidates. Four candidates have nullity two and fourteen
have nullity one. Every isolation inequality and disk-separation comparison is
positive, and every degree-count residual is zero. Root refinement differences
are at most `4.495e-50` in inverse-time reference units. Physical and determinant
pullback residuals are literal zeros. All new metadata groups are present.

For the physical reference and both ends, ten candidates match the existing
fixed-frequency Fourier chart and eight are on its other branch. The first-axis
mutation also computes those counts while making the old axial quotient
identically singular; the physical solve and alternate field charts remain
available. These are candidate records, not counts of open flux channels.
The algebraic PIT packet retains four unresolved sheet labels.

The physical-input right-equation residual maxima by restored `[L,T,M]` are
`3.8292e-15` at `[-2,-2,1]`, `2.1593e-15` at `[-3,-1,1]`, and `2.1453e-15`
at `[-1,-2,1]`. Left-equation maxima are `3.5467e-15` at `[-1,0,0]` and
`2.6853e-15` at `[0,0,0]`. The separate PIT packet has a maximum left residual
`7.994e-15` at `[0,0,0]` and projector residual `4.443e-14` at `[0,0,0]`.
These numerical residuals do not replace the exact polynomial root-count
certificate or establish a global error bound on modal fields.

## Four-case evidence and preservation

The full run exited 0 with empty stderr in `3963.933` seconds (66.1 minutes),
peak RSS `1724396` KiB. Source hashes before and after the run match. The
[95,506,147-byte main transcript](/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S11c_d_mixing_scattering_sympy_audit.out)
has SHA-256 `1db1d902fb35e00002e375b33a6e04db350896fe5f26a2aa6435463e0f5aca77`,
33,436 unique tags, one completion marker, no metadata gaps, and three empty
dimensional-constraint records. Its emitted and recomputed input digests agree.

All 24 native reference/left/right packets have degree-11 physical polynomials,
nine distinct radical roots and 18 distinct `(k,q)` candidates: 432 records in
total. Each root-isolation certificate accounts for the full polynomial degree;
all disks are isolated and mutually disjoint. The denominator, normal-threshold
and radical-branch gcd degrees are zero in these packets. The largest 50-to-80
digit root-refinement difference is `5.796e-50` in inverse-time reference units.

| Input family | Packets | Candidates | Nullity 1 / 2 | Sheet true / false / unresolved |
|---|---:|---:|---:|---:|
| Explicit physical input | 12 | 216 | 168 / 48 | 120 / 96 / 0 |
| Algebraic PIT sample 0 | 12 | 216 | 168 / 48 | 72 / 96 / 48 |

Every computed algebraic/geometric multiplicity difference is zero. Each input
family has 168 thickness-like and 48 transverse-like candidate classifications.
The sheet counts use the existing fixed-frequency chart and do not count open
flux channels. The 48 unresolved native PIT labels remain explicit.

All physical/quotient matrix and determinant pullback residuals are zero. The
right-equation maxima by restored `[L,T,M]` are `4.997e-15` at `[-2,-2,1]`,
`3.228e-15` at `[-3,-1,1]`, and `7.274e-15` at `[-1,-2,1]`. Left-equation
maxima are `3.250e-15` at `[-1,0,0]` and `7.994e-15` at `[0,0,0]`.
Projector-idempotence maxima are `4.443e-14` at `[0,0,0]`, `5.565e-15` at
`[1,0,0]`, and `1.449e-15` at `[-1,0,0]`. Radical residuals are at most
`1.813e-48` for physical inputs and `8.408e-52` for PIT, both at `[0,-2,0]`.

All 12 earlier full-sector symbol payloads and all 24 legacy mode packet
payloads (528 candidate records) match `c5af3181` literally. On the eight common
end/PIT0 packets, the native roots and earlier regular roots each number 18;
the largest nearest-root difference is `2.456e-16` in the declared numerical
frame. This comparison excludes the legacy chart-limited candidates.

The unchanged inverse-Fourier inventory finds all 96 carrier inverses and their
288 dimensioned residual entries zero; the 96 source-image residuals, 288
remainder entries and 24 branch residuals are also zero. All projections of
the 222 integral and eight row residual fingerprints are zero. The five-slot
census and metadata coverage remain complete. Among pre-existing engine
definitions, only `run` changed; the new constructor is `EndSpectrumCoverage`.

The [run record](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_runs.json)
retains commands, source/cache provenance, resource measurements and output
hashes. The [full inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_inventory.json),
[focused inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_focused_inventory.json)
and [inverse-Fourier preservation inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_end_spectrum_inverse_inventory.json)
retain the computed summaries. All six successful transcripts are published
under `scripts/out/` (98,218,875 bytes total); atomic replacement preserved the
old annex payload. The engine and both new instruments compile. The directive's
`reduction/derived_or_declared.py` and `reduction/engine_output_checks.py` remain
absent; no claim is made to have run them.

## Remaining scope

The broader full-spectrum TODO remains open. This construction establishes
regular finite algebraic coverage at bound inputs; generic sheet continuation,
cut-bank/continuum treatment, threshold/denominator intersections and
classification across mixed degeneracies, and generalized modes at defective
roots still need treatment before a complete physical end-channel space is available.
The completed regular jets and all previous reduction constructions are retained.
No upstream repair has been indicated by these checks.

Full nonlocal current/flux normalization, complete two-ended scattering,
poles/Riesz/overlap, survival, flux bookkeeping, weak coefficients, Section 5
controls and the own-row export remain uncomputed. The transform and coordinate
checks above are not Section 5 physical profile-FORM controls. Section 1 inputs
remain supplied and unfalsifiable here; the separate shear-normalization and
c2 cross-engine operand/sign debts remain open. No review leg, comparator,
Wolfram engine, downstream stage or commit is part of this build.
All ten live TODO items remain. `scripts/S11c_d_exports.py` remains absent
because its required constructions are unfinished.
