# Mathematica c2 pressure-trace repair

The user approved the repair after audit checkpoint `946acc84`. Native c2 now
extracts the reference-pressure evaluation map from its actual face law and
solves its ordered inverse before filling pressure and normal-jet slots.
The original physical c1 response remains a separate operand. The two-leg
conversion and the height map acting on the three-leg first-shape response
keep output, input and intermediate momenta distinct through the full eta/sigma
rectangle. Both weak contractions and the native closed-slot construction use
the computed reference images. Four new trace tags expose physical/reference
responses, reconstructed physical pressure and the literal residual.

The [focused checkpoint](S11c_wolfram_pressure_trace_focused_checkpoint.json)
records **218 zero residual scalars** across the generic ordered matrix joins,
separately targeted zero-frequency solves, 16 native anchoring/density/route/face
row cases and four exact shifted-wave regressions. Each of six missing/doubled/
reversed map mutations gives three nonzero matrix components. Identifying the
intermediate momentum with the input changes four components on each face.
The independent boundary-amplitude matrix solve does not call the native
reference-pressure recurrence. The native conormal/pressure kernels are shared
inputs to that ordering check, not independently re-derived physics. The shifted
wave regression independently reconstructs pressure and its normal continuation.

The complete focused run took 20.51 seconds at 398,136 KiB peak RSS, with exit
zero, empty stderr and stable source hashes. Its 1,642,937-byte
[transcript](../mathematica/out/S11c_wolfram_pressure_trace_repaired.out) retains
both operands, raw/retained/discarded parts and restored dimensions/grades.
The first development attempt is preserved: repeated general series expansion
was replaced in the checker by direct coefficient extraction when the denominator
is independent of eta/sigma. This changes no physical producer operand. The
rerun completed in 20.80 seconds with 158 zero development residual scalars.

The five existing ablation-harness selector literals still occur exactly once;
no historical ablation result is being treated as a new repaired run. The old
main c2 output remains preserved in annex history. The native
development case completed in 1,004.16 seconds at 588,792 KiB peak RSS with empty
stderr and all source/prior-output pins stable. Its transcript has 126,423,192
bytes. The [development validation](S11c_wolfram_pressure_trace_native_development_checkpoint.json)
finds zero nonzero numerators across nine required residual families, with no
zero sample denominators. All 144 new trace components have complete family,
face, grade and restored-dimension joins.

The [four-case native checkpoint](S11c_wolfram_pressure_trace_native_checkpoint.json)
records 44 tags and **9,968,256 zero numerator evaluations** across the nine
required residual families, with no zero checked sample denominators. Each case
uses 17 recorded real-axis cells, three primes and eight draws (408 sample rows).
All 576 new trace components have complete family/face/grade/dimension joins.
The CAS run took 3,196.55 seconds at 1,062,960 KiB peak RSS; inert validation
took 601.61 seconds. Both exited zero with empty stderr and stable source pins.
The validated 496,254,757-byte [main transcript](../mathematica/out/S11c_c2_N6_mathematica_audit.out)
has SHA-256 `08bfb8be8097e97b711918d32bc7ad2127298fb0305e25c0c1ff2faf3ea9c254`.
Both completion watchers delivered their events without model polling.

The [emission inventory](S11c_wolfram_pressure_trace_emission_delta.json)
finds 26 unchanged assignments, 14 changed assignments, four added trace tags
and none removed. Unchanged assignments include the carrier/source operands,
covariance families and slot guards. Pressure-contracted operands, channel
representations, closure guards and associated metadata change. Changed hashes
do not imply a nonzero residual; the residual validation above is separate.

The first validator incorrectly required raw `R_N6` equality. The adopted
[Reading B](S11c_c2_N6_RESOLVED.md) requires operator covariance and channel
reconstruction; raw coefficients in different field frames need not coincide.
Both the old and repaired development transcripts have the same 18 nonzero raw
entries (20,736 numerator samples). The failed validation and source snapshot
are retained. The validator now records those diagnostics and checks covariance,
its increment, carrier/cross channels, the split, slot/closure guards and physical
trace reconstruction. This corrects validation scope; no producer rerun or new
physical repair is indicated. The full run retains 20,736 nonzero raw numerator
samples in each of LAB_HELD/RHO4_CONSTANT, LAB_HELD/RHOBR_CONSTANT and
MATERIAL_ADVECTED/RHOBR_CONSTANT; MATERIAL_ADVECTED/RHO4_CONSTANT has zero.
These are representation diagnostics under Reading B. Cross-engine operand
agreement remains an open debt.

b/c1 sources and outputs, all SymPy producers/exports and accepted LEFT/RIGHT
normalization, authorities, S10/Lean and the retained solver/export contract are
unchanged. Matching-denominator and threshold exclusions remain explicit. This
is a repair of the existing c2 N6 instrument, not a full c2 self-energy or d
engine. No review, comparator, recurring scheduler or push has run. This closes
the authorized pressure-trace repair on its validated domain. REFERENCE
normalization has resumed after publication commit `5acfdf30` and full-hash
verification. The main transcript is a locked annex symlink with key
`MD5E-s496254757--7a429434ad5c5157b0c5a75287db6784.out`.
