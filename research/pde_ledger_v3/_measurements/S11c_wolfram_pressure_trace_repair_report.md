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
main c2 output remains historical until native regeneration completes. The native
development case completed in 1,004.16 seconds at 588,792 KiB peak RSS with empty
stderr and all source/prior-output pins stable. Its transcript has 126,423,192
bytes. The [development validation](S11c_wolfram_pressure_trace_native_development_checkpoint.json)
finds zero nonzero numerators across nine required residual families, with no
zero sample denominators. All 144 new trace components have complete family,
face, grade and restored-dimension joins. The four-case regeneration is next.

The first validator incorrectly required raw `R_N6` equality. The adopted
[Reading B](S11c_c2_N6_RESOLVED.md) requires operator covariance and channel
reconstruction; raw coefficients in different field frames need not coincide.
Both the old and repaired development transcripts have the same 18 nonzero raw
entries (20,736 numerator samples). The failed validation and source snapshot
are retained. The validator now records those diagnostics and checks covariance,
its increment, carrier/cross channels, the split, slot/closure guards and physical
trace reconstruction. This corrects validation scope; no producer rerun or new
physical repair is indicated. Cross-engine operand agreement remains an open debt.

b/c1 sources and outputs, all SymPy producers/exports and accepted LEFT/RIGHT
normalization, authorities, S10/Lean and the retained solver/export contract are
unchanged. Matching-denominator and threshold exclusions remain explicit. This
is a repair of the existing c2 N6 instrument, not a full c2 self-energy or d
engine. No review, comparator, scheduler or push has run. REFERENCE normalization
resumes after the native regeneration checkpoint closes this repair.
