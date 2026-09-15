# Parallel preflight: standalone transcript recovery

The four prefix workers and four complete underresolved workers exited zero with
empty stderr. The original preflight compared all numerical prefixes, native
terms and complete actions before reaching emission. The second standalone
transcript failed its decoder because the preflight cleared tag state but reused
the first transcript's shared-payload encoder. Both original transcripts and all
numerical packets remain preserved in the validation run.

Reset the encoder before each independent preflight emission. Prove the entire
preflight AST rejoins its frozen version after removing exactly those two resets.
Keep the physics engine, codec, serial constructor, emitter and numerical worker
implementation byte-identical. Use a separate recovery validator to:

1. Pin every original artifact and source before reading saved results. Verify
   every child exit/log, source/provenance join and partial/record packet hash.
2. Recompare all four 65,536-node prefixes and replay only the native frequency
   census over their node positions. Reconstruct complete actions from the saved
   integrals and compare every native term and held contribution. Do not repeat
   source/profile or momentum integration.
3. Account explicitly for workspace telemetry: a reused serial process reports
   the cumulative peak; an isolated worker reports its own peak. Require the
   serial value to equal the running maximum of the independently saved worker
   values. Every remaining numerical record operand must agree exactly.
4. Require the original second transcript to fail standalone decoding and to
   decode exactly when read after its original first stream. Re-emit both saved
   results with independent encoders, replay every metadata path, and compare
   every decoded payload with both original streams. No decoder exemption.
5. Verify original and copied packets remain byte-identical after emission. The
   failed run did not save its serial monotonic timer: record the file-time phase
   interval and partial interval separately from measured parallel elapsed time;
   any speedup inferred from file times is approximate, not an exact benchmark.

Commit acceptance before launching the already prepared four-worker production
continuation. These deliberately underresolved finite quadratures establish
scheduling equivalence only. No physical convergence, domain limit or scattering
result follows. Use repository scratch and the existing silent completion watcher.
