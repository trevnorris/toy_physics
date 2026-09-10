# S11c-d SymPy build — COMPLETION mandate (resume pass; extend the existing prototype)

⭐ **This is a RESUME of an incomplete build, not a new build.** Your authority is unchanged: the cleared build
directive `directives/S11c_d_sympy_build_directive.md` (governs the build-mechanical layer) and, through it, the
cleared physics spec `directives/S11c_d_SHARED_PHYSICS.md` (v10, `399a8516`, the physics authority). Read both.
Model: `gpt-6-astra`. All of the directive's obligations bind — the three script clauses, input-driven construction,
the §3 census (all five payload slots), leak discipline, and the builder lane (build → run → report → **stop**; ⛔ no
legs, comparator, WL engine, downstream steps, or commits).

## What already exists (⛔ READ and EXTEND — do NOT rewrite from scratch, do NOT discard the reduction)
`scripts/S11c_d_mixing_scattering_sympy_audit.py` (a running prototype) + its `.out` +
`_measurements/S11c_d_sympy_builder_report.md` (its own account). The prototype implements, and you must **preserve**:
- the frozen import wiring (the 3-parent `load_model` fold; `IMPORT_KEYS`);
- the §3 3-D→1-D **Fourier reduction** (`EdgeReduction`, `hat()` at every carrier argument, the action-integral
  reduction with tangential-constraint solve + Jacobian + normal integrals, ONM processed at all three sites, the
  five-slot census). ⭐ Keep its correct discipline: `hat()` computes the engine's own reduction and does **NOT**
  read `FOURIER_PROFILE_BINDINGS` as the answer; `profile_definition()` processes that slot as the **(i) 3-D operand**
  only, ⛔ never substituted as the reduced answer (directive §3 F1 classification).

## COMPLETE these remaining obligations (the report's "Unimplemented obligations", in build order)
Implement each per the directive/spec (⛔ compute, never type; ⛔ no empty-set / zero / baseline substituted for an
unexecuted construction — the prototype correctly refused to, keep that discipline):
1. **Finish the reduction layer:** per-object restored `[L,T,M]` dimensions + able-to-fail dimensional checks,
   computed `(ε,η,σ_W)`/`λ` multigrades, the both-row **reconstruction round-trips**, and the zero-jet
   unequal-asymptote distributional prescription (spec §1c).
2. **§2:** the reduced full pencil `𝓛` from the reduced closed rows, the canonical off-diagonal extraction vs the
   reduced kernel (emit the residual as a **surfaced finding**, ⛔ do not silently override), computed `K₀,K₋,K₊`
   and reference/end baselines, the two full asymptotic block pencils, modes + classifiers.
3. **§3a:** modal energy current (S11b bilinear on the reduced operator), nonlinear-pencil normalization, the
   **complete two-ended channel S-matrix** (both incident ends, every open channel), conversion amplitude, continuum
   `T→H` flux functional.
4. **§3b:** pole set + normalized Riesz residues/projectors + sheet/normalizability/width/all-channel-closure tests +
   spectral overlap (⛔ no capture probability), and the transverse survival functional.
5. **§3c/§3d:** amplitude components + multigrade, physical flux baseline/interference/quadratic slots, total
   conversion fraction, induced-field quadratic form, N12 baseline/interference operands; weak Taylor coefficients +
   the named strong-edge handoff.
6. **§5 controls:** 5a coordinate-covariance regression + shape-sensitivity mutations (⛔ RHO4 absence is a computed
   structural absence, not an `A−A`), 5b the three uniform regressions, 5c profile-FORM ablation + edge-vs-bump +
   modulus discriminants (two separate operands, ⛔ not the identity as one payload), 5d the computed
   flux-normalized dimensionless conversion FORM (⛔ no typed expected shape). Each = the object **and** its literal
   residual (both operands).
7. **The export:** write `scripts/S11c_d_exports.py` as the own-rows delta (directive §5 membership: S-matrix /
   conversion amplitude / continuum flux + survival / **bound pole set + Riesz + spectral overlap** / §3d weak
   coefficients + recursive closure), with the bind-closure `assert_delta_is_minimal` guard and the casewise
   compact-vs-expanded semantic check. Everything else is EMIT-only (§5).

## Heavy-CAS discipline (governing — the S-matrix + pole solve + full symbolic construction can blow up)
- **Measure the process that runs** (the prototype's reduction ran ~152 s at ~1.66 GB — cheap; the construction on
  top may not be). ⛔ Never run two memory-heavy CAS jobs concurrently.
- If a symbolic object is impractical full-symbolic (e.g. the pole solve, a full symbolic S-matrix over 4 cases),
  use a **carrier-first + numeric-PIT fingerprint** representation (a compact numeric residual/fingerprint + a SHA
  digest of the symbolic object), ⛔ not a multi-hundred-MB full-symbolic dump. Store EXPORTS transparently factored
  (⛔ not `sp.expand`ed, ⛔ not opaque hold). Defer a genuinely heavy control in-band→out-of-band
  (`DEFERRED_HEAVY_RUNS.md`) and NAME it — ⛔ never silently drop it.

## Run + report + stop
Re-run the completed script (detached-safe; emit-before-guard). Produce the full `.out`, and UPDATE
`_measurements/S11c_d_sympy_builder_report.md` (frozen import/reduction map + per-object which-line-computed-it +
literal §5 residuals — ⛔ no verdicts). ⚠ **If you still cannot complete ALL of it in this pass, STOP at a clean
checkpoint** (a running script) and clearly report exactly what remains, so it can be resumed — ⛔ do NOT fake,
empty-set-substitute, or type any unimplemented object. Then STOP: ⛔ no review legs, no downstream, no commits.
