**Verdict: CLEAR FOR THIS BOUNDED INGREDIENT**

This covers only the source-bound 5×5 inverse-symbol ingredient and its regular spectral density. It does not clear the outgoing kernel, the FORM roots, or A11/A12, and a clean exit from the run would not be scientific acceptance. I read only files in the packet and ran nothing.

## Substantive blocking findings

None. The checks below found no transposition, sign, Fourier, unit, join, persistence or resource error that invalidates this ingredient.

- **Transposed cofactor** (`S11c_d_reference_kernel.py:127-131`): entry `(i,j)` removes row `j` and column `i` and applies `(-1)^(i+j)`. That is the correct adjugate. The permutation sign at `:115-117` counts inversions correctly.
- **Identities** (`:229-236`): the left form is `Σ_k P_ik adj_kj − det·δ_ij` and the right form is `Σ_k adj_ik P_kj − det·δ_ij`. Both are checked with generic placeholder symbols for the matrix entries, then the exact source entries are substituted. Recording the identity domain as `Ne(det,0)` (`:239`) is honest.
- **Orientation and units**: `strong_matrix` (source-definitions.md, `ConstantEndPencil.strong_matrix`) builds column `j` from field `j`. `strong_matrix_units` gives entry `(i,j)` = row(i) − field(j). The worker matches this at `:191-201`:
  - it derives the row unit from column 0;
  - it then checks row − field − entry = 0 for all 25 entries;
  - inverse `(i,j)` = field(i) − row(j) (`:202`) is the correct unit for a map from rows to fields;
  - the kernel unit adds the `dk_n` measure (`:252-254`).
- **Fourier convention** (`:260-263`): the symbol is the coefficient of `e^{+i k_n z}` acting on the field. So `P(∂_z)·(1/2π)∫e^{ik(z−z')}P⁻¹ = δ(z−z')` needs the phase `e^{+ik_n(z−z_p)}`, with `z` the observation point and `z_p` the source. That is what the worker uses. `fourier_mass` (`EdgeReduction.__init__`) is `∫dq∫du e^{−au²+iqu} = 2π`, so dividing by it is the right normalization. The time character `e^{−iωt}` matches `wave_phase`.
- **Source joins** (`:144-173`): all routes are hash-checked before any data is loaded. Both the uniform packet and the reduction packet are joined through the checkpoint. The worker requires a resolved 5×5 strong matrix (`:175-178`) and binding by exact name for every live symbol (`:279-282`), and fails closed otherwise.
- **Independent binding** (`:286-331`): the scalar inverse comes from exact substitution into adj/det. The comparison inverse is an mpmath partial-pivot LU solve on a separate numeric evaluation of the source. The probes are fixed and a failure is not rescanned.
- **Persistence before guards**: every `Journal.op` writes its arguments and start receipt before running, and its result before any later `require`. The consumed-route records, unit operands, certificate, probe operands and posthashes are all written before the checks that use them.
- **Resources**: the worker checks cgroup, affinity, nice and thread limits before importing sympy (`:54-75`), then sets RLIMIT_AS to 2 GiB and disables core dumps. This is consistent with `s11c_guarded_run.py`.

## The outgoing prescription question

Stopping here is honest and within the approved first stage. The supplied conventions fix what a prescription would have to be, but this worker cannot build it from its approved inputs without either a new zero census or inputs it does not consume. The reasoning:

1. **What the conventions determine.** With time dependence `e^{−iωt}` and the bulk radiation branch `q_out` (c1 SHARED_PHYSICS:109-118 and :499-500: continue the branch already selected, never re-select it), the outgoing kernel is the limit `δ→0⁺` of `P0(k_n; ω+iδ, k_∥)⁻¹`. The continuation must carry `q_out` forward from its selected branch. Re-evaluating the real-axis Piecewise is not a continuation, as the implementation doc (`:60-66`) already says.
2. **At the approved point.** By the plan's hand substitution (plan:55-60), `q² = −1/25 − k_n²`. So the bulk radical has no branch point on the real `k_n` axis; its branch points are at `±i/5`. Any real-axis singularities left would be real zeros of `det P0`, which are open slab channels. `construct_end` requires open incoming channels (`I` has shape `(5,2)`), and the reference records carry `EXACT_REAL_NORMAL`, so such zeros are plausible. For a simple zero, the side on which the contour passes follows from the sign of `∂_ω det / ∂_{k_n} det`. The source ties this to the current through the identity for `∂_{k_n}𝓛` and `∂_ω𝓛` (d SHARED_PHYSICS:454).
3. **Why it cannot be done here.** That needs three things:
   - the real zero set of `det P0` at the point, with each zero shown to be simple;
   - left and right null vectors, or the signed-current records;
   - for the lossy memory-kernel pencil, a premise that no zero crosses the real axis along the `ω+iδ` path.

   The first two live in the saved REFERENCE modal records, which are not among this worker's inputs. Computing them fresh would be the root search the plan excludes (plan:94). The third is not supplied. Even with all three, the result would be a prescription at one numerical `(ω,k_∥)`, whereas this artifact keeps both symbolic.

**Justified improvement (record only, no new computation):** make the `reason` in `prescription-status` (`:342`) name this specific dependency instead of "separate coupled outgoing-contour adjudication". It should say:
- the limit is `ω+i0` with `q_out` continued from the radiation-selected branch;
- the remaining input is the real `det P0` zeros, taken from the saved REFERENCE modal records (`EXACT_REAL_NORMAL`, `SIGNED_CURRENT`) together with the ∂k/∂ω/current identity;
- a premise that no zero crosses the axis is still needed for the lossy pencil.

## Optional observations (not blocking)

1. **Shadowed alarm.** The guard's 900 s wall clock (`s11c_guarded_run.py:127`) starts before the supervisor and the worker, so it fires before the worker's own `alarm(900)` (`:74`). The worker then gets SIGTERM and writes no `failure.json`. Per-operation receipts survive, and the guard records `wall-time limit`. A slightly shorter native alarm, or a SIGTERM handler, would restore the worker-side failure record.
2. **Decoder.** `SavedCodec` (`:78-83`) allows all of `builtins` and every `sympy.*`/`numpy.*` module. That covers `eval`/`exec`-capable callables, so it prevents importing producer modules but is not a sandbox. It is acceptable for hash-pinned local artifacts and matches the existing A9 reader. If you want it tighter, use a builtins allowlist, but don't describe it as a security boundary.
3. **Probe norms mix units.** The off-diagonal entries of `PR−I`, `RP−I`, the difference norm and the condition estimate combine entries with different dimensions. So the `1e12` and `1e-40` thresholds depend on the chosen reference frame, as the implementation doc says. With 80 digits there is ample margin; the only risk is a spurious condition stop.
4. **Branch evidence.** The plan asks for "source-derived outgoing branch/phase evidence" (plan:170). Cheap evidence would be to evaluate the saved radical and Piecewise nodes at the three probes and save the active branch and the sign of `Im q`. At present only the symbolic context is saved.
5. **Consistency joins.** Adding these would be low-cost:
   - require `fourier == constantFourierMass`, or save their difference;
   - require `spec['physicalInput']` and `spec['endCheckpoint']` to be in `by_path`, as is already done for the pickles (the checks currently hold only because the manifest lists them);
   - if the checkpoint records the hash of the physical input, join it.
6. **Mutation control.** The mutation response (`:316`, `:333`) equals `−2a_ij·e_i·(row j of R)`, so it can only fail if the chosen entry evaluates to a near-zero number that is not exactly zero. It shows the residual instrument reacts to the source. It is not independent evidence.
7. **Probes near real poles.** The probes are not checked against the saved open-mode `k` values. If a probe lands near one, the condition guard stops the job with no rescan. That is correct under the stopping rule, but it is a known way the run could fail.
