**NEEDS REVISION**

I read the method, the evidence guide and the first 764 lines of the 980-line packet index. I did not open the last 216 index lines. I inspected the native c2 routes, the reference census, the final-slot routing, the pressure-omission record, the THETA/U1 pressure censuses, the source grade split, the physical-domain and speed-reuse records, and the reference/consumer worker sources. I did not open the `direct/*` certificates, the `continuation/*` control records, the closure-operand details, or the remaining `raw/` and `consumer/` files, so statements about those rest on the workers' own summaries. I executed nothing.

## What checks out

- **Native epsilon:** The pressure-slot row coefficients carry `epsilon_shape` once (`native-pressure-consumer-omission.json`). The sources are per unit amplitude (`consumer-worker.py:306-317`). Counting epsilon once in the final increment is consistent.
- **Normal multiplier:**
  - In `build_face` (`native-c2.py:566-572`) the normal jet is `d/dN(exp(i f qo (N-ref)) * reference_matrix[.,.])`.
  - The reference matrices do not depend on N, so the jet is `i f q(l)` times the reference entry for every entry, with `qo` the output depth.
  - `jetKernels` in `retained-response-census.json` matches this for flat, height, slope and mixed, on both faces (`+i qo` plus, `-i qo` minus).
- **Flat delta:**
  - The native diagonal `reference_matrix[0,0]` uses `qo` in prefactor and pole. The saved pressure flat is the equivalent `qi` form. The method's `delta(l-k)` identification at `q(l)` is correct.
  - `p0` is the already delta-reduced single-momentum route (`native-c2.py:468-470`). The method's warning against a second delta is right.
- **Fourier convention:** `native-fourier-contract.json` gives profile forward power -3, source inverse power -3 and source forward power 0. Working this through, the plain `dl dk` composition with a normalized transform and no extra 2π is correct, and a constant multiplier reduces by its delta alone.
- **Grade triples:** The counts 9, 3, 3, 1 (16 in total) are right.
- **Consumer and source coverage:**
  - Only THETA_BALANCE and E_W_BALANCE carry the four pressure slots. The pressure-slot coefficient is grade (0,0) and the normal-slot coefficient is grade (1,0) (`η·w1_profile`). U0/U1/U2 have none.
  - The saved `plus-source-grade-split.json` source has denominator `(1+η w1_profile)`, which is finite and nonzero at η=0. It contains `η·w1_profile·e_W_d1d1`, so a source10 spatial jet exists and the `p→k` control is plausible.
- **Frequency:** The omission record is bound at ω=3 (the 109 denominators). `effective-speed-reuse-domain.json` reports `sourceSymbolsCs:false` with saved and new frequency both 3.
- **Control arithmetic:** At cs=√(10/7) and ω=3, `q(3/2)=2`, `q(2)=3/2` and `q(30/13)=25/26` all check, including the 1/25+1/100 transverse term.
- **(2,1) remainder:** The old `selectedIncrement ∝ η²σ·D` and `mixedPerSource=ablatedRowPerD=0` match the method's account that the normal consumer begins at grade (1,0).
- **Scope:** The explicit UNRESOLVED treatment of momenta outside [-3,3] is honest. I found no hidden cutoff or global claim.

## Blocker (scientific / claim)

**B1. The direct (1,1) addend is not in the native c2 route. The method counts it as one of the "retained" native pieces, and no packet observation routes it natively.**

- **Where:**
  - Method §4, "Direct F_(1,1) … enters once", and §6, "direct(1,1)/source00/consumer00 must match …".
  - Native side: `native-c2.py:409-417` builds `z_three` with a hard 0 at [0,2] and takes `second` only from the iterated product.
  - Saved side: `reference-worker.py:235,239,256,297` (`F[0,2].subs({D:0,…})`, `noDirectInNativeSecond`, `directNotNativeSecond:true`).
- **Evidence gap:**
  - The final-slot join (`*-final-native-slot-routing.json`, `reference-worker.py:276-281`) is tested only on columns 1 and 2 with D=0. It records `sourceCompositionPerformed:false`.
  - `jetKernels.mixed` (`reference-worker.py:298`) multiplies `normal*(directTag+…)` by hand.
  - Native `kernel_apply` would send any [0,2] entry through `p2`, which integrates a spurious middle variable. The method correctly forbids that, but the replacement `p1`-type route for the direct is therefore not native.
- **Consequence:** Anyone implementing the inventory could read "native retained pressure part" as covering the direct. The two-face direct block, its `i f q(l)` jet and its immunity to trace subtraction would then be treated as natively verified when they are inherited or hand-formed.
- **Minimum correction:**
  1. State in §4 that the composed object is native c2 plus an inherited, non-native direct addend. Give its provenance chain: the raw-closure `raw_direct` tag, then `rawKernelPlus` in `raw/*-retained-increment.json`, then the restored closed density in `direct/closed-density.json`.
  2. Make the direct a separate off-diagonal address type with plain `dl dk` and no middle integral.
  3. Require a formal join, with an independent placeholder D, of the actual `build_face` `reference_pressure`/`normal_jet` assignment text for both faces. It must show that the direct's pressure part is `R(qo) D R(qi)`, unchanged by the `T01` subtraction, and that its jet part is `i f qo` times that.
  4. Label the result an algebraic routing check, not native execution.

## Non-blocking corrections (tooling and documentation)

- **Pressure census:** The saved censuses count only four exact names. The rerun must scan every top-level child for any atom containing `delta_p` or `d_w_`, so the "report and stop" rule can actually fire.
- **Regulated frequency:** The response objects run at Ω=3+iδ, while the source, consumer and `-iω` factors are saved at real ω=3. State whether any δ>0 route continues them analytically. They are rational in ω. Say explicitly that none is claimed here.
- **Source grade split:** The saved split is zero-grade versus lumped higher grades. The per-grade (1,0), (0,1), (1,1) pieces and the discarded (2,0)/(0,2) must be computed at runtime from the saved `full` expression. This is feasible without replay.
- **Lower face:** The velocity-normalization join exists only for the upper face (`consumer-worker.py:308`). The pressure-omission record removes only `delta_p_plus`. Lower-face source and omission equality are runtime obligations, as the method already says.
- **Control applicability:** The source10 and consumer10 coefficients for the nonzero-row requirement must be read from the runtime grade split. Record exact absence if one is missing.

## What a successful bounded instrument would establish

- A complete, hash-addressed inventory of every ordered triple (explicit zeros included) for both faces and slots, with epsilon count one.
- Formal Fourier routes with distinct p, k, l, r and the flat support rule.
- Argument-bearing whole-tag signatures H, Jwhole and Dwhole with the saved hashes.
- Exact reconstruction of the affine rows and of the saved total-mixed tags as sum checks.
- Control sensitivities that are labelled formal where they pass through tags.

It would not establish a value, a convergence result, a grazing limit or any finite-solver readiness.

## Unresolved before a finite near-unity defect pilot

- Global-momentum domain and test-space for the tanh-profile transforms times PV/contact responses, with no cutoff.
- Closure of the composed (1,1) expressions outside [-3,3].
- The reflected-depth sheet for `q(l-t)` away from the certified domain.
- Evaluation and quadrature of the contact, PV and J, D integrals.
- Lower-face source equality.
- The regulated-frequency convention.

## Runtime evidence needed

- Pinned constructor-text restoration of the sources, the slab rows and the density join, with a full pressure-atom scan.
- Per-grade source and consumer splits, with finiteness certificates at η=σ=0.
- Exact reconstruction of the affine rows.
- The B1 direct-routing join.
- Saved-hash joins for the response pieces: flat, height, slope, mixed pieces and the total.
- Applicability records for each control, with exact movements.