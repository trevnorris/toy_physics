**CLEAR FOR THIS BOUNDED RAW-INCREMENT BUILD**

This is a source review only. I ran nothing, and I can't prove runtime restoration or cost from source. I found no finding that changes the result or blocks execution. The tooling items below are fail-closed risks, not wrong-answer risks.

## What I checked against source

- **Lower boundary** (`worker.py:352-376`).
  - I re-derived the four amplitude equations by hand, and they reproduce the saved `A00`, `A10`, `A01`, `A11` and the saved `mixed`.
  - The face factors cancel identically in the code (`:359-362`). The lower face enters only through the native normal derivative of `-1`, `:345-350`.
  - The native minus trace gives `∂trace/∂jet = -W₀·η·w1/2` (`native-selected.json:385`). The native background normal is `(-σ·w1_d/2, …, -1)` (`:260`).
  - The native `dtn_first_kernel` (`native-c1.py:606`) is face-agnostic. It reproduces the height, slope and flat coefficients (`:385-391`).
- **Factorization** (`:393-410`). I verified by hand that `H·(qi-qh)(qi+qh)/(…)` and `H·(qs-qo)(qs+qo)/(…)` reproduce `N` exactly. So `C = H·B` holds modulo both dispersion identities, and the contact term is `N(H=0, qh=qi, qs=qo) = 0`.
- **Saved selected integrand** (`:426-427`). By hand, at `k=0` and `Wωρ=0.3` the raw kernel gives `-5√3/16`, matching the saved integrand.
- **Closure.**
  - `∂/∂D (I+aZ)⁻¹Z |₀₂ = M⁻¹₀₀M⁻¹₂₂`, which equals `Rprod` (`:508`).
  - The reference factor `Rprod/trace0` matches native `reference_pressure_kernels` and `build_face` (`native-c2.py:478-592`): output-leg trace, with the normal jet equal to `i·f·qo` times the reference.
  - It also matches the saved upper factor and jet at `qpoint`. I checked the algebra.
- **Slot census and grades.**
  - Native `selectedPressureJetChildren` are per-slot `Mul` terms.
  - The jet-slot terms are `η²σ`, so `diff(·, η, σ)` at 0 drops them (`:557`). The retained increment is carried only by the δp slots.
  - Both faces enter before grade extraction. THETA and E_W both have nonzero `δp_minus` coefficients, so the lower omission and doubling controls respond.
- **Fourier contract** (`fourier_contract`). It matches native `c2` lines 129-131, 461-470 and 850-861. The edge reduction `(2π)⁻³·(2π)² = (2π)⁻¹` matches the saved normalized hat `A(0) = 1/2π`. The reduced-kernel dimension join gives `(-2,-1,1)` with `ρ_m = M/L⁴`.
- **Execution safety.**
  - Containment comes before the SymPy import.
  - Inputs and returns are written with `'x'` mode.
  - There are no duplicate emit names.
  - The gate joins worker, manifest, guard and supervisor hashes.
  - There are no timers or retries.

## Math, method and claim notes (none blocking)

1. **`:378` mirror check is near-tautological.** The independence of the lower mixed term rests on the native normal and trace joins (`:342-351`) and the native linear joins. The `(1,1)` term itself has no native join. Say so in downstream claims.
2. **Jet and normal-jet factors are emitted, not consumed.**
   - The `ht0` term in `trace0` is zero at zero grade (`:514-515`).
   - The jet slots drop at `(1,1)`.
   - `lower-reference-jet-sign` (`:597`) is correctly only a trace control. Don't cite it as an exercised jet consumer.
3. **Frequency is not required to match.** `saved['binding']['frequency']` is only emitted (`:307`). The bound and combined-source joins (`:494`, `:551`) enforce it implicitly. Optionally add `require(saved['binding']['frequency'] == freq)`.

## Tooling findings (non-blocking)

- **`:632` posthashes cover only `manifest['sourcePins']`.** Add the manifest (`args.inputs`), the gate, `g['buildReviewRecord']` and the copied `saved-operands/*`.
- **`Journal.zero` (`:232`) normalizes with `cancel(together())` only.** Upstream used `factor(cancel(simplify()))` (`trace-worker.py:117`). Joins with `√3` or `I` and different rationalizations (`:324`, `:526-527`, `:563`) could leave a non-literal-zero residual, which is fail-closed. Smallest fix: fall back to `simplify` or `radsimp` before `require`. The raw residual is already persisted.
- **`J.nonzero` (`:241`, `:579-581`, `:594`, `:596`) needs `is_zero is False` on sums of four unrelated radicals.** It may return `None` and abort after a long run. With one authorized run and no retry, use `abs(sp.N(movement, 50)) > 0` instead.
- **Heavy `cancel` calls run before their operands are persisted.** Those at `:416`, `:520` and `:523` happen before any emit. Emit the operands first.
- **`:409` emits `'lowerContact': 0` as a literal.** `:408` derives it, so optionally emit the computed value.
- **Native row files are pinned under `build-review/packet/native-rows/`** (`input-manifest.json:141-145`). They are hash-pinned and joined to `expandedConstructorSha256` (`:544`), so this is sound. Make sure that directory survives until the gated run.

## Not provable from source

- SymPy `decode` and `one_symbol` uniqueness at runtime. Any assumption-variant duplicate of a name, such as `W_0` plain versus positive, fails closed.
- The cost of `cancel` over `Piecewise`/`sinh` expressions (`:523`, `:563`) under the 4 GiB limit. A `MemoryError` would be caught and preserved by `main`.
- The Fourier and edge reduction is checked only syntactically, plus my hand arithmetic. No integral is executed.

These are source findings only, not an executed ablation, runtime restoration proof, or result acceptance.