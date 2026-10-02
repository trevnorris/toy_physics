**NEEDS REVISION**

I read the method, guide, worker, manifest and launcher/guard/supervisor/hook sources, plus the saved response, source, direct, native-route and census files that the worker cites. I did not run anything. I could not reach the saved files the worker never loads or only partly joins. These runtime items are left to the future run: the 4 GiB cgroup memory ceiling on the large `J.zero` residuals, `bind()` hitting an ambiguous symbol assumption, and the key-schema match of the saved `raw/*` files.

## Scientific / validation blocker

**B1. The source-side chemical, live-density and raw-to-bound joins are not performed (`worker.py:342-351`, `:245`).**
- **What the worker checks:** it joins only the velocity amplitude: `bind(nativeVelocity/eps)` against `velocityAmplitude` (`:345`). It also checks the parts sum and the coefficient×amplitude products, and that `combined` equals `saved['full']`.
- **What it never loads:** `consumer/native-chemical-amplitude.json`, `chemical-amplitude-domain.json`, `inherited-source-normalization.json` and `selected-fourier-contraction.json`.
- **What it never uses:** `inp['raw']`. That is the only expression still carrying `s11cc1_mu_*`, `s11cc1_V_*` and `rho_br_bg_rho4_constant`. No `J.zero` takes it through the stage-2 identifications plus the live-density map to `combined`.
- **What stands in for the density join:** line 350 only checks that `rho_br_bg_rho4_constant` is absent from the bound result. Line 245 asserts `densityRestored: True` as a literal with no evidence behind it.
- **Consequence:** the chemical amplitude `M` and the live-density binding reach every address only through the saved bound `S_f`. For both faces, method §3 ("M = native chemical source/epsilon", density-field binding, per-face joins) is therefore not evidenced by this worker. Both-face source equality (`:369`) compares two saved objects, not two joins.
- **Minimum correction:** for each face, add `J.zero` joins for the following.
  - `chemicalAmplitude` against the native chemical amplitude file.
  - `bind(raw.subs(stage2 identifications, density map))` against `combined`.
  - `inherited-source-normalization` against `velocityCoefficient`.
  
  Load them through `load()` so they appear in `used`. Alternatively, record them as inherited-only and drop the `densityRestored: True` wording.

## Weaker validation points (not blockers)

- **V1. Control baselines leave out the consumer scalar (`:872-908`).** The slope and direct controls form `baseline = response × source` and never multiply in the consumer field or its transform tag. `omit-consumer10` and `omit-source10` therefore have `movement = ∓baseline`, which only shows the baseline is nonzero. They are not sensitive to the addressed consumer coefficient. This does not contradict the stated formal-tag scope. Optional fix: carry the consumer transform tag coefficient as a named formal factor in the baseline.
- **V2. Hard-coded literals.** The flat check `3/10/(3/2+beta3)` (`:884`) and the `10/7` control point are not read from saved data. I checked `mu = ω/10` against the saved flat coefficient and the three q values against `q² = 25/4 − p²`, and both are correct. Optional fix: derive them from the saved objects.
- **V3. Zero-profile reductions of the whole tags are declarations.** Setting the `Hwhole`, `Jwhole` and `Dwhole` tags to zero (`:626-631`) is a formal profile-degree map, labelled as such. That is adequate for this inventory. The absent nonzero constant-height reduction is correctly recorded as `NOT_COMPUTED`.

## Tooling

- **T1. The launch command is not pinned to the gate's guard and supervisor paths.** `launcher.verify()` (`launcher.py:34-39`) compares `gate['command']` to a literal list. `worker.verify_gate` hashes `gate['sharedGuard']` and `gate['supervisor']`, but nothing asserts those paths equal `scripts/s11c_guarded_run.py` and `S11c_d_end_normalization_run.py` in the command. Fix: assert the two path equalities.
- **T2. `except BaseException` (`:947`)** also records `SystemExit` and `KeyboardInterrupt` as `FAILED_PRESERVED`. That is acceptable, since it still exits 1 and writes `failure.json`.

## What checked out

- **Typed factors:** `Rprod`, `E` and `Dwhole` are kept as three separate objects. There is no density-minus-factor comparison and no extra resolvent on `Dtag`.
- **Isolated-factor and full-density pairs:** the isolated-factor pairs (`:397-411`) and the full-density original-input/canonical-return pairs (`:388-396`, `:484-493`) are joined separately and not replayed. I checked `beta(ω=3) = (30+9i)/109` against the saved factor.
- **Independent-D adapter:** its argument order matches the native `kernel_apply` signature. The native `reference`, `extension`, `jet_transfer`, `normal_jet` and `reference_pressure` assignments are exec'd in dependency order. I confirmed they yield `i·f·qo·P` and a retained `P`. The excluded (2,1) term equals `−H·i·f·qo·P`, with `H` a symbol bound to `eta_bg·hat(l−k)`.
- **Trace-matrix action:** `T·Δ = Δ` is equivalent to `T⁻¹Δ = Δ`, and the new `Δ` leaves the [0,1]/[1,2] routing unchanged.
- **Row census and ablation:** broad `delta_p`/`d_w_` text-count against literal-AST-hit cross-check; unsupported constructor spellings are refused. The four slot ablations, affine reconstruction and child-hash join against the saved census are present.
- **Quotient grades:** regular quotient recurrence from the full cancelled fraction, with the excluded remainder preserved.
- **Address coverage and row substitution:** 320 coverage entries (5 rows × 2 faces × 2 slots × 16 triples), zeros retained. The placeholder substitution into the actual affine rows is keyed identically to the address tokens. Comparing it with the address sum does test the grade expansion, though the comparison is partly tautological for the consumer and source factors.
- **Control support:** the off-delta eligibility certificate is an exact tanh-polynomial test, with no Fourier-value claim. The control support points satisfy the stated q values at `cs² = 10/7`. Positive-branch checks and the `H` bound variable are present.
- **Scope and containment:** historical (2,1) remainders, no-deadline containment, hook-first handshake, write-before-guard journaling and posthash records are sound. No old producer, triangular solve or response function is called.

## Coverage limits and runtime evidence still needed

**Limits**
- Source coefficient dimensions after binding are required units only. The native jet-dimension rule is read, not executed.
- Every control movement is a formal tag coefficient, not a field value.
- Global momentum, test-space and grazing composition stay UNRESOLVED outside `[-3,3]`.

**Runtime evidence needed**
- The run must complete the many large `J.zero` residuals within the 4 GiB limit. Any fail-closed refusal (e.g. an unexpected jet name, a symbolic `den0`, or a symbol assumption mismatch in a map) must be inspected.
- The post-hashes must come back clean.
- The fixed B1 joins must also be present in the saved output.

This is a source-only verdict and does not accept any runtime restoration or result.