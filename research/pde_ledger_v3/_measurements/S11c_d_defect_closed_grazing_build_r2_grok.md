I'll review this closed-grazing build from the fixed source packet only, starting with the packet index, plan, worker, and inputs, then the method and the saved operands those joins depend on.The packet index is a source-only both-face grazing certificate with no peer reports. Next I’ll read the plan, method, inputs, and worker before checking the saved operands those joins name.The worker’s certificates are explicit enough to check against saved operands. I’ll read the helper surface, launcher, and the JSON fields those joins actually use.**CLEAR FOR THIS BOUNDED CLOSED-GRAZING BUILD**

The worker restores the saved both-face nongrazing objects and derives only the closed grazing certificates. I found no substantive mathematics or validation blocker, and no concrete runtime or tooling fault. Method support and this reading do not accept a runtime result. No READY gate or scientific run is present.

## Symbol transport and complex-frequency contact

`worker.py` lines 126–160 take the actual symbols out of the saved coefficient and factor with `one_symbol`, then replace them simultaneously. The old frequency is `Symbol('omega', real=True)`. The old depths are the saved `Symbol('q_i'|'q_h'|'q_s'|'q_o', nonzero=True)` objects from `raw/physical-factorization-input.json` and `raw/physical-factorization-return.json`. They are sent to fresh unrestricted symbols `grazing_qi`, `grazing_qh`, `grazing_qs`, `grazing_qo`, and `grazing_unrestricted_frequency`. Line 143 requires the new frequency to have `is_real is None` and the old frequency to be real. Argument identity is kept: `H` and `increment_transfer` both become the new real transfer, and `increment_difference` becomes `l-k`, with `k` left as the saved real momentum.

The complex-frequency proof uses those transported saved expressions. It clears `C-tB` with the two dispersion rules `qh**2 = qi**2-t(t+2k)` and `qs**2 = qo**2+t(2l-t)`, requires powers 0, 1, or 2, and saves the numerator before the zero guard. Contact is the direct substitution `t=0`, `qh=qi`, `qs=qo` on each transported face coefficient. The lower face is the unique saved mode grade `(1,1)` from `raw/lower-boundary-operands.json` zipped with `raw/lower-boundary-return.json` (lines 88–94 and 154–158), which is the same mixed coefficient as the upper factor. No boundary constructor is called.

## Closed density, native law, and the qs route

Lines 161–198 build

`Bc = -i ω ρ_m/[(qi+β)(qo+β)] [k(2l-t)/(qs+qo) + k(t+2k) qi/(qh(qh+qi)) + qi²/qh]`

and `G = (WL/(4i)) A(t) A(l-k-t) Bc`, with `β = Λ_A0 ω/(ρ_m(1-iωτ_A))`. `Bu*R` is exactly `Bc`. Each face joins its saved mixed coefficient at frequency 3, its raw kernel to `(pref·B)` at frequency 3, its closure factor to `R`, its zero-grade trace to 1, and its normal jet to `i f qo`. The saved reference

`qi qo (60+91i) / [qi qo(60+91i) + (qi+qo)(9+30i) + 9i]`

equals `qi qo/[(qi+β₃)(qo+β₃)]` because `(60+91i)β₃ = 9+30i` and `(60+91i)β₃² = 9i` at `β₃=(30+9i)/109`. The plus jet is `i qo` times that reference; the minus jet in `raw/minus-closed-before-cancel.json` is the minus sign. Height signs are `+η w/2` and `-η w/2`. The ordered weight join at lines 168–170 matches `A(t)=(5/2) t/sinh(5π t)` and `jet = 5 A(Q-t)` against `raw/raw-ordered-before-cancel.json`, with plain middle measure true in `raw/edge-delta-reduction.json`.

Native roots at lines 199–218 are the four saved `Piecewise` sheets in `raw/physical-sheet.json`, transported onto `k`, `k+t`, `l-t`, and `l`. Positive and decaying branches match `sqrt(9/cs²-1/20-p²)` and `i sqrt(-(that))`, strict inequalities match the radicand, and the zero branch remains `nan` with condition `true`.

`qs` stays the reflected route `q(l-t)`. Lines 258–259 prove `(l-t)²-(k+t)²=(l+k)(l-k-2t)`, hence the squares agree at `l=-k`. The collision record then treats `qs` and `qh` as the same function while the integrable envelope remains a sum. Same-momentum endpoints are `[0,-2k,2k,0]`; opposite-momentum endpoints are `[0,-2k,0,-2k]`.

## Bounds and the L1 argument

The checked constants match the stated domain `cs∈[1,2]`, `|k|,|l|≤3`, `δ∈[0,1/10]`:

- `β_min=3000/11101` is the minimum real part at `D_max=11101/10000`.
- `κ_min²=879/400` is the exact corner `cs=2`, `δ=1/10`.
- `|Ω|≤4` and `|qi|,|qo|≤5` are the checked loose comparisons at lines 240.
- The bracket bound `(18+3|t|)/|qs|+(43+3|t|)/|qh|` follows from `|k|≤3`, `|l|≤3`, `|qi|≤5`, and `|u+v|≥max(|u|,|v|)` on the closed first quadrant.
- For `|t|≥12` and endpoints in `[-6,6]`, `61+6|t|≤12|t|`, distance at least `|t|/2`, and `|Q-t|∈[|t|/2,3|t|/2]`. The constant `C=(WL/4)(3L²/2)(4ρ_m/β_min²)(24√2/√κ_min)` and the displayed primitive of `t³ e^{-10π t}` are the tail in lines 241–250.
- `1-e^{-10π}>1/2` follows from `π>3` and `e^x≥1+x`. The local integral of `|t-a|^{-1/2}` over a set of measure `m` is at most `2√(2m)`.

Those identities, the uniform majorant, and the written absolute-continuity-plus-tail argument are sufficient for L1 convergence of this kernel on that domain, including `δ→0` and `qi` or `qo→0`, with the contact term absent because `C=tB` and the `t=0` value are zero before the limit. This is an analytic argument attached to exact identities. The symbols `grazing_effective_speed` and `grazing_delta` are only positive and nonnegative; `cs≤2` and `δ≤1/10` are the written domain of those gap identities.

Real-frequency source and row coefficients stay at frequency 3. Only `β(Ω)`, the four depths, and the explicit `ω` in `Bc` are continued.

## Formal jets, normalization, and Fourier scope

The retained D-slot formula in `source/raw-worker.py` lines 563–600 is `η σ D · (reference or jet) · source`, then the mixed `η,σ` derivative at zero. Here `η` and `σ` are `eta_bg` and `sigma_W`, because the saved zero-grade source is that substitution and no longer contains those symbols. Native jet children in the THETA and `E_W` censuses already carry one `eta_bg`, so that grade derivative vanishes. The surviving pieces are the `delta_p` children.

For THETA that factor is `-i Λ_A0 ε / (ω ρ_m τ_A + i ρ_m)` at the bound numbers. For `E_W`, `Lambda_X_0` is 0 in `physical-input.json`, and the surviving factor is `W_0 ε/2`. Dividing by the saved reference and `epsilon_shape` (`worker.py` lines 278–297) leaves a t-independent linear form in the actual source jets `e_W`, `e_W_t`, `e_W_d*`, `theta`, `theta_d*`, and `u_*_d*`, with finite constant coefficients. Depth, transfer, and old frequency symbols are rejected. Both raw D-slot kernels are joined to the saved ordered kernel, and the raw row density is joined to `ε (m₊+m₋) G` at frequency 3. U0–U2 are the saved zeros required by an empty pressure census.

`fourier-convention.json` is not a runtime operand. Lines 261–263 and 294–297 record that no general or lower-face Fourier map is certified. The two modal matches are dispersion checks only. This is enough for an L1 statement about the kernel multiplied by those fixed finite formal-jet coefficients. It does not bind a physical incoming mode.

## Controls

The three new controls are computed from the new expressions:

- Missing closure: the residue of `qi B` at `qi=0` is `-i ω ρ_m k(2l-t)/(qo(qs+qo))`. At the left match `cs=√(3/2)`, `k=√(595)/10`, `l=0`, `t=1/10`, it reduces to `i · 3/100/(qs+qo)` with positive real quotient.
- Wrong sheet: the saved slope and output `Piecewise` values are joined to `√(κ²-t²)` and `κ`, the slope root is negated, and the closed limiting density moves by a nonzero exact fraction while the negated root fails the nonnegative test.
- Lower jet: both saved jet/reference ratios are joined to `±i qo` before the plus ratio is substituted for the minus ratio. Profile factors `A(1/10)` and `A(κ+1/10)` are positive, and numerator and denominator each need a finite signed nonzero component.

Old `finish/` responses are not loaded or rerun.

## Execution, gate, and launcher

Field names, grades, epsilon, D-slot names, traces, and Piecewise branches match the saved JSON. `decode` rebuilds `srepr` with the helper’s `Str` hook; booleans such as `plainReducedMiddleMeasure` stay Python `True`. Helpers are the ten named functions from `source/raw-worker.py`, exec’d only after the gate, with SymPy imported only after `containment()` (`worker.py` lines 367–374; `tooling-tests.py` lines 54–59). Substitutions that must not capture `k` inside `l-k` are simultaneous. `Journal.zero` writes the raw residual before `require`. Failures land in `failure.json`, posthashes cover sources and byte copies, and `scientificAcceptance` stays false.

`verify_gate` (lines 41–59) requires the method record’s joint clearance, a build record with both literal verdicts `CLEAR FOR THIS BOUNDED CLOSED-GRAZING BUILD`, matching worker, manifest, guard, and supervisor hashes, authority with retry off, one run, and `durationLimits is None`. `launch.py` lines 41–98 arms the completion hook before the coordinator starts the guarded command. The command is the 4 GiB pooled guard, zero swap, one CPU, 32 tasks, 4 GiB host reserve, around `runtime-source/supervisor.py`, with `RuntimeMaxUSec=infinity`. The worker does not call the old boundary, factorization, profile, or row constructors.

The 30-second `select` and the five-second watcher poll are startup arming checks. The scientific child has no wall-clock, CPU, or inactivity deadline. Low-memory stop remains the guard’s host reserve.

## Optional wording

`Journal.sinh_zero` records “intrinsic odd symmetry,” while the helper only expands arguments inside `sinh`. The equality test is the stricter one. The emitted flags `tIndependent` and `csIndependent` sit beside the actual symbol requires; those requires are the control.

## Inspection limits and the strongest future claim

I did not execute SymPy, restore science, or parse the megabyte `THETA_BALANCE` and `E_W_BALANCE` expression strings byte by byte. Those joins were assessed from the producer in `source/raw-worker.py`, the executed censuses, the U-row schema, and the closure, reference, and source-zero-grade files. Tooling tests were read, not run. The method-record JSON is pinned by hash and is not in this packet.

If a later contained run saves zero residuals and the failure flag stays unset, its output can support exact both-face joins of this direct `(1,1)` kernel, the unrestricted-frequency contact identity, and an L1 limit for that kernel times the certified real-frequency-3 formal `(0,0)` jet coefficients, on strict-rest `LAB_HELD/RHO4_CONSTANT`, saved edges, `cs∈[1,2]`, `k,l∈[-3,3]`, and finite positive permeability and memory. That future output still would not be an evaluated integral, a finite solve, a full operator, a calibrated or draining model, a defect sweep, or a loss result.
