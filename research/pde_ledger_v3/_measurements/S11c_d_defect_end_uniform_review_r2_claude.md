## Assessment: corrected weak ends vs the selected uniform subspace

**Packet coverage actually read.**
- **Fully read:** `method.md`, `evidence-guide.md`, `review-prompt.md`, `packet-index.json` and `ends/method.md`.
- **Partly read:**
  - `source/uniform-worker.py`: lines 205–359 (`selected_lift`, `EndBinding`, `native_leg`, `retained`, `on_wave`, `end_restriction`) and 925–1020 (main driver, speed schedule, result scope).
  - `source/engine.py`: lines 3270–3340 (field ansatz, leg maps, power maps, `CLOSED_PENCIL_LEGS`).
  - `source/frequency-source.py`: lines 160–210 (`end_sources`).
  - `uniform/opaque-interface-receipts.json`: first 120 lines.
  - `ends/symbols.json`: cell 0 in full, plus a key and row/field census (2,200 row/field/side tags, so 200 cells with 4 tags per cell appears consistent).
  - `uniform/selected-result.json`: status grep only.
  - Row and field label definitions in the weak-end workers.
- **Not read:** `controls/*` (including the 28-record evidence chain), the other `ends/*.json`, `native/*`, `uniform/plan.md`, `physical-input.json`, `current-source.py`, `pairing-check.py`, `c1.py`, `c2.py`, the control-continuation sources, most of `engine.py`, and most of `symbols.json`.
- **Limitation:** I decoded none of the opaque blobs. Statements about the old operands rest on what the producing source requires. I did not check any hash.

### Blocking defects
None found. The method's gates are mathematically consistent, and its stop conditions (`SOURCE_MAP_UNRESOLVED`, unresolved grade attribution, unresolved grazing) cover the failure modes I could construct. I checked the following.

- **Identity `R−A=I_old`.** It follows algebraically from the definitions of A and R.
- **Remainder identity `A=ΣηᵃσᵇJ_ab−H_phys`.**
  - It needs only the exact definitions of H and J.
  - It needs `P_raw(origin)=P_old`.
  - It needs a grade-independent L. No Taylor convergence is required, because H is the exact remainder.
- **Grade-independence of L is checkable without calling a producer.** `selected_lift` (`uniform-worker.py:234-238`) binds the restored `uniform['curl']` only from `params` and the normal momentum `k`. Any `eta_bg` or `sigma_W` symbol there would raise a `KeyError`. The worker can verify this from the free symbols of the restored curl.
- **The retained grade set matches the old convention.** `B.retained` (`:305-310`) is multilinear, a,b∈{0,1}, which matches G.
- **Branch is consistent.**
  - `CLOSED_PENCIL_LEGS[0]` carries `+ω` and `+k`, with `qright` and `exp(+i q·depth)`. This is the positive-carrier, outgoing-branch leg (`engine.py:3284-3321`, `3330`).
  - Leg 1 flips the signs of `ω`, `k` and `q`. The method's ban on interchanging the legs is correct.
  - At real ω=3 the first-quadrant limit of the weak-end root gives positive real q or positive imaginary q, as the method states.
  - `q+β` never vanishes on the closed outgoing branch. For q≥0 real, Im(q+β)=9/109. For q=it with t≥0, Re(q+β)=30/109. The worker should record this one-line certificate rather than assert it.
- **Inherited correspondence targets are grazing points.**
  - LEFT cs²=3/2, p²=595/100 gives 9/cs²−p²−1/20 = 6−5.95−0.05 = 0.
  - RIGHT cs²=150/101, p²=601/100 gives 6.06−6.01−0.05 = 0.
  - So the method is right that these targets go through the closed `R0` route and the four limit operands, not through raw subtraction.

### 1. Joins (row, field, unit, epsilon, Fourier, memory, depth)
The method is adequate. Concretely:
- **Labels.** Weak rows are `U0,U1,U2,THETA_BALANCE,E_W_BALANCE` and fields are `u_1,u_2,u_3,theta,e_W` (`full-weak-worker.py:23-24`, `ends-worker.py:19-20`). The old lift field order is `u1,u2,u3,theta,eW` (`uniform-worker.py:223`). Matching by position is only a hypothesis. The method correctly requires a source-derived join, and `U0` versus `u_1` is exactly where it could go wrong.
- **Scalar factors.**
  - `real_field=ε(plus·phase+minus·phase)/2` (`engine.py:3279-3282`) puts a 1/2 in the trial normalization. Both trial and residual fields carry ε there, so some old operands are ε² objects.
  - The E cells are the coefficient of one ε (`epsilonPower` 1, or 0 for exact zeros).
  - A stray scalar such as 2, 1/2 or ε-power from this is the most likely legitimate map. It must be derived from the native definitions and fixed before looking at any residual, not fitted.
  - The method already forbids fitting and requires `SOURCE_MAP_UNRESOLVED`. Naming these candidates in the worker spec would help.
- **Carrier and pairing.** Derivatives stay on the trial, the pairing is `2π∫v̂(−p)ᵀE û(p)`, and the carrier is exp(ipz). Then E is the plateau limit of the strong operator symbol, so no integration by parts is needed. The method's refusal to use `plus_power_map`, `weak_matrix` or the Euclidean Gram is correct. `plus_power_map` is a ∂²-of-the-source-coefficient object (`engine.py:3315-3318`).
- **Grade-variable identity.**
  - The method requires the coordinate maps to be grade-independent.
  - It does not separately name the identity of the expansion parameters themselves. E's η is the height amplitude with `H_-=0, H_+=W/2`. The old `eta_bg` and `sigma_W` are bound through `profileEndpoints` and `origin`.
  - A constant rescale of η would show up as a spurious `J_10` and also change the origin value.
  - The method's instruction to join endpoints and origin bindings, and its warning about W/2 versus W, mostly covers this. I recommend making it an explicit gate.

### 2. Assembly, attribution and outcomes
The two-axis outcome design is sound (retained, finite, and unavailable are kept separate).

**Non-blocking, but real.**
- **`I_old` is zero only on the wave.** `end_restriction` saves both `rawResidual` and `onWaveResidual`, and only `on_wave(...)` is required to be zero (`:353-356`). The raw `P L − L D` need not be identically zero.
  - The consistency join `R−A−I_old` should therefore be saved in raw form using the restored raw `I_old`.
  - The ZERO/NONZERO statuses should be on-wave reductions.
  - If the worker assumes raw `I_old=0`, it can report a spurious inconsistency.
- **A=R=0 is a hard target.**
  - H contains the exact excluded grades. For a rational `P_raw` it is generically nonzero.
  - The likeliest honest results are therefore `FINITE_TRUNCATION_DIFFERENCE` or `RETAINED_MISMATCH`.
  - Section 4's "if A=R=0" branch is correctly conditional. The report should treat J_ab as the primary axis and not as a fallback.
- **Rational regularity.** Regularity at η=σ=0 is checked on the symbolic denominator D(0,0). That expression can still vanish on a subvariety of (p,q,cs). Each excluded locus should be reported per coefficient, as the method already says.

### 3. Wave-surface certificates
- **Wave-surface reduction.** It is sound: `a+bq` with a,b q-free vanishes as a function on the surface iff a=b=0. That only holds if q is not in the coefficient field, so the certificate is for the symbolic surface.
- **NONZERO_CERTIFIED.** It needs a witness at a concrete point with the actual outgoing root, because at a rational point with a perfect-square radicand `a+bq` can vanish with a,b≠0. The method already says this.
- **Symbolic versus point/path coverage.**
  - The old worker holds `cs` and `k` as free symbols (`:930-932`), so symbolic (p,cs) certificates are plausible.
  - Whether the saved operand is symbolic or point-bound can be decided only after hash-checked decoding.
  - Anything other than a symbolic (p,cs) operand must be demoted to point coverage.
- **`R0` and `D_old` at grazing.**
  - The restriction check raises if D contains q (`:358`), so D is a q-free rational function of (p,cs). `R0` is evaluable when its denominators do not vanish at the grazing point.
  - `A` at grazing needs the joined limits, as the method says.
- **Grazing attribution.** The grade split at grazing needs `H_lim` from the same limit operands. That should be stated as unavailable otherwise (a non-blocking suggestion).

### 4. Scope of a successful correspondence
- **What it covers.** The claims listed in `method.md` Section 4 are appropriately limited to the selected real-ω=3 equation.
- **Coverage gap.**
  - L has support only in rows 0–2 (`uniform-worker.py:239`, `vstack(curl[:,1:3], zeros(2,2))`).
  - A and R therefore test only E's `u_1..u_3` columns (all five rows). The `theta` and `e_W` columns are untouched, and they are the pressure-bearing channel for control 3.
  - The motivating gap was nonuniform and direct-term applicability. The cell table should record which nonzero local and pressure cells actually act through `E·L`, with their grade ancestry.
  - If the pressure-bearing and direct-term cells are zero-weighted by L, a pass does not discharge the earlier conditional applicability. The claim text should say so.

### 5. Controls and handling of a substantive mismatch
- **Control 1.** Omitting a term of E that acts through L is sound, provided the movement is recorded before the guard.
- **Control 2.**
  - The method's fallback for a symmetry-silent case is correct.
  - It must not be presented as a direction proof, and the method already says so.
- **Control 3.** It tests sheet sensitivity only on the pressure column and is not a mode check. This is acceptable as an addressing control.
- **Restored controls and failure handling.**
  - The three completed controls are restored and not re-run.
  - The late JSON failure is correctly treated as tooling, outside method acceptance.
  - The stop-and-report behavior for a nonzero residual is correct.

### Non-blocking suggestions
1. Name the ε/2 and ε²-power candidates (see Section 1) among the pre-fixed scalar candidates.
2. Make the η/σ ↔ `eta_bg`/`sigma_W` identity (including origin values and endpoint normalization) an explicit gate.
3. Save raw and on-wave `I_old`.
4. Provide a per-cell "acts through L" coverage table.
5. Record the `q+β≠0` certificate.

CLEAR FOR THIS BOUNDED WEAK-END UNIFORM-COMPARISON METHOD