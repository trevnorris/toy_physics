I'll assess the rest-bulk near-unity defect method from the staged sources only, starting with the method document and the packet index.The method and index are in. Next I’ll read the review request and the saved operands that the selected-source, boundary, and finite-solve claims actually rest on.The selected source restrictions and the three nonzero reconstruction remainders are in hand. I’ll now check the closure, pressure consumers, and the native insertion points those claims depend on.The assessment still needs exact addresses for the two-transfer rule, the source/consumer limits, and the lower-face premises. I am reading those passages in the staged input before writing the three verdicts.The selected source and consumer evidence does not protect a finite-contrast defect calculation. The missing scalar-sector direct increment still has to be built. The proposed one-dimensional route is the right bounded task, and it is not yet specific enough to implement.

- SELECTED EVIDENCE: SUPPORTED WITH STATED LIMITS
- FINITE-PILOT REUSE: NOT ESTABLISHED
- NEXT INCREMENT METHOD: NEEDS REVISION

No producer was run, no saved integral was recomputed, and no future worker is cleared. A completed process and a zero guard are not scientific acceptance. Every worker in this packet sets `scientificAcceptance` false.

## Selected evidence

The upper-face mixed term is a nonzero unintegrated integrand at one rest point. `bare/profile-action-return.json` has `integrationNotPerformed: true`. The integrand is

`-5*sqrt(3)*t*(10*t - 1)/(16*q_t*sinh(5*pi*t)*sinh(5*pi*t - pi/2))`,

with branch points `t = ±1/5`. At the midpoint `t = 1/20`, `bare/responsive-controls.json` records the sample `-sqrt(5)/(32*sinh(pi/4)^2)`, and the wrong-sheet and slope-omission images differ from that sample. This is the coefficient of independent `eta_bg*sigma_W` on the upper face at profile momenta `0 → 1/10`, `omega = 3`, `c_s0 = 10`. The development file `physical-input.json` still has `omega = 1`; the workers override the frequency. It is not an `eta^2` correction and not a finished integral.

The generic upper coefficient that must be reused is `bare/boundary-coefficient-return.json` field `mixed`:

`I*omega*rho_m*(-H*q_i**2 + k*q_h*q_i + k*q_h*q_o - k*q_h*q_s - k*q_i**2)/(q_h*q_i*q_o)`.

The field `selected`, `-I*H*omega*q_i*rho_m/(q_h*q_o)`, is only the `k = 0` reduction. `q_s` and the slope transfer `S` are still live in `mixed`.

The one-face closure multiplies that inherited whole kernel by a nonzero reference factor. `trace/physical-factors.json` keeps `integralEvaluated`, `sourceIsTransverseMode`, and `fullSlabContraction` false, and its source channel is per unit normal velocity. `source/consumer-worker.py` `selected_increment` (lines 148–154) writes the minus pressure and jet slots to 0. `consumer/selected-fourier-contraction.json` records `lowerFaceCorrectionConstructed: false`. The symbol `inherited_whole_bare_mixed_kernel` is a comparison object at this point.

The transverse zero is the flat zero-grade source on a curl. `consumer/plus-source-restrictions.json` gives `restrictions.TRANSVERSE = 0`, while `LONGITUDINAL`, `THETA`, and `E_W` are nonzero. The flat `source00` uses velocity only through `u_1_d1`, `u_2_d2`, and `u_3_d3`, together with `theta`, `e_W`, their Laplacians, and `e_W_t`. `source/native-c2.py` `wave_jet` (lines 139–143) builds that transverse trial as the curl of compact-support potentials, without a Fourier projector. The grade string on the restriction file limits the test to `kernel(1,1) * source(0,0) * consumer(0,0)`. The full chemical amplitude still contains undifferentiated velocity times profile slopes; those terms sit above this grade. `consumer/selected-fourier-contraction.json` is an off-shell witness at `kin = (0, 1/5, 1/10)`, `kout = (1/10, 1/5, 1/10)`, `omega = 3`: the `u_1` channel is 0, and the `u_2`, `u_3`, `theta`, and `e_W` channels are not. `incidentTransverseMode` is false.

Pressure consumers split by row. `consumer/weak-directions.json` counts zero `delta_p` and `d_w_delta_p` occurrences on all three `U_BODY_BALANCE` rows, and `FACE_GENERALIZED_FORCE_ROWS` `U` is `Tuple(0,0,0)`. The same file’s note says the source restriction precedes kernel action. `THETA_BALANCE` and `E_W_BALANCE` row factors in the Fourier contraction are nonzero and carry `epsilon_shape`; the normalized source does not. `consumer/zero-grade-jet-consumers.json` is 0 on every tested row, including the two scalar rows, so that jet zero does not describe the pressure-value coupling. Removing `delta_p_plus` in `continued/native-pressure-consumer-omission.json` sends the selected row-per-`D` factor to 0 while the recorded movement is nonzero. The probe there is `amplitude_e_W = 1`. `continued/routing-control.json` has actual vector and curl 0; the misrouted scalar-to-vector image is nonzero.

The representation correction checks out as coefficient extraction. The continuation tail in `continued/resumed-source-tail.json` calls `operator_form` only on the divergence-omission amplitudes of `THETA_BALANCE` and `E_W_BALANCE`. Residuals `000`–`002` have raw remainder `Integer(0)`. Residuals `003`–`005` have structurally nonzero raw remainders in `independentJet2` and `independentJet3`, and cancelled remainder `Integer(0)`. The two denominators are

`D1 = 60*sqrt(3) + 210 + 482*I + 500*sqrt(3)*I`,

`D10 = 600*sqrt(3) + 2100 + 4820*I + 5000*sqrt(3)*I`.

Term by term, `D10` is `10*D1`. `D1` is the `E_W_BALANCE` denominator in the Fourier contraction. `source/consumer-continue-worker.py` `resumed_operator_form` (lines 214–235) cancels that residual before the zero decision. The collected curl-jet coefficients stay nonzero, which is what omitting one divergence piece must do. The full transverse restriction remains the separate integer 0. This identity does not prove the curl annihilation, and it is not an operand defect.

`retained_shape` in `source/native-c2.py` (lines 664–668) keeps grades `(0,0)`, `(1,0)`, `(0,1)`, and `(1,1)`. `shape_coefficients` pushes `Integral` through that rectangle (lines 716–718). `grades` (lines 807–811) reports a shape-dependent denominator by literal degree and does not assign that report to an unexpanded inverse. `source/finite.py` `matrices` and `solve` (lines 57–93) assemble and factor a full `5*N` matrix, with both-end trace rows, and label the result `RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED`.

Both faces are restricted in `source/consumer-worker.py` (lines 295–357). `plus-source-restrictions.json` and `minus-source-restrictions.json` have the same packet hash `a3f08120…`, while the source-input hashes differ. The saved flat `source00` carries no face label. The plus and minus raw operators do (`s11cc1_mu_theta_lab_held_plus` versus `_minus`). The identical restriction files are the saved output of that face loop. They were not recomputed here.

`bare/dimensions.json` records residual 0 for the reduced kernel `[-2,-1,1]`, equal to impedance `[-3,-1,1]` plus one length, with two edge deltas factored. `DIMENSION_SCHEMA` gives `rho_m` the declaration `[-4,0,1]`. That residual is a saved worker result.

The conditional statement in `method.md` lines 47–63 is limited correctly. A change confined to `L11`, with transverse `psi00`, unchanged boundary data, and a unique reference inverse, would give `delta_psi11 = 0` only after the reference solution, both-face source and consumer maps, regularity, and the exterior prescription are joined. Those joins are absent. Proving that formal line would not remove the untruncated `5*N` calculation: the scalar and longitudinal sources are nonzero, and the solve inverts the coupled matrix. No separate proof job is warranted.

`uniform-context.json` opens with `scientificStatus: SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE`. That file concerns constant ends. It was not re-audited here. The effective bulk speed stays a later declared parameter; the rest dispersion used above does not contain it. The calibrated ratio and the drain stay outside this increment.

## Finite-pilot reuse

No retained-solution exemption removes the scalar-sector increment. The transverse column of flat `source00` is silent. `THETA_BALANCE` and `E_W_BALANCE` are not, the longitudinal and scalar restrictions are not, and `finite.py` does not truncate the inverse to the rectangular grade set. An induced scalar field can enter the block the direct mixed term has not yet filled. Poles and the grazing threshold are outside the conditional grade identity. The old benchmark stays unlabeled; `finite.py` `control_report` leaves `analogLightCalibration` open, and the numeral cited in `method.md` line 62 was not in the finite source that was read.

These objects can be carried in unchanged once the new coefficient is actually joined:

- The generic upper `mixed` coefficient, the native height and slope joins, and the saved boundary residuals for that upper-face derivation.
- The tanh transforms and the distributional rule in `bare/profile-transform-return.json`: forward `exp(-i k y)/(2 pi)`, inverse `exp(+i k y)`, height `delta(t)/2 + PV[L/(4 i sinh(pi L t/2))]`, with `t*delta(t) = 0` and `t*PV(1/t) = 1`.
- Source, density, memory, and `epsilon_shape`, including `m1_profile` and `mu_R_bg`.
- The existing first-shape path. `kernel_bridge` in `source/native-c2.py` (lines 409–411) puts a literal 0 in `z_three[0,2]`. `trace/kernel-measure.json` says the first-shape products are formal composition integrands. Adding a new direct entry once in that zero slot leaves the first-shape iteration once.
- The finite assembly and LU/SVD driver, as machinery, after the new entries exist.
- Constant-end trace data only where a fresh zero-jet reduction of the new both-face coefficient reproduces the end operator already used. The saved check `(h*s*mixed)` at `s = 0` is the upper-face slope reduction in `source/bare-worker.py` line 234.

These cannot be called corrected:

- Any `5*N` matrix or solve, until the new direct scalar entries and their contact, principal-value, and branch controls are in the joins.
- `inherited_whole_bare_mixed_kernel`, and every row built by `selected_increment` with the minus slots set to 0.
- The `k_in = 0` integrand, as a stand-in for general profile momentum.

## Next-increment method

The scope in `method.md` lines 65–74 is the right task: one-dimensional tanh, both faces, live profile momenta, fixed edge momenta `(1/5, 1/10)`, rest bulk, raw integrand, stop at the incremental operator. Three method sentences do not yet name the objects an implementer must build.

**Blocker 1, equation and method.** `method.md` lines 85–90 say “Insert both assignments of height and slope to the two transfers” and then say that those assignments must not double a coefficient already generated by the ordered expansion. The saved action is one ordered product. `source/bare-worker.py` lines 323–328 set `H = t` inside `selected` and multiply by the weighted height at `t` and the slope jet at `Q - t`. Lines 374–376 say that this convolution already includes both assignments, and that the symmetrized form integrates identically under `t → Q - t` only with the factor `1/2`. At general `k`, `mixed` is not symmetric in `(H, S)`: `H` appears as `-H*q_i**2`, while `S` appears through `-k*q_h*q_s`. Adding the swapped assignment, or averaging with `1/2`, builds a different operator from the ordered coefficient.

Smallest correction: the integrand is the single product `C(k; H=t, S=Q-t)` times the height transform at `t` times the slope transform at `Q-t`, with depths of `k`, `k+(t,0,0)`, `k+(Q-t,0,0)`, and `k+(Q,0,0)`. The `k_in = 0` expression with the factor `1/2` may be used only as an identity check of that already-built special integral.

**Blocker 2, equation and method.** The same lines treat the contact and principal-value terms as something to rederive, without naming which pieces survive. At `k_in = 0` every surviving term carries `H = t`, and `source/bare-worker.py` lines 306–318 cancel `t` against the height transform before `t = 0` is used. Every term in `mixed` that is independent of `H` multiplies the full height transform, including `delta(t)/2`. That contact is absent from the saved integrand. `kernel_apply` (lines 455–472) integrates its `second` argument over all three middle momenta and over output momentum, input momentum, and the source coordinate, with measure `1/(2*pi)**3`. The reduced profile convention is the other measure: one profile integral and two factored edge deltas. A reduced integrand that has already been integrated over the transfer is a different object from that `second` argument. Putting it there integrates it again. The direct slot that can receive one new unintegrated coefficient is `z_three[0,2]`, which the native constructor currently sets to 0. The first-shape iteration stays in the off-diagonal transfers.

Smallest correction: write the contact from the `H`-independent part of `C` against the saved height transform, in the reduced convention, and keep the principal-value integral separate. Pass an unintegrated coefficient into the direct slot only. Leave the inherited whole `D` as a one-point comparison.

**Blocker 3, equation and method.** `method.md` lines 77–85 require a lower-face derivation from the native normal, location, pressure convention, and outgoing depth sign, and one physical outgoing dispersion, without writing either rule. `reference_pressure_kernels` (lines 504–506) already sets `reference = face*W_0/2` and `extension = exp(I*face*qo*(NORMAL-reference))`. The graph normal in `bare/native-graph-normal.json` uses `face` in `normal_exact`; the saved upper normal is the `face = 1` evaluation, whose linear tangential slope is `-1`. The trace worker’s `reference` (lines 284–286) drops `face` because that witness is upper-face only. The bare worker hardcodes `face = 1` at line 189. Copying that witness and also changing the sign of `q` applies the face sign twice. `outgoing_spectral` (lines 439–446) takes `sp.solve(...)[-1]` and multiplies by `sign(omega)` only on the propagating piece. That is a solver-order representative, not a branch rule. `trace/middle-domain.json` excludes `q = 0`. The `k_in = 0` endpoint check does not cover a depth that vanishes at a modal match.

Smallest correction: for `omega > 0`, use one root on the input momentum, both transfer momenta, and the output momentum,

`q = +sqrt(omega**2/c_s0**2 - |k|**2)` when the radicand is positive, and `q = +I*sqrt(|k|**2 - omega**2/c_s0**2)` when it is negative,

with `|k|**2` including the edge momenta `(1/5)**2 + (1/10)**2`. Stop before a defect run if any radicand is zero or any depth in a denominator vanishes. Build the lower face from `normal_exact` at `face = -1`, `reference = -W_0/2`, `exp(-I*q*(normal-reference))`, and the minus slots `delta_p_minus` and `d_w_delta_p_minus` of `reference_pressure_kernels`. Keep this same positive `q`.

Step 4 of the method then does the double-counting control that the closure source supports: one new direct contribution per face, the existing first-shape iteration left once, full slab rows contracted before the grade restriction, pressure and normal jet kept as separate factors. Memory and permeability stay in the reused source. A source-silent `U` row or a silent transverse column is a selected-subspace fact; the exterior inverse in `finite.py` still couples all five fields through both end trace maps.

No tooling defect blocks this. The continuation’s exact-cancel predicate is the representation fix, and the cancelled residuals above are zero. No further diagnostic of the conditional `L11` identity is useful.

Uncertainty that the packet leaves open: the face-independence of the flat source is a saved identity of the two restriction files, not a fresh expansion of both chemical operators. Consumer dimension totals were checked against the targets declared in `source/consumer-worker.py` lines 544–573 (`THETA_BALANCE` `[-3,-1,1]`, `E_W_BALANCE` `[-1,-2,1]`); each emitted total in `continued/consumer-unit-joins.json` was not re-added by hand. The uniform grazing joins were not re-read past the status and the opening speed records.

## Next action

Revise only those three method sentences, then construct the raw increment. Use the saved generic upper `mixed` coefficient as one ordered convolution at general profile momentum. Rederive its contact and principal value in the reduced Fourier convention. Derive the lower-face coefficient from the native normal, `face*W_0/2`, and `exp(I*face*q*(normal-reference))`, on the single positive outgoing root above, and stop if a depth vanishes. Emit raw middle-momentum integrands and separate contacts. Add that direct coefficient once per face beside the existing first-shape iteration, and contract the full slab rows, including `THETA_BALANCE` and `E_W_BALANCE`, before any grade cut. Do not evaluate the integral, do not rerun the `k_in = 0` action or the uniform schedule, and do not change the finite matrices until every contact, principal-value, and branch piece has an explicit integrability statement.