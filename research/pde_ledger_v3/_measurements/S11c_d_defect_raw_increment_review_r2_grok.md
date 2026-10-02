I'll review the staged both-face raw-increment build from the packet only: the worker, plan, method, and the indexed inputs and native source. Scientific execution stays off; the verdict will rest on what the source actually does.The review is source-only: I'll read the staged worker, plan, method, and packet index, then check the native bindings those files actually use.The packet index pins the worker and its restored operands. Next I’ll read the instrument itself and the method it claims to implement, then check those bindings in the staged sources.The worker is loaded. I’ll now check the manifest pins and the native formulas the new lower-face, closure, and Fourier joins actually call.The saved upper coefficient and the theta row are in hand. Next I’ll compare the lower-face normal and closure with the native formulas, and check whether the selected-integrand residual is exact.The review is still open: the lower-face resolvent text, the native kernel and Fourier bodies, and the row census have not been read yet. I will finish those source addresses before giving a verdict.The resolvent definitions match on both faces. Next I will check the saved integrand join, the native normal and kernel assignments, and the row census that the worker actually executes.NEEDS REVISION

The bounded worker stops at the saved selected-integrand join. That check is required before the both-face closure, the expanded-row census, and the nonzero-input controls, so this source review does not treat those later stages as executed. No scientific worker, ablation, or computer-algebra run was performed.

## Math / method

**Saved selected-integrand residual cannot be recognized as zero.** `input/worker.py` lines 416–427 build the ordered density and compare it with `input/saved/action.json` lines 2–4. `Journal.zero` (`input/worker.py` lines 228–233) emits the raw residual, then accepts the check only when `cancel(together(left-right))` is the integer 0.

The saved integrand is

`-5*sqrt(3)*t*(10*t-1)/(16*q_t*sinh(5*pi*t)*sinh(5*pi*t-pi/2))`.

Its `srepr` has the second hyperbolic factor `sinh(5*pi*t - pi/2)`. That form is the stored return of `sp.factor(...)` in `input/source/bare-worker.py` line 325. The new density is built from `A(x)=L*x/(4*sinh(pi*L*x/2))` with `L=10` and `Q=1/10`, so the slope factor is `sinh(5*pi*(1/10-t))`.

Those two expressions are the same rational density after the exact identities `sinh(5*pi*(1/10-t)) = -sinh(5*pi*t-pi/2)` and `(10*t-1) = -10*(1/10-t)`. `cancel` and `together` leave the two `sinh` atoms unchanged, so the structural residual stays nonzero and `require` raises `saved-selected-integrand`. The failure is preserved (`failure.json`, `automaticRetry` false, raw residual already written), and the run does not reach closure or the row controls.

Smallest correction: in that one comparison, rewrite the constructed slope factor by those two identities, emit the unreduced residual and the rewritten residual, and apply the existing `cancel(together)` test to the rewritten pair. Keep the ordered density formula as it is.

Uncertainty: SymPy was not executed here. The conclusion is the expression-tree mismatch between the saved `srepr` and the constructed `sinh(pi*L*(Q-t)/2)` under the normalizer actually called. If a future run shows `cancel(together)` collapsing this particular pair, the residual file is the evidence that would retire the finding. Cost of the later triangular closure and the full-row census is not known from source, because those stages sit after this raise.

## Source that agrees with this bounded scope

These were read as source and were not rerun.

- Saved JSON copies are pinned by path and sha256 in `input/input-manifest.json`, and `scientific_work` copies them before use. `omega` is joined to the saved integer 3; the development string `"1"` stays in `originalParameters`.
- The lower boundary uses face `-1`, lab slope `face*s`, extension factor `I*face*depth*height_face`, and the native `normal_exact` assignment. With the positive outgoing normal, the flat, height, and slope channels reproduce `rho*omega/qi`, `I*omega*rho*(qh-qi)/qh`, and `k*omega*rho/(qi*qs)`. `dtn_first_kernel` (`input/source/native-c1.py` lines 606–627) matches those after the `c_s^{-2}` substitution. The outward mirror is a later residual.
- Both dispersion substitutions in `factorization` match `qh^2-qi^2=-H*(H+2*k)` and `qo^2-qs^2=-H*(2*k+2*Q-H)`. The direct slot is the symbol `D` inside the executed `three_inverse` assignment; the native `z_three[0,2]` literal remains 0, and `kernel_apply` is not called.
- Both `LAB_HELD` / `RHO4_CONSTANT` resolvent definitions are `I + c*Z` with the same `c = Lambda_A_0*rho_m**-2*(1-I*omega*tau_A)**-1` and the face’s own DTN symbol (`input/native-selected.json` lines 88 and 112). The affine reconstruction is checked before `R_f` is used. `uniformGrazingLimitEstablished` is false and `exactMatch` is `UNRESOLVED`.
- `THETA_BALANCE` and `E_W_BALANCE` pressure-slot coefficients are free of `eta_bg` and are nonzero; the jet-slot coefficients already contain `eta_bg`. Both slots are substituted before the mixed derivative. The lower omission and doubling controls rebuild that substituted row. The separate jet-sign check is only a nonzero trace factor.
- `profile_bindings` and `kernel_apply` (`input/source/native-c2.py` lines 849–861 and 454–471) match the executable Fourier contract: profile forward power `-3`, source inverse power `-3`, unnormalized source forward, and two edge coordinates. The contract records that no integral is evaluated.
- Schema units used by the dimension join are `rho_m = [-4,0,1]`, `omega = [0,-1,0]`, and `W_0 = [1,0,0]` (`input/source/native-c2.py` line 59). With the reduced hats and `dt` of dimension `(-1,0,0)`, the sum is the saved `(-2,-1,1)`.
- Containment is checked before the SymPy import. The guard verifies `RuntimeMaxUSec=infinity` and `Restart=no`. There is no scientific retry loop. Failed journal inputs stay exclusive-create files, and posthashes run in `finally`.

## Tooling

`input/runtime-source/supervisor.py` lines 72–78 and `input/runtime-source/shared-guard.py` lines 414–429 record stderr size and capture stdout in separate files. Neither compares that stdout with `checks.json` nor requires stderr to be empty. `input/plan.md` line 77 asks for stderr evidence at runtime acceptance. This does not change the density. Smallest correction, when a run is accepted: hash-compare the captured stdout with `checks.json` and require stderr length 0, with no retry.

The recorded stdlib test (`input/tooling-checks.json`) is `physicsValidated: false`. It does not cover this integrand identity.

## Optional

`input/worker.py` line 420 stores the uncancelled height/PV law as a string. The coefficient contact at `H=0` is a real residual. Emitting `(W0/4)*C(k;0,Q)*j_hat(Q)` as an expression would make that distributional operand the same kind of object as the ordinary density.