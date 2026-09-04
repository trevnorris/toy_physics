# S11c-c1 SymPy engine — repair directive (fold of the 2-leg build review)

## Authority and boundary
Repair `research/pde_ledger_v3/scripts/S11c_c1_bulk_closure_sympy_audit.py` (committed reviewed baseline
`65afa1cd`). `CLAUDE.md` binds. `directives/S11c_c1_SHARED_PHYSICS.md` (the migrated spec) is the sole physics
authority and wins every conflict; `directives/S11c_c1_sympy_build_directive.md` governs the export topology
(fold over the frozen base, the exact 44-root `IMPORT_KEYS`, the own-rows delta, the §D3 guard). This directive
adds **no physics** — every construction it names is already in the spec; it points at the spec section and fixes
an implementation defect the two review legs found. Add no expected value or acceptance criterion (rule 5): every
residual below is a computed measurement adjudicated after the run, ⛔ never asserted, ⛔ never stated to be zero.

⛔ **Do NOT touch the confirmed-sound core.** The two review legs (independent derivations) verified the
two-momentum DtN kernel (both branch legs live; the η/σ_W multigrade split reconstructs the two-leg `k·k'` form —
`_measurements/c1_tangential_kernel_adjudication.py`), the operator-inverse face response `[I+(Λ_A/ρ_m²)Z]⁻¹`,
`Λ_X`-only-in-traction, opaque `μ_θ`, and the exported 44-row own-rows delta. ⛔ Leave `dtn_first_kernel`,
`closed_coefficients`, `response_operator_case`, the `IMPORT_KEYS` set, and the export/guard topology unchanged
except where a fix below requires it. The 5 exported model-level rows must be **byte-identical** after the repair
(all 5 fixes are in emit-only controls) — the builder confirms this in the §report.

## The five repairs (each folds one verified finding; ⛔ fold once, do not re-scope)

**R1 — the energy-balance second route must be genuinely INDEPENDENT (spec §3b obj. 3).** Current defect: the
"far-field flux" operand is `z_energy/(iρω)` multiplied back by `iρω`, i.e. `δp` **at the face** — the spec
explicitly forbids this ("⛔ not δp at the face"; `½Re(δp·V*)` equals the bulk flux by the acoustic identity and
never sees `t_s`), so `ENERGY_RESIDUAL[BASELINE]` is a structural `A−A` zero that a branch/impedance/`q_out`-sign
error cannot move. Fix per §3b obj. 3: construct the bulk field `φ` (the outgoing half-space solution driven by
the face normal velocity under the §1b radiation condition) and compute the **outgoing far-field Poynting flux**
from `φ` on a control surface at `|w|→∞` — a function of the outgoing-wave amplitude, ⛔ not of the face `δp`.
Keep the face operand the true-area **traction** pairing `½Re Σ_s ∫ a_s^0 t_s·v_face,s*` built from the §3b
`t_s` object. Emit both operands and their residual. A one-sided reversal of the `t_s` sign must move **only** the
face operand (residual becomes nonzero); a `q_out`/branch-sign error in `φ` must move **only** the bulk operand.
⛔ The two operands must not share a common `z_energy` factor — that is the defect.

**R2 — the representation-invariance control must be a genuine SECOND CONSTRUCTION (spec §5a).** Current defect:
the `EULERIAN` and `HANZAWA` branches of `dtn_first_kernel` return the **same** rational function (`route=` only
regroups the algebra); `EULERIAN − HANZAWA ≡ 0` identically, so `REP_INVARIANCE_RESIDUAL` cannot detect a
disagreement. Fix per §5a: build the `HANZAWA` route as an **independent** radiation-preserving boundary
flattening — a cutoff/Hanzawa change of coordinates equal to the face map near the boundary and the identity at
infinity (⛔ not the global scaling `w′=[w−ζ_c]/[W_bg+δW]`, secular at infinity), OR a boundary-integral /
layer-potential construction whose outgoing kernel is the §1b radiation-selected one — with its transformed
radiation condition stated. The two routes must be structurally different constructions of the **same** DtN, so a
one-sided corruption of the Eulerian route alone drives `REP_INVARIANCE_RESIDUAL` nonzero (verify by one-sided
corruption, §5a). ⛔ Emit the residual `sp.simplify`-reduced enough that an identical pair prints `0` (the current
`object_difference` leaves an `evaluate=False` `Add` that never prints `0`).

**R3 — the on-shell reduction must actually fire (dispersion substitution, not an `Add`-key `xreplace`).** Current
defect: `on_shell` keys on `Add(k_out_1²,k_out_2²,k_out_3²)`, but `sp.factor` distributes `c_s0²` across the
`k_i²` so the bare `Add` never matches — `DTN_RIGID_SHIFT_RESIDUAL` is emitted **off-shell** (it is `0` only under
a direct `q_out²` substitution), and the §5d `ZERO_JET` operand carries the same off-shell term, reading as a
spurious thickness-dependence finding. Fix: reduce on shell by substituting the dispersion on `q_out²` directly
(`q_out(ω,k)² → ω²/c_s0² − |k|²` on each leg) — or `collect`/`together` so the reduction is representation-robust
— everywhere the rigid-shift cancellation (§3a) and the §5d zero-jet regression are formed. The rigid-shift
cancellation must be **emitted as the computed residual** (§3a), ⛔ not asserted. ⛔ Do not change the kernel; only
the on-shell reduction of the residual/regression operands.

**R4 — the dissipation Hermitian objects must be fed the FULL bare DtN `Z_0+Z_1` (spec §3b obj. 1–2).** Current
defect: `DTN_HERMITIAN_PART` (`hermitian_kernel(direct_kernels[...])`) and `port_matrix` (`p_v = … z1 …`) are fed
the first-shape kernel `z1` only, so both vanish at `(η,σ_W)=0` and cannot audit the leading bulk radiation
`H_a[Z_0]` (the S11b object, §3b obj. 1) nor the closed port map at leading order (§3b obj. 2). Fix: feed the
**full bare DtN** `Z = Z_0 + Z_1` (the flat symbol `dtn_flat_symbol` on the appropriate leg plus the first-shape
kernel) to `DTN_HERMITIAN_PART` and to the port map, so at `(η,σ_W)=0` the Hermitian part reduces to `H_a[Z_0]`
and the port Hermitian form to the leading closed-port map. Keep the true-area boundary pairing and the
`s→−s` adjoint (`hermitian_kernel`'s swap) intact. ⛔ Restrict any sign-definiteness claim to the subspace where
the zeroth-order Hermitian form is nondegenerate (§3b), emitting `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER` on its
nullspace — ⛔ state no sign in prose.

**R5 — `assert_delta_is_minimal` must check the INTENDED export closure, not the delta itself.** Current defect:
`own_closure = resolved_keys` is the delta's own key-set, so the call asserts `delta == delta` and an accidental
re-accumulation would pass. Fix: compute `own_closure` as the bind-closure of the **five intended model-level
export roots** (`dtn_flat_symbol`, `dtn_operator`, `dtn_kernel`, `s11c_c1_face_response`,
`s11c_c1_face_response_coeffs`) over the NEW `s11cc1_` symbols they introduce — derived from the roots, ⛔ never
read back from the written delta — and pass that as `own_bind_closure`. The check must then RAISE if a stray
inherited or non-closure row is added to the delta (verify with a one-row injection). The lookup smoke-test and
`check_consumer` are already genuine — leave them.

## Method and obligations
The three script clauses are exact (spec §6): a script prints computed CAS objects and states no conclusion; for
every comparison emit operand A, operand B, and `A−B` before any guard; a physics disagreement emits and
continues (exit 0), nonzero exit only for operational failure. Every new operand (the far-field Poynting `φ`, the
Hanzawa route) re-enters the chain at the **ansatz / bulk operator / radiation condition / supplied law**, ⛔
never at a result. ⛔ No emission is conditional on a payload's value. Re-run the engine to regenerate
`S11c_c1_exports.py` and re-run the §D3 guard; the export must remain the 44-row own-rows delta with the 5
model-level rows byte-identical.

## Builder report (§report, under 25 lines)
Files written, the diff scope (which functions changed), confirmation that `dtn_first_kernel` /
`closed_coefficients` / `IMPORT_KEYS` / the export topology are unchanged and the 5 exported model-level rows are
byte-identical, the guard outcome, runtime, and any object it could not compute. ⛔ Report no computed value and
no conclusion about the physics.
