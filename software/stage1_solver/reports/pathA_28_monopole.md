MONOPOLE_DIPOLE_RETURN_CONDITIONAL

## Scope Caveat (read first)

This is a VERIFIED CONSTRAINT-SPECIFICATION, not a falsifiable test of no-monopole-radiation. What is verified-solid: the raw ell=0/1/2 outgoing amplitudes and their radiation orders, and the exact moments/orders at which any brane<->bulk return must cancel. What this report does NOT establish: that monopole/dipole radiation is suppressed or unavoidable.

- The cancellation condition (`R0=-M0`, `R1=-D1`) is the algebraic bookkeeping of what the return must cancel. Its symbolic "derivation" is an identity (`x - x`), not a deep load-bearing result.
- The `UNAVOIDABLE` rung is NOT decidable in-scope: `cancellation_possible` is a parameter (a literal flag), not a computed quantity. Deciding suppression-vs-unavoidable requires the track-3 brane<->bulk return, which is out of scope here.
- The real falsification lives in track 3: whether an admissible return can actually deliver `R0=-M0`, `R1=-D1`. This report only specifies the target it must hit.

## Computed Verdict

The raw ell=0 and ell=1 outgoing c_s-wave amplitudes are present.  Suppression is therefore not automatic; it is conditional on a brane<->bulk return that cancels the specific source moments below.  Solving that return is track 3 and is out of scope here.

## Raw DtN Amplitudes

- `ell0`: `A_raw = I*M0*a*omega/cS`, `p_raw = 1`
- `ell1`: `A_raw = I*D1*a**3*omega**3/(2*cS**3)`, `p_raw = 3`
- `ell2`: `A_raw = I*Q2*a**5*omega**5/(27*cS**5)`, `p_raw = 5`

The normalized outgoing DtN admittances have first radiation-phase terms `i*a*k`, `i*a^3*k^3/2`, and `i*a^5*k^5/27` for ell=0,1,2 respectively, with `k=omega/c_s`.

## Cancellation Condition

Headline condition for pde_ledger open-item #9:

- `ell0`: track 3 must deliver `R0(omega)=-M0(omega)`, where `M0=int_brane S_leak d^3x`.  This cancels the raw `O(omega^1)` outgoing coefficient.
- `ell1`: track 3 must deliver `R1_i(omega)=-D1_i(omega)`, where `D1_i=int_brane x_i S_leak d^3x + int_brane j_i d^3x`, including any carried odd wake.  This cancels the raw `O(omega^3)` outgoing coefficient.

These are conditions, not closure.  The report does not use `S_leak=0`, strict recovery, projection-locking, derivative-vertex suppression, or the statement that leaked medium enters the bulk as a kill.

## Raw Vs Vertex

The raw branch uses the projected-continuity source moments directly.  The derivative outlet vertex `g_W0(omega)=eta*omega` is reported only as a `branch_assumption`; in the two-vertex kernel bookkeeping it adds two powers of `omega` and is not used for the verdict.

## Controls

These controls do NOT test whether suppression occurs. What they confirm is narrow: the source moments `M0`/`D1` are kept live (no `S_leak=0`, no strict-recovery basis, no projection-locking that would zero them out by construction). Beyond keeping the source live, the controls pass by construction — they are not able-to-fail probes of the physical question. Treat them as guards against the obvious tautologies, not as evidence of suppression-vs-unavoidable.

- `raw_monopole_present`: same_pipeline=`True`, fired=`True`
- `steady_no_radiation`: same_pipeline=`True`, fired=`True`
- `quadrupole_survives`: same_pipeline=`True`, fired=`True`
- `return_necessity`: same_pipeline=`True`, fired=`True`
- `anti_tautology_no_S_leak_zero`: same_pipeline=`True`, fired=`True`
- `anti_tautology_no_strict_recovery_basis`: same_pipeline=`True`, fired=`True`
- `anti_tautology_no_projection_locking`: same_pipeline=`True`, fired=`True`
- `anti_tautology_derivative_vertex_not_basis`: same_pipeline=`True`, fired=`True`
- `anti_tautology_no_track3_bulk_kill`: same_pipeline=`True`, fired=`True`

## Engine Status

- SymPy: `timeout 600 python3 software/stage1_solver/tools/pathA_28_monopole_sympy.py exited 0`.
- Mathematica: `timeout 600 math -script software/stage1_solver/tools/pathA_28_monopole.wl exited 0 and asserted PASS`.
- Engine agreement: `PASS`. The `26` shared-expression count includes a few trivial entries (e.g. integer radiation orders and bare coefficients); the load-bearing matches are the ell=0/1/2 amplitudes, kernels, and radiation orders.

## Provenance

- DtN/outgoing expansions reuse `research/4d_2_5pn` spherical-Hankel machinery for ell=0,1,2.
- The source moments `M0`/`D1` are this script's OWN constructions, built to be consistent with `research/pde_ledger` Part VIII projected open-system continuity (and the Stage 243/244 leakage lane checks). They are NOT verbatim-cited Part VIII objects.
- `G=c=c_s=1`, `K_eos=1/500`, and `(a*,L*)=(4731/2500,18121/10000)` are checked in the frozen slice.
