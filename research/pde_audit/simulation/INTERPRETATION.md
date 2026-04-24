# Simulation Interpretation

The current simulation layer is a target-blind reduced FEM surrogate.  It is a
branch-search and handoff test, not yet a validation of the full nonlinear PDE
claim.

The bundle also contains a small nonlinear manufactured exporter lane.  That
lane verifies the nonlinear solver mechanics and emits frozen V2-22B-compatible
packets from the pre-target protocol, but it is still a manufactured readiness
model rather than the final physical moving-throat PDE exporter.

The strict physical nonlinear moving-throat exporter is explicitly blocked in
this bundle.  The ledger currently supports an effective wall/interface closure
and reduced BdG/Maxwell blocks, but it does not yet freeze the promoted
parent-level `S_eta`/`S_Sigma` action or the full coupled physical residual
equations.  `verify_physical_model.py` records that boundary and confirms that
no physical branch packets are emitted.  The report also hashes the cited
source notes and checks the supporting phrases used for each blocker.

The default `operator_v1` protocol freezes 192 V2-22B-compatible solver packets
before target evaluation.  The evaluator then classifies those frozen packets
through the existing V2-22B -> V2-22A -> V2-21 chain.

Current referee-run result:

- The nonlinear readiness checks pass 5/5.
- The physical model guard passes with strict physical export still blocked.
- The nonlinear manufactured exporter freezes 3 V2-22B-compatible packets.
- 3/3 nonlinear packets validate against the V2-22B handoff schema.
- 3/3 nonlinear packets pass the open and stability gates.
- 0/3 nonlinear packets pass the target packet.
- The nonlinear manufactured packets require at least a `60.81133778677986`
  multiplier on `C` or `D0` at fixed `A` to reach the one-pole surface.
- 192/192 frozen packets validate against the V2-22B handoff schema.
- 189/192 packets pass the open and stability gates.
- 0/192 packets pass the target packet.
- The reduced open+stable packets require at least a `7.387352901601946`
  multiplier on `C` or `D0` at fixed `A`; the median required multiplier is
  `16.30132163440465`.
- The reduced open+stable median `P0` normalization multiplier is
  `171.65261223353198`.
- The best post-hoc candidate has score `3.2433255817471953`.
- Among the 189 open+stable candidates, the one-pole ratio
  `D0*C/(3 A^2)` ranges from `0.0033775383274364888` to
  `0.1353664855760648`, with median `0.06134471930726503`.
- 0 open+stable candidates lie within `|ratio - 1| < 0.1`.

This is evidence against the reduced surrogate family and the current
manufactured nonlinear readiness family realizing the target.  It is not
evidence against every possible nonlinear moving-throat PDE branch.  The next
claim-relevant step is to replace the manufactured forcing/export path with the
real nonlinear PDE/continuation branch equations while keeping the same freeze,
target-blind guard, evaluator, obstruction diagnostics, and required-deformation
map.

`mechanism_gap_report.json` summarizes this state as a large one-pole support
deficit: the current frozen families under-supply the `C`/`D0` channel relative
to `A`, local continuation is not recommended, and strict physical export
remains blocked until the parent wall dynamics and coupled residual equations
are frozen.
