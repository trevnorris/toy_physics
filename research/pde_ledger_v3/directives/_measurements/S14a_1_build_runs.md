# S14a part 1 build runs — 2026-10-21

Deliverable: two independently constructed CAS engines, their final-code guarded transcripts, and this run record. Contract: `directives/S14a_1_decision_list.md`, revision 3, at HEAD `e9eed6145b4b170719d0c852e2199d44e9a86da8`. No commits were made.

The Wolfram source was written before the SymPy source. Neither engine reads the other engine or its output. No cross-engine comparison or scientific interpretation was performed. Formal integrals retain the live fields; velocity objects carry their nonzero-density domain. All contract OPEN items remain OPEN.

The requested `_scratch/s14a/guard` scratch path is a symlink to `_scratch/s11c/s14a`, satisfying the unchanged guard's log-root restriction. Runs were sequential, without time limits; each used pool `s14a`, a 50 GiB cap, zero swap and the default 32-task cap.

**SymPy**

- Repair scope: B1 only, against reviewed baseline `0ee93e9b9b437a11c49b0cc247ad0a2ef418c275`. The bulk-direction derivative integral supplies endpoint contributions to the interior and closed-projection face terms; dependent channel and shell objects are recomputed. Field content, domains and the spherical specialization are unchanged.
- Script: `research/pde_ledger_v3/scripts/S14a_1_projection.py`; SHA-256 `0a1c5f6d2aef7416881738a213794d2770235fc3e8fd9b9190d6b0ed1f281811`.
- Transcript: `research/pde_ledger_v3/scripts/out/S14a_1_projection.out`; SHA-256 `5ae0c78ae44934378e65cc851523119e18e41b08d46c3aee73ac230446884f80`.
- Guard outcome: exit `0`; child exit `0`; `limitsVerified = true`; `guardReason = null`; stderr `0` bytes.
- Runtime: child `2.0667957560217474` s; guard launcher `2.1435441979847383` s. Peak cgroup memory: `62517248` bytes (whole guarded cgroup, not process RSS).
- Requested objects not produced: none. The transcript is byte-identical to its captured stdout; the requested object tags and the projection-integration operands and residuals are present.
- No commit, review round or cross-engine comparison was performed in this repair; contract OPEN items remain OPEN.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s14a --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s14a/guard/sympy_B1_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_projection.py
```

**Wolfram**

- Replacement scope: B1, B2 and B3, against the engine preserved at `0ee93e9b`. This engine was written anew from the supplied equations and U1 encoding; it reads neither the SymPy engine nor either prior transcript. There was no cross-engine comparison.
- Construction: pull back the spacetime current to a fixed transverse interval, integrate its transverse derivative, and separately differentiate the regional `HeavisideTheta` currents. Each face-concentrated coefficient enters the closed projection with full weight. The supplied face factor is printed only. The native fields remain live; the spherical specialization keeps radius, `w` and time.
- Script: `research/pde_ledger_v3/mathematica/S14a_1_projection.wl`; SHA-256 `809e981244dbb47953b4bff8ba4ab99ac1c3ef6c73467883ddc2ecafa67ce6e4`.
- Transcript: `research/pde_ledger_v3/mathematica/out/S14a_1_projection.out`; SHA-256 `d578f221c1fc344866070fafaa781e27c3b0021905c68ba11805e4ac2729a9c3`. It is byte-identical to the final-code guarded stdout.
- Guard outcome: exit `0`; child exit `0`; `limitsVerified = true`; `guardReason = null`; stderr `0` bytes. Pool `s14a`; memory cap `53687091200` bytes; zero swap; 32 tasks; both Wolfram thread controls installed; no time limit.
- Runtime: child `5.0703665019827895` s; guard launcher `5.156231011002092` s. Peak cgroup memory: `229928960` bytes (whole guarded cgroup, not process RSS).
- Requested objects not produced: none. Formal integrals retain the unrestricted fields. Projected velocity and its volume flux, including their dimension objects, carry the nonzero projected-density condition. Every requested shell channel, storage, balance, radius dependence and dimensional object is printed.
- Validation: the final transcript contains the coordinate-divergence and closed-integration operands and residuals for each balance in each coordinate setting, as well as the concentrated coefficients and all requested object tags. No Wolfram message lines were found. No commit or review round was performed; all contract OPEN items remain OPEN.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s14a --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s14a/guard/wolfram_rewrite_03 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S14a_1_projection.wl
```

Validation: Python source parsed; Wolfram delimiters checked; both final sources ran through the guard. Each transcript is nonempty, includes the requested object tags and dimensional objects, and is byte-identical to its own captured stdout. No error-token matches were found. The repository-referenced `derived_or_declared.py` and `engine_output_checks.py` tools are absent from this checkout; no substitute review was launched.
