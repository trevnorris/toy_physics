# Measurements — S11c-b #88 blast-radius build directive (premises only)

This file records the commands that back the **premises** the directive states about the
engine. It deliberately does **not** state the #88 result (the per-row disturbance) — that
is the withheld measurement the build reproduces independently; it lives in
`~/.s11_build/S11c_b_88_probes/` (outside the repo) until the post-build RESULT record.

## P1 — 15 candidates/source, engine selects 8, omits 7 (frozen quotient)
Reproduced from the engine's own enumeration (imports the audit module, calls
`enumerate_new_candidates` + `basis_euler_signatures` + `quotient_independent_indices`):

    python3 ~/.s11_build/S11c_b_88_probes/probe_omitted_enumeration.py
    → N_CANDIDATES: 15
      N_SELECTED: 8  indices(0-based): (0, 3, 4, 5, 6, 8, 9, 12)
      N_OMITTED:  7  indices(0-based): (1, 2, 7, 10, 11, 13, 14)
      OMITTED_LABELS: FIRST_JET_CONTRACTION_{02,03,08,11,12,14,15}

Matches the committed #86 record's frozen 8-selection {01,04,05,06,07,09,10,13}.

## P2 — the frozen quotient never differentiates the spurion (source of the 8-vs-15)
`basis_euler_signatures` (`scripts/S11c_b_brane_operator_sympy_audit.py:936-969`) builds
`derivative_maps` from `fields` = `basis_fields` (`:1025-1028`, DOF only: `bu`, `btheta`,
`be` + their jets). The abstract spurion `bg` (`:1014`) is **not** in `fields`, so
`basis_dx` (`:947`) never lifts it ⇒ the second background derivative (Hessian) is absent
from the quotient. Confirmed structurally in
`~/.s11_build/S11c_b_88_probes/probe_correct_el.out`: the omitted candidates 03/12/15 have
**identically-zero frozen signatures** (they look like null Lagrangians only because the
Hessian is frozen); adding the spurion to the derivative map restores a nonzero EL equal
to the pure Hessian term (`CORRECT−FROZEN = −Hess·DOF`).

## P3 — strong rows freeze the Hessian; coupling/admissibility retain it
- Global `dx` (`:616-621`) differentiates only `DERIVATIVE_MAP` keys; the map carries only
  the **first** background jet (`:611-613`: `DERIVATIVE_MAP[i][W_bg]=grad_W[i]`,
  `[mu_R_bg]=grad_mu[i]`). No `grad_W[j]→…` key ⇒ second derivative frozen.
- `operator_from_density` (`:1459-1511`) builds the emitted strong rows
  (`U_BODY_BALANCE`/`THETA_BALANCE`/`E_W_BALANCE`) with that global `dx`.
- The committed routines that **do** install the second jet, in a local map copy:
  `operator_dx` (`:1850-1874`, `derivative_map[grad_W[jet_direction]] =
  sigma_W*w_profile_second/L_W`, `:1864-1868`) and `background_dx` (`:2122-2137`,
  `:2127-2132`). The directive's Hessian-retaining `dx` reuses this pattern verbatim.

    grep -n 'DERIVATIVE_MAP\[i\]\[W_bg\]\|derivative_map\[grad_W\|second_w\[' \
      research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py

## P4 — retained-grade rule "η,σ_W ≤ 1"; the Hessian is σ_W¹ (retained)
`first_shape_series` (`:713-725`): subs `PROFILE_GRADE_SUBS`, `series(…, eta_bg, 0, 2)`,
keep terms with `eta_bg power ≤1 AND sigma_W power ≤1`. The second background derivative
enters as `sigma_W*w1_profile_d{i}d{j}/L_W` (`:1865`, `:2117`) — a **single** `sigma_W`
(σ_W counts background-factor count, not derivative order) ⇒ the Hessian survives the
truncation. Confirmed: the harness run reports `grade_ok=True` for every row (all
surviving terms σ_W≤1).

## P5 — DOF→family-row map
Committed comparator `scripts/S11c_b_cross_engine_comparator.py` `extract_slab` (~L760-799):
`U_BODY_BALANCE→(U_MOMENTUM_ROWS, U)`, `THETA_BALANCE→(MU_THETA, THETA)`,
`E_W_BALANCE→(THICKNESS_ROW, E_W)`, `ADVECTIVE_MASS_OPERAND→(MASS_EVOLUTION_ROW, MASS)`.
No transverse/longitudinal split at the strong-row level (that split lives only in
`COUPLING_KERNEL`, `...sympy_audit.py:1724-1754`).

## P6 — the completed basis is 15/source, all independent (the #86 result, supplied)
`research/pde_ledger_v3/directives/_measurements/S11c_b_86_reference_result.md`: corrected
§3a basis = 40 = 10 uniform + 15 `∂W_bg` + 15 `∂μ_R,bg`; per-source 15 exact (nullity 0);
reduces to the engine's committed frozen 26 (anchor
`grep -a S11CB_ENERGY_BASIS_COUNT scripts/out/S11c_b_brane_operator_sympy_audit.out`
→ `Integer(26)`, both anchorings).

## P7 — controls are satisfiable harness invariants
The reconstruction identity `RESIDUAL − (UNIFORM_HESSIAN + SELECTED_SPURION_HESSIAN +
OMITTED_CORRECT) == 0` holds exactly by linearity (verified in
`~/.s11_build/S11c_b_88_probes/probe_decomp.out`, which found the missing cross-term). The
freeze+drop-omitted reduction to the engine density is reproducible
(`~/.s11_build/S11c_b_88_probes/probe_blast_radius.out` → `Control A ... PASS`), but that
freeze-vs-freeze form is TAUTOLOGICAL (both decision legs, CAS) and was replaced by
`CONTROL_ENGINE` (instrument EL vs the engine's committed `operator_from_density`/`EXPANDED`
rows, ε-stripped, inertia-free).

## P8 — engine EL normalization + acceptance-tool existence (folded from the decision legs)
- The emitted engine rows carry `epsilon` and an inertia term:
  `sed -n '1465,1491p' scripts/S11c_b_brane_operator_sympy_audit.py` →
  `u_local = tuple(epsilon*diff(density,u[a])…)`, `… − (epsilon*rhobr*u_tt[a] if include_kinetic…)`,
  `theta_expanded = epsilon*mu_theta_amplitude`, `… − (epsilon*mu_W*W_bg**2*e_tt …)`. The
  inertia carries no background spurion ⇒ cancels from every residual; the directive works with
  the ε-free, inertia-free amplitude and reconciles in `CONTROL_ENGINE`.
- `find . -name derived_or_declared.py` → **absent** (repo has `reduction/s11_*` only). The
  original acceptance item citing it was unsatisfiable; removed.
- The abstract spurion symbol is `M.bg` (`scripts/…_sympy_audit.py:1014`); enumeration is on
  `M.bg`, live substitution per source is applied only after index selection.

## Provenance — 2 decision legs run BEFORE any builder (rule 7); folded once (rule 7)
Orchestrator-written directive → Codex + Grok. Prompt
`directives/_legs/S11c_b_88_blast_radius_directive_review.md`; raw transcripts
`~/.s11_build/S11c_b_88_decision_{codex,grok}.txt`. Both verified the premises + Hessian-σ_W¹
retention in CAS (no answer leak found) and converged on the material fixes now folded:
three-disjoint-driver decomposition + reconstruction identity; constant-coefficient
absorbability span (DOF-only rank invalid); non-tautological `CONTROL_ENGINE`; ε/inertia;
`enumerate(M.bg)` + index convention; sentinel/Jacobian controls; LAB_HELD/PY-witness scope.
Folded once, then build (not iterated to green).
