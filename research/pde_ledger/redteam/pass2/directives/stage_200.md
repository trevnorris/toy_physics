---
unit_id: 200
batch: V.3
created_at: 2026-06-09T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-09T13:51:17-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
orchestrator_authored: true
---

# Codex directive — unit 200 (CHECKPOINT)

This directive was authored by the orchestrator after the pass-2 ground-truth `.wl`-vs-`.py`
read overturned the audit agent's "borderline → accept" disposition on the checkpoint. The
user has authorized the re-author (re-author-vs-accept is a user-level call).

**You DESIGN and WRITE the new Mathematica route.** This directive states the REQUIREMENT and
the ACCEPTANCE CRITERIA only — it deliberately does not prescribe the specific algorithm. Pick
an independent method, justify it in your `## Applied` note, and make the script run clean.

Apply F1. F2 is informational (no Codex action). After applying, append an `## Applied: F1`
block with `files_changed`, `summary` (one sentence), `route_chosen` (the independent method you
used for each section), and `deviation` (or "none").

Do NOT touch the SymPy script `scripts/moving_throat_pde_stage200_reference_free_home_stretch_theorem_sympy_audit.py`
— it is the reference. Do NOT touch paper.tex, notes/, or any prose document. Edit only the `.wl`.

After editing, RUN `timeout 600 math -script <the .wl>` and iterate until it exits 0 with every
`expectZero`/`pass` check passing. A timeout (exit 124) is a FAILURE — reformulate the math, never
raise the cap.

## F1 — mathematica_transliteration (load-bearing; re-author the `.wl`)

**Severity:** high (checkpoint; dual-engine independence is the whole point of the second engine)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl`

**What's wrong:** the `.wl` is a near-line-by-line port of the `.py`. Every load-bearing object
is produced by the SAME operation in both engines, on the same intermediate quantities:

- **§I — the load-bearing compiler matrix `M_*`.** Both engines build the three physical
  monomials, apply the SAME `ratioSubs` dict, the SAME `ratioToLogs` (ratio→`Exp[Δ]`) dict, form
  `qPair = Log[ratio]`, and then obtain `M_*` by **symbolic autodiff of those same logs** —
  SymPy `q_pair.jacobian(Dvec)` (`...sympy_audit.py:154`) vs Mathematica
  `Table[D[qPair[[i]], Dvec[[j]]], …]` (`...mathematica_audit.wl:157`), both checked against the
  SAME hard-typed `Mexpected` literal (`py:155-161` / `wl:158-162`). `.jacobian` vs `Table[D[…]]`
  is the same differentiate-the-log operation — not an independent derivation.
- **§III — the dependent-triple orbit.** `KetaOrbit`/`TOrbit`/`muOrbit` are **posited as the same
  hand-written closed forms** in both engines (`py:286-298` / `wl:223-231`); neither solves the
  orbit system. (Contrast verified siblings 198/199, whose `.wl` actually `LinearSolve`/`Solve`s
  the orbit where the `.py` posits it.)
- **§V — the Packet-A linearization.** Both linearize `chiFromDef` by `Series[…,{eps,0,1}]`
  (`py:409-420` / `wl:281-292`) — same series operation, same `…Expected` literal.

§II (witness invariance) and §IV (cocycle) are downstream consequences of §I's `M_*` and the
shared monomial definitions; once §I/§III/§V derive their objects independently, §II/§IV inherit
the independence so long as they are built from the independently-derived `M_*` and monomials.

The shared physical MONOMIAL DEFINITIONS (`ctrMonomial`/`cntMonomial`/`epsEtaMonomial`) are the
premise and may stay — they encode the physics, not a method choice. The defect is that
everything DOWNSTREAM of the monomials (how `M_*`, the orbit, and the linearization are extracted)
mirrors the `.py` operation-for-operation.

**Why this matters:** a checkpoint's second engine must independently verify the result, not re-run
the first engine's algebra in a second CAS. The first pass already "de-transliterated" this `.wl`
once (hand-collapsed ratios → the current Jacobian route) and it remained a port — so the bar is an
operation genuinely different from the `.py`, not merely a different simplifier.

**Required change (you design the route):** re-author the `.wl` so each load-bearing object is
derived by a method genuinely INDEPENDENT of the `.py`'s route:

1. **`M_*` (§I)** must NOT be obtained by an autodiff Jacobian (`Table[D[…, Dvec[[j]]]]` / `D[…]`)
   of the same `Log[ratio /. ratioToLogs]` quantities the `.py` differentiates with `.jacobian`.
   Extract the matrix by a different operation. (Feasible routes you may choose among — design the
   actual code yourself: read each row's exponents directly off the expanded/`PowerExpand`'d logs
   via `Coefficient`/`CoefficientList` over the `Δ` symbols; or assemble the rows from the
   primitive exponent-weight vectors as verified sibling 199's `.wl` does. The V.2-185 re-author is
   the precedent for a `Coefficient`-based exponent read.)
2. **The orbit (§III)** must be SOLVED, not posited: obtain `KetaOrbit`/`TOrbit`/`muOrbit` from the
   orbit/closure equations (e.g. `Solve`/`LinearSolve` on the log-linear system), the way 198/199's
   `.wl` do — do not retype the `.py`'s closed forms.
3. **The Packet-A linearization (§V)** should be obtained by an operation other than `Series` on the
   same `chiFromDef` (e.g. an analytic base-point derivative `D[…,eps]/.eps->0`), so it is not the
   `.py`'s Series re-typed.

**Preserve EXACTLY (these are the checkpoint deliverables — they must not move):**
- the carried Stage-192 matrix `Mexpected` (`wl:158-162`) and the assertion `M_* − Mexpected == 0`;
- the witness-invariance zeros (§II): `Log[Ctr/Cnt/epsEta(witness) / (*)] == 0`;
- the mismatch chart `q = ((1+chi0_*) Log[m_T], Log[m_mu] − Log[m_K] − F_* Log[m_T], −Log[m_K])`
  (§III) and `M_* · Δx_mis − q == 0`;
- the cocycle `q^(31) − q^(32) − q^(21) == 0` (§IV);
- the Packet-A linear law `Delta_Q^lin = eps(5 eps_beta + dSigma0/(3 S) + 9 dSigma5/S)` (§V).

No checkpoint constant changes; every printed deliverable and every `expectZero` target stays the
same value. Only the METHOD that produces them changes.

**Verification (how the verifier confirms the fix landed):**
- the `.wl` no longer contains an autodiff-Jacobian of the `ratioToLogs` logs feeding `M_*`
  (no `Table[D[qPair…, Dvec…]]` / equivalent `D[]` of those logs);
- `M_*` is produced by an operation distinct from the `.py`'s `.jacobian` (Coefficient-read /
  weight-vector assembly / etc.);
- §III's orbit triple is produced by a `Solve`/`LinearSolve`, not retyped from the `.py`;
- the script exits 0 with all checks passing; every preserved deliverable above still reads `= 0`
  / `PASS` in the refreshed committed output.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage200_reference_free_home_stretch_theorem_mathematica_audit.wl`
- summary: Re-authored the Stage 200 Mathematica audit so the compiler matrix, dependent orbit, and Packet-A linearization are derived by independent methods while preserving the checkpoint checks.
- route_chosen: §I assembles `M_*` from primitive monomial exponent weights; §III solves the dependent log-invariant system with `LinearSolve`; §V computes the base-point `eps` derivative instead of using `Series`.
- deviation: none

## F2 — stale_output (informational; no Codex action)

Both committed `.txt` outputs predate their scripts and carry the stale pre-renumber banner
`STAGE 183`. This is resolved by the orchestrator's independent re-run after your re-author (which
regenerates both `.txt` with the canonical `STAGE 200` banner). Do not hand-edit the `.txt` files.
