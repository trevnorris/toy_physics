# Independent BUILD review — S11c-c2 N6 companion diagnostic (a SCRIPT)

## Artifacts
- **Companion (primary):** `/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py`
- **Audit edit (secondary):** the change in `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` — inspect via
  `git -C /var/projects/toy_physics diff a24acfb6~1 a24acfb6 -- research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py`
  (it must replace the obsolete cross-anchoring `REP_INVARIANCE_*` loop with `ANCHORING_L_MINUS_M = increment[LAB] −
  increment[MATERIAL]`, raw, no `representation_pullback`, keeping the density loop).

Astra(Codex)-written diagnostic that computes the c2 **N6 representation-invariance control**. Working dir
`/var/projects/toy_physics`. This is a **SCRIPT review** — derive independently, then ABLATE. ⛔ You are NOT handed,
and there is NOT, any expected value for the N6 residual — §5c/the spec fix no target.

## The governing physics (the script must implement THIS — derive from it yourself, ⛔ not from the script's output)
- `directives/S11c_c2_SHARED_PHYSICS.md` §5c (the object), §3a/§3c/§1d/§5d/§6.
- The **cleared route-2 construction spec** `research/pde_ledger_v3/_measurements/S11c_c2_N6_route2_spec_astra.md`
  (fresh-Claude + Grok CLEAR) — the authoritative construction. Key settled points:
  - route 2 = `extract(close(SLAB_M) − SLAB_M)` with **NO `T` on the increment** (native S11c-a MATERIAL face sources
    folded into the same δp symbols, differenced directly against route 1 — the parent `task_rep_invariance` pattern);
  - **material μ_θ** bound at BOTH the face substrate AND the c1 response source (route 1 = Eulerian μ_θ);
  - the **reverse-u channel is grade-suppressed** — the script must REPORT reverse-block sensitivity and PERMIT
    computed absence, ⛔ NOT force a reverse channel, ⛔ never credit `extract`'s U-row curl as an N4 witness;
  - residual + guards + controls tested by **exact finite-field PIT**; residual-zero is ⛔ never an exit/assert.

## Required method — SCRIPT branch (derive, then ablate; a prose re-derivation is worth nothing)
1. **Derive independently.** Write your OWN derivation/probe script BEFORE trusting the artifact; save it + its literal
   stdout to named /tmp paths and report those paths. ⛔ Without them your derivation claims are discarded.
2. ⛔⛔ **FORM ABLATION IS MANDATORY.** Change the STRUCTURE of a load-bearing object (flip a sign AND an off-diagonal;
   collapse two independent symbols into one; drop a δp slot) in a /tmp COPY, re-run, and report the LITERAL diff. A
   COEFFICIENT rescale tests arithmetic; only a FORM change tests physics. If ablating a load-bearing check produces
   BYTE-IDENTICAL output, that object is not computed — report it.
3. **One-sided route corruption (the independence test).** Corrupt ONE route only (e.g. perturb route 2's material
   source, or route 1's Eulerian carrier) and report which operands move. If breaking route A moves route B's operand,
   the two routes were never independent and `R_N6` is decoration. Verify the two routes are structurally distinct
   constructions (route 1 imported-Eulerian; route 2 native-material) — ⛔ not one derived from the other.
4. **Probe for the classic defects:** a value verified with the predicate that produced it; a conclusion emitted as an
   unconditional literal; a hand-typed CAS payload with no data-dependence (delete the derivation → does the emit
   move?); an answer-bearing tag NAME; a suppressed-identical payload; an `assert` BEFORE the emit it guards (report
   any). For every emitted object ask **WHICH LINE COMPUTED THIS** — give the line or report it uncomputed.

## N6-specific checks (catch a WRONG or VACUOUS control)
- **No `T` / no annihilation:** route 2 must be the native material face-fold (not a θ-shift / `material_pullback` on
  rows / a re-pullback of route 1); `I_{M→E}` must be a nonzero linear carrier, not annihilated to ≈0.
- **μ_θ = material at both sites**, δp symbols shared with route 1 (so `I_E − I_{M→E}` is meaningful), c1 response
  opaque (no θ-sub / Jacobian into `DELTA_P`/`Z`/resolvent).
- **Finite-field PIT soundness:** several primes; ≥ several generic draws/prime/cell shared across routes; joint
  singular rejection; branch-cell selection BEFORE modular reduction; on-shell `q` NOT sampled independently of `k`;
  degree bounds; an honest false-negative bound (bad-prime condition handled — ⛔ not "multiple primes ⇒ done"); a
  certified nonzero is decisive, all-zero is only a bounded conditional pass; residual-zero ⛔ never an exit.
- **Controls:** advection = tag `t` on the material bulk-θ-advection `u·∇ρ₄/ρ₄` inside material μ_θ, bound at both
  sites, one route only; **RHO4_CONSTANT** structurally absent (emit the computed absence, ⛔ not A−A), **RHOBR_CONSTANT**
  live. Tilt = an injectable Eulerian carrier factory PIT-checked vs the imported carrier (⛔ not the `FLIP_FACE_SLOPE`
  c1-DtN-jet override). Each control must MOVE its own route while leaving the other unchanged.
- **Emit / dimensions:** operands PRINTED before residual; `[L,T,M]` + able-to-fail dimensional consistency +
  `(ε,η,σ_W)` multigrade on every object; ⛔ no giant symbolic dumps. ⚠ `reduction/derived_or_declared.py` could NOT
  classify the emit JSONL format during the build — assess whether the emit schema is correct per §3c/the emit
  contract, or a real defect.
- **Audit edit:** raw `LAB − MATERIAL`, no `representation_pullback`, density loop kept.

## Ablation sandbox / ops
⛔ Copy the artifact to /tmp and ablate the COPY; ⛔ never modify the working tree. Pure SymPy (no Mathematica). ⚠ The
finite-field PIT is per-case compute-heavy — run **one case at a time**, wrap each run in `timeout 600`, ⛔ never run a
full-symbolic zero-test over all grades (that is the F/G wall the design avoids). Save every ablation script + its
literal stdout to named /tmp paths and report them.

## Physics filter
Report a finding only if it catches a way the N6 physics could be **wrong, vacuous, intractable, or answer-leaking** —
not "wrong on a different input."

## Output
Findings each with the line, the ablation + its literal diff, why it matters, minimal fix. If nothing outstanding
changes what is computed or may be claimed, say **BUILD CLEAR**. Evidence-first, brief.
