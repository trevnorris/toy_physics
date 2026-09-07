# Decision review — the S11c-c2 N6 RECONCILE build directive (physics-bearing pre-builder directive)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_directive.md` — an
orchestrator-written directive for a diagnostic SCRIPT (astra will author `scripts/S11c_c2_N6_reconcile_sympy.py`)
that adjudicates whether the per-engine N6 residual `R_N6 = I_E − I_{M→E}` **vanishes modulo the fixed
material↔Eulerian defining relations** (⇒ representation invariance holds) or leaves an unexplained retained-order
discrepancy (⇒ a real finding). Working dir `/var/projects/toy_physics`; paths under `research/pde_ledger_v3/`.

## Role + what this review is (and is NOT)
This is a **pre-builder directive review of physics-bearing content**. Review the requested DECISIONS + their
supporting evidence + the physics. ⛔ Do NOT ablate a fictional script (the deliverable does not exist yet); ⛔ do NOT
run CAS or a build. Executable script-control tests (FORM ablation, one-sided corruption, PIT soundness on the actual
script) are the **build legs'** job, ⛔ not this review's — a directive review never pays the build tax. Your job: is
the directive's physics right and complete, is anything under-specified or leak-prone, before astra builds against it.

## Method — DOCUMENT branch (read the SOURCES first, then the directive)
Read, and form your own view, BEFORE reading the directive's claims:
- The finding + object: `scripts/S11c_c2_N6_diagnostic_sympy.py` — esp. `run` (:788-899), `build_increment`
  (:476-513, note signature 0's bare `−p` term :491-506), `residual` (:576), `constitutive` (:318-360), the shared
  sampling in `pit` (:668-773).
- The cleared route-2 construction `_measurements/S11c_c2_N6_route2_spec_astra.md` — esp. §2 pullback (:39-83), §3
  native face fold + pressure identities (:87-128), §4 what survives pressure projection (:130-150), §5 close +
  combined source `b_{r,s}` + carrier `S_{P,r}=Σ C_{r,p}p` + increment `J_r` (:152-178).
- `directives/S11c_c2_SHARED_PHYSICS.md` §5c (:303-395).
- The E1 question-vet that produced the corrected question: `_legs/S11c_c2_N6_reconcile_question_vet.md` and the
  Codex-sol response it was answered with (you may derive the same conclusions independently).

## What to check (substantiate each independently against the sources)
1. **The corrected question.** Is "does `R_N6` vanish modulo the fixed defining relations, at retained order
   `(η^{≤1},σ_W^{≤1})`" the correct operationalization of N6 — and is the directive right that there is **no nonzero
   "representation offset J" to subtract away** (both operands already in the common Eulerian basis), so an `R_N6 = J`
   framing would be author-freedom leakage? If the correct test is different, state it.
2. **The frozen defining relations** (directive "frozen defining bridge relations"). Are `g_i`, `a_ρ`, `h_α`, the
   field maps `θ↦θ+a_ρ`, `e_W↦e_W+h_α`, their prolongations, the Jacobian (quadratic-projection-only), the
   covector map, and the RHO4/RHOBR laws COMPLETE and faithful to route-2 spec §2? Is anything a material↔Eulerian
   identity that belongs in the frozen set but is missing (which would let a real discrepancy be mislabeled), or an
   over-inclusion (a relation NOT sanctioned, which would wrongly explain away a discrepancy)? Confirm the fixed-
   anchoring restriction (⛔ never bridging `LAB_HELD ↔ MATERIAL_ADVECTED`).
3. **The two discriminating bridges.** (a) Carrier `C_E = C_M` tested at the **coefficient** level (row/face/slot/
   grade) — is the directive right that `Increment(C_E−C_M,S_E)=0` alone is weaker (kernel/source nullspace)? (b)
   Combined-source `b_{E,s}=b_{M,s}` with `b_{r,s}=(1+Λ_V/ρ_m)V̄_{r,s}+(Λ_A/(ρ_m ρ_br^bg))μ̄_r` — is this the right
   discriminating object and formula (route-2 spec :157-159), and correctly kept OUTSIDE the opaque c1 response
   `Z[b]`? Are these two bridges together sufficient to localize a discrepancy (carrier vs source vs cross-term)?
4. **The exact affine split.** Verify `build_increment(C,S)` is affine (not linear) in `S` via the signature-0 bare
   `−C·p` (:491-506); confirm `R_N6 = I(C_E−C_M,S_E) + B(C_M,S_E−S_M)` with **term 2 a signatures-{6,9,12}-only
   contraction** (no re-added `−C_M·p`), and that calling `build_increment(C_M, S_E−S_M)` would be WRONG. Is the
   `SPLIT_CHECK = TERM1+TERM2−R_N6` an adequate structural guard on the split identity?
5. **Emit contract + leakage.** Any supplied expected residual value, any residual-zero exit/assert/loop predicate,
   any tag NAME that leaks a value/sign, any A−A control, any place emission is conditioned on a payload's value? Is
   the shared-sample PIT + one-sided-certificate semantics (nonzero = discrepancy certificate; all-zero = "no nonzero
   found" at conditional δ, ⛔ never "certified zero") correctly stated?
6. **Under-specification for astra.** Anything astra would have to guess that could change what is computed — the
   term-2 contraction path, the `C_{M→E}` carrier construction (native material factory, already covector-mapped),
   the bridge keying, the "independent of `R_N6`" construction, or the able-to-fail control structure?

## Physics filter
Report a finding only if it catches a way the reconcile could be **wrong, vacuous, leak-prone, mislabel a real
failure as invariance (or vice versa), or intractable** — not "wrong on a different input."

## Output
Findings each with `file:line`, why it matters for what is computed or may be claimed, and the minimal fix. End with
**DIRECTIVE SOUND — CLEAR TO BUILD** or the exact fold list. Brief, evidence-first, grounded in the cited lines.
