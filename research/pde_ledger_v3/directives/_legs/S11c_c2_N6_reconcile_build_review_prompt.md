# Independent BUILD review — S11c-c2 N6 reconcile instrument (a SCRIPT)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py` (astra-written, 300 lines).
Working dir `/var/projects/toy_physics`. This is a **SCRIPT review** — derive independently, then ABLATE.

## Role + the ONE boundary you must respect
The instrument adjudicates whether the per-engine N6 residual `R_N6 = I_E − I_{M→E}` **vanishes modulo the fixed
material↔Eulerian defining relations**. ⛔⛔ **Your job is to verify the instrument COMPUTES the reconcile correctly —
NOT to decide whether `R_N6` vanishes.** The disposition (does the residual vanish / is there an unexplained
remainder) is the orchestrator's adjudication AFTER your review. ⛔ You are NOT handed, and there is NOT, any expected
value for `R_N6` or any bridge/channel residual. Report a defect only if the instrument could compute the reconcile
**wrongly, vacuously, or leak-prone** — ⛔ never "the residual came out nonzero/zero" (that is not yours to judge).

## The governing construction (implement THIS — derive from it, ⛔ not from the script's output)
- Directive: `directives/S11c_c2_N6_reconcile_directive.md` (CLEAR TO BUILD). Cleared route-2 construction:
  `_measurements/S11c_c2_N6_route2_spec_astra.md`. Object + finding: §5c of `directives/S11c_c2_SHARED_PHYSICS.md`
  and the N6 diagnostic `scripts/S11c_c2_N6_diagnostic_sympy.py` (the instrument imports its machinery).
- Settled points to verify the script implements:
  - **The corrected question:** N6 requires `[I_E − I_{M→E}]_retained = 0` **in the quotient by the fixed defining
    relations** — ⛔ NOT `R_N6 = J` (both operands already in the common Eulerian basis; a subtractable "justified J"
    is leakage). The evidence is the residual factoring through the defining relations, tested via two bridges.
  - **Carrier bridge** `CARRIER_BRIDGE_RESIDUAL = C_E − C_M` at the **coefficient** level (imported `e_coeff` :797 vs
    material `m_coeff` :809, row/face/slot/grade) — ⛔ not `Increment(ΔC,S)`.
  - **Source bridge** `SOURCE_BRIDGE_RESIDUAL = es − ms` where `es`/`ms` **are the diagnostic's `source_terms`
    circuits** (:378-389, :841-843) — ⛔ not a re-coded `Λ` expression (that double-counts ε via `V̄=V/ε` or misses
    imported slot factors). The `b_{r,s}` formula is a slot identification only.
  - **The three-way affine split** (`build_increment` is affine in S: signature 0 = bare `−C·p` :491-506):
    `R_N6 = CARRIER + SOURCE + CROSS` with `CARRIER=I(ΔC, ms)` (via `build_increment` on the carrier diff at the
    MATERIAL source), `SOURCE=B(C_M, ΔS)`, `CROSS=B(ΔC, ΔS)` — SOURCE/CROSS are **signatures-{6,9,12}-only
    contractions with NO bare term** (⛔ `build_increment(·, ΔS)` would re-add `−C·p`). `SPLIT_CHECK =
    CARRIER+SOURCE+CROSS − R_N6`.
  - **Frozen defining relations** predeclared + emitted, including the jet-vocabulary bridge `a.grad_theta[i]
    (theta_d1) ↔ b.grad_theta[i] (grad_theta_1)` (route-2 :128; diagnostic `N6_JET_BRIDGE` :351-356); fixed-anchoring
    (⛔ never bridging `LAB_HELD ↔ MATERIAL_ADVECTED`).
  - **PIT:** one in-process joint `pit()` per case over ALL reconcile objects (shared samples); ⛔ NOT a cross-run
    join of stored diagnostic PIT tables. Nonzero modular numerator = one-sided certificate; all-zero = "no nonzero
    found" at conditional δ (⛔ never "certified zero"); residual-zero is ⛔ never an exit/assert.

## Required method — SCRIPT branch (derive, then ablate; a prose re-derivation is worth nothing)
1. **Derive independently.** Write your OWN derivation/probe script BEFORE trusting the artifact; save it + its literal
   stdout to named /tmp paths and report those paths. ⛔ Without them your derivation claims are discarded. Confirm the
   affine identity `I(C,S) = −C·p + B(C,S)` and that `R_N6 = I(ΔC,S_M) + B(C_M,ΔS) + B(ΔC,ΔS)` on the actual objects.
2. ⛔⛔ **FORM ABLATION IS MANDATORY.** In a /tmp COPY, change the STRUCTURE of a load-bearing object (e.g. make SOURCE
   include signature 0 by calling `build_increment` instead of the {6,9,12} contraction; swap `ms`→`es` in CARRIER;
   collapse `C_E` and `C_M`), re-run, report the LITERAL diff. A byte-identical output for a load-bearing change means
   that object is not computed — report it. **In particular: corrupt ONE channel and confirm `SPLIT_CHECK` moves** (the
   split identity must be load-bearing, not a tautology). A COEFFICIENT rescale tests arithmetic; only a FORM change
   tests physics.
3. **The one-sided bites (independence + able-to-fail), each on a /tmp COPY:**
   - **Carrier bite:** corrupt the MATERIAL covector/normal-map feeding `C_M` (`material_inverse_transpose` / material
     normal, route-2 :113) with `ms` held ⇒ must move `CARRIER_BRIDGE_RESIDUAL`/`CARRIER_CHANNEL`, leave `C_E`, the
     Eulerian operand, and `SOURCE_BRIDGE_RESIDUAL` unchanged.
   - **Source bite (`RHOBR_CONSTANT` only):** `a_ρ→0` on one material constitutive route ⇒ must move
     `SOURCE_BRIDGE_RESIDUAL`/`SOURCE_CHANNEL`, leave `CARRIER_BRIDGE_RESIDUAL` and the Eulerian operand unchanged.
     For `RHO4_CONSTANT` (`g_i=0⇒a_ρ=0`) confirm the instrument emits computed ABSENCE, ⛔ not an A−A.
   If a bite does NOT move its own residual, or moves the OTHER route, the routes were never independent — report it.
4. **Probe the classic defects:** a value verified with the predicate that produced it; a conclusion emitted as an
   unconditional literal; a hand-typed CAS payload with no data-dependence (delete the derivation → does the emit
   move?); an answer-bearing tag NAME; a suppressed-identical payload; an `assert` BEFORE the emit it guards. For every
   emitted object ask **WHICH LINE COMPUTED THIS** — give the line or report it uncomputed.
5. **PIT soundness:** shared samples across all objects; several primes; ≥ several draws/prime/cell; joint singular
   rejection; on-shell `q` derived from `k` (⛔ not sampled independently); an honest false-negative bound (bad-prime
   handled — ⛔ not "multiple primes ⇒ done"); residual-zero ⛔ never an exit. `SPLIT_CHECK` samplewise-zero (⛔ not a
   `number(0)` node, ⛔ not a structural node-identity).

## Ablation sandbox / ops
⛔ Copy the artifact to /tmp and ablate the COPY; ⛔ never modify the working tree. Pure SymPy (no Mathematica). ⚠ The
PIT is per-case compute-heavy — run ONE case at a time, wrap each run in `timeout 900`, ⛔ never a full-symbolic
zero-test over all grades. astra's baseline `.out` are at `/tmp/S11c_c2_N6_reconcile_sympy.<ANCHORING>.<DENSITY>.out`
(you generate your own from the copy). Save every ablation script + its literal stdout to named /tmp paths and report
them.

## Physics filter
Report a finding only if it catches a way the reconcile INSTRUMENT could be wrong, vacuous, intractable, or
answer-leaking — ⛔ not "wrong on a different input", ⛔ not the disposition of `R_N6`.

## Output
Findings each with the line, the ablation + its literal diff, why it matters, minimal fix. If nothing outstanding
changes what the instrument computes or may be claimed, say **BUILD CLEAR**. Evidence-first, brief.
