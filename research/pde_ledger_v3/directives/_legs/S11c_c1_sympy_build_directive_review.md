# Independent review — S11c-c1 SymPy BUILD DIRECTIVE (orchestrator-written §-wiring + an export-narrowing change)

You are one of two independent review legs (the other is a different engine). This is a **document review** of
an orchestrator-written build directive. Derive your own view from the source specs and rules FIRST, then read
the directive and check it. Quote both sides for every finding. ⛔ A prose "looks fine" is worthless — for any
claim about what an object contains or what the next step needs, cite the exact spec line / LEDGER row / rule
that grounds it.

## Artifact under review
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c1_sympy_build_directive.md`
(committed `d1fe1bf0`). It instructs Codex to build the SymPy engine
`scripts/S11c_c1_bulk_closure_sympy_audit.py` and its export `scripts/S11c_c1_exports.py`. It claims to **add no
physics** — to be pure chain-wiring over the spec. Its most consequential clause is an **export-narrowing**
(the "EMIT vs EXPORT" section): it exports to the plain-git LEDGER only a strict subset of what the engine
emits.

## Source of truth — form YOUR view from these BEFORE reading the directive
- `research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md` — the physics authority (the spec the directive
  wires). Esp. **§0** (what is reserved for S11c-c2 vs c1), **§3b** (the closing sentence naming what c2
  consumes), **§4** (the full emit set), **§7** (chain output — "the flat plain-git LEDGER of the emitted
  `PY_S11CC1_*` objects", BUILD_INPUT_DIGESTS, F9).
- `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` — S11c-b's slab operator with its **symbolic
  `δp_s`/`Λ_I(ω)` face slots**; and what S11c-c2 (the self-energy fold) will substitute into those slots.
- `research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` — the export-chain rules **F1–F9**; and
  `research/pde_ledger_v3/S9_REWRITE_PLAN.md` **D1** ("export what the primary derivations emit").
- The LEDGER reality: `research/pde_ledger_v3/scripts/{S11b,S11c_a,S11c_b}_exports.py` (what is already in the
  chain and can be bound vs must be freshly exported).
- The comparator: `research/pde_ledger_v3/scripts/S11c_b_cross_engine_comparator.py` (what channel it reads).

## What to check — findings only where the directive would produce a BROKEN engine, a BROKEN chain, a WRONG export, or a LEAKED answer

1. **UNDER-EXPORT (the dangerous one — probe hardest).** The directive's export keep-set is
   `{dtn_flat_symbol, dtn_operator, dtn_kernel, face_response, face_response_coeffs}` + their free symbols +
   the carried-forward inherited LEDGER. Independently determine, from the c1 spec §3b + §0 and the **S11c-b
   slab operator's symbolic face slots** (what S11c-c2's fold actually substitutes), the EXACT set of c1 objects
   S11c-c2 must import. Is anything c2 needs **missing** from the keep-set? Consider specifically: does c2 need
   the noninvertibility condition, any `dtn_by_regime_pair`/`dtn_by_parity` structural view, any dissipation or
   energy object, or the Hermitian/reactive split — or does the spec have c2 **re-derive** those on the coupled
   operator (§0/§3b)? Cite the spec line that settles each. A single missing consumed object is a must-fix.
2. **FREE-SYMBOL CLOSURE.** Do the kept objects (`dtn_kernel`, `dtn_operator`, `face_response`,
   `face_response_coeffs`, `dtn_flat_symbol`) reference any free symbol that is **neither** already inherited in
   `S11c_b_exports.py` **nor** in the directive's listed new-symbol set (`q_out, kappa_out, zeta_out, k_out_*,
   k_in_*, V_±, Lambda_*, tau_*`, the profiles, density representatives)? If a referenced symbol would be
   exported by neither route, the `_restore` reviver / c2 import breaks (the engine's D3 round-trip should catch
   it, but flag any gap you can see structurally).
3. **OVER-EXPORT / POLICY CONSISTENCY.** Is the directive's "refine D1 → export what the next step CONSUMES"
   consistent with the committed export-chain rules (F1–F9, D1) and spec §7's "LEDGER of the emitted objects"?
   §7 reads broad, §3b reads narrow — does the directive's resolution contradict §7 or any F-rule, or is it a
   legitimate refinement? If it contradicts a committed rule, that is a finding (the rule wins unless the spec
   overrides it).
4. **EMIT-vs-EXPORT SOUNDNESS.** The directive claims narrowing the export loses nothing for review because the
   cross-engine comparator and the WL engine work off the **stdout `.out`** tag streams, not the exports LEDGER.
   Verify: does `S11c_b_cross_engine_comparator.py` read the `.out` streams (not `_exports.py`)? Does the WL
   engine emit any exports LEDGER at all? If the comparator or any downstream consumer actually needs the
   exported rows, narrowing is unsafe — flag it.
5. **NO ADDED PHYSICS / NO LEAKED ANSWER.** Does the directive state any expected value, parity, sign, regime
   outcome, coefficient, or residual — or a prohibition that reveals the answer (rule 5 / the "a prohibition
   leaks" rule)? It should be pure wiring. (Two known-benign items: `v_bulk_normal_0=0` is the supplied grazing
   *domain limit* the directive marks step-record-only; "passive/lossless" appears only in the spec-verbatim
   prohibition against forcing that form. Flag anything else.)
6. **CHAIN BINDINGS.** `μ_θ` carried as the **opaque** imported `mu_theta_operator` (not expanded — that is c2);
   inherited constants **bound** to LEDGER rows (F9b corroboration of a re-derived constant — same symbol/key,
   accreting `corroborated_steps` — is correct, not re-origination); `BUILD_INPUT_DIGESTS` pins the **five**
   inputs (own audit + S11b/S11c_a/S11c_b exports + spec) per §7; F9c prefix `s11c_c1_` non-colliding with
   `s11_/s11b_/s11c_a_/s11c_b_`; the audit-script product name; F6 publish gated on §4 completion. Any binding
   the directive gets wrong (or an inherited object it tells the engine to re-originate) is a finding.
7. **AMBIGUITY THE BUILDER WOULD FILL WITH AN INVENTED DECISION.** Any place the directive is under-specified
   such that Codex would have to invent physics, an expected result, or a self-review mechanism to proceed.

## Method + filter
- Read the sources above, write down your own answer to "what does c2 consume from c1" and "does the keep-set's
  free-symbol closure hold" **before** trusting the directive's own claims. Cite lines.
- ⛔ Do not modify any file. This is a read-only document review.
- Report a finding only if it catches a concrete way the directive is wrong (broken engine/chain, wrong export,
  leaked answer, unbindable object). Not style, not "could be clearer." If you find nothing on a check, say so
  and name what you verified.

## Output
Findings ranked most-severe first (must-fix vs nit), each quoting (a) the spec line / rule / LEDGER row you
derived your expectation from and (b) the directive line it contradicts or satisfies. End with a one-line
verdict: is the directive (especially the export-narrowing) safe to build against as-is, or are there must-fix
items first?
