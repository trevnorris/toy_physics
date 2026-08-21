# S11b SymPy engine — fix round 1 (emission fidelity)

## Authority and boundary

Repair the existing engine `research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py`
(round-1 baseline committed `7dd89076`) in place; its products stay the same two — the flushed stdout tag
stream and `research/pde_ledger_v3/scripts/S11b_exports.py`.

⛔⛔ **The physics is already correct and ⛔ must NOT be re-derived or re-valued.** Two independent
derivations (the round-1 script legs) reproduced every load-bearing object. `directives/S11b_SHARED_PHYSICS.md`
(`1a2395a3`) remains the sole physics authority and `directives/S11b_sympy_build_directive.md` (`9bd2f184`)
the wiring authority. This directive changes **only two emission classes** named below. `CLAUDE.md` binds:
the three script clauses and the anti-tautology corollaries (`.claude/skills/build/SKILL.md`) apply — ⛔ point
at them, do not restate. Add no expected value (rule 5).

⛔ **Everything not named in Fix 1 / Fix 2 stays byte-identical in behaviour.** ⛔ Do not change any
physics-object **value**, and ⛔ do not touch the export wiring (the S11-LEDGER import and carry-forward, the
object bindings, `F9`, `BUILD_INPUT_DIGESTS`, the `D3` round-trip, `_RELATIONALS`, the `MappingProxyType`
freeze). The orchestrator re-verifies carry-forward, wiring, and every physics-object value against `7dd89076`
after the build; a change there is a regression. ⚠ Fix 1 recomputes the impedance but the computed object
**equals** the current typed one, so ⛔ no impedance-dependent value may change — a change surfaces as a
finding.

## The decisive test used throughout — ONE-SIDED CORRUPTION, not "the row moves"

For every claim below that an object is genuinely computed / a check genuinely fires, the test the
orchestrator runs (and the property the repair must satisfy) is **one-sided corruption** (`.claude/skills/build/SKILL.md`
corollary 3): corrupt **only** the route the object/check claims to depend on, and the object's **value** —
or the check's **residual/status**, i.e. the *tautological element itself* — must move, while its independent
reference route stays unchanged. ⛔⛔ **"The emitted row's bytes move under some ablation" is NOT sufficient**
— a tautological residual with a moving ornament in the same payload passes that weaker test (measured this
round on `convention_check_conservative` and `kernel_propagation_residuals`). ⛔ Two routes that share a
constructor (e.g. both call `derive_model()`) are **not** independent.

## Fix 1 — every impedance instance must DESCEND from ONE acoustic-field-computed object

**The object:** **every** impedance any downstream tag consumes — the general impedance, **each regime form**
(propagating / evanescent / grazing), and the **breathing-slice** form.

**The defect (round-1):** the general impedance is a typed literal `Z_GENERAL` (~line 138); the bulk-response
**solve** that computes the same object from the outgoing acoustic field (`z_plus`/`z_minus`, ~612–620,
emitted as `Z_IMPERMEABLE`) is **discarded**; and the regime and slice forms are **independently re-typed**,
so they do not descend from any solve at all: the propagating/evanescent forms at ~714–720
(`PERMEABLE_DISSIPATIVE_BY_REGIME_AND_PARITY`) and ~990–1012 (`BULK_OPERATOR_BY_REGIME`,
`MASS_INTERPRETATION_VALID_WHERE`), and the `k=0` slice `z_value=rho_m*c_s0` at ~1335–1342
(`BREATHING_SLICE_DISPERSION`). Spec §2 (`S11b_SHARED_PHYSICS.md:175-176`) forbids typing `Z` **and its three
regime forms**.

**Required property (⛔ not a value):** there is **one** impedance object, **reached by computation from the
bulk acoustic field** (the construction that produces `z_plus`), with `q_out` symbolic; **every** regime and
slice form is obtained from it **by substitution/limit**, ⛔ never re-typed. The decisive test (one-sided
corruption): a **FORM** change to the bulk-solve construction must **move** every impedance-dependent tag —
`ADDED_MASS`, `Z_BY_REGIME`, `Z_BY_PARITY`, `FACE_RESPONSE`, `PERMEABLE_DISSIPATIVE_BY_REGIME_AND_PARITY`,
`BULK_OPERATOR_BY_REGIME`, `MASS_INTERPRETATION_VALID_WHERE`, `BREATHING_SLICE_DISPERSION`, and the
dispersion/roots. If any stays byte-identical, an impedance is still typed. ⛔ Do not state the impedance's
closed form; compute it. ⚠ Value preservation (above) still binds — the computed object equals the current
one.

## Fix 2 — the export must carry no tautological check

Round-1 exports many ledger rows that are self-checks whose residual is `0`/`−1` **by construction** and
whose **residual is inert under a one-sided corruption of the object they claim to police** — they present as
corroboration to a downstream consumer but verify nothing (`CLAUDE.md` rule 2; corollary 3; the `x === x`
prohibition). ⛔ There is **no** genuine exemplar to copy: `convention_check_conservative` (below) is itself
one of the defective rows — both its operands come from the stored energy `U_LONG`, never from the assembled
equation, and its "movement" under an energy ablation is a moving ornament (`omega_squared`) beside an inert
residual.

**Per-row disposition — apply exactly one:**
- **(G) make it genuine:** build the check's two operands from **independent routes**, at least one of which
  is the **live assembled object** it is specified to police (named per row below), so that a one-sided
  corruption of that object moves the **residual/status**. ⛔ The reference route must not share a constructor
  with the tested route.
- **(L) localise / delete:** only where the row is a **pure ornament** (a hardcoded literal with no object
  behind it) or a genuinely engine-local, **non-mandated** diagnostic — make it `§10`-LOCAL (⛔ not an
  exported ledger row) or delete. ⛔ `§10` `_LOCAL_` is reserved for objects that cannot exist in both engines
  (`S11b_SHARED_PHYSICS.md:938-943`); it is **not** an escape hatch for a required-but-badly-implemented check.
- **mixed row** (real content + a tautological element): keep the content, **remove the fake residual** (or
  make that residual `§10`-local); ⛔ do not re-derive the content.

⛔⛔ **(G) is REQUIRED — (L)/delete is forbidden — for every check §6/§B7 mandate:** the convention checks
(`S11b_SHARED_PHYSICS.md:420`), the three energy discriminators (`:452-489`), kernel extraction "from the
equations actually used", and B7 termwise homogeneity with a failing-corruption demonstration (`:843-850`).

**The worklist (closed; each row is EXPORTED at round-1):**

| exported row | line(s) | tautological element | disposition + reference route |
|---|---|---|---|
| `convention_check_inplane` | ~937–941 | fresh `σ−μ`, `∂/∂μ=−1`; never reads `MODEL["inplane"]` | **(G)** vs live `MODEL["inplane"]` (emitted at ~923). §6(a) mandated |
| `convention_check_conservative` | ~943–948 | both operands from `U_LONG`; residual is a chain-rule identity; `omega_squared` is a moving ornament | **(G)** equation-side operand from the live assembled `MODEL["thickness"]` (~503) vs the reduced stored energy. §6(b) mandated |
| `pressure_work_sign_check` | ~885–912,~967 | `∓Re(Z)|V|²` from the **same** pair ⇒ 0 | **(G)** L/R from the live EOM pairing vs the bulk flux. §6 mandated |
| `full_two_port_balance_check` | ~862–913,~968 | placeholder products regrouped | **(G)** from the assembled model. §6 mandated |
| `energy_sinks`, `energy_sources` | ~851–968 | `±½Re(...)` from free symbols, never see `MODEL` | **(G)** `d/dt(T+U)` from the live equations. §6 mandated |
| `unattributed_sink_terms`, `unattributed_exchange_terms` | ~965–966 | typed empty `Tuple()` | **(G)** derive the residue from the live accounting (an honestly-computed empty set is fine; a typed empty is not). §6 mandated |
| `kernel_orientation_identities` | ~808–826 | `∂(Λ_I·coord)/∂coord − Λ_I ≡ 0` on a typed copy | **(G)** extract `K_{A,V,X}` from the equations **used** (the closure in `face_solution`); one-sided corrupt that closure. §6 mandated |
| `kernel_propagation_residuals` | ~828–843 | both routes call the same `derive_model()` (`MODEL = derive_model()` at ~2035) ⇒ residual inert to a form change inside it | **(G)** make the reference route genuinely independent of the tested constructor, or keep only the moving content and drop the residual-as-check |
| `conservative_positivity_inequality` | ~948–955 | `positive_energy` is a **typed** coefficient inequality | **(G)** compute the positivity condition from the live reduced-energy Hessian; **or** keep the computed `conservative_inequality` and **(L)** the authored implication |
| `two_port_power_identity` | ~1137–1148 | both sides from the same fresh `V,mu_s,p,J` ⇒ residual 0 | **mixed:** keep the **identity object** (B2d wants the identity, `:704-711`); **(L)** the residual-as-corroboration. ⛔ do NOT re-derive it from the EOM (that duplicates the §6 discriminator) |
| `parity_even_jw`, `parity_odd_jw` | ~581–590 | hardcoded even/odd algebra on placeholders; never use the live projection source (~576–577) | **(G)** evaluate against the **live projection source** (⛔ not `MODEL`). B0a (`:597-603`) mandated |
| `homogeneity_*` family + `homogeneity_ablation_demo` | ~1453–1472 | stamped dimension tuples compared to themselves; the demo corrupts the **stamp**, not a live term | **(G)** derive each term's dimension from the actual equation term so the check fails on a corrupted **live** term, and make the ablation demo corrupt a **live** term. ⭐ keep the `DIM_<name>` **values** (~1450). B7 (`:843-850`) mandated |
| `dim_route_kind_*` | ~1451 | `route_kind`/`route_operands` typed in the `dimensions` dict | decide explicitly: compute the route labels, or if they are irreducibly authored metadata state that and keep — ⛔ do not leave the disposition ambiguous |
| `parity_interval` (last field) | ~598 | unconditional literal `true` | **(L)** the ornament field (keep the real interval content) |
| `sheet_of_each_root` (one conjunct) | ~1306 | `Eq(−q_out,−q_out)` | **(L)** the tautological conjunct (keep the real roots) |
| `branch_degenerate_point` (rank) | ~671 | `Matrix([[1,0],[1,0]]).rank()`; not from the bulk ODE | **(G)** the degenerate-locus rank from the actual solve, **or (L)** the typed rank |

⚠ `causality_check` (~935) copies `causality["orientation"]` from the `kernel_orientation`/propagation block —
fixing those carries it; it is not a separate row.

**Two hard constraints:**
- ⛔⛔ **Do NOT add new checks, and ⛔ do not add lines that carry no physics** to make a check "pass"
  (measured S9: a hardening pass added 183 no-physics lines, 19/23 still unfailable). The repair **removes**
  tautologies — it does not proliferate.
- ⛔ **After the repair, no EXPORTED check-row's residual/status may be inert under a one-sided corruption of
  the object it polices.** The orchestrator enforces this by one-sided corruption, ⛔ not by a self-report and
  ⛔ not by mere row-byte movement.

## Out of scope — ⛔ do not touch

- The **energy-basis construction and its count** — a **withheld physics result** (`S11b_SHARED_PHYSICS.md:272-278`);
  the round-1 legs disagreed on it, which the comparator will **report, not resolve**. ⛔ Do not change the
  energy-basis construction and ⛔ do not act on the disagreement here.
- The **longitudinal dispersion / mode fate** (B5) — emitted honestly as an unevaluated determinant + symbolic
  roots; resolved downstream (comparator + step record), ⛔ not here.
- Every physics-object **value**, and the whole **export wiring** — unchanged (above).

## Conflicts

Use the spec's `§10` (and the `§13` report) for anything unclassifiable, unavailable, or that you cannot emit
under one-tag-per-object. ⛔ Do not fill a gap with new physics, an expected result, or a self-review
mechanism. If making a mandated check genuine would require a premise the model does not carry, ⛔ do **not**
localise it — report the obstruction in the `§13` report so the orchestrator can decide; an honest
"cannot police this in-engine" is a finding, ⛔ never a `§10` hide.
