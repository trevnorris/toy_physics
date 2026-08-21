# S11b SymPy engine — fix round 2 (make the tail checks genuine)

## Authority and boundary

Repair `research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py` (baseline `6d57b27e`,
the round-1 repair result) in place; products stay the same two — the flushed stdout tag stream and
`research/pde_ledger_v3/scripts/S11b_exports.py`.

⛔⛔ **The physics is correct and ⛔ must NOT be re-derived or re-valued.** `directives/S11b_SHARED_PHYSICS.md`
(`1a2395a3`) is the physics authority; `CLAUDE.md` binds. Add no expected value (rule 5). This directive
changes **only the four check-genuineness fixes** named below.

**The three script clauses — non-negotiable, carried verbatim:**
> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit`/`Print` payload must
> be a CAS object — an expression, a solved root, a boolean from a symbolic test. ⛔ Never prose describing a
> result.
> **2. PRINT the residual; do NOT assert it.** `assert residual == 0` **is the builder writing down the
> expected output**, and it turns an informative value into a binary crash. Compute → emit → *then* assert.
> **3. Interpretation belongs to the STEP RECORD.** ⛔ The script does not editorialise.

> **The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE ACTION and the
> ANSATZ. Every other expression involving them must be REACHED BY COMPUTATION. Every control re-enters the
> chain at the ACTION, ⛔ never at a result.**

**Preservation boundary (value/wiring).** ⛔ Do not change any **non-target** physics-object value, the
impedance construction (already computed), the primary §6 checks (already genuine —
`convention_check_inplane`/`_conservative`, `pressure_work_sign_check`, `kernel_propagation`, and the
membership/wiring), or the export wiring (S11-LEDGER import + carry-forward, bindings, `F9`,
`BUILD_INPUT_DIGESTS` topology/logic, `D3`, `_RELATIONALS`, freeze). ⭐ The **four target diagnostics** below
**must** change (that is the repair): the parity payloads, the `onsager_reciprocity` cross-check, the
`control_no_reciprocal_traction` residuals, and `kernel_orientation_identities` — together with their copied
containers (`causality_check`, growth/decay diagnostics) and the regenerated audit-file `sha256` digest. The
orchestrator re-verifies **non-target** values + carry-forward + wiring against `6d57b27e`; a change **there**
is a regression, a change in the four targets + their copies + the digest is expected.

## The decisive test — ONE-SIDED CORRUPTION of the RESIDUAL, ⛔ not "the row moves"

For each check below, corrupt **only** the live object it claims to police, re-run, and require the check's
**residual/status** — the tautological element itself — to move. ⛔ "The row's bytes moved" is not
sufficient (a moving ornament passes it). ⛔ The two routes of a check must not share a constructor.

## The four fixes

**1 · `parity_even_jw` AND `parity_odd_jw` (both exported) — rebuild BOTH branches from the live source.**
Round-2 baseline: both branches reduce the live projection to its **outer** integral coefficient and then use
typed placeholder algebra — `even_pair = integral_coefficient*(Omp*je + (-Omp)*je)` and
`odd_pair = integral_coefficient*(Omp*jo + (-Omp)*(-jo))` (~669–670); neither evaluates the live integrand, so
a corruption **inside** the integral body is invisible to both. B0a (`S11b_SHARED_PHYSICS.md:595–603`) requires
evaluating the **live source** `source_finite` (~654, its `Integral` atom and integrand) for both even and odd
`j^w` on the stated interval. Rebuild **both** even and odd from `source_finite` under B0a's parity/interval
assumptions.
⭐ **Value-free acceptance (⛔ state no outcome):** a one-sided corruption of the **even** channel of the live
source must move the even residual, and a one-sided corruption of the **odd** channel must move the odd
residual. ⚠ The genuine object is the live source under B0a's assumptions — ⛔ do NOT engineer an artificially
sensitive residual (B0a asks for the parity evaluation on a stated interval, not integrand-oddness
sensitivity on a symmetric one).

**2 · `control_no_reciprocal_traction` (⚠ `export=False`, a B8 control — still must fail) — genuinely test the
dropped-channel trap, per residual.**
Round-2 baseline: `control_f_model = derive_model(Λ_X0:0)` vs `MODEL["thickness"].subs(Λ_X0,0)` (~1919–1942) —
**both** routes drop `Λ_X`, so `operator_residual` is 0 by construction; the power residual is a typed
identity on fresh control symbols. B8 F names the **dropped-channel trap**. Rebuild with genuinely independent
routes and give **each residual** its own acceptance:
- the cut recomputation vs the **full B2a operator** after symbolic `Λ_X⁰ = 0` — a corruption confined to the
  cut recomputation must move this residual;
- the **full-minus-cut traction effect** checked against the independent placeholder constructor already in
  the engine (~901–938, which still carries `ℓ_X`), ⛔ not a second `derive_model` call — deletion/corruption
  of the live full traction channel must move this residual;
- the power comparison against the **actual emitted B2d power identity** (~1332–1343), ⛔ not fresh
  control-only symbols — a corruption confined to the recomputed power route must move this residual.

**3 · `onsager_reciprocity` cross-check (exported) — two conversions derived SEPARATELY, then compared.**
Round-2 baseline: `reciprocity_crosscheck = flux_res − force_res` with `force_res ≡ −flux_res` by construction
(~1370–1385), so the cross-check is structurally `2·flux_res` and carries no route-agreement information. The
spec **mandates** it (`S11b_SHARED_PHYSICS.md:733–742`: the all-flux and all-force conversions "must give the
**same** relation — report a disagreement rather than choosing"). ⛔ Do not drop it. Derive the two
conversions **independently** from the **live** mixed law (~1363) — the all-flux form (solve the second row,
formal `d = ε ≠ 0`, clear powers of `ε`, take `ε → 0`) and the all-force form (solve the first row where
`a ≠ 0`, clear `a`, extend to `a = 0`) — canonically normalise **each** relation up to a nonzero scalar/sign,
then compare. ⛔ Do **not** merely normalise the two typed dual matrices (they are algebraically dual by
construction). Acceptance: a one-sided corruption of **one** conversion route moves the cross-check.

**4 · `kernel_orientation_identities` (exported) — compare against a SEPARATELY constructed retarded
reference.**
Round-2 baseline: `extracted` is spec-correct coefficient extraction, but the reference `supplied` is the
**module-global** `Λ_A/V/X` (~942–957, defined ~154–156), so a global retarded→advanced transpose moves both
sides and the residual stays 0. Build the reference retarded kernel **separately** from `Λ_I⁰, τ_I, ω` (the
supplied retarded form `Λ_I(ω) = Λ_I⁰/(1 − iωτ_I)`, `S11b_SHARED_PHYSICS.md:216`) and compare it to the `K_I`
coefficients extracted from the **live** equations (spec-decisive, `:381–406`). Acceptance: a retarded→advanced
transpose of the kernel in the live equations moves the residual whenever `τ_I > 0`; the spec's
**indistinguishable** cases (`τ_I = 0`, `Λ_I⁰ = 0`, `:404–406`) must be reported as indistinguishable, ⛔ not
as checked passes (the engine's `indistinguishable` at ~956 already does this — keep it).
⭐ `causality_check` (~1123–1127) and the growth/decay diagnostics (~1518–1522) carry the **same** payload —
they repair automatically with this fix; ⛔ they are not separate work.

## ⛔ Leave unchanged (out of scope)

- `homogeneity_{thickness,mass,affinity}` and the `HOMOGENEITY_*` family — **already able to fail**:
  `trace_dimension` (~1598–1624) recurses into every nested `Add`/`Mul`/`Pow` (the collapsed fraction's
  numerator/denominator terms), and the spec asks for **one** demonstrated corruption (inplane, supplied
  ~1783–1789). ⛔ Do not rebuild these from the uneliminated equations — that would weaken coverage of
  face-elimination errors.
- `energy_sources`/`energy_sinks` and `unattributed_*` — mandated signed/semantic views of the live
  accounting, ⛔ not independent-corroboration tautologies. Leave.
- `two_port_power_identity`, the impedance, the primary §6 checks, every non-target physics value, and the
  whole export wiring — unchanged.

## Hard constraints

- ⛔⛔ **Do NOT add new checks or no-physics lines** (the S9 hardening trap). Make the four existing checks
  genuine; ⛔ do not proliferate.
- ⛔ After the repair, **no targeted check-row — exported or not, including the `export=False` control and
  each of its separate residuals** — may have a residual/status that is inert under a one-sided corruption of
  the object it polices, and no check's two routes may share a constructor. The orchestrator enforces this by
  corruption, ⛔ not by a self-report.

## Conflicts

Use the spec's `§10` (and the `§13` report) for anything unclassifiable or that you cannot emit under
one-tag-per-object. ⛔ Do not fill a gap with new physics, an expected result, or a self-review mechanism. If a
genuine second route for a mandated check is genuinely unavailable, ⛔ do not localise it — report the
obstruction in the `§13` report so the orchestrator can decide.
