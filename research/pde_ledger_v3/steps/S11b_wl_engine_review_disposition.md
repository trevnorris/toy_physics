# S11b blind Wolfram engine — two-leg review, adjudicated (orchestrator fold, rule 7)

Engine: `mathematica/S11b_interface_coupling_law_mathematica_audit.wl` (Codex-built, 178 tags, blind —
no `Import/Get/<<`, no VERDICT). Legs (byte-identical prompt `_legs/S11b_wl_engine_review_prompt.md`,
serialized, FORM ablation mandatory): **Grok** `ba3ry3rnm`; **fresh Agent** `a869b2dcd4e01c9c2`.
Each leg wrote + ran its own derivation before opening the artifact; both saved scripts+stdout under `/tmp`.

⭐ **Physics VALUES are correct and genuinely computed** (both legs, moved under FORM ablation, matched
independent derivations, matched the SymPy engine): impedance / regimes, added mass, face response,
`ZPERM_SLICE`, transverse B6 (`{0,0}`, `ω²=μ_R k²/ρ_br`, `Im ω=0`), breathing B4 (`K₀ = B_ρ³−2CW₀+k_W W₀²`,
radiation-damped dispersion). The defects below are **one localized sign bug + a decoupled verification
layer**, NOT wrong physics in the derived responses.

## CONFIRMED findings (each verified by the orchestrator reading the code, rule 13)

### F-WL-1 — `ZPERM_SLICE_MAP` emits the withheld G13 coefficient with the WRONG SIGN  (PHYSICS)
- Engine (L856–867): `rawPressureCoefficient = Coefficient[equationZeroForm[flux_eq], facePressure]`,
  then `Λ_p := rawPressureCoefficient`. `equationZeroForm[a==b] := a−b` (L367), so it extracts the
  coefficient of `facePressure` in the **residual** `faceFlux − RHS`, which is **−(coefficient in the flux
  RHS)**. Result: `Λ_p = +Λ_A/ρ_m`.
- Correct: on the `μ_s=0` slice `J = Λ_A(−δp/ρ_m) + Λ_V V = −(Λ_A/ρ_m)δp + Λ_V V`, and the raw-pressure
  closure `J = Λ_p δp + Λ_V V` (spec §4: "raw-pressure-driven **with coefficient Λ_p⁰**") forces
  **`Λ_p⁰ = −Λ_A⁰/ρ_m`**.
- Independent confirmations of `−Λ_A/ρ_m`: SymPy engine G13 (`sympy_audit.py:1329–1334`, both SymPy legs);
  Grok leg; orchestrator SymPy derivation (`steps/_measurements/s11b_zperm_slicemap_sign.py`,
  `Λ_p + Λ_A/ρ_m = 0`); first principles.
- ⚠ **Leg disagreement (resolved):** the fresh Agent leg reported "AGREE" with the engine's `+Λ_A/ρ_m`,
  but its script (`/tmp/s11b_wl_leg/indep_derivation.py:57`) defined the closure as `J = −Λ_p δp + Λ_V V`
  (a built-in minus), a non-standard convention. Under the spec's stated convention (`Λ_p⁰` is *the*
  coefficient of the raw-pressure closure), `−Λ_A/ρ_m` is correct and the engine is wrong. Grok is right.
- Scope: `ZPERM_SLICE_MAP` is a terminal emitted tag; it does **not** feed downstream (`ZPERM_SLICE` is
  computed from the affinity directly and is correct). The bug is isolated to this tag's extraction.
- ⛔ This is a genuine cross-engine disagreement (SymPy `−` vs WL `+`) resolving in SymPy's favour — a code
  bug in the WL extraction, ⛔ NOT a physics disagreement to preserve. The T7 comparator will also flag it.

### F-WL-2 — the §6 energy-accounting layer is TAUTOLOGICAL / decoupled from the assembled EOMs  (FIDELITY)
- `PRESSURE_WORK_SIGN_CHECK` (L1168–1181): `slabPressurePower = −½Re[X]`, `outgoingBulkPower = +½Re[X]`
  with the **same** `X` → `RESIDUAL = −½Re[X]+½Re[X] ≡ 0`, `TEST_OBJECT → True` for any input. The two
  "operands" are literally ±the same object — operand theatre (build-skill corollary 3).
- `FULL_TWO_PORT_BALANCE_CHECK`, `TWO_PORT_POWER_IDENTITY`, `ENERGY_SINKS`, `ENERGY_SOURCES`,
  `UNATTRIBUTED_*` run off **hand-typed** `energyEquationRules`/`energyKinematicRules` (L1116–1135), NOT
  the assembled `fullSystem`. A traction-sign FORM ablation that moved 29 dispersion/root/response tags
  left every one of these byte-identical (both legs). A sign error that flips `Im ω` (B5's deliverable)
  passes all of them undetected. Informationless.

### F-WL-3 — kernel-orientation / causality / grazing checks police SHADOW laws, not the assembled objects  (FIDELITY)
- `KERNEL_ORIENTATION_IDENTITIES` (L1020–1047), `CAUSALITY_CHECK` (Grok: `Cases` on nested `Association`s
  → `{}`, `And[] = True`, unconditionally true), `GRAZING_MODE_CLASSIFICATION` (Grok: `Limit[q·MATRIX,q→0]`
  = zero matrix, rank 0, cannot distinguish bound / threshold / BIC / resonance). All extract from separate
  hand-written laws (`closureLawForExtraction` L1020, `requiredRealAxisObject` L523), not `solveFace` /
  `qContinued` / `fullSystem`; byte-identical under FORM ablation (both legs). Same decoupling as F-WL-2.

## CROSS-ENGINE DISAGREEMENT for the comparator (not a WL bug on its own)
### X-1 — energy-basis count: SymPy **11** vs WL **10**
- WL emits `ENERGY_BASIS_COUNT = 10`, omissions `{∇θ·∇e_W, |∇θ|², e_W·∇·u, θ·∇·u, (∇·u)²}`. **Both legs
  independently derived 10** (Agent: 2+3+3+2; Grok: 11 before the total-divergence quotient, 10 after,
  since `(∇·u)² ≡ tr((∇u)²)` modulo a total divergence). SymPy emitted **11** (counts `(∇·u)²` separately).
- ⇒ a genuine convention-level disagreement (mod-divergence quotient taken or not). The T7 comparator must
  report it; resolution (which count the ledger adopts) is a step-record decision, ⛔ not silently dropped.

## Non-findings (measured null, kept as controls)
- `HOMOGENEITY_ABLATION_DEMO` can fail (corrupt→False, restore→True). `CONVENTION_CHECK_INPLANE` /
  `CONVENTION_CHECK_CONSERVATIVE` carry real operands (residual 0 with a live `−∇(δU/δθ)` / energy-route =
  equation-route). Transverse coupling `{0,0}` is a **computed** vanishing (moved under FORM-B, 50 tags).
  No `Abort` precedes a physics emit. No duplicate/verdict tags.

## Disposition — [⚠ HISTORICAL / EXECUTED: comparator DONE 17fe32c8, WL repair DONE bd598ae7; SymPy X-1 = current NEXT. Original proposal below.]
1. **Pin the reviewed baseline**: commit `.wl` + captured `.out` before any repair (rule 6).
2. **Comparator FIRST** (honours "the disagreement is the measurement"): author + freeze + commit the T7
   comparator, run it against the *current* engines so it MEASURES F-WL-1 (sign) and X-1 (count) plus any
   DIMENSION-tag mismatches the legs never checked — then fold ALL findings (2 legs + comparator) into ONE
   WL repair directive (rule 7, fold once). ⛔ Do not repair the WL sign before the comparator records it.
3. **Repair scope** (one round): fix F-WL-1 extraction (coefficient of the flux RHS, or a consistently
   signed zero-form); make F-WL-2/F-WL-3 checks genuine (wire to `fullSystem`) OR emit the objects honestly
   with no tautological residual ("no independent second route here"). X-1 count resolution decided in the
   step record, applied to whichever engine the decision goes against.
4. **Codex-consult the disposition + repair directive** before the builder (SymPy round-2 precedent: my
   disposition was wrong there; Codex corrected it).

---

## Codex consult — corrections ADOPTED (each re-verified by the orchestrator, rule 13)
Consult prompt `directives/_legs/S11b_wl_disposition_codex_consult.md`; Codex measurements under
`steps/_measurements/s11b_codex_*`. The consult confirmed F-WL-1 and comparator-first, and corrected the
disposition in three ways + found two misses. Each adopted below was re-derived/re-read by me.

- **F-WL-1 (sign) — CONFIRMED**, and **EXTENDED**: the emitted `ZPERM_SLICE_MAP` payload (`.out:52`) is
  `RAW_PRESSURE_COEFFICIENT -> lambdaA0/(rhoM(1 - I ω τ))` — i.e. **also the wrong LEVEL** (the dynamic
  `Λ_A(ω)/ρ_m`), where B2c asks for the **static** `Λ_p⁰ ↔ Λ_A⁰`. Repair: emit static `Λ_p⁰ = −Λ_A⁰/ρ_m`
  (optionally a separate dynamic relation).
- **F-WL-3 — SPLIT (I over-bundled).** Verified in code: grazing consumes `fullSystem["MATRIX"]` (L1312)
  and kernel-propagation touches `fullSystem` (L1048) — these are NOT shadows. Three distinct defects:
  - **3a shadow extraction** — `KERNEL_ORIENTATION_IDENTITIES` extracts from `closureLawForExtraction`
    (L1020), not the assembled closure ⇒ extract from the actual closure equations.
  - **3b broken aggregator** — `CAUSALITY_CHECK`: `Cases` on nested `Association`s → `{}`, `And@@{} = True`
    unconditionally ⇒ explicit association lookups.
  - **3c degenerate grazing** — `Limit[q·fullSystem["MATRIX"], q→0]` = zero matrix, rank 0 ⇒ stratum-aware
    Laurent / nullspace classification.
  (F-WL-2 stands: `PRESSURE_WORK_SIGN_CHECK` tautological; energy/two-port are detached hand
  transcriptions L1116–1144 / L1184–1228 ⇒ contract the actual EOM.)
- **X-1 — NOT open; the SPEC decides it, and it is a SYMPY defect.** §5 (L286-287) mandates equivalence
  **modulo total divergences**. Verified: `st_squared = (2/3)(∇·u)² + (1/2)curl²` modulo a total divergence
  (from `(∇·u)² − tr((∇u)²) = ∂_i(u_i ∂_j u_j − u_j ∂_j u_i)`), so among SymPy's `{curl²[0], st_squared[1],
  (∇·u)²[2]}` only TWO are independent. **WL's 10 is correct; SymPy's 11 over-counts** (`sympy_audit.py:
  470-510` judges independence over the full 9-component gradient, not modulo divergence; `st_squared`'s
  coefficient is degenerate with `B_div`/`mu_R`). ⛔ This **REOPENS the committed SymPy engine** for X-1
  (basis count + redundant coefficient); the SymPy physics EOM is unaffected (the surplus term is a total
  divergence). ⇒ X-1 is a **SymPy** repair scope, ⛔ NOT part of the WL repair.
- **Q5-2 (missed) — thickness-response DIMENSION difference is normalization, not physics.** WL uses
  displacement/pressure (L763-792); SymPy uses displacement/force-per-x-volume (`sympy_audit.py:1203-1212`).
  ⇒ the T7 comparator must **normalize by `W₀`** on that row, ⛔ not report it as a physics disagreement.

## Path — CONFIRMED (user 2026-08-21): comparator FIRST, then SPLIT repairs
> ⚠ **STATUS 2026-08-22: comparator DONE+FROZEN+RUN (`17fe32c8`, X-1 confirmed); WL repair DONE (`bd598ae7`).
> This path is mostly executed — the CURRENT NEXT is the SymPy X-1 repair. See `STATUS.md` + the memory
> frontmatter for authoritative current state.**
1. [✅ DONE 17fe32c8] Author + FREEZE + commit the T7 comparator (join by emitted name, residual with symbol transliteration,
   ⛔ reject native boolean as a residual operand, three-valued undecided, COMPARE dimension tags with the
   `W₀` normalization above, synthetic pass/fail self-tests + repoint ablation) BEFORE it sees the engines'
   real output (rule 5). Run it against the pinned baseline (`S11b_exports.py` rows ↔ WL `.out` tag stream);
   the live F-WL-1 sign defect and X-1 count are useful failure cases.
2. [WL repair ✅ DONE bd598ae7; SymPy X-1 repair = current NEXT] THEN two SEPARATE 2-legged repair directives: **WL** (F-WL-1 sign+level, F-WL-2, F-WL-3a/b/c) and
   **SymPy** (X-1 basis count + redundant coefficient). Each: directive → 2 legs (rule 7) → build → 2 legs.
