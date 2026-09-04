# S11c-c1 SymPy engine — build directive

## Authority and boundary

Write `research/pde_ledger_v3/scripts/S11c_c1_bulk_closure_sympy_audit.py` in full. Its products are the flushed
stdout tag stream and `research/pde_ledger_v3/scripts/S11c_c1_exports.py`. Those two products are its only
writes; every other file, including `S11b_exports.py`, `S11c_a_exports.py`, and `S11c_b_exports.py`, is
read-only.

`CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md` (first committed `db5cbf88`;
§1/§7 amended by the export-architecture migration to the fold/own-rows-delta topology — its exact content sha
is pinned in `BUILD_INPUT_DIGESTS`) is the sole physics authority and the sole physics input: every equation,
inherited object, ansatz,
premise, supplied law, face map, radiation condition, branch object, derivation route, control and tag
obligation enters from it alone, and it wins every conflict with this directive. Implement its §§0–8 directly.
Point to rather than duplicate: §0 (scope, and what is reserved for S11c-c2/S11c-d/S11c-e); §1 (the inherited
S11b/S11c-a/S11c-b setup — the DOF and two face variables, the bulk acoustics and radiation condition and branch
object `q_out`, the flat B0b/B0c reference used only as the §5c regression operand, and the face closure laws
with the `Λ_X`-in-traction placement); §2 (the two disconnected curved half-spaces and the driven boundary
data, the standing rest-frame limit and its non-uniform validity domain, the multigrade bookkeeping); §3a (the
two-momentum DtN/impedance **operator**) and §3b (the permeable curved face response, its Fredholm/finite-dim
loci, and the three-object dissipation audit); §4 (the objects to compute and emit); §5 (the
representation-invariance / one-sided independence / form-ablation / uniform-limit / zero-jet / branch-momentum
liveness controls); §6 (method, dimensions, the locus protocol, the three script clauses, no VERDICT); §7 (names,
F9, chain output, tag grammar); §8 (supplied vs computed). Add no expected value or acceptance criterion
(`CLAUDE.md` rule 5) — S11c-c1 withholds every computed value (the operator/kernel structure, the parity, the
regime-pair behaviour, the loci, the dissipation signs, the residuals), and this §-wiring adds no physics.

Do not add parallel machinery: no per-cell completion registry beyond §4/§5's task structure, and no
directive-local exit policy. §6's script clauses and the emission-not-conditional-on-value rule apply as stated,
without named exceptions.

## The conventions this directive pins (both engines must anchor them identically — naming/binding, not physics)

- **`μ_θ` is the reserved imported operand (spec §1c/§1d/§3b/§4).** Bind `S11CB_MU_THETA_OPERATOR` to the
  imported row `LEDGER['mu_theta_operator']` (`S11c_b_exports.py`) and carry it as the opaque named operator
  feeding `μ_s = μ_θ / ρ_br,bg⁰` (branched `μ_s^α ≡ μ_θ^α / ρ_br,bg^{0,α}` on MATERIAL_ADVECTED, anchored at `χ`
  as every other `Q_bg`). ⛔ Do **not** re-derive `μ_θ`, ⛔ do **not** expand it into slab DOFs, and ⛔ do **not**
  construct it from the uniform energy — its composition with the slab variables (and the `Σ_E`/`μ_R,bg`
  channels living inside it) is **S11c-c2's** work, reserved (§0, §5a, §5b).
- **The §5c/§5d regression operand B is the UNMODIFIED imported S11b flat object.** Bind the flat B0b/B0c
  reference from the LEDGER — `z_impermeable`, `z_by_regime`, `z_by_parity` (if present), `added_mass`,
  `grazing_behaviour`, `face_response`, `face_response_coeffs`, `permeable_dissipative_by_regime_and_parity` —
  as the §5c/§5d regression targets. ⛔ §5d operand B is the **bare** `z_impermeable` with **no**
  `W_0→W̄₀(1+η)` substitution and **no** two-face re-solve (spec §5d: re-solving a finite-gap cavity reproduces
  the very `O(η)` cavity error the control exists to catch). Where the spec says "imported/re-derived" (§1c), the
  imported row is authoritative for the regression operand; the engine's own zeroth-shape-order flat symbol
  `S11CC1_DTN_FLAT_SYMBOL` is its **independent** re-derivation and the two are compared, ⛔ not silently
  identified.
- **The two-momentum kernel leg convention is fixed once (spec §3a/§5e): in `Z_s(ω;k,k′)` with `Ŵ_bg(k−k′)`,
  `k` is the OUTPUT leg and `k′` the INPUT leg.** Anchor it identically so the §5e sign-flip and momentum-freeze
  controls corrupt the **same** named leg in both engines. The spec runs **both** one-leg freezes to remove the
  residual convention ambiguity — implement both.
- **`v_bulk_normal_0` binds to the imported drain row `LEDGER['v_bulk_normal_0']`; ⛔ it is never aliased to
  `v_0`** (spec §7/§2b). The engine computes only the **rest-frame** bulk object (§1b/§2b, `N11a`); ⛔ it does
  **not** construct the convective operator, and the non-uniform validity domain and "grazing = strict
  `v_bulk_normal_0=0`" (§2b) are recorded by the **step record**, ⛔ not carried as an operator term or emitted
  as a computed value by the engine.

These are naming/binding/anchoring conventions, ⛔ not new physics; every one is already in the spec, and the
spec wins any conflict.

## Delta export — generate over the frozen base (`directives/export_ledger_bind_closure_design.md` §D1–§D3)

⭐ **Topology (§D2): import the fold, write only c1's own-rows DELTA — ⛔ NOT the accumulated model.** Obtain the
inherited model with `from ledger_fold import load_model, check_consumer, assert_lookups_equal_manifest,
assert_delta_is_minimal` then `fold, fold_audit = load_model('research/pde_ledger_v3/scripts/S11c_b_exports.py')`.
`S11c_b_exports.py` is the **single atomic frozen base** (c1 is the first delta-producer, so there are no prior
deltas); it already carries the whole F9-resolved S11b + S11c-a + S11c-b model flat (2441 rows), so every object
c1 binds is present in `fold`. ⛔ Do **not** `import LEDGER from S11c_b_exports.py`; ⛔ do **not** re-apply F9 on
the fold (`load_model` is last-wins on exact key and never re-routes). **Throughout this directive `LEDGER[k]`
denotes `fold[k]`.** Write `S11c_c1_exports.py` as c1's **own-rows delta only** — the rows c1 defines/re-derives
that a later step binds (§D1) — ⛔ never the previous LEDGER carried forward.

⭐ **Declare c1's `IMPORT_KEYS` explicitly and exactly (§D3.1; ⛔ not "closure will find them")** — the exact list
of F9 write-keys c1 binds directly (every key looked up in `fold` as a construction or control/regression
operand), grounded present-and-closing against the real frozen base (both review legs folded;
`directives/_measurements/c1_migration_import_keys_scan.py`: **44 roots → 193-row recursive closure, 0 ambiguity**):

- **μ_θ operand + bulk/closure constants, frequency, and grade bookkeepers:** `mu_theta_operator`, `c_s0`,
  `rho_m`, `rho_br`, `W_0`, `e_W`, `v_bulk_normal_0`, `q_out`, `omega`, `epsilon_shape`, `Lambda_A_0`,
  `Lambda_V_0`, `Lambda_X_0`, `tau_A`, `tau_V`, `tau_X` (⭐ `omega` = the frequency ω in `q_out(ω,k)`/`Λ_I(ω)`/the
  dissipation; `epsilon_shape` = the wave-amplitude bookkeeper `ε` grading every object, the third of `(ε,η,σ_W)`
  with `eta_bg`/`sigma_W` — both are direct binds a faithful build looks up, not "closure will find them");
- **the nine S11c-a T-substrate rows:** `face_normal`, `conormal_deriv`, `face_measure_shape_deriv`,
  `face_velocity`, `relative_flux`, `kinematic_balance`, `traction`, `face_shift`, `closure_shape_deriv`;
- **the S11c-a varying profiles and the two ρ_br,bg⁰ density representatives** (⭐ density enters c1 ONLY through
  `μ_s = μ_θ/ρ_br,bg⁰`; the bulk is constant `ρ_m`, §1b): `W_bg`, `w1_profile`, `L_W`, `sigma_W`, `eta_bg`,
  `rho_br_bg_rho4_constant` (the RHO4-const rep — `rho_br` above is the RHOBR-const rep);
- **the S11b flat B0b/B0c §5c/§5d regression rows (spec §1c):** `z_impermeable`, `z_by_regime`, `z_by_parity`,
  `added_mass`, `grazing_behaviour`, the **bare** `face_response`, `face_response_coeffs`,
  `permeable_dissipative_by_regime_and_parity`, and the five `degenerate_loci_{equations,solution,
  identically_satisfied,inconsistent,real_admissible}` (the design names the loci a c1 regression operand,
  `export_ledger_bind_closure_design.md`:49; §5d uses only `z_impermeable`, the rest are §5c's flat package).

⛔ **NOT bound by c1 — reserved for c2 or structurally absent** (both legs, verified against the real base: none
is referenced by any c1 root's closure): `mu_R_bg` and `m1_profile` (the modulus channel `μ_R,bg`, inside opaque
`μ_θ`, reserved for c2 — spec §5a/§5b/§8); `mu_R` (the bare constant, reachable only through `μ_θ`'s closure —
⛔ never a direct lookup); `rho_4D_bg_rho4_constant`/`rho_4D_bg_rhobr_constant` (slab-side 4D-density reps — the
bulk background is the constant `ρ_m`); `e_W_bg` (the background contrast, carried by `η` = `eta_bg`).

⚠ **The manifest is grounded, but the build's bidirectional smoke-test is the FINAL arbiter.** If the engine's
actual construction/regression lookups differ from this list in either direction — an undeclared lookup
(under-export) or a declared-but-unused key (stale manifest) — `assert_lookups_equal_manifest` raises. Treat any
such mismatch as an **operational** failure: report the exact offending keys under §8 for a one-key manifest fix,
⛔ never resolve it by adding a spurious lookup, by re-originating an inherited row (severs chain identity), or by
silently dropping a genuinely-bound key.

Bind the inherited objects to the imported rows per §1/§2/§3 of the spec, ⛔ never re-declaring them under a new
identity: the constants `c_s0` → `LEDGER['c_s0']`, `ρ_br⁰`/`W̄₀`/`e_W` →
`LEDGER['rho_br']`/`LEDGER['W_0']`/`LEDGER['e_W']`, `ρ_m` → `LEDGER['rho_m']` (an S11b `KNOB`; ⛔ do not
re-originate it — a fresh declaration severs chain identity in every bulk/flux/balance/traction/closure object
that uses it), `v_dr` → `LEDGER['v_bulk_normal_0']`, the branch object `q_out` → `LEDGER['q_out']`, the frequency
`ω` → `LEDGER['omega']`, the amplitude bookkeeper `ε` → `LEDGER['epsilon_shape']`; ⚠ `μ_R` is **not** bound
directly (it lives inside the opaque `μ_θ` and is reserved for c2 — it reaches c1 only through `μ_θ`'s closure,
⛔ never a direct `LEDGER['mu_R']` lookup). Bind the S11c-a shape-derivative substrate the curved face objects
are built on (T-a `face_normal`, T-a′ `conormal_deriv`, T-a″ `face_measure_shape_deriv`, T-b `face_velocity`,
T-c `relative_flux`, T-c′ `kinematic_balance`, T-d `traction`, T-e `face_shift`, T-i `closure_shape_deriv`) and
the S11c-a varying profiles the bulk closure consumes — `W_bg`, `w1_profile`, `L_W`, `sigma_W`, `eta_bg`, and the
two ρ_br,bg⁰ density representatives `rho_br`/`rho_br_bg_rho4_constant` (⛔ **not** `mu_R_bg`/`m1_profile` — the
c2-reserved modulus channel; ⛔ **not** the slab-side `rho_4D_bg_*` 4D-density reps or `e_W_bg` — the bulk is
constant `ρ_m` and the contrast is `η`=`eta_bg`) → their `S11c_a`-originated rows; the S11c-b `mu_theta_operator`
→ its row; the S11b flat B0b/B0c reference rows named above → their rows. ⛔ Never re-originate any of these. Only S11c-c1's **own new objects** — the curved
DtN operator, its flat symbol, the two-momentum kernel, the permeable curved face response, the formal
noninvertibility condition, and any new permeability/memory-kernel constant §3 emits — originate in this step
and have no upstream row; ⛔ do not manufacture one, and ⛔ do not let any of them reuse a constant key (§7's
reservations bind: the new operators/kernels/constants take fresh names, never `W_0`/`mu_R`/`rho_br`/`e_W`/`v_0`,
and `v_bulk_normal_0` is never aliased to `v_0`).

**⭐ EMIT vs EXPORT are SEPARATE channels** (the survivor of the superseded `F10`, retained in
`directives/export_ledger_bind_closure_design.md` §D1: the LEDGER carries the bind-closure model-level rows;
step-level diagnostics are emit-only) **and spec §7** (the §3b consume-set). The two channels:

- **Emit (stdout tag stream → the annexed `out/*.out`):** the engine emits **every** §4 primary object **and**
  every §5/§6 control operand and residual, unchanged and complete. That stream is what the cross-engine
  comparator and the review legs consume, and it is where the large diagnostic objects live — the `.out` is
  git-annex+GIN backed, so it carries **no** size cap. ⛔ Nothing is dropped from the emit (rule 11:
  correctness is king; the dissipation audit, loci, and structural views are all still computed, printed,
  cross-engine-compared, and reviewed).
- **Export (the plain-git `S11c_c1_exports.py` delta → folded by S11c-c2):** c1's **own-rows delta only** — the
  rows c1 defines/re-derives that a later step **binds** (§D1), plus the **new symbols** those rows introduce.
  ⛔ **Not** merged onto a carried-forward LEDGER: the delta holds only c1's own rows; c2 obtains everything else
  from the frozen base through `load_model(base ⊔ deltas)`.

**The EXPORTED delta population is EXACTLY** (§3b: *"the response objects `(δp_s, J_s, t_s)`, the DtN operator,
and the flat symbol are the c1 exports S11c-c2 consumes"* — that sentence is the authority; nothing else):

- the curved-face **DtN operator** and its **flat symbol** (§3a): `dtn_operator` (composition), `dtn_flat_symbol`,
  `dtn_kernel` (two-momentum `k,k′`) — fresh names, F9-written **bare** (no collision in the frozen base);
- the **permeable curved face response** (§3b): its object name `face_response`/`face_response_coeffs`
  **collides** with the inherited S11b flat rows, so F9c writes c1's under the **producer-prefixed** keys
  `s11c_c1_face_response`, `s11c_c1_face_response_coeffs` (the bare predecessors stay the S11b rows c1 imported as
  §5c/§5d regression operands — ⛔ never overwrite them). Report the F9c write under §8;
- the **new symbols** those five rows introduce (any new permeability/memory-kernel constant §3 emits, the new
  `k`/`k′` legs, `V_±`, …) — ⛔ **not** the inherited symbols they reference (`q_out`, the `Lambda_*`/`tau_*`
  knobs, the profiles, the density representatives, …): those already live in the frozen base and are resolved by
  the fold, ⛔ not re-emitted into the delta.

⛔ **NOT exported — emit-to-stdout ONLY, because these are step-level DIAGNOSTICS/EVIDENCE, not model-level rows a
later step binds** (per §D1's recovery clause, if a later step turns out to bind one it is promoted then via an
explicit **producer promotion-delta** — `promotion_delta` — the object being still present in the annexed
`.out`; ⛔ a bare `IMPORT_KEYS` addition does **not** materialize it): the dissipation
audit (`permeable_dissipation_vs_omega_tau`, `permeable_port_hermitian`, the `energy_*` operands and residual),
the DtN Hermitian/reactive split (`dtn_hermitian_part`, `dtn_reactive_part`, `dtn_inertial_loading`), the DtN
structural views and traces (`dtn_by_regime_pair`, `dtn_by_parity`, `dtn_grazing_behaviour`, `dtn_term_origins`,
`dtn_rigid_shift_*`), the loci and noninvertibility condition (`noninvertibility_condition`,
c1's own `degenerate_loci_*`), and **every §5 control (5a–5e) and §6 homogeneity operand** (the S11b-`B8`
analog). §7's engine-local tags are engine-local, ⛔ not exports. ⛔ **A missing manifest is a HARD ERROR, ⛔
never a fallback to "export every primary"** (§D3.1 — that fallback is what put 17 MiB of `term_origins` in plain
git). If an object's chain role is genuinely unclear, default it to emit-only and report it under §8 — but ⛔ do
not drop any object §3b names as consumed. A row's `class` follows from what its object **is** — a declared
symbol carries its declared class, ⛔ never a table lookup by row name. Use the class vocabulary at
`research/pde_ledger_v3/S9_REWRITE_PLAN.md:41` and `research/pde_ledger_v3/REBUILD_HANDOFF.md:40`; report anything
unclassifiable under §7/§8 rather than guess.

**⭐ Validate the delta before publication (§D3, from `ledger_fold.py`) — a required real-artifact acceptance
step, not a synthetic check:**

- `check_consumer(fold, IMPORT_KEYS)` — the manifest's recursive `srepr`-identity + `dimension_key` closure
  resolves with no `AmbiguousSymbolError`/`ClosureError` (a missing declared key raises `ManifestError`);
- `assert_lookups_equal_manifest(build_fn, fold, IMPORT_KEYS)` — the **bidirectional smoke-test**: route **every**
  construction/regression `fold` lookup through the access-recording proxy this passes to `build_fn`, so the set
  of observed `__getitem__` keys **equals `IMPORT_KEYS` exactly** (an undeclared lookup = under-export risk; a
  declared-but-unused key = stale manifest — both raise);
- `assert_delta_is_minimal(delta_ledger, own_bind_closure, infra_keys)` — the written delta's key set equals c1's
  own bind-closure (the five model-level rows above + their new-symbol closure) plus any explicitly named
  infrastructure; an accidental re-accumulation fails loudly.

Emit `check_consumer`'s `closure`/`symbol_edges`/`dimension_edges` counts and each assertion's outcome as
computed objects in the stdout stream (rule-2 operands; ⛔ not a bare "OK"). A guard failure is an **operational**
failure (nonzero exit), distinct from a physics disagreement (§6: emit-and-continue, exit 0).

`research/pde_ledger_v3/directives/S11_export_chain_decisions_v2.md` **F1** (keys stay FLAT; `D5` unchanged) is
the storage rule, and **F9** alone decides which key form a collision writes. **Apply both whole; ⛔ do not
restate them, specify a tag-to-key transform, or substitute an earlier export proposal** — §7 already directs
the engine to apply the inherited F9 collision check. ⭐ **F9's collision rule is three-valued and its
comparison is TOTAL over the imported row shapes** (elementwise into containers; a residual that raises on a
tuple/boolean/string/relation makes `F9b` silently never fire, `S11_export_chain_decisions_v2.md:171-186`):
apply it whole. Emit as computed objects each candidate row's engine-derived key and the imported `LEDGER` keys
having those spellings (rule-2 operands, ⛔ not a claim or count; add no guard). Apply **F3** (a re-derived row
carries its own evidence, in the row). Two S11c-c1-specific points pointing cannot supply: **(1)** on an `F9c`
write the fresh step prefix is **`s11c_c1_`** — `s11_`, `s11b_`, `s11c_a_`, `s11c_b_` are already taken in the
imported `LEDGER`; verify non-collision against the imported keys, and report any `F9c` write in the §8 report;
**(2)** S11c-c1 **has** a §6 locus protocol (the flat/finite-dimensional `S11CC1_DEGENERATE_LOCI_*`), but that
protocol governs §3b **physics** loci and its typed CAS booleans are objects consumed by the downstream T7
comparator (which rejects native booleans) — it is **not** an F9 override; F9's proved-different-vs-undecided
row-collision reporting still uses the **inherited `T5` default**, with no spec-specific override.

`BUILD_INPUT_DIGESTS` is a `MappingProxyType` sha256-pinning **exactly the six inputs the c1 spec §7 requires**
— the full ancestry (differing from the earlier sub-steps, which pinned only the immediate parent) **plus the
shared fold module** (a 6th executable input, `export_ledger_bind_closure_design.md` §D3 / build-plan step 1):
`{ 'S11c_c1_bulk_closure_sympy_audit.py' (this engine's own source), 'S11b_exports.py', 'S11c_a_exports.py',
'S11c_b_exports.py', 'S11c_c1_SHARED_PHYSICS.md', 'ledger_fold.py' }`. All six are read to compute their sha256
even though `load_model` reads only `S11c_b_exports.py` for data — the `S11b`/`S11c_a` pins are the byte-provenance
of the frozen base's ancestry (each base file already pins its parent's content sha), retained unchanged. Freeze
the export as S11b/S11c-a/S11c-b do — outer `MappingProxyType`, per-record `MappingProxyType`, then `del _LEDGER`.

Choose **F6**'s first branch: publish `S11c_c1_exports.py` only if every **primary** task of §4 (the §3a DtN
operator objects and the §3b response objects) completed. Use §4/§5's observed run record; add no completeness
field.

## Three deviations from copying S11 verbatim (inherited from S11b `G8`)

- **(a) The cross-engine comparator is a SEPARATE downstream artifact — ⛔ NOT built here.** It is authored
  after both S11c-c1 engines build, to the frozen **`T7`** contract
  (`S11_C17_C18_spec_repair_decisions_v2.md:53-60`) and S11b's **`G8(a)`** (`S11b_unified_decisions.md:82-87`) —
  ⭐ applied **whole**, ⛔ **not restated here** (a paraphrase drops `G8(a)`'s **repoint ablation** and weakens
  `T7`'s undecided from a `T5` coverage-bounding finding to an ordinary third state — the same defect class as an
  F9 paraphrase). ⚠ c1's §6 locus protocol emits typed CAS booleans (the `STATUS_TOKEN` fields), so the T7
  contract's "reject a native boolean as a residual operand" is on the critical path for this join. This engine
  emits no comparator and no cross-engine residual.
- **(b) Restore the `D3` in-run exec-and-residual round-trip** of the finished export module
  (`S9_REWRITE_PLAN.md:217`: reconstruct what was written and compare against the live object, with a
  human-readable rendering alongside).
- **(c) Include the `_RELATIONALS` reviver block** in `S11c_c1_exports.py` (S11b/S11c-a/S11c-b carry it): if any
  emitted row's `value` is a `sympy` relational, `_restore` must revive it. ⚠ c1's noninvertibility condition
  and some loci are relational, so the block is likely used here; ⛔ do not drop it, and ⛔ do not manufacture a
  relational to justify it.

## Mechanical precedent

S11c-b's exporter (`scripts/S11c_b_brane_operator_sympy_audit.py`'s export-emission and
`scripts/S11c_b_exports.py`) is **mechanical precedent, not authority** — the row schema (`display`, `value` as
`_restore("<srepr>")`, `value_kind`, `class`, `step`, optional `dimension_key`, `f9_operands`,
`corroborated_steps`), the `_restore`/`_RELATIONALS` reviver, the self-referential digest, and the freeze footer
are the working shape. ⛔ Non-verbatim: any S10/S11/S11b/S11c zero-collision-residual assertion is not reused —
`F9b` carries **no** residual (the operand pair is the object). ⛔ Do not import or transcribe any sibling
computed physics beyond the chain bindings named above — S11c-c1's own §3–§5 objects are re-derived from its
spec.

## Memory / stdout hygiene (⛔ NOT a change to what is computed — every object and case is still emitted)

⚠ This operator/kernel build is memory-heavy (the S11c-b build OOM'd a 30 GB box). Flush stdout per emitted tag,
and prefer per-case streaming (per anchoring / regime-pair / parity / direction / mutation) over holding the
fully expanded two-momentum kernel across all `(k,k′)` regime pairs in memory simultaneously. This is stdout and
memory hygiene only — ⛔ it drops no object, narrows no case, and never conditions emission on a payload's value
(rule 11: correctness is king; every §4/§5 case is still computed and emitted).

## Conflicts

`F9` is not an acceptance harness; ⛔ do not encode any review checklist as guards, registries or invariants in
the engine, and ⛔ do not make a check audit the input from which it was built. Use the spec's §7/§8 (and the §8
report) for any conflict, ambiguity, unavailable construction, or object you cannot emit under
one-tag-per-object. ⛔ Do not fill such a gap with new physics, an expected result, or a self-review mechanism.
⛔ Do not run the full task set as a completeness loop that hides a per-task failure — §6's exit rule
(emit-and-continue on a physics disagreement; nonzero only on operational failure) is exact. ⛔ A growing or
decaying root, a real or imaginary impedance/kernel part, is emitted as the object the computation returns —
never suppressed, re-branched, or reclassified to reach a passive or lossless form (spec §1b/§6).
