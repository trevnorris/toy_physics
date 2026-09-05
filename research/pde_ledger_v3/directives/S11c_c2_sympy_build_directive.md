# S11c-c2 — SymPy build directive (the self-energy fold; SymPy engine)

⭐ **THIN directive.** All physics is the cleared spec `directives/S11c_c2_SHARED_PHYSICS.md` (v2, committed
`16849fc6`; 2 decision legs + fold). This directive ⛔ does **not** restate or re-derive the physics (a re-wording is
weaker and drifts) — it POINTS at the spec and fixes only the **build-mechanical** layer: the import wiring, the fold
symbol map, the emit tags, and the three script clauses. Model: **CODE build → `gpt-6-astra` high**.

- **Deliverable script:** `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` (`PY_S11CC2_*`).
- **Deliverable export:** `scripts/S11c_c2_exports.py` (own-rows delta; §5).
- **Physics authority (SUPPLIED, unfalsifiable in this build):** the whole c2 spec `16849fc6`. This directive governs
  the **SymPy** engine only; the blind Wolfram engine + the T7 comparator are separate downstream artifacts.

---

## 0 · The three mandatory script clauses (verbatim — `.claude/skills/build/SKILL.md`, non-negotiable)
> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit`/`Print` payload
> must be a CAS object — an expression, a solved root, a boolean from a symbolic test. ⛔ Never prose
> describing a result.
> **2. PRINT the residual; do NOT assert it.** `assert residual == 0` **is the builder writing down the
> expected output**, and it turns an informative value into a binary crash. Compute → emit → *then* assert.
> **3. Interpretation belongs to the STEP RECORD.** ⛔ The script does not editorialise.

**Structural rule (verbatim):** *The ONLY place the physical symbols may be combined by hand is in CONSTRUCTING THE
ACTION and the ANSATZ. Every other expression involving them must be REACHED BY COMPUTATION. Every control re-enters
the chain at the ACTION / the imported operands, ⛔ never at a result.*

**Corollaries (all bind):** a hand-typed CAS object is still hand-typed (delete `Solve`/`subs` and the emit must
move); a tag **name** names the object, ⛔ never its value/sign/order/parity/grade; **no tautological residual** — the
§3c increment is an **EXPORT REPRESENTATION**, ⛔ not a check (emit both operands and the increment; ⛔ do not dress a
structural zero as a check); the §3b adjointness residual and the §3d.4 pairing residual are emitted **only if a
genuine independent second route exists** (⛔ else emit the objects and say there is none); emission is **never**
conditional on a payload's VALUE (only on which package/quantity). ⛔ Every anti-example in this directive uses
**placeholder** symbols. ⛔ Run `reduction/derived_or_declared.py` on the deliverable after the build (triage, not a
verdict).

---

## 1 · What to build — POINTERS to the cleared spec (⛔ do not restate the physics)
Build the c2 self-energy fold exactly per `S11c_c2_SHARED_PHYSICS.md`:
- **§0** scope; **§1** inherited setup (SUPPLIED); **§1b** the c1-import AGREE/UNDECIDED disposition — ⛔ do **not**
  treat any UNDECIDED item (density freeze, `t_s`, DtN whole-form, flat-resolvent leg-labeling, ENERGY) as
  cross-engine-closed; **§1c** the real θ-row + the fold operation; **§1d** the Λ-channel routing.
- **§2** the non-commuting ordering — build **`extract(close(SLAB))`** (close FIRST, then re-extract); the three
  distinct named objects (`S11CC2_CLOSED_COUPLING_KERNEL`, `S11CC2_ORDERING_COMMUTATOR`,
  `S11CC2_SELF_ENERGY_INCREMENT`).
- **§3a** close; **§3b** re-extract (adjointness residual only if independent); **§3c** the increment as an EXPORT
  REPRESENTATION that **SURFACES** (⛔ not cancels) S11c-b's two cross-engine-UNVALIDATED sign conventions
  (face-generalized-force PY `+diff` / WL `−linearVirtualVariation`; #90 closure-fold); **§3d** the SIX
  re-adjudications.
- **§5** the controls (5a ordering ablation; 5b routing ablation; 5c the **N6 two-anchoring** rep-invariance route;
  5d density field-vs-field; 5e the **three DISTINCT** limits `Z→0` / `Λ_A⁰=0` / impermeable + the `μ_R,bg` **FORM**
  ablation).
- **§6** method/dimensions/script obligations; **§7** export schema; **§8** supplied-vs-computed + builder report.

⛔ **HELD PHYSICS the fold MUST carry (spec; the v1 decision gate caught each of these — do not regress):**
close-then-extract ordering; **substitute the closed `δp_s`(V_s,μ_θ)+w-jets, ⛔ NOT a closed `J_s`** (there is no
`J_s` slot — the θ-row flux is already `Λ_A𝒜_s+Λ_V V_s`; adding a `J_s` double-counts); **`Λ_X` in the traction
`t_s` ONLY** (⛔ not in the θ-row flux, ⛔ not in `closure_shape_deriv`); the increment **SURFACES** the two sign
conventions (⛔ does not cancel them) and the **§3d.4 traction-slab mechanical-power pairing ADJUDICATES the
face-force sign**; the three mandatory rule-17/UNDECIDED carry-ins — **§3d.1** density **field-vs-field** (bind
`rho_br_bg_rho4_constant` to `background_density_map` **before** the fold, ⛔ never a bare constant; `∇ρ→0` is **not**
an accepted corruption), **§3d.2** `t_s` in its **native covector form** (⛔ not a pre-collapsed scalar), **§3d.3**
DtN **kernel-vs-whole-form** (the load-bearing bulk object is the AGREE'd two-momentum `dtn_kernel`; emit whole-form
dependence separately so a leg can ablate it).

---

## 2 · Import wiring — VERIFIED against the real two-parent fold (build-mechanical)

Read the inherited model via the **positional** call (spec §7; signature `load_model(base_path, *delta_paths)`,
`scripts/ledger_fold.py:102` — ⛔ **NOT** `load_model(base=…, deltas=[…])`, which TypeErrors):

```python
fold, audit = load_model("scripts/S11c_b_exports.py", "scripts/S11c_c1_exports.py")
```

Verified on the real files (reproduce; do not hard-code): base **2441** rows + c1 delta **44** rows → fold **2485**,
exact-key intersection **empty**, `overwrites == []`, `check_consumer` closure resolves (355 keys, no ambiguity).

**`IMPORT_KEYS` — the rule (⛔ not a fixed number to copy blindly).** The fold guard
`assert_lookups_equal_manifest` requires `IMPORT_KEYS` to equal **exactly the set of keys the build looks up by
`fold[key]`** — an **undeclared** lookup **and** a **declared-but-unused** key each fail. Therefore fix `IMPORT_KEYS`
= the build's actual direct-lookup set, subject to:
- **MUST include** these **18 closure-covering object roots** (verified: their `check_consumer` closure covers the
  entire §7 consume-set — all 38 keys — so every constant, resolvent, momentum leg, and `w1`-FT object is
  **closure-pulled** and need not be an explicit root):
  ```
  S11c-b (13): slab_operator, slab_operator_term_origins, mu_theta_operator, closure_shape_deriv,
               face_normal, conormal_deriv, face_measure_shape_deriv, face_velocity, relative_flux,
               kinematic_balance, traction, face_shift, background_density_map
  S11c-c1 (5): s11c_c1_face_response, s11c_c1_face_response_coeffs, dtn_operator, dtn_flat_symbol, dtn_kernel
  ```
- **MAY additionally declare** any constant / sub-object row the build reads **directly by key** (for its symbol or
  `dimension_key`) — e.g. the constants `rho_m, rho_br, W_0, Lambda_A_0/V_0/X_0, tau_A/V/X, c_s0, omega,
  epsilon_shape, W_bg, w1_profile, L_W, sigma_W, eta_bg, rho_br_bg_rho4_constant`, the 8
  `s11cc1_response_resolvent_*`, the momentum legs `s11cc1_q_out_{input,output}` / `s11cc1_k_{input,output}_{1,2,3}`,
  `s11cc1_w1_profile_hat_transfer`, `s11cc1_w1_profile_jet_hat_{1,2,3}`, `s11cc1_dtn_operator_{α}_{s}` — **iff** the
  build looks them up by key (c1's own `IMPORT_KEYS` declared **44** roots this way — a valid non-minimal template).
- **MUST additionally declare `coupling_kernel` as the mandatory 19th root** — a §5a ORDERING-ABLATION **regression**
  operand ONLY. It is neither among the 18 nor closure-pulled (verified: `COUPLING_KERNEL_IN_18_ROOT_CLOSURE=False`),
  and the bind-closure design **promotes** a declared regression operand into `IMPORT_KEYS`
  (`export_ledger_bind_closure_design.md:148`). ⛔ It is **never** a §3a/§3c construction operand (c2 re-extracts the
  closed kernel from the closed operator).
- **MUST NOT declare** `v_bulk_normal_0` (in the fold but **not** referenced by c2's consume-set → declared-but-unused
  fails) or any other key the build does not look up.
⇒ The build reports its final `IMPORT_KEYS`; the two decision legs verify it against `assert_lookups_equal_manifest`
on the **real** two-parent fold (⛔ not by reading — run it).

**⛔ Provenance the guard will NOT catch — the directive + legs must:** bind the **`s11c_c1_`-prefixed**
`s11c_c1_face_response` / `s11c_c1_face_response_coeffs` (step `S11c-c1`), ⛔ **NEVER** the bare `face_response` /
`face_response_coeffs` (step `S11b`, the open/flat regression operands c1 imported). Both keys exist and their values
differ; `check_consumer` and `assert_lookups_equal_manifest` pass on **either** (key-existence only), so a swap is
silent — it folds the wrong physics (the N14/F9 false-equal). ⛔ `coupling_kernel` is bindable **only** as the §5a
ORDERING-ABLATION regression operand (`export_ledger_bind_closure_design.md:148-153`), ⛔ **never** as a §3a/§3c
construction operand (c2 re-extracts the closed kernel from the closed operator).

---

## 3 · The fold symbol map — VERIFIED against the real rows (build-mechanical; spec §3a "the build directive freezes the exact symbol map")

The two engines use **disjoint names** for the same objects (verified: the slab δp_s slots do not appear in
`s11c_c1_face_response`), so the fold is an **explicit, per-case `subs`**, ⛔ never a silent name identification. The
four slots are **bare `Symbol`s** (class `COORDINATE`, step `S11c-b`) present in **`slab_operator`,
`closure_shape_deriv` (θ-row flux), AND `traction` (mechanical `t_s`)**:
```
delta_p_plus , delta_p_minus , d_w_delta_p_plus , d_w_delta_p_minus
```

**3.1 · The closed-pressure source — `DELTA_P` in `s11c_c1_face_response` (⛔ NOT `PRESSURE` in `_coeffs`).** The
closed operator-valued `δp_s` per case is `s11c_c1_face_response → CASES → (anchoring∈{LAB_HELD,MATERIAL_ADVECTED},
face∈{+1,−1}, density∈{RHO4_CONSTANT,RHOBR_CONSTANT}) → VALUE → `**`DELTA_P`** (quantity keys
`['RESOLVENT','RESOLVENT_DEFINITION','MU_S','DELTA_P','J','T']`). ⛔ **NOT** the `PRESSURE` object on
`s11c_c1_face_response_coeffs` — that is the **scalar flat/coefficient** regression object (its
`RESOLVENT_SYMBOL_OUTPUT` is a scalar `1/(Λ_A⁰ω/(ρ_m q_out(1−iωτ_A))+1)`, no nonlocal `Z`, no
`s11cc1_V_*`/`s11cc1_mu_theta_*`); binding it is exactly the "scalar division" this directive forbids and the
close-the-extracted-kernel route the design bars. `DELTA_P` carries the c1 symbols (`s11cc1_V_*`, `s11cc1_mu_theta_*`,
the resolvent, `s11cc1_dtn_operator_*`) — the genuine nonlocal closed pressure.

**3.2 · The w-jets are COMPUTED, not imported.** c1 exports **no** `d_w_delta_p_*` object (verified:
`C1_RESPONSE_EXPLICIT_D_W_DELTA_P_REFS=[]`). `d_w_delta_p_±` are genuine **normal jets** (the same shape-derivative
that produced the OPEN slots, S11c-a `S11c_a_interface_geometry_sympy_audit.py:150`). ⇒ **compute** `d_w_delta_p_±` as
that shape-derivative `d_w` applied to the closed `DELTA_P` of §3.1, then substitute **both** `(delta_p_±,
d_w_delta_p_±)` slots. ⛔ Do not leave the `d_w_delta_p_±` slots unsubstituted, and ⛔ do not point at a nonexistent c1
export for them.

**3.3 · ε-normalization — reconcile ONE ε (⛔ do not double-count to O(ε²)).** c1 builds its inputs as `ε·symbol`
(`S11c_c1_bulk_closure_sympy_audit.py:834`), so the closed `DELTA_P` is **ε-scaled**; the slab θ-row slot carrier is
**itself** multiplied by an outer ε (the real row `(1/4)·ε·[−4I·Λ_A(−delta_p_plus/ρ_m + μ_θ/ρ_br)/(ωτ_A+I) + …]`). A
**direct** substitution therefore double-counts ε → **O(ε²)** (verified: coeff·`DELTA_P` factors an
`epsilon_shape**2`), whereas the transverse↔thickness coupling is the inherited **O(εη)** (§4/N12). ⇒ **reconcile the
ε-convention so the closed operator carries the correct `(ε,η,σ_W)` order** — substitute the ε-consistent closed
pressure (e.g. `DELTA_P/epsilon_shape` against the ε-bearing slot, or the equivalent single-ε normalization the build
freezes and **states**) — and **EMIT the `(ε,η,σ_W)` multigrade** of the closed operator and the increment as the
order-check (⛔ PRINT the order; ⛔ do not assert an expected value).

**3.4 · Stage-2 identifications — the FOUR c1 identifier classes → slab fields (per case).** After the §3.1
substitution the slab rows carry c1's opaque symbols; identify each explicitly (the build freezes the map, states
which line applies it; the legs verify each against the real rows):
| c1 symbol (per case) | slab field / bridge | ⚠ note |
|---|---|---|
| `s11cc1_mu_theta_{lab_held,material_advected}_{plus,minus}` | slab `mu_theta_L` (lab_held) / `mu_theta_M` (material_advected) | **ASYMMETRY (verified OK)**: slab μ_θ is per-anchoring only; c1 μ_θ is per anchoring×face — both faces' c1 μ_θ map to the single slab `mu_theta_{L/M}` (μ_θ is the face-independent held-fixed operand, §1d; resolve via `mu_theta_operator` per `(anchoring,density)`). |
| `s11cc1_V_{lab_held,material_advected}_{plus,minus}` | the slab **`face_velocity`** VALUE for the matching `(α, s, {DELTA_W\|ZETA_C})` case | ⛔ **NOT `kinematic_balance`** (that is the `delta_v_bulk_*` OPERAND_A/B/RESIDUAL identity — a different object). `face_velocity` has 8 cases (`DELTA_W` e.g. `W_0·e_{W,t}·ε/2`, `ZETA_C` e.g. `ε·ζ_{c,t}`, per α×face); ⛔ **pin the representation** — the build states whether it closes the `DELTA_W` or `ZETA_C` case and uses it consistently. |
| `s11cc1_dtn_operator_{α}_{s}` (whole-form `Z`) **and** `s11cc1_response_resolvent_{α}_{s}_{ρ}_constant` = `[I+(Λ_A/ρ_m²)Z]⁻¹` | **bridge each to the AGREE'd two-momentum kernel `dtn_kernel`** per case (its `q_out(k),q_out(k′)` legs + on-shell dispersion) | **§3d.3 — the load-bearing bulk object is the KERNEL, ⛔ not the whole-form.** ⛔ Do NOT let the raw whole-form `dtn_operator` become the silent construction operand; freeze the per-case kernel representation and emit the whole-form dependence **separately** as `S11CC2_DTN_WHOLEFORM_DEPENDENCE` for ablation. ⛔ Do not assume whole-form AGREE. |

⛔ **Order of operations:** bind the RHO4-live density **before** §3.1 (§3d.1) — `rho_br_bg_rho4_constant →
W_bg·ρ_br/W_0` from `background_density_map`, ⛔ never a bare constant. Then §3.1 substitute `DELTA_P`, §3.2 the
computed w-jets, §3.3 the ε-consistent normalization, §3.4 the four identifications. ⛔ The build **reports the fully
frozen map** (source tag, ε-normalization, each of the four identifications, per case) and states which line applies
each; the legs verify against the real rows.

---

## 4 · Objects to emit (spec §4) + provenance
Per anchoring `α∈{L,M}` and density representative `ρ∈{ρ_4D,ρ_br}`, each object carrying its computed `(ε,η,σ_W)`
multigrade and restored `[L,T,M]` dimension (⛔ no object reported without both):
- `S11CC2_CLOSED_SLAB_OPERATOR` (assembled two-face), `…_TERM_ORIGINS` (per face), `…_PARITY_BLOCKS` — §3a.
- `S11CC2_CLOSED_COUPLING_KERNEL` (both off-diagonal blocks), `…_TERM_ORIGINS`; `S11CC2_SELF_ENERGY_ADJOINTNESS_RESIDUAL` (only if an independent route exists) — §3b.
- `S11CC2_SELF_ENERGY_INCREMENT` **and its two same-extract operands** (`extract(close(SLAB))`, `extract(SLAB)`) — §3c.
- The six §3d objects: `S11CC2_DENSITY_LIVE_MINUS_FROZEN`, `S11CC2_TRACTION_MECHANICAL_CONTRIB`,
  `S11CC2_DTN_WHOLEFORM_DEPENDENCE`, `S11CC2_TRACTION_SLAB_POWER_PAIRING` (+ its residual),
  `S11CC2_FLAT_SYMBOL_USAGE`, and the `μ_R,bg`-form ablation output.
- The §5 controls: `S11CC2_ORDERING_COMMUTATOR`, `S11CC2_REP_INVARIANCE_RESIDUAL`, and each limit/routing/form
  residual — **each emitted as the object and its literal residual**.

**F9 / N14:** every new object gets a **fresh injective** `mechanical_lower_camel` write-key; ⛔ never reuse an
imported key (`slab_operator`, `coupling_kernel`, `closure_shape_deriv`, `dtn_operator`, `dtn_kernel`,
`face_response`, `traction`, any constant). Per emitted object, the builder report states **which line computed it**.

---

## 5 · Export schema (spec §7 / `export_ledger_bind_closure_design.md` §D1–D3)
Write `scripts/S11c_c2_exports.py` as an **own-rows delta** (⛔ not the accumulated whole-model file). Membership =
the bind-closure (D1); `assert_delta_is_minimal` requires the delta's key-set = c2's own bind-closure ∪ infra.
`BUILD_INPUT_DIGESTS` pins `{this SymPy audit, scripts/S11c_b_exports.py, scripts/S11c_c1_exports.py, the c2 spec,
scripts/ledger_fold.py}` (§D3). ⛔ Never `git add -f` a big `.out`; ⛔ never annex an `*_exports.py`.

---

## 6 · Supplied / withheld (leak discipline)
Everything in spec §1 is **SUPPLIED and unfalsifiable in this build** — state so in the builder report so a passing
build does not read as if it verified it. The **outputs** are §§3–5; ⛔ this directive and the script state **no**
component value, sign, order, parity, or grade, and carry **no acceptance criterion referencing an expected value**
(the diff happens on OUR side; a genuine cross-engine disagreement is a **finding**, ⛔ not a build failure). The six
§3d dispositions and every §5 residual are read in the **step record**, ⛔ not asserted by the script.
⚠ **Before launch, leak-gate this directive**: `rg` for co-occurring step symbols (e.g. `delta_p` with a typed
coefficient) and read every hit; a residual "must vanish"/"must equal" phrasing is a leak — remove it.

---

## 7 · Run discipline
Detached launch (`setsid` + a completion marker + `Monitor`; the harness reaps `run_in_background` ~87 s). c2's fold
threads the nonlocal `Z`, so the full increment may be **heavy** — measure the process that runs; defer heavy
controls in-band→out-of-band (`DEFERRED_HEAVY_RUNS.md`); ⛔ never run two memory-heavy CAS jobs concurrently. Astra
needs `--sandbox danger-full-access` to run SymPy. Verify the **deliverable** (script + export exist, non-empty,
plausible token count), ⛔ not the exit status; then the two re-review legs (Codex-written → **fresh Claude agent +
Grok**) launch on sight.
