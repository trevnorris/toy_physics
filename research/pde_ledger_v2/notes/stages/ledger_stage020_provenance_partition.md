# ledger_stage020 — the `54/5=2·27/5` provenance partition + the CALIBRATED verdict label (Check II-G4c)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 radiative sector). The THIRD leg of a 4-way split of
`pathA_33`: this stage carries the **`54/5=2·27/5` provenance-partition + the CALIBRATED-verdict-label component (3/4) of the
joint `QUAD_CALIBRATED`** — the leg that LANDS the CALIBRATED headline. The outgoing ℓ=2 DtN Hankel fingerprint + χ_Q sign
was **stage 018** (II-G4a, done); the squared-denominator prefactor algebra was **stage 019** (II-G4b, done); the μ̂₀-free
`[P₀^phys]=1` dimensional closure is **stage 021** (II-G4d, the COMPLETING leg).

**Verdict.** `QUAD_CALIBRATED` (JOINT, 4-stage) — landed here as a **PARTIAL** component (018/019 done, 020 EARNED; 021
PENDING). Ledger earned-label: `PROVENANCE_PARTITION_CALIBRATED`. ⭐ This is the leg where the CALIBRATED *classification* is
DETERMINED — the honest verdict is `QUAD_CALIBRATED` (not `QUAD_PASS`) **because** the assembled magnitude `54/5` and Newton's
`G` are `external_bridge_input` by provenance.

**Status.** Exact closed-form / symbolic / rational (units-bearing but float-free for the earned content): the a⁻⁵ scaling
(exact `a`-power via 018's frozen slot), the Γ5/χ_Q equivalence (`compact(expr)==0` residuals), the `54/5=2·27/5` SymPy
rational identity (BOUND to `target_rhs`/`v5_slot`), and the 4-way provenance partition (exact class comparison via tag
dominance); `expect_zero`/`expect_bool`/`expect_fail` asserts, no floats/tolerances in the earned content. Dual-engine SymPy
**74 PASS** / Mathematica **82 PASS**, exit 0, CWD-independent.

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_33_quadrupole_normalization_{sympy.py,.wl}` +
> `reports/pathA_33_quadrupole_normalization.md` (the 020 slice = report :3, :5, :29–30, :36, :39–40, :49) + the original
> directive `directives/pathA_33_quadrupole_normalization.md`. The report/directive are cited for provenance only; the
> derivation below is self-contained. ⚠ The source `build_partition` emitted `54/5 = 2*27/5` as a **typed string**; this
> reshape makes it a SymPy-VERIFIED identity bound to 020's own assembled expressions (§1.2).

---

## 1. What this stage earns

The ℓ=2 quadrupole normalization's assembled magnitude is a CALIBRATED bridge to General Relativity — and this stage makes
the earned-vs-calibrated split precise, and PROVES the verdict is `QUAD_CALIBRATED` (not `QUAD_PASS`) for the right reason
(provenance), immune to a G-invariance trap.

### 1.1 The assembled magnitude and the Γ5/χ_Q equivalence
The ℓ=2 quadrupole radiative normalization the gravity sector must match is the Burke–Thorne form
```
    m̂₀²P₀ = target_rhs = 54·G·c_s⁵/(5·a⁵·c⁵)   ⟺   γ_quad^eff = 2·G/(5·c⁵)   (Burke–Thorne)
```
The equivalence is verified as a self-contained algebraic bridge. With `P0_physical=(c_s/a)²·(N0/D0)` (019's prefactor) and
018's outgoing `χ_Q=+1`, the forward/reverse residuals
```
    forward = target_rhs·χ_Q·a⁵/(27·c_s⁵) − 2G/(5c⁵) = 2G(χ_Q−1)/(5c⁵)
```
vanish **iff χ_Q=1 and 54/27=2** (the `27` is 018's radiative slot `1/v₅ᶻ`). So the outgoing sign 018 derived is exactly the
sign General Relativity's positive quadrupole demands; a wrong `2/5→3/5` breaks the bridge (`3e` → `FAIL_EQUIVALENCE`).

### 1.2 The `54/5 = 2·27/5` decomposition (the headline — a bound SymPy identity)
The assembled coefficient factors as
```
    54/5 = 2 · 27/5 ,   with   the 27 = derived_in_gate (018's fingerprint 1/v₅ᶻ)
                                the 2/5 + G = external_bridge_input (GR Burke–Thorne 2G/5c⁵, G = GENUINE_BLOCKED)
```
⭐ This is a **SymPy-VERIFIED rational identity BOUND to 020's own assembled expressions** — NOT the source's typed string
`"54/5 = 2*27/5"`, and NOT a bare-literal `Rational(54,5) − 2·Integer(27)/5` arithmetic tautology. The magnitudes are
EXTRACTED: `mag = compact(target_rhs / (G·c_s⁵/(a⁵·c⁵)))` → `54/5` (from the assembled magnitude), `27_from_slot =
compact(a⁵/(c_s⁵·v5_slot))` → `27` (from 018's frozen radiative slot `v5_slot = a⁵/27c_s⁵`), then
`compact(mag − 2·27_from_slot/5) == 0`. Mutating `target_rhs` (54→55) or the identity (`2·26/5`) fires the assert — the
decomposition tracks the actual assembled magnitude, not a self-evident literal.

### 1.3 The a⁻⁵ target scaling (strengthened via 018's frozen slot)
The assembled magnitude carries `a⁻⁵` as a DERIVED consequence of a-cancellation, not a typed target power: the Burke–Thorne
coefficient `gamma_target=2G/5c⁵` is a-free (`a_power=0`), and 018's radiative slot `v5_slot=a⁵/(27c_s⁵)` supplies the `a⁵`,
so `derived_power = a_power(gamma_target) − a_power(v5_slot) = 0 − 5 = −5`, and `target_from_bridge = gamma_target/(χ·v5_slot)
= target_rhs` (a-power −5). The `a⁵` is typed ONCE (in `v5_slot`). The wrong-scaling probe (`a⁵→a⁴` on the assembled target)
fires at the assembled residual (`3c` → `FAIL_SCALING`). ⚠ HONEST: the ACTUAL branch a-scaling (from `N0/D0`) is DEFERRED
Gate-6 (report :49); 020 checks the TARGET rhs a-power.

### 1.4 The 4-way provenance partition and the PROVENANCE-driven CALIBRATED verdict
Every quantity is classified from its provenance TAGS by `classify_provenance`, a tag-DOMINANCE rule
`deferred > external > derived > convention`:
```
    fingerprint_27 {dtn_radiative_slot}                       → derived_in_gate
    PN_2_over_5     {external_pn_bridge}                       → external_bridge_input
    G              {external_gr_constant}                      → external_bridge_input
    assembled_54/5 {external_pn_bridge, dtn_radiative_slot}    → external_bridge_input   (external DOMINATES)
```
The assembled `54/5` mixes an EARNED tag (the `27`, from the DtN radiative slot) with an EXTERNAL tag (the `2/5`, from the
PN bridge); **external dominates**, so the assembled magnitude is `external_bridge_input`. The 020-local verdict then reads
the classes with the source-faithful rule
```
    if g_class == derived_in_gate AND mag_class == derived_in_gate:  QUAD_PASS
    else:                                                            QUAD_CALIBRATED
```
Since `G` and the assembled `54/5` are both `external_bridge_input`, the verdict lands **`QUAD_CALIBRATED`** — the honest
landing: the PDE delivers the FORM/branch (the `27` fingerprint), NOT Newton's `G` (`G=GENUINE_BLOCKED`).

### 1.5 ⭐ The verdict is PROVENANCE-driven, NOT G→λG-invariance-driven (the invariance-only trap)
A separate, non-verdict-driving diagnostic makes the point sharp: `54/5` is a pure number, so it is **G-invariant**
(`.subs(G, λ_G·G)` leaves it unchanged) — yet it is `external_bridge_input` by provenance. So an invariance-only test would
MISLABEL the calibrated `54/5` as earned (`invariance_only_trap_catches_54_over_5 = True`). The classification is set by the
provenance tags, NOT by G-invariance. This is why the partition (not an invariance check) drives the verdict.

### 1.6 Units-bearing scope (no dimensional-homogeneity gate)
⭐ Unlike stage 019 (units-free abstract algebra), 020 IS units-bearing — `{c_s, a, c, G}` live in the scaling/equivalence
algebra (the `54Gc_s⁵/5a⁵c⁵` target, the `2G/5c⁵` bridge). BUT 020 does ALGEBRA + PROVENANCE, **not** the dimensional
homogeneity gate `[P₀^phys]=1` (that is stage 021's). The 020/021 cut is by OPERATION (algebra+provenance = 020; dimension =
021), enforced by a runtime free-symbol NAME allowlist (`{a,c_s,c,G,N0,D0,chi_Q,lambda_G}`; no `μ̂₀`/`mtilde0` symbol — the
source's legacy `gamma_quad_eff` μ̂₀ string is DROPPED), and a structural no-dimensional-dependency cut (no `dim_of`/`[·]`
construct; a structural assertion that the 020-local verdict's signature/bytecode carries no dim parameter/dependency — the
`.wl` checks its definition arity/text).

---

## 2. The able-to-fail battery (020-owned)

The verdict runs a SCOPED 020-local gate chain (`scaling_ok ∧ equivalence_ok ∧ provenance_ok` + the provenance-driven
CALIBRATED/PASS determination); 018/019/021's fingerprint/prefactor/dimensional gates are NOT computed here. The 020 teeth:

| tooth | mutation → verdict | notes |
|---|---|---|
| a⁻⁵ scaling (`3c`) | `a⁵→a⁴` on the assembled target → the bridge residual ≠ 0 → `FAIL_SCALING` | the `−5` = `a_power(gamma_target)−a_power(v5_slot)` from 018's frozen slot (`a⁵` typed once) |
| Γ5/χ_Q equivalence (`3e`) | `2/5→3/5` → `wrong_reverse ≠ 0` → `FAIL_EQUIVALENCE`; corrupt χ_Q → `forward ≠ 0` | self-contained bridge; `forward=2G(χ−1)/5c⁵` |
| `54/5=2·27/5` bound identity | mutate `target_rhs` (54→55) or `2·27/5` → `compact ≠ 0` fires | `mag` from `target_rhs`, `27` from `v5_slot` — not a bare literal |
| classifier dominance (truth-table + key-class + tag-mutation) | wrong dominance rule / wrong key class / strip a tag → the truth-table/class assert fires | ⚠ NOT `partition_ok` (trivially True since emitted defaults to computed) |
| `3f_partition_mislabel` (BOTH directions) | `{"G":"derived_in_gate"}` AND `{"fingerprint_27":"external_bridge_input"}` → `FAIL_PROVENANCE_PARTITION` | each with a DYNAMIC 020-local self-ablation |
| CALIBRATED verdict (QUAD_PASS + REQUIRED MIXED controls) | constant-CALIBRATED → the QUAD_PASS control fires; the inverted `PASS-unless-both-external` rule → the MIXED control fires | proves the verdict READS provenance with the source-faithful default |
| g-invariance trap | force a G-invariance-based classifier → the `54/5` mislabel is caught by the partition | the diagnostic is SEPARATE / non-verdict-driving |
| free-symbol allowlist / structural dim-cut | inject a `mu_hat0` symbol → allowlist fires; wire a dim param into the verdict → structural cut fires | enforces units-bearing-but-no-dim-gate (not a source grep) |

Adversarial per-tooth ablation (`ADVERSARIAL_ISSUES` → remediated → `REVERIFY_CLEAN`): ~65 driven `.py` mutations + 7 native
`.wl` mutants. All four key genuineness risks — the bound `54/5=2·27/5` identity, the MIXED control catching the inverted
rule, the proven classifier, and the genuinely-authored `.wl` — confirmed CLEAN. Three LOW-severity vacuous/subsumed teeth
were flagged and remediated make-genuine (no de-counts): a near-tautological `dimensional_ok`-independence `f(x)==f(x)` →
a structural signature/bytecode assertion (fires if a dim input is wired into the verdict); a subsumed `P0=N0/D0`-disjoint
firewall → a `{N0,D0}⊆Gamma5 ∧ ∉residual` before/after that runs before the residual-form asserts; a subsumed tag-mutation
`≠EXTERNAL` → an independently-classified baseline-external "from" side complementing #18's "to" side. Fresh-agent re-verify
confirmed each fires at its own named assert and goes vacuous when the fix is neutered, both engines, no regression.

---

## 3. Honest scope

- **EARNED classification / CALIBRATED magnitude.** 020 EARNS the a⁻⁵ scaling, the Γ5/χ_Q equivalence, the `54/5=2·27/5`
  decomposition, and the 4-way provenance partition (the classification machinery). The `54/5` MAGNITUDE and `G` are
  CALIBRATED (`G=GENUINE_BLOCKED` — the PDE delivers the branch/form, not Newton's `G`).
- **`c` is a GR/PN units bridge, not EM propagation.** The `c⁵` in `54Gc_s⁵/5a⁵c⁵` and `2G/5c⁵` is the GR-matching / λγ
  units bridge (`P₀ ∝ c_s⁵/c⁵ = 1/λγ⁵`); `c` = the light cone `c_γ` in its GR-bridge role, a benchmark carrier cited from
  the PN corpus — not a fresh model knob, not EM propagation. (The gravity radiation rides the medium's own density ripple
  at `c_s`.)
- **The actual branch data is deferred.** The ACTUAL branch a-scaling and the numerical `(D_n, N_n)` port scalars remain
  Gate-6 sim-deferred (report :49); 020 checks the TARGET rhs's structure, with the numerical branch realization deferred.
- **Units-bearing, not dimension-closing.** 020 carries `c_s`/`a`/`c`/`G` in the algebra (a reversal from 019); the
  dimensional homogeneity closure `[P₀^phys]=1` is stage 021's.

---

## 4. Consumed / exported

- **Consumed — PROVENANCE (cite, no theatrical dual-site; §1c of the directive).** 018's χ_Q=+1 and the `27` (`=1/v₅ᶻ`)
  genuinely ENTER 020's self-contained equivalence bridge (a corrupt cited value breaks `forward` — that is the equivalence
  tooth itself, not a manufactured cross-stage guard); 019's `P0=N0/D0` appears ONLY in the `Gamma5` definition, absent from
  the equivalence `ok` residual (a `{N0,D0}⊆Gamma5 ∧ ∉residual` firewall guards this); 017's ℓ=2 port-kernel D-lanes are the
  provenance of the abstract port scalars. NO cross-stage dual-site (the earned content is 020's own self-contained bridge +
  partition; 018/019/017 supply the provenance of why χ=1, why 27, what `P0` is).
- **Register.** ZERO new counted knobs (a classification slice, like 018/019). New structural edge **R39** (the
  `54/5=2·27/5` provenance-partition + the CALIBRATED-verdict provenance: the assembled `54/5` and `G` are
  `external_bridge_input`, so the verdict lands `QUAD_CALIBRATED` not `QUAD_PASS`; the invariance-only trap; discharges
  nothing — a classification landing, not a reduction). `G=GENUINE_BLOCKED`/external (already registered); `c` = the `c_γ`
  GR-units-bridge (cited benchmark); the `2/5` = GR bridge; the `27` = 018's derived fingerprint. Part-II counted CALIB set
  UNCHANGED at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) + `{T_Ω, β₂}` (017) = 6.
- **Exported.** The `54/5=2·27/5` provenance partition + the CALIBRATED verdict label + the Γ5/χ_Q equivalence chain + the
  a⁻⁵ target scaling → stage 021 (the μ̂₀-free `[P₀^phys]=1` dim closure completes the joint) + stage 022 (Gate-4
  non-regression) + stage 027 (pathA_43 closure slot).

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O. ⭐ Unlike stages 018/019 (which kept the source `.wl`'s
already-native fingerprint/prefactor blocks and severed only the YAML), the source `.wl` had **ZERO 020 content** — so the
020 `.wl` is a **genuinely fresh independent authoring**: native `Exponent`/`Together`/`Cancel`/`FullSimplify` for the
scaling + the Γ5 bridge residual, native `Rational`/`Simplify` for the `54/5=2·27/5` identity bound to its own
`targetRHS`/`v5Slot`, a native `Association`-based provenance classifier using a rank-based `MaximalBy` (a genuinely
different mechanism from the `.py` if-chain, same dominance order), and `expr /. G -> lambdaG G` for the g-invariance trap —
NO `Series` (the 020 slice has none), NO `.py` mirroring. Agreement is transcript-level: both engines emit the same `a⁻⁵`,
`forward=0`/`reverse=0`, `54/5=2·27/5` TRUE (bound), the provenance classes (`G`:external, `27`:derived, assembled
`54/5`:external), the CALIBRATED verdict with the QUAD_PASS + MIXED controls, and both `3f` directions firing. The stage-007
unevaluated-leakage failure mode is actively guarded (arity + leakage self-check over the nine authored helpers).

**Directive review** used the Codex→Grok→Codex bookend: Codex `DIRECTIVE_ISSUES` (4 BLOCKING genuineness gaps + 3 nits —
the strengthened-3c duplicated-literal risk, the trivially-true `partition_ok`, the missing verdict-only positive control,
the unsound no-μ̂₀ mechanism incl. the false "ZERO μ̂₀ refs" claim vs the source L525 string) → all folded → Codex confirm
`DIRECTIVE_CLEAN`; ⭐ **Grok-4.5 compute-verify then caught the CALIBRATED verdict-rule INVERSION** (the source default is
`else → QUAD_CALIBRATED`; the directive had inverted it to `PASS unless both external` — a shippable rig, since both
endpoint tests still pass) + 3 nits → folded (added a REQUIRED MIXED control) → a Codex confirm caught 2 consistency-sweep
gaps (stale bare-literal forms + the mixed control unpropagated) → swept → final Codex confirm `DIRECTIVE_CLEAN`.

**Tri-review** on fresh agents: `FIDELITY_CLEAN` (an independent read hand-re-derived the Γ5 bridge `2G(χ−1)/5c⁵`, the bound
`54/5=2·27/5`, the classify_provenance dominance, and the g-invariance trap; confirmed the `.wl`'s genuinely-different
rank-based classifier and no 018/019/021 leakage) + `ADVERSARIAL_ISSUES` (per-tooth ablation: 4 key risks CLEAN, 3
LOW-severity vacuous/subsumed teeth) → Codex remediated all 3 make-genuine (no de-counts) → fresh-agent `REVERIFY_CLEAN`
(the coupling meta-test on each remediated tooth: fires at its own named assert under mutation, goes vacuous when neutered,
both engines, no regression). Tallies unchanged: SymPy 74 / Mathematica 82.
