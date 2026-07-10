# ledger_stage021 — the μ̂₀-free `[P₀^phys]=1` dimensional closure (Check II-G4d)

**Part / anchor.** Part II — Gravity (the frozen-throat ℓ=2 radiative sector). The FOURTH, COMPLETING leg of a 4-way split of
`pathA_33`: this stage carries the **μ̂₀-free `[P₀^phys]=1` dimensional-closure component (4/4) of the joint
`QUAD_CALIBRATED`** — the leg that lands the joint fold COMPLETE. The outgoing ℓ=2 DtN Hankel fingerprint + χ_Q sign was
**stage 018** (II-G4a, done); the squared-denominator prefactor algebra was **stage 019** (II-G4b, done); the `54/5=2·27/5`
provenance partition + the CALIBRATED verdict label was **stage 020** (II-G4c, done).

**Verdict.** `QUAD_CALIBRATED` (JOINT, 4-stage) — this leg lands the joint **COMPLETE** (018 ∧ 019 ∧ 020 ∧ 021). Ledger
earned-label: `DIMENSIONAL_CLOSURE_CERTIFIED`. ⭐ **COMPLETE ≠ PASS:** the joint token STAYS `QUAD_CALIBRATED` (calibrated,
not `QUAD_PASS`) — stage 020 already determined the CALIBRATED classification (the assembled `54/5` and `G` are
`external_bridge_input`, `G=GENUINE_BLOCKED`); 021 supplies the dimensional closure that makes the assembled magnitude
dimensionally sound and completes the fold, but does NOT flip the verdict.

**Status.** Exact symbolic / dimensional-vector `(L,M,T)`-triple algebra, float-free: the μ̂₀-free gate
(`dim_of(P₀^phys) == ZERO_DIM`), the natural-units-trap catch (dropping the frequency normalization leaves `T²≠1`), the
corrupt-dim scope truth-table, and the anti-v1 back-solve-tautology demotion; `expect_zero`/`expect_bool`/`expect_fail`
asserts, no floats/tolerances. Dual-engine SymPy **42 PASS** / Mathematica **50 PASS**, exit 0, CWD-independent.

> **Provenance.** Reshaped from `software/stage1_solver/tools/pathA_33_quadrupole_normalization_{sympy.py,.wl}` +
> `reports/pathA_33_quadrupole_normalization.md` (the 021 slice = report :3, :5, :22–25, :37–38, :49) + the original
> directive `directives/pathA_33_quadrupole_normalization.md`. The report/directive are cited for provenance only; the
> derivation below is self-contained. ⚠ The source `.py` kept the μ̂₀ homogeneity back-solve as a verdict-bearing check in an
> earlier form; this reshape DEMOTES it to a labelled non-verdict diagnostic and makes the verdict the μ̂₀-FREE
> `[P₀^phys]=1` gate (§1.4).

---

## 1. What this stage earns

The ℓ=2 quadrupole normalization's assembled magnitude `m̂₀²P₀ = 54Gc_s⁵/(5a⁵c⁵)` must be dimensionally consistent — and
this stage certifies it with a **μ̂₀-FREE** dimensional gate that catches the natural-units trap and reads exactly the
SOURCED port dimensions, immune to the v1 back-solve tautology.

### 1.1 The μ̂₀-free `[P₀^phys]=1` gate
With the SOURCED port dimensions `[N₀]=L⁻¹M` (the density/continuity-port numerator) and `[D₀]=L⁻¹T⁻²M` (the carried reduced
static conservative denominator `D₀=K−B₀−Z₀`), and the frequency normalization `(c_s/a)²`:
```
    [P₀_raw]   = [N₀/D₀]           = T²
    [(c_s/a)²]                     = T⁻²
    [P₀^phys]  = [(c_s/a)²·N₀/D₀]  = 1   ⟹   dimensional_ok = True
```
The verdict-bearing check is `dimensional_ok = (dim_of(P₀^phys) == ZERO_DIM)` — a COMPUTED `(L,M,T)`-triple equality reading
the sourced dims, NOT a typed `"1"`. ⚠ **μ̂₀ (`mu_hat0`) does NOT enter the verdict** — it is absent from `raw_symbol_dims`
(the gate's dimension map), and the verdict's computed read-set excludes `mu_hat0`/`mu_dim`/`homogeneity_pass`.

### 1.2 The natural-units-trap catch (the dimensional milestone)
The handoff's `P₀=N₀/D₀` silently DROPS the `(c_s/a)²` frequency-normalization factor (a natural-units artifact where `c_s/a`
was set to 1). The gate CATCHES it: `[P₀_raw]=T²≠1`, so the `3d` drop-normalization probe fires `FAIL_DIMENSIONAL`. The
frequency normalization is REQUIRED for `[P₀^phys]=1`.

### 1.3 The corrupt-dim scope truth-table (the gate is able-to-fail on the dims that enter `P₀^phys`)
```
    corrupt [N₀]→0  ⟹  [P₀_raw]=(1,−1,2)=L T² M⁻¹ ,  [P₀^phys]=(1,−1,0)=L M⁻¹ ≠ 1  →  FAIL_DIMENSIONAL   (3d′)
    corrupt [D₀]→0  ⟹  [P₀^phys]=(−1,1,−2)=L⁻¹T⁻²M ≠ 1                              →  FAIL_DIMENSIONAL
    corrupt [c_s]→0 ⟹  [P₀^phys]=(−2,0,2)=L⁻²T² ≠ 1                                 →  FAIL_DIMENSIONAL
    corrupt [G]→0   ⟹  [P₀^phys]=1  STILL  (G ∉ free_symbols(P₀^phys))              →  NO_FAIL  (scope diagnostic)
    correct         ⟹  [P₀^phys]=1                                                 →  NO_FAIL
```
⚠ Note corrupt-`[N₀]` gives `[P₀^phys]=(1,−1,0)` — the `(c_s/a)²` factor REMAINS (it is NOT `−[D₀]=(1,−1,2)`, which is
`[P₀_raw]`). The gate fires **precisely on the dims that enter `P₀^phys`** (`{N₀,D₀,c_s}`) and is G-independent (`[G]` is a
dependency-scope diagnostic, asserted via the runtime test `G ∉ free_symbols(P₀^phys)`).

### 1.4 ⭐ The μ̂₀ homogeneity diagnostic is NON-verdict — the anti-v1 discriminator
The full closure `m̂₀²P₀ = 54Gc_s⁵/(5a⁵c⁵)` can be made dimensionally homogeneous by back-solving the free carrier
`[μ̂₀]=L⁻¹T⁻¹M⁻¹ᐟ²` from `(rhs_dim−p0_dim)/2`, giving `[lhs]=[rhs]=L⁻²T⁻²M⁻¹`, `homogeneity_pass=True`. ⚠ **This is a
TAUTOLOGY** — re-solving `[μ̂₀]` after ANY corruption of the dimensions keeps `homogeneity_pass=True` (the back-solve absorbs
any change), so it can NEVER fail. That is exactly the v1-rejected rig, and 021 correctly DEMOTES it to a labelled
diagnostic (`participates_in_verdict=False`, computed from the verdict read-set). The decisive anti-v1 tooth is TWO-part:
(a) a computed read-set assertion that the verdict excludes `mu_hat0`/`mu_dim`/`homogeneity_pass` (and `mu_hat0 ∉
raw_symbol_dims`), and (b) a deliberately-wired MUTANT verdict = the back-solved `homogeneity_pass`, RE-RUN over all four
corruptions and shown to stay `NO_FAIL`/True on ALL of them — proving the back-solve is vacuous, so the real verdict must be
the μ̂₀-free `[P₀^phys]==ZERO` gate (which fires on 3 of the 4). So the μ̂₀-free gate's `[N₀]/[D₀]/[c_s]` FAIL rows (not the
`[G]` row) are what reject the all-`NO_FAIL` back-solved gate.

### 1.5 The `Yhat` dimensionless check
The outgoing physical expansion `Ŷ = 1 + u₂ω² + u₄ω⁴ + i·v₅ω⁵` with 018's frozen fingerprint `u₂=a²/9c_s²`, `u₄=4a⁴/81c_s⁴`,
`v₅=a⁵/27c_s⁵` (a provenance-frozen LOCAL fixture, not the 018 machinery) is dimensionless (each term `T²·T⁻²=1`, etc.), so
`[Ŷ]=1`. The check is wrapped in a structured `try_dim_of`/`Catch` so a corrupted ω-power (`u₂ω³`) fires the NAMED `Yhat
dimensionless` assert (a caught mismatch), not an uncaught `DimError` that would bypass the tooth.

### 1.6 The COMPLETING landing
021 lands the joint `QUAD_CALIBRATED` **COMPLETE** (018 fingerprint ∧ 019 prefactor ∧ 020 provenance-partition+CALIBRATED ∧
021 dim closure). ⭐ The joint token stays `QUAD_CALIBRATED` (calibrated, not PASS) — the PDE delivers the FORM/branch (the
`27` fingerprint + the dimensionally-sound `[P₀^phys]=1` closure), NOT Newton's `G` (`G=GENUINE_BLOCKED`). 021 does not
re-present 018/019/020's accounting; it cites them DONE and declares the 4-leg fold complete.

---

## 2. The able-to-fail battery (021-owned)

The verdict runs an 021-LOCAL gate over `dimensional_ok` (the μ̂₀-free `[P₀^phys]=1` gate); 018/019/020's gates are NOT
computed here. The 021 teeth:

| tooth | mutation → verdict | notes |
|---|---|---|
| μ̂₀-free gate + natural-units trap (`3d`) | drop `(c_s/a)²` (`P₀^phys→P₀_raw`) → `[P₀_raw]=T²≠1` → `FAIL_DIMENSIONAL` | reads the COMPUTED `dim_of(P₀^phys)` vs `ZERO_DIM`, not a typed string |
| corrupt-`[N₀]` (`3d′`) | `[N₀]→0` → `[P₀^phys]=(1,−1,0)≠1` → `FAIL_DIMENSIONAL` | the `(c_s/a)²` factor remains — not `−[D₀]` |
| corrupt-`[D₀]` / corrupt-`[c_s]` | `→0` → `[P₀^phys]≠1` → `FAIL_DIMENSIONAL` | the middle cases that pin the μ̂₀-free default (a back-solve would NOT fire) |
| corrupt-`[G]` scope control | inject `G` into `P₀^phys` → `G ∈ free_symbols` → fires; the corrupt-`[G]`-dim row is correctly `NO_FAIL` | the runtime `free_symbols` G-independence carrier |
| anti-v1 read-set exclusion | wire `mu_dim`/`homogeneity_pass` into the verdict read-set → fires | the verdict's read-set is computed/traced, excludes the forbidden names |
| anti-v1 wired back-solve mutant | the back-solved `homogeneity_pass` mutant is RE-RUN per corruption → all-`NO_FAIL` (vacuous) → the real gate must be `[P₀^phys]==ZERO` | genuinely re-run, not stamped |
| μ̂₀ `participates_in_verdict` (make-genuine) | add a forbidden name to the read-set → `participates_in_verdict=True` → fires | derived from the computed read-set |
| DYNAMIC 021-local self-ablation (`rerun_gate_logic`, make-genuine) | neuter the self-ablation to a constant (with==without / single call) → `rerun_gate_logic=False` → fires | a two-verdict re-run over the 021-local verdict, NOT the joint `base_verdict`, NOT a constant |
| `Yhat` structured-catch | corrupt a term's ω-power (`u₂ω³`) → the NAMED `Yhat dimensionless` assert fires | `try_dim_of`/`Catch`, not an uncaught raise |

Adversarial per-tooth ablation: `ADVERSARIAL_ISSUES` → the two KEY teeth (the anti-v1 read-set/wired-mutant, and the
corrupt-`[G]` scope control) confirmed GENUINE; 5 LOW-severity stamped/subsumed teeth (a stamped `rerun_gate_logic`, a stamped
`participates_in_verdict`, a tautological QUAD-landing comparison, and two subsumed aggregate summaries) → remediated (2
make-genuine: `rerun_gate_logic` derived from the actual re-run, `participates_in_verdict` derived from the computed read-set;
3 de-count: the QUAD-landing tautology + the two subsumed summaries retained as labeled prints) → fresh-agent `REVERIFY_CLEAN`
(coupling meta-test: each make-genuine tooth fires under a mutation of its object and goes vacuous when the fix is neutered;
the de-counts do not lose per-row/per-mutant coverage; no KEY-tooth/earned-logic regression). Tallies SymPy 45→42 /
Mathematica 53→50 (−3/engine from the honest de-counts).

---

## 3. Honest scope

- **EARNED closure / CALIBRATED magnitude.** 021 EARNS the μ̂₀-free `[P₀^phys]=1` dimensional closure, the natural-units-trap
  catch, and the corrupt-`[N₀]`/`[G]` discriminator (the gate reads the sourced port dims, not a back-solved μ̂₀). The `54/5`
  MAGNITUDE and `G` are CALIBRATED (`G=GENUINE_BLOCKED`); the dim closure makes the FORM dimensionally sound, it does NOT
  deliver `G`.
- **Units-bearing, dimension-checking.** 021 is the OTHER half of the 020/021 operation-level cut — 020 did algebra +
  provenance; 021 does the `[·]` dimensional-homogeneity gate. It is the only pathA_33 leg with a dimensional gate.
- **The sourced port dims are provenance.** `[N₀]=L⁻¹M` (pathA_43 density-port numerator) and `[D₀]=L⁻¹T⁻²M` (the carried
  reduced static conservative denominator `D₀=K−B₀−Z₀`) are SOURCED dimensional inputs (Cluster C stages 024/027 + pathA_34
  not yet built in build order); they genuinely enter the gate (so the corrupt-`[N₀]` tooth is genuine), but 021 does not
  re-derive them.
- **The branch data is deferred.** The actual branch a-scaling and the numerical `(D_n,N_n)` port scalars remain Gate-6
  sim-deferred (report :49).

---

## 4. Consumed / exported

- **Consumed — PROVENANCE (cite, no theatrical dual-site).** The SOURCED port dims `[N₀]=L⁻¹M` / `[D₀]=L⁻¹T⁻²M` (pathA_43
  density-port numerator + the carried reduced static conservative denominator), 018's fingerprint `u₂/u₄/v₅` (a LOCAL frozen
  fixture for the `Yhat` check, not `build_fingerprint`), 019's `P0=N0/D0` (the prefactor structure entering `P0_raw`), and
  020's assembled `target_rhs` (the μ̂₀ DIAGNOSTIC's rhs only). The sourced port dims genuinely ENTER the μ̂₀-free gate (the
  corrupt-`[N₀]` probe is the genuine able-to-fail tooth); no cross-stage dual-site.
- **Register.** ZERO new counted knobs (a dimensional-closure slice, like 018/019/020). `μ̂₀` is a free-carrier
  non-verdict diagnostic (NOT a counted knob); the SOURCED port dims are dimensional provenance; `c`/`G` are already
  registered. New structural edge **R40** (the μ̂₀-free `[P₀^phys]=1` dimensional closure: dimensionless from the sourced port
  dims, μ̂₀ non-verdict; the natural-units trap CAUGHT; corrupt-`[N₀]`/`[D₀]`/`[c_s]` fire while corrupt-`[G]` does not; the
  back-solve tautology demoted — discharges nothing, like R37/R38/R39; records that 021 COMPLETES the joint `QUAD_CALIBRATED`,
  which stays CALIBRATED not PASS). Part-II counted CALIB set UNCHANGED at `{μ_η, T_w, β}` (013) + `{Vp0/ℓ_c}` (015) +
  `{T_Ω, β₂}` (017) = 6.
- **Exported.** The μ̂₀-free `[P₀^phys]=1` dimensional closure + the COMPLETED joint `QUAD_CALIBRATED` (all four Gate-4 legs)
  → stage 022 (Gate-4 non-regression: pathA_34 cross-ℓ cites the completed quadrupole normalization) + stage 027 (pathA_43
  port checks + the `P0_phys=(c_s/a)²N0/D0` closure slot — the same dim closure, marked shared) + the Part-VII whole-system
  dimensional-firewall check.

---

## 5. Dual-engine and verification

Both engines are standalone, print-only, assert-zero, ZERO file I/O. ⭐ Like stages 018/019 (and unlike 020's fresh
authoring), the source `.wl` HAS the 021 dimensional block, so the 021 `.wl` KEEPS its native, already-independent Wolfram
dimensional-vector engine (a `Which`-based `dimOf` over `Times`/`Power`/`Plus`, native `<|…|>` `rawDims`, `Series`-free
dimension algebra) and severs only the YAML machinery; the 018 fingerprint (`j2`/`y2`) and 019 prefactor blocks are removed
and the `Yhat` coefficients become local `u2Sourced`/`u4Sourced`/`v5Sourced` fixtures. It ADDS (the source `.wl` lacked
both): the corrupt-`[G]`-does-not-fire scope control and the DYNAMIC 021-local self-ablation. The `.wl`'s local pass-guard
uses a split-out `yhatOk` and excludes `homogeneityPass` (which is a print-only diagnostic). Agreement is transcript-level:
both engines emit `[P₀_raw]=T²`, `[(c_s/a)²]=T⁻²`, `[P₀^phys]=1`, `dimensional_ok=True`, the `3d`/`3d′` `FAIL_DIMENSIONAL`,
the corrupt-`[G]` `NO_FAIL`, and `homogeneity_pass=True` (diagnostic). The stage-007 unevaluated-leakage failure mode is
actively guarded (def/call arity scan + leakage transcript scan over the `.wl` helpers).

**Directive review** used the Codex→Grok→Codex bookend: Codex `DIRECTIVE_ISSUES` (6 BLOCKING — most notably that the
originally-proposed corrupt-`[G]` anti-v1 discriminator was BACKWARDS: a back-solved μ̂₀ gate re-solves `mu_dim` after every
mutation, so it fires on nothing, and it is the μ̂₀-free gate's `[N₀]/[D₀]/[c_s]` FAIL rows that reject it; plus the
corrupt-`[N₀]` dim `(1,−1,0)` not `−[D₀]`, the `homogeneityPass`-in-the-`.wl`-guard, the self-containment fixture, the `[D₀]`
provenance, and the `Yhat` structured-catch) → all folded → Codex confirm (1 residual dim-label fix) → **Grok-4.5
compute-verify `DIRECTIVE_CLEAN`** (SymPy-confirmed the μ̂₀-free gate dims, the corrupt-dim truth-table, the
back-solve-is-a-tautology crux, and the `Yhat`/structured-catch; independently reproduced the adversarial contrast that the
original framing was valid only for a *pinned* μ̂₀) → final Codex confirm `DIRECTIVE_CLEAN`.

**Tri-review** on fresh agents: `FIDELITY_CLEAN` (an independent read hand-re-derived the μ̂₀-free gate dims, the corrupt-dim
truth-table, and the back-solve-tautology table, and confirmed no 018/019/020 leakage, a genuine independent `.wl`, and the
COMPLETE≠PASS landing) + `ADVERSARIAL_ISSUES` (per-tooth ablation: the two KEY teeth genuine; 5 LOW-severity stamped/subsumed
teeth) → Codex remediated (2 make-genuine, 3 de-count) → fresh-agent `REVERIFY_CLEAN`. Tallies SymPy 42 / Mathematica 50,
CWD-independent.
