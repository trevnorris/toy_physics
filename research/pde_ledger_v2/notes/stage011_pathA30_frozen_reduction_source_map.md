# II-G1a (ledger_stage011) source map — pathA_30 frozen-reduction certificate

> Running-start prep captured 2026-07-08 (post stage010) from a focused pathA_30 investigation so a fresh session can
> author the reshape directive without re-discovery. Verify line refs against the sources before finalizing.
> Companion: `part2_gravity_atomic_split.md` (rows 011/012 + the Cluster-A reshape-cost bullet + the pathA_30 trip-ups).
> Build-order id **011**, Part II. **pathA_30 SPLITS into 011 (this = the frozen-reduction certificate) + 012 (the
> DtN pole ladder + Robin falsifier)** — same two-stage pattern as pathA_29 → 009/010.
> **Source top-line: `DN_UNITTEST_BC_DEPENDENT`** (the joint pathA_30 verdict; 011 carries the reduction-certificate
> component, 012 carries the D/N-ladder + the `BC_DEPENDENT` landing).

## ⭐ The headline difference from 009/010 (READ FIRST)
009/010 were **pure bridge-strips** (the physics was already computed, just de-plumbed). **Stage 011 is a bridge-strip
PLUS a de-rig.** Three 011 headline booleans in the source are currently **X≡X / literal tautologies** (the stage006
B2 rig class): the reshape must make them genuinely computed, i.e. *the reduction must be PRODUCED, not typed*:
1. **`operator_is_helmholtz`** (`.py` L537–539): compares two **identical hand-typed** Helmholtz expressions
   (`reduced_operator` L537 ≡ `ideal_operator` L538) → tautologically True. **Fix:** assemble `L_s` from the Helmholtz
   core + the Σ of computed intruding terms (each from `build_reduction_certificate`, asserted 0-or-deferred), so it
   can FAIL if a term survives. (Directive §4 firewall: "reduction first, comparison second; `L_s` is *produced* by
   the reduction, the idealized operator is the target.")
2. **`speed_is_cs = True`** literal (L540): **Fix:** compute the ψ-coefficient of the assembled `L_s`, equate to
   `(ω/c_S)²` with `c_S²=5Kρ*⁴/m`; show it is the bulk sound speed, not wall/healing-renormalized.
3. **`domain_is_L0 = expr_equal(L0−0, L0)`** (L541) = `L0==L0`: **Fix:** derive the cap endpoint from the regularity
   `R0(L0)=0` pinch-off, not `L0−0==L0`.
This is the sharpest fidelity finding and the reshape's central burden. The adversarial leg will hunt exactly this.

## File inventory
- **Report:** `software/stage1_solver/reports/pathA_30_dn_unit_test.md` (5.5KB) — 011 content: `## Reduced Operator`
  + `## Reduction Certificate` (report :3–24); background :19; `mutation_fires`/`fail_suppressed` :80.
- **`.py`:** `software/stage1_solver/tools/pathA_30_dn_unit_test_sympy.py` (39KB). **`.wl`:** `…/pathA_30_dn_unit_test.wl`
  (10KB). Directive `…/directives/pathA_30_dn_unit_test.md` (21KB); `pathA_30_results.yaml`.

## §1 The 011 / 012 split (`.py` line ranges)
**011 slice (frozen reduction → Helmholtz operator + c_S² + projection measure + certificate):**
- symbol setup + speed law + Helmholtz operator OBJECT: L408–425 (`k=ω/cS` L423; `cs_squared_from_eos=5Kρ*⁴/m`
  **L424**; `ode = ψ''(s)+k²ψ(s)=0` **L425** — 011 owns the operator; the `dsolve` at L426 opens 012 — clean cut at
  the `dsolve` call).
- `build_reduction_certificate` — the 011 CORE — **L333–404** (per-term fate ledger: projection measure, ρ0-grad,
  c_s-grad, δV_conf, ∇Q, BdG-k⁴ ratio + the `kξ≪1` deferral, `unsuppressed_operator_intrusion` L403).
- reduced-operator/speed/domain checks L537–541 (the three rigs above); the 011 payload keys L615–621.
- 011 verdict rungs in `compute_verdict` L719–738: `operator_is_helmholtz`→`FAIL_OPERATOR_INTRUSION` (L724,726–727);
  `speed_is_cs`→`FAIL_WRONG_SPEED` (L728–729); `domain_is_L0`→`FAIL_WRONG_DOMAIN` (L730–731).
- `build_dimensional_check` L200–330 is **SHARED**: the c_S²=5Kρ*⁴/m walk + `[K]=[P]−5[ρ]` chain (L217–219, 273–287)
  + the **corrupt-[K] probe** (L241–251) are **011**; the `tan_argument`/`Z00_prefactor` walk (L233–234, 244–245) are
  012. (⚠ split-note L97 files "the K-dim mutation probe" under 012's trip-ups — RESOLVE: c_S² leg → 011, tan-arg leg
  → 012.)

**012 slice (NOT this stage — draw the boundary):** `dsolve` L426 → D/N matrix/det/bc L435–447 → DtN
`−(ω/c_S)tan(L0ω/c_S)` L448–454 → pole ladder L456–468 → static series L470–474 → Robin block L476–529 → round-trip
L531–535 → provenance L542–543,571–578 → all `mma_exports` L580–587. Verdict rungs `FAIL_POLE_LADDER`,
`FAIL_COUNTERFACTUAL`, `BC_DEPENDENT` L732–737.

## §2 The 011 claim-set (derive + assert; report quotes)
- **Frozen-wall reduction:** `R=R0, η=0; matter perturbation stays dynamic e^{−iωt}` (report :5); background `straight
  finite throat, η=0, s∈[0,L0], R0(L0)=0, ρ0(s)=ρ*, A_M0=0` (:19).
- **Helmholtz operator:** `L_s ψ̂ = (c_S²ψ̂'' + ω²ψ̂)/c_S²` (:9) → `ψ''(s)+(ω/c_S)²ψ(s)=0` (:13).
- **`c_S²=5Kρ*⁴/m`** (`build_reduction_certificate` L343–345; dim-rewalked L217–219) — **THIS IS Part I edge R1**
  (see §4): CITE, don't re-derive; 011 only evaluates R1 at the throat reference density ρ*. `c_S` is the bulk sound
  speed, NOT wall/healing-renormalized (:15/:22).
- **Projection measure:** `√γ0(s)=A_perp0` constant → `d/ds(log√γ0)=0`, no first-derivative measure term (L344,354;
  report :20).
- **The "certificate"** = an emitted **per-term fate ledger** (NOT one assert-zero): each intruding term computed with
  its vanishing condition — ρ0'/ρ0=0, c_s'=0, δV_conf=0 (frozen η=0), ∇Q=0, BdG `k⁴`-term `ħ²k⁴/(4m²)` deferred only
  under `kξ≪1, ξ=ħ/(m c_s)` (L348–398; report :24). The reshape ASSEMBLES these into `L_s` (see the ⭐ header).

## §3 Reshape cost (the bridge to sever)
Scratch-yaml payload-mirror (Cluster-A variant). Strip: `write_sympy_exports` L661–670 (writes
`_scratch/pathA_30_sympy_results.yaml` + `_scratch/pathA_30_sympy_exprs.wl`); `load_engine_agreement` L673–716 (reads
the MMA scratch yaml + `digest_matches`/`engine_agreement`); `digest_mapping`/`yaml_read`/`yaml_write` L63–82;
`mma_exports`/`expression_digest` L580–588 (all 6 exported exprs are 012 — none are 011); the report/feed/yaml
writers L768–931. **Zero file I/O** in the reshaped stage scripts.
- **⭐ The `.wl` is ALREADY a GENUINE independent route** (transfer-matrix native: `transferMatrix` L25–28 →
  `dtnTransfer`/`robin*`/`staticSeries` L30–48; dim block rebuilt natively L104–135). It imports the `sympy*` values
  **only as cross-check RHS** in equality tests (L50–57), never as computation inputs. **Reshape = delete L18–19
  `Get[sympyExprFile]` + the L50–57 `sympy*`-comparison lines; keep the native transfer-matrix + dim derivations and
  print them.** (This corrects an initial "the `.wl` is a mirror" read — it is not.)
- **⚠ 011 `.wl` GAP:** the `.wl` has **no reduction-certificate content** (projection measure, ρ0/c_s/Q fates, BdG-k⁴
  are SymPy-only). For dual-engine, the 011 `.wl` must independently re-derive what it can: `c_S²=5Kρ*⁴/m` (already
  there L120), `d/ds(log√γ0)=0`, `d/ds(ρ*)=0`, the BdG ratio `ħ²k²/(4m²c_S²)` — else the certificate is single-engine.

## §4 Consumed / exported
- **Consumes (cite, don't re-derive):** `c_S²=5Kρ⁴/m` = **Part I edge R1** (register; home **stage005**
  `ledger_stage005_sound_speed_light_ratio`) — ⚠ caveat: the source evaluates it at the **throat centerline density
  ρ*** (not asymptotic ρ0) and names the symbol `c_S`; cite as "R1 evaluated at the straight-reference density ρ*".
  EOS exponent-5 closure stays `IMPOSED`. Also cites **stage004's `{L,T,M}` dictionary + the `[K]=[P]−5[ρ]` chain**
  (don't re-derive the dim foundation); `ħ, m` = 004 ACTION primitives; `ξ=ħ/(m c_s)` DERIVED (R2 family).
- **Exports** (split-note cross-stage flow): `011/012 export the frozen throat packet + D/N provenance + Helmholtz
  operator → 013 (β) + 017 (calibration input)`. 011's share = the frozen const-coeff Helmholtz operator `L_s`, the
  domain `[0,L0]` (cap `R0(L0)=0`), `c_S`, and the reduction certificate (validity window `{ρ0'/ρ0=0, kξ≪1,
  δV_conf=0, √γ0 const}`). (Distinct from 009/010's *bulk* Helmholtz mode → 024/026 — do NOT conflate the
  frozen-throat mode with the bulk mode.)

## §5 Teeth candidates
Keep/assign to 011: **projection-measure tooth** (mutate `A_perp0 → A_perp0(s)` so `d/ds(log√γ0)≠0` → induces
`f(s)d/ds` intrusion → `FAIL_OPERATOR_INTRUSION`); **BdG-k⁴ tooth** (drop `bdg_k4` WITHOUT the `kξ≪1` condition →
`FAIL_OPERATOR_INTRUSION`; wire `unsuppressed_operator_intrusion` from the fate ledger, currently hardwired False);
**K-dim mutation probe** (L241–251 + self-ablation L751–765: corrupt `[K]` by one power of L → c_S² dim flips →
`DN_UNITTEST_FAIL_DIMENSIONAL`; report :80 `mutation_fires:True`). PLUS the three de-rig fixes (§⭐) each become
genuine teeth (a surviving intrusion term / wrong speed / wrong domain now FAILs). ⚠ **`.wl` arity scan + unevaluated-
leakage transcript scan** (stage007/008/009 silent-skip lesson). Robin guard + tan-arg = 012, not 011.

## §6 Register expectation
**Zero new counted knobs.** Cites `K, ρ0(→ρ*), m, ħ` (Part I ACTION primitives, 004) + `c_S` (DERIVED via R1,
stage005). `L0` (throat depth) = POSTULATED straight-reference geometry — `ACTION`-geometry (like stage009's `d`,
"NOT a medium constant"); first enters at 011 (domain `[0,L0]`, cap `R0(L0)=0`); register once (011 or jointly
011/012). `ℓ_c` (confinement length) = INERT in the frozen η=0 test (`δV_conf=0`) → tracked, not counted.
`ξ=ħ/(m c_s)` DERIVED (R2 family), not new. **New edges:** no *reduction* edge (speed = R1, banked). Two
**structural** edge candidates (like R22/R25 — NOT reductions): (a) the frozen-reduction validity record (`L_s` =
const-coeff Helmholtz conditional on `{ρ0'/ρ0=0, kξ≪1, δV_conf=0, √γ0 const}`); (b) the **`ξ ≠ ℓ_c` firewall**
(healing length vs confinement length — directive §1 warns not to conflate; analogous to R22's `μ_R≠μ_R⁽⁴⁾`).

## Verdict tokens + honest scope
011 component of **`DN_UNITTEST_BC_DEPENDENT`**: the frozen reduction genuinely yields the const-coeff Helmholtz
resonator operator (EARNED, once the X≡X rigs are fixed) with a computed validity certificate; the BC-dependence and
the `BC_DEPENDENT`/banked-calibration landing are **012's** (the D/N derivation is `bc_derivation_emitted=False` — a
banked CALIBRATION input, earning it → PASS is an optional later upgrade). Nothing here selects the throat depth `L0`
(postulated geometry). CITED: R1 (speed law, stage005). Sim-deferred: none new (the deep throat solve stays off the
critical path).

## Process (unchanged, calibrated)
Author reshape directive (⭐ header de-rig + §1/§2 faithful cover + §3 bridge-strip + the `.wl` certificate-content
requirement) → Codex xhigh design-review → fold to `DIRECTIVE_CLEAN` (no GLM on Parts I–VI) → **pre-exec USER GATE**
→ Codex builds the two scripts (`--sandbox danger-full-access`, background, `< /dev/null`, absolute paths, xhigh) →
dual-engine exit 0 (repo root + foreign CWD) → arbiter re-run → tri-review (fidelity + adversarial-with-ablation,
incl. the X≡X-derig verification + arity scan + tally spot-check) → remediate → fresh-agent re-verify → registration
10→11 + parameter register + Codex-verify → note/card/`\input{stages/stage_011}` → PDF → commit + docs/memory sync.
Target stem: `ledger_stage011_frozen_reduction_certificate` (confirm slug at authoring).
