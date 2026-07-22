# III-2 / III-3 (ledger_stage030 + 031) source map — pathA_35 gateL couple-stress no-go

> **⚠ RECLASSIFIED 2026-07-22 → FAILURES-PAPER SEED (no longer a ledger source map).** Under the "ledger = surviving
> solution only" standing rule, the couple-stress no-go is a post-mortem of the retired `P` and was DROPPED from the
> ledger (stages 030/031 cancelled). This distillation is PRESERVED here as the seed for the future "approaches that
> didn't work" paper — see `ledger_exclusions_failures_paper_backlog.md` Exclusion 1. It is accurate source material;
> it is simply not ledger content. (The rest of this doc is the original running-start, kept intact for the paper.)

> Running-start prep for the two new Part-III stages, captured 2026-07-21 from a source distillation so a fresh
> session can author the reshape directives without re-discovery. Verify against the cited sources before finalizing
> each directive. Companion: `part3_light_atomic_split.md` (rows 030, 031). Build-order ids **030** (III-2, the
> live-`P` horn + gapped control) and **031** (III-3, the slaved escape + C5 + aggregate). This ONE map covers BOTH
> stages (the couple-stress no-go is split at the §2 seam).
> **Reshape of an EXISTING pair** (NOT fresh-authored), FULL tier (heaviest `pathA_35` pair):
> `software/stage1_solver/tools/pathA_35_gateL_sympy.py` (~1182 lines) + `pathA_35_gateL.wl` (~987 lines).
> **Headline verdict (verbatim, report line 11):** `FAIL_COUPLE_STRESS_NOGO`, provenance tag `CONDITIONAL_ON(both)`.
> **⚠ Decision-16 wrinkle (READ FIRST):** this no-go is the analysis of the **now-retired** polar field `P`. It
> references the **HISTORICAL** `P`-bearing `S_G0` tier (`T0_SHEAR_FROZEN(d9520d3819c3)`, T0 block `8fa41ac51e88`
> byte-embedded), NOT the post-D16 operative subset. Frame it as *"the historical `P`-sourcing attempt and why it
> failed"* — the **in-ledger home of Decision-16's "light" payoff-failure evidence** (see §9). This is the light-sector
> *reason* `P` was retired; NO new physics vs Part I.

## File inventory
- **Report (authoritative source):** `software/stage1_solver/reports/pathA_35_gateL_light.md` (~8.9 KB) — verdict
  `FAIL_COUPLE_STRESS_NOGO`; Section 2.6 frames the "four-way no-go"; the dimensional firewall PASS/ENGINE_AGREE.
  Machine ledger: `reports/pathA_35_gateL_light_results.yaml`.
- **Scripts to reshape:** `tools/pathA_35_gateL_sympy.py` (~1182 lines) + `tools/pathA_35_gateL.wl` (~987 lines)
  (see §5 for anatomy). Current agreement is via `argparse --compare` + JSON scratch payloads (`_scratch/
  pathA_35_gateL_{sympy,mathematica}.json`) — a **payload mirror**.
- **Provenance-only (already folded, do NOT re-derive):** `reports/pathA_35_G0_freeze.md` = the frozen constitutive
  action, folded into Part I `ledger_stage007`. The freeze report contains **NO gateL content** (scope line 9:
  "no Gate L verdict computed"). The couple-stress no-go lives ONLY in the gateL report/scripts.
- **Directive (methodology to carry):** `directives/pathA_35_shear_surface_brane_gates.md` (§1 binding methodology,
  the Gate-L protocol, §10 dimensional firewall).
- **Decision-16 doc (cite as the downstream landing):** `software/stage1_solver/decisions/16_retire_brane_polar_field.md`
  (evidence item 2 "Light" + item 4 "Active harm").

## §1 The physical question + the frozen ingredients consumed
**Question:** can light's transverse shear stiffness `μ_R` be **SOURCED from the T0 polar field `P`** acting as a
Cosserat/couple-stress micro-rotation reservoir, instead of being a bare postulated modulus? (The frozen action
explicitly reuses `P` for this: "The Cosserat micro-rotation reservoir reuses the T0 polar field `Pⁱ`. No independent
micro-rotation field is added.") It is NOT a test of whether light exists — `L_Mac` with `μ_R` postulated already
carries the transverse wave (that's stage003). It tests whether the `P`–`u` coupling `L_Pu` can supply the surface
traction so `μ_R` is *derived provenance*, WITHOUT (a) leaking hidden propagating `P` modes, (b) unbounded Hamiltonian,
(c) broken angular-momentum closure, (d) a longitudinal zero mode.

**Frozen ingredients consumed (from `S_G0`, historical tier):**
```
L_Mac = ½ρ_br(∂_t uᵃ)² − ½μ_R Ω_uᵃΩ_uᵃ ,   Ω_uᵃ := (∇_∥×u)ᵃ ,   μ_R>0 postulated modulus
L_Pu  = −λ_Pu ϖ_a Ω_uᵃ = −λ_Pu ŵ·(P_∥×(∇_∥×u)) ,   ϖ_a := (ŵ×P_∥)_a   (parity-EVEN; direct P_∥·Ω_u is parity-ODD, excluded)
L_uw  = ½ρ_br(∂_t u_w)² − ½ρ_br Ω_w² u_w² ,   Ω_w>0 bare gap
S_brane = ∫dt d⁴X g_ℓ(w)[L_Mac + L_Pu + L_uw] ,   g_ℓ(w) = exp(−(w/ℓ_g)²)/(√π ℓ_g)
```
**DERIVED in-gate (re-postulated, NOT imported)** — surface projection of the inherited T0 `P` inertia/stiffness through
`∫dw g_ℓ = 1`, mode weight `∫dw ℓ_g g_ℓ = ℓ_g`:
```
J_P     = m ρ a² ℓ_g            (P surface inertia)
Γ_P     = m ρ c_s² a² ℓ_g       (P surface Frank stiffness)
M_gap_P = J_P Ω_P²              (gapped-control mass)
```
**Load-bearing object** — the transverse `(u_T, ϖ_T)` Hessian (second variation of the derived quadratic energy —
NOT typed):
```
H = [[ μ_R k² ,  λ_Pu k ],
     [ λ_Pu k ,  Γ_P k² ]]
```
**IMPORTED, not earned here:** `c_γ² = μ_R/ρ_br` — a definitional identity of the frozen `L_Mac` (transverse
`ρ_br ω² = μ_R k²`), earned in stage003; gateL only re-records + dim-re-checks it as the "MacCullagh feed-forward
speed for later gates." **Do NOT present `c_γ²` as a gateL-earned result.**

## §2 The three-horn structure + the split seam (030 | 031)
The report is a "four-way no-go" (Section 2.6): three `P`-sourcing branches + a shared C5/`φ` horn. The stage split:

- **`stage030` (III-2) = the live-`P` horn + the gapped control** (branches i + ii).
- **`stage031` (III-3) = the slaved-rigid escape + the C5/`φ` obstruction + the four-way aggregation** (branch iii +
  the shared horn + the `FAIL_COUPLE_STRESS_NOGO` roll-up + the `FREE_LIGHT_OK_CONDITIONAL` positive control).

**Seam rationale:** the live-`P` instability and the slaved escape are genuinely different derivations reaching
**different config verdicts** (`FAIL_COUPLE_STRESS_NOGO` for Config A vs `FAIL_C5_LONGITUDINAL_ZERO_MODE` for Config B).
The L(c) Magnus-leak and L(d) `u_w`-gap sub-hurdles PASS in both configs (orthogonal guardrails) and ride with 031.
**Do NOT split per-sub-hurdle** — L(a-i)…L(d) are cross-cutting across both configs and the verdict is inherently the
*aggregation*; splitting per-hurdle fragments the linkage logic that is the whole point.

## §3 Derivation content per stage

### stage030 (III-2) — the live-`P` horn (Config A) + gapped control
**(a-i) Traction provenance (Config A `ARROWS_SUPPLY_TRACTION`).** Virtual work from `L_Pu` gives a **rank-2 cut-traction
matrix** (MacCullagh-like) — establishing live `P` genuinely IS the traction reservoir under test (report:64,
`.py:608–664`); Frank-only reference (rank-0) is the discriminator (`FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION`). This
positive result is load-bearing (without it the horn shows `P` harmful without establishing it was the candidate
reservoir). *(The Config-B downgrade of L(a-i) to `POSTULATED_SURFACE_ELASTICITY` is stage031's.)*

**(i) Live/massless `P` (Config A baseline: massless `Pⁱ`, parity-repaired `P–u`, no `φ`, gapped `u_w`).**
- 7×7 principal symbol (degree 7 in `ω²`); at `λ_Pu=0` factorizes:
  `(−ω²ρ_br)(μ_R k² − ω²ρ_br)²(Γ_P k² − J_P ω²)³(ρ_br(Ω_w² − ω²))`.
- **Mechanism 1 — hidden modes:** 5 gapless positive branches at `λ_Pu=0` ⇒ **3 extra `P` spin waves**; nonzero
  `λ_Pu` gives 2 low-`k` negative `ω²` branches. → `FAIL_HIDDEN_PROPAGATING_MODE`.
- **Mechanism 2 — unbounded energy:** transverse Hamiltonian principal minor `k²(Γ_P μ_R k² − λ_Pu²)`, **negative on
  `0 < k² < λ_Pu²/(Γ_P μ_R)`** ("an energy-matrix failure, not a dispersion-only claim"). `λ_Pu→0` restores the minor
  to `Γ_P μ_R k⁴` (bounded) — proving the coupling itself sources the instability. → `FAIL_NOT_BOUNDED_BELOW`.
- Config A verdict: `FAIL_COUPLE_STRESS_NOGO`.

**(ii) Gapped `P` (a control leg, target L(b)).**
- Derived gapped response `P = −λ_Pu Ω_u/(M_gap_P + Γ_P k²)`; couple-stress divergence `Γ_P k² P`.
- Low-`k` the divergence → 0 and the **closure residual remains `2 μ_R Ω_u`** — angular-momentum closure lost (the
  reservoir can't absorb the antisymmetric stress). → `FAIL_GYROSTAT_NO_CLOSURE` (same residual as the
  `omit_couple_stress_reservoir` control).

**stage030 earned residue (carry forward):** `c_γ² = μ_R/ρ_br` re-recorded; `μ_R>0` an independent postulated modulus.

### stage031 (III-3) — slaved-rigid escape + C5 obstruction + aggregation
**(iii) Slaved-rigid `P` (Config B: `P_∥ = ŵ×Ω_u`, no independent `P` spin modes).**
- Substitution ⇒ `ϖ = ŵ×P_∥ = −Ω_u`; retaining the bilinear ⇒ `μ_eff = μ_R − 2λ_Pu`,
  `ω² = ((μ_R−2λ_Pu)k² + Γ_P k⁴)/ρ_br`. Slaved determinant (degree 4):
  `(−ω²ρ_br)(Γ_P k⁴ + k²(−2λ_Pu + μ_R) − ω²ρ_br)²(ρ_br(Ω_w² − ω²))`.
- **This branch ESCAPES the live-`P` horn:** L(a-ii) PASS (0 extra gapless `P` branches); L(b) PASS conditional on
  `μ_R − 2λ_Pu > 0`; closure residual `0`.
- **But it fails on provenance + inherited C5:** (a) `k⁴` dispersion (exact nondispersive equality false at finite `k`;
  clean only for `k² ≪ (μ_R−2λ_Pu)/Γ_P`); (b) L(a-i) traction provenance drops `ARROWS_SUPPLY_TRACTION →
  POSTULATED_SURFACE_ELASTICITY` (the clean `k²` light traction IS the postulated surface modulus — `μ_R` back to a
  bare knob, defeating the point); (c) L(a-iii) still `FAIL_C5_LONGITUDINAL_ZERO_MODE`. Config B verdict:
  `FAIL_C5_LONGITUDINAL_ZERO_MODE`.

**Shared 4th horn — the C5/`φ` obstruction (kills BOTH configs, orthogonal to couple-stress):**
- No scalar-potential `φ` frozen ⇒ kinetic longitudinal zero mode remains (`FAIL_C5_LONGITUDINAL_ZERO_MODE`).
- `φ = u_w` rescue "collides with the required u_w gap" (the descendant would need to be massless while L(d) requires
  `Ω_w>0`). An independent variational `φ` "is a fresh G0, not a Gate-L pass" (`FIRED_PASS_FIXTURE_ONLY`).

**The aggregation (`.py:938–949`) — represent verbatim, NOT a collapse onto one hurdle:**
```
FAIL_COUPLE_STRESS_NOGO  fires iff
  FAIL_HIDDEN_PROPAGATING_MODE  AND  (FAIL_NOT_BOUNDED_BELOW OR FAIL_GYROSTAT_NO_CLOSURE)
  AND  FAIL_C5_LONGITUDINAL_ZERO_MODE  AND  gapped-control-fires
```
Report: "the overall label is `FAIL_COUPLE_STRESS_NOGO`, not a collapse onto L(b) alone."

**Guardrail PASSes (ride with 031):** L(c) leak — flat wave `T_na=0` (bent interface ⇒ `FAIL_LEAK_BREAKS_MAGNUS`);
L(d) `u_w` gap — `Ω_w>0` protects (ungapped ⇒ `FAIL_BENDING_MASSLESS_FIFTH_FORCE`).

**stage031 earned residue:** `μ_R` = an honestly-postulated modulus (the no-go rules out only *deriving* it from `P`;
light stands on the bare modulus, with `pathA_36`/stage003 getting photons `P`-free).

## §4 Able-to-fail teeth / controls (the `controls_summary`, `.py:998–1013` / `.wl:906–920`)
| control mutation | target | fires → | stage |
|---|---|---|---|
| Frank-only reference (rank-0 `u` traction) | L(a-i) | `FAIL_FRANK_TORQUE_NOT_MACCULLAGH_TRACTION` | **030** (the discriminator for Config-A `ARROWS_SUPPLY_TRACTION`; Grok MINOR-3 — this row supersedes the earlier "031") |
| Cauchy reference (3 propagating elastic modes) | L(a-ii) | `FAIL_CAUCHY_STRAY_LONGITUDINAL` | 030 |
| no-`φ` C5 branch | L(a-iii) | `FAIL_C5_LONGITUDINAL_ZERO_MODE` | 031 |
| independent variational `φ` fixture (removes zero mode) | L(a-iii) | `FIRED_PASS_FIXTURE_ONLY` (fresh-G0, not a pass) | 031 |
| raw `div u=0` projector (no variational provenance) | L(a-iii) | `FAIL_C5_LONGITUDINAL_ZERO_MODE` | 031 |
| omit couple-stress reservoir | L(b) | `FAIL_GYROSTAT_NO_CLOSURE` (residual `2 μ_R Ω_u`) | 030 |
| large-gap `P` leg | L(b) | `FAIL_GYROSTAT_NO_CLOSURE` (residual `2 μ_R Ω_u`) | 030 |
| bent interface (nonzero direct `T_na`) | L(c) direct | `FAIL_LEAK_BREAKS_MAGNUS` (flat ⇒ `T_na=0`) | 031 |
| Frank-only indirect channel (zero induced `P`) | L(c) indirect | `FIRED_ZERO_SOURCE` | 031 |
| nonplanar indirect (nonzero advective curl) | L(c) indirect | `FIRED_NONZERO_SOURCE` | 031 |
| ungapped `u_w` (`Ω_w→0`) | L(d) | `FAIL_BENDING_MASSLESS_FIFTH_FORCE` | 031 |
| `λ_Pu→0` (Hessian sanity) | L(b) | minor restores to `Γ_P μ_R k⁴`, bounded | 030 |
| **`good_structure_fixture` (all-sub-hurdles-pass)** | aggregator | **`FREE_LIGHT_OK_CONDITIONAL`** — the able-to-**PASS** tooth, proving the verdict machinery is NOT hardwired to fail | 031 |

**⭐ The `FREE_LIGHT_OK_CONDITIONAL` able-to-pass fixture is load-bearing** (anti-tautology positive control) and MUST be
preserved through the reshape — the adversarial leg will ablate it.

## §5 The script pair — FULL reshape (what to strip/keep)
- **STRIP:** `argparse --compare` (`.py:1164–1166`); the JSON scratch payload writes
  (`_scratch/pathA_35_gateL_{sympy,mathematica}.json`); the `--compare` byte-equality `agreement_payload` diff; the YAML
  ledger write (`reports/pathA_35_gateL_light_results.yaml`). → two standalone print-only / assert-zero engines,
  independent in-process, both exit 0.
- **⚠ The `.wl` is a payload-mirror AND hardcodes branch counts the `.py` genuinely computes** via `sp.roots`:
  `gaplessMasslessLam0=5`, `gaplessGappedLam0=2`, `gaplessSlaved=2`, `positiveMasslessLam0=6`, … (`.wl:403–411`), then
  only `assertZero`-checks the determinant factorizations (`.wl:386–401`). **The independent `.wl` re-author must
  re-derive those counts from its OWN root/rank analysis** (native `Solve`/`Roots`/`Eigenvalues`/`CountRoots`), NOT
  literal-stamp them. This is the sharpest reshape tooth (worse mirror than stage007's).
- **KEEP (the load-bearing checks):**
  - Freeze-fidelity: re-extract the `freeze-action` fence by parsing (not fixed offset) + SHA-256 == `d9520d3819c3…`
    (+ T0 `8fa41ac51e88…` byte-embedded) + the 5 required-line needles (`.py:110–119`).
  - `build_dimensions`: ~42 `check.check` + 2 `expect_fail` ablations (matches the report's "42 expressions").
  - `build_derived_quantities`: ~25 raise-on-mismatch checks (`g_ℓ` norm, mode weight, `J_P`, `Γ_P`, radial,
    Hessian-match, minor, `λ_Pu→0` sanity, low-`k` sign, `μ_eff`, `k⁴`, closure residual; degrees `(7,7,4)`; gapless
    `(5,2,2)`; positive `(6,6,3)`; zero `(1,1,1)` — **⚠ the THIRD entry of each tuple (degree 4 / gapless 2 / positive 3
    / zero 1) is the SLAVED config = stage031; stage030 owns the live+gapped PAIRS `(7,7)`/`(5,2)`/`(6,6)`/`(1,1)`
    only** (Codex r1 BLOCKING-1); low-`k` count 3; gapped closure `2 μ_R Ω_u`; direct `T_na`
    flat/bent; indirect flat curl / Frank-only; nonplanar curl).
  - Six sub-hurdle audits × 2 configs (rank/det raise-guards), the aggregation, the `good_structure_fixture`
    (able-to-pass), the `hypothetical_pass_fixture`, and the 13-row `controls_summary`.
  - **Total ≈ 80+ able-to-fail assertions** (≈2× the stage007 pair). Split across 030/031 at the §2 seam.

## §6 Dimensional targets (units-restored, able-to-fail; report firewall PASS/ENGINE_AGREE, 42 expressions)
`Pⁱ` dimensionless · `uᵃ, u_w` = L · `ρ_br` = M L⁻³ · `μ_R, λ_Pu` = M L⁻¹T⁻² · `Ω_w, Ω_P` = T⁻¹ · `a` = L · `c_s` = LT⁻¹ ·
`g_ℓ` = L⁻¹ · `ℓ_g` = L · `J_P` = M L⁻¹ (surface inertia) · `Γ_P` = M LT⁻² (surface Frank stiffness) ·
`c_γ² = μ_R/ρ_br` = L²T⁻² · `T_wa = mρ v_w v_a` = M L⁻²T⁻².
**Able-to-fail dim ablations (both FIRED):** `drop_m_from_T_wa` (`ρ v_w v_a` as `L⁻²T⁻²` not `M L⁻²T⁻²`);
`P_u_without_curl_or_cut_gradient` (`λ_Pu P u` as `M T⁻²` not `M L⁻¹T⁻²`). *(Confirm `J_P`/`Γ_P` dim triples against the
script at authoring — restated here from the surface-projection integrals.)*

## §7 Inputs / provenance (strings, not runtime deps)
gateL CONSUMES the whole G0 freeze (hash-verified: G0 `d9520d3819c3` + T0 `8fa41ac51e88` byte-embedded, 5 needles).
Depends on: the frozen `S_G0` action (`L_Mac`, `L_Pu`, `L_uw`); the `P` field definition + its T0 coefficients
(`mρa²`, `mρc_s²a²`, `mρc_s²`); `{ρ_br, μ_R, λ_Pu, Ω_w, g_ℓ/ℓ_g, a, c_s, Ω_P}`. **Historical tier** (P-bearing), per
the Decision-16 wrinkle.

## §8 Reshape trip-ups (pin in each directive)
1. **Independent `.wl` must re-derive the branch/mode counts** (§5) — the single sharpest tooth; the source `.wl`
   literal-stamps `5/2/2/6`. A faithful independent route computes them from its own `Roots`/`Eigenvalues`.
2. **Payload-mirror → two standalone engines** — strip `--compare`/JSON/YAML; each asserts in-process, both exit 0.
3. **Hash fidelity by fence-parsing**, not fixed byte offsets (revision-sensitive); or embed the frozen bytes.
4. **Preserve the `FREE_LIGHT_OK_CONDITIONAL` able-to-pass fixture** (§4) — the anti-tautology positive control; the
   aggregator must NOT be hardwired to fail.
5. **`c_γ²=μ_R/ρ_br` is IMPORTED, not earned** (§1) — re-recorded/dim-checked only; do not claim it as a gateL result.
6. **Historical-tier framing** (Decision-16 wrinkle) — reference the `P`-bearing `S_G0`, present as the retired-`P`
   sourcing attempt; cross-ref Decision 16 + stage007's amendment. NO new physics vs Part I.
7. **The verdict is an AGGREGATION** (`.py:938–949`), not a collapse onto one hurdle — 031 owns the roll-up logic;
   keep the linkage `AND`/`OR` structure intact and able-to-fail.
8. **`μ_R` notational collision** (as in stage007): the 3D brane modulus `μ_R` (M L⁻¹T⁻²) ≠ stage006's 4D
   `μ_R⁽⁴⁾` (M L⁻²T⁻²); register edge R17 (PENDING). Do not conflate.

## §9 Decision-16 linkage (the downstream landing — cite in both stage notes)
Decision 16 cites the couple-stress no-go as one of `P`'s failed payoffs, twice:
- Evidence item **2. Light** (verbatim): "Light: pathA_35 `FAIL_COUPLE_STRESS_NOGO` (light's stiffness cannot be
  derived from `P` substructure); pathA_36 derives both photons with no `P` sector present."
- Evidence item **4. Active harm**: cites the pre-registered G0.3 `FAIL_NOT_BOUNDED_BELOW` exposure firing —
  `det = k²(2A_P μ_R k − λ_Pu²) < 0 for k < λ_Pu²/(2A_P μ_R)` — the same transverse Hessian minor as gateL L(b) (later
  confirmed structural by U1 Phase A: `INSTABILITY_CONFIRMED_STRUCTURAL`).
So `stage030`/`stage031` are the in-ledger home of the light-payoff-failure evidence for the P-retirement we landed in
Part I (stage006/007 amendment). Frame accordingly: the brane DOF budget was already spoken for (inflow=gravity,
shear=light, swirl=gravitomagnetism); every computed `P` payoff failed independently, and the couple-stress no-go is
the one showing light's stiffness could not be derived from `P`.

## §10 Downstream consumers / parameter-register additions
- **→ Part V (magnetism):** `c_γ² = μ_R/ρ_br` + the transverse `u_T` sector forward (`S_{T+move}` reuses this kinetic +
  curl); stated in 031's downstream-consumers section.
- **stage003 "Next step" staleness fix (owed, note-only):** amend stage003's "Next step" to point at 030/031 and drop
  the superseded `c_L²=B_eff/ρ_br` cone-pair (model_map §6). No executable-core change, no re-review.
- **Parameter register:** `μ_R`, `ρ_br`, `λ_Pu`, `Ω_w`, `g_ℓ`/`ℓ_g` rows already exist (stage007). **030/031 ADD:**
  `J_P` [M L⁻¹], `Γ_P` [M LT⁻²], `M_gap_P`, `Ω_P` [T⁻¹] as DERIVED-in-gate manifestations (surface projections of the
  T0 `P` coefficients — track, do NOT count as new independent knobs); a structural edge recording the couple-stress
  no-go (`μ_R` NOT reducible via `P` — a CLOSED-NEG reduction route: the `P`-sourcing of `μ_R` is closed negative,
  reinforcing `μ_R` as postulated). Confirm dim triples for `J_P`/`Γ_P` against the script at authoring.

## Verdict tokens + honest scope
Headline `FAIL_COUPLE_STRESS_NOGO` (`CONDITIONAL_ON(both)`). EARNED (carried): `c_γ²=μ_R/ρ_br` feed-forward (imported),
`μ_R>0` a legitimate postulated modulus, `POSTULATED_SURFACE_ELASTICITY` (Config B). CHARACTERIZED-DEPARTURE: `μ_R`
cannot be micro-sourced from `P` (hidden modes / unbounded / no closure / provenance downgrade / C5 zero mode). The
no-go rules out only *deriving* `μ_R` from `P`; light stands on the bare postulated modulus (`pathA_36`/stage003 gets
photons `P`-free). This is Decision-16's light-payoff-failure evidence, formalized in-ledger.

## Process (reshape variant — same calibrated per-stage flow, ×2 for 030 then 031)
Per stage: author the reshape directive (contract + the §8 pins + faithful cover of §3/§5 blocks) → Codex xhigh
design-review → fold to `DIRECTIVE_CLEAN` → **Grok-4.5 headless design-review** → assess/fold → Codex confirm-pass
(no GLM on Parts I–VI) → **pre-exec USER GATE** → Codex builds the two scripts (`--sandbox danger-full-access`,
background, `< /dev/null`, xhigh) → dual-engine both exit 0 → arbiter re-run → full tri-review on fresh agents (arbiter
+ fidelity + adversarial-scoped-to-reshape-integrity, per-tooth ablation) → remediate → bump counts → update +
Codex-verify the parameter register → note/card/`\input{stages/stage_030}`(then 031)/registration → PDF → commit +
docs/memory sync. Target stems (confirm on authoring): `ledger_stage030_couple_stress_live_p_horn`,
`ledger_stage031_couple_stress_nogo_aggregate`. Build 030 first (031 consumes its live-`P`/gapped results as
provenance for the aggregation).
