# pathA_43 ↔ 2.5PN match-back — consolidation note (Phase A3)

> **Status: CONSOLIDATION NOTE, not a new derivation.** 2026-07-06. This note records a consistency
> match that was *already computed* (dual-engine) inside `pathA_43`; it introduces no new math and no new
> verdict. It DOES add a dedicated dual-engine verification artifact for the match-back itself (§2a), so the
> claim is a self-contained, able-to-fail unit rather than a pointer to another gate's overlay. It is the
> last step of **Phase A** (finish gravity) of the ledger-v2 rebuild (`notes/ledger_v2_rebuild_plan.md` §4,
> step A3). It will be carried into v2 **Part II (Gravity)** during Phase B; until then it lives here.

> **Citation shorthand** — every `path:line` below is either a full repo-root path or one of these declared
> abbreviations (line numbers index into the named file):
> - `_sympy.py` = `software/stage1_solver/tools/pathA_43_density_quadrupole_port_sympy.py`
> - `.wl` = `software/stage1_solver/tools/pathA_43_density_quadrupole_port.wl`
> - `results.yaml` = `software/stage1_solver/reports/pathA_43_density_quadrupole_port_results.yaml`
> - `report.md` = `software/stage1_solver/reports/pathA_43_density_quadrupole_port.md`
> - `directive.md` = `software/stage1_solver/directives/pathA_43_density_quadrupole_port.md`
> - `4d_2_5pn.tex` = `research/4d_2_5pn/paper/4d_2_5pn.tex`
> - `matchback.py` = `software/stage1_solver/tools/pathA_2_5pn_matchback_sympy.py`  *(this note's own artifact, §2a)*
> - `matchback.wl` = `software/stage1_solver/tools/pathA_2_5pn_matchback.wl`

---

## 0. What this records (one paragraph)

The separate audited PN corpus `research/4d_2_5pn` proved a **conditional** 2.5PN theorem: the full
point-particle Burke–Thorne radiation-reaction sector follows **if** a single scalar normalization holds on
the moving-throat branch — `m̂₀²Γ₅ = 2G/(5c⁵)`. `pathA_43` (Phase A1, `DENSITY_PORT_HOSTED`) built the ℓ=2
quadrupole radiative port **natively on the density/`c_s` mode** (retiring the old `A_w`/`U,W` vector
scaffold). Its closure overlay carries the calibrated quadrupole moment coefficients and **computes,
dual-engine, that they satisfy exactly that normalization** — `Γ̄₅ − 2G/(5c⁵) = 0` — together with the
moment consistency relation `K̄₄ − 4K̄₂²/K̄₀ = 0`. So the density-mode port meets `research/4d_2_5pn`'s single
open item **at reduced-closure**. The magnitude is **calibrated** (`G` is a genuine gap); this is a
consistency check, not a first-principles derivation of the Burke–Thorne coefficient (see §4).

---

## 1. The `research/4d_2_5pn` open item this addresses

`research/4d_2_5pn` is titled *"A Conditional Derivation of the Point-Particle 2.5PN Sector from the Unified
4D Toy Model"* (`4d_2_5pn.tex:26`). Its own statement of the remaining problem:

- `4d_2_5pn.tex:57-60` — *"The remaining serious open problem is no longer generic dissipative physics, but
  the final passive/outgoing quadrupole normalization, `m̂₀²Γ₅ = 2G/(5c⁵)`, on the actual moving-throat
  branch."*
- Boxed target, `4d_2_5pn.tex:819-824` — *"Matching to Burke–Thorne is therefore equivalent to the scalar
  equality"* `m̂₀²Γ₅ = 2G/(5c⁵)` (called *"the sharpest normalization target used in the present paper,"*
  `4d_2_5pn.tex:825`).
- The corpus also states the same content as the moment-invariant pair, `4d_2_5pn.tex:469-473`:
  `K̄₄ K̄₀ = 4 K̄₂²` and `Γ̄₅ = 9 K̄₂^{5/2}/K̄₀^{3/2}`.
- Headline landing, `4d_2_5pn.tex:333-355` / `4d_2_5pn.tex:491` — *"the 2.5PN program is now a normalization
  problem rather than a whole-model rescue problem."*

So the corpus reduced 2.5PN to one scalar normalization and left it open. That is the item this note records
`pathA_43` matching at reduced-closure.

## 2. What `pathA_43` computes (dual-engine, on record)

`pathA_43`'s closure overlay (SymPy `closure_overlay` `_sympy.py:701-722`; Mathematica `closureOverlay`
`.wl:353-364`) carries these **calibrated** ℓ=2 moment coefficients (hardcoded closure inputs — see §4 for
why they are calibrated, not derived):

```
K̄₀ = 54 G c_s⁵ / (5 a⁵ c⁵)      (= m̂₀² P₀)
K̄₂ =  6 G c_s³ / (5 a³ c⁵)
K̄₄ =  8 G c_s  / (15 a c⁵)
```
`_sympy.py:705-707`; `.wl:355-357`.

From these it forms `Γ̄₅ ≡ K̄₀ · a⁵ / (27 c_s⁵)` (equivalently `m̂₀² P₀ · a⁵/(27 c_s⁵)`, since `K̄₀ = m̂₀² P₀`;
the `27` is the outgoing ℓ=2 density-Hankel fingerprint, **earned** in pathA_43) and checks **both** residuals
to zero, in **both engines independently**:

| Invariant | SymPy (`sp.simplify`) | Mathematica (`FullSimplify`) | Emitted `= 0` |
|---|---|---|---|
| `K̄₄ − 4K̄₂²/K̄₀ = 0` | `_sympy.py:708` | `.wl:358` | `results.yaml:192` |
| `Γ̄₅ − 2G/(5c⁵) = 0` | `_sympy.py:709-710` | `.wl:359-360` | `results.yaml:193` |

Neither residual is hardcoded to `0`: each is built as a symbolic expression and reduced
(`sp.factor∘sp.cancel∘sp.simplify`, resp. `FullSimplify`). The gate requires **both** residuals to vanish —
`closure["ok"] = (K̄₄_residual == 0) ∧ (Γ̄₅_residual == 0)` (`_sympy.py:721`; `.wl:363`) — which is a
necessary condition of the `DENSITY_PORT_HOSTED` verdict (`_sympy.py:631`; `.wl:423`); each engine's
self-test asserts that verdict (`_sympy.py:775`; `.wl:466`, so `closure["ok"]` — both residuals — is
entailed) and additionally asserts the `K̄₄` residual directly (`_sympy.py:783`; `.wl:473`). Dual-engine
independence is on record — `results.yaml:28-30` (`sympy_consumes_mathematica: false`,
`mathematica_consumes_sympy: false`), `results.yaml:31-32` (`agreement.verdict: true`), and `exit_status: 0`
for both engines (`results.yaml:22`, `results.yaml:27`); report prose `report.md:11-14` (*"neither engine
consumes the other's numerator or booleans"*). Closure/residual prose: `report.md:78-80`.

**Reading (and its limit).** pathA_43's density-mode closure coefficients hit the exact Burke–Thorne
normalization `research/4d_2_5pn` singled out, and both engines verify it symbolically rather than by
assertion. But be precise about what the residuals *test*: given the hardcoded `K̄₀ = 54Gc_s⁵/(5a⁵c⁵)`, the
check `Γ̄₅ − 2G/(5c⁵) = 0` reduces to the arithmetic identity `54/27 = 2`, and `K̄₄ − 4K̄₂²/K̄₀ = 0` to
`8/15 = 8/15`. They confirm the calibrated-plus-earned coefficients are **mutually consistent** and land on
Burke–Thorne — **not** that any physics was independently derived *here*. The physics content is upstream:
the `54/5` is calibrated and the `27` is the earned density-Hankel fingerprint (§3); these residuals verify
those pieces are arithmetically arranged to reproduce `2G/(5c⁵)`.

## 2a. This note's own dual-engine match-back artifact

pathA_43's closure overlay only checks its own form and the `K̄₄` invariant. To make *this* match-back a
self-contained, committed, able-to-fail artifact (rather than leaning on another gate's overlay plus
in-context arithmetic), it has its own independent dual-engine check: `matchback.py` (SymPy) + `matchback.wl`
(Mathematica). Each engine builds and reduces the residual set below to exact `0` by its own route
(`matchback.py:157-188`; `matchback.wl:136-164`) — both exit 0, agree per residual, and emit distinct
canonical surface forms for the nonzero mutation residuals (so the `.wl` is a genuine independent route, not
a payload mirror):

- **INV1** `K̄₄·K̄₀ − 4·K̄₂² = 0` (corpus moment invariant).
- **INV2** `K̄₀·a⁵/(27·c_s⁵) − 2G/(5c⁵) = 0` (pathA_43 form → Burke–Thorne).
- **INV3** `9·K̄₂^{5/2}/K̄₀^{3/2} − 2G/(5c⁵) = 0` (the corpus's *own* form → Burke–Thorne — the tie-out
  pathA_43 does **not** carry).
- **INV4** `K̄₀·a⁵/(27·c_s⁵) − 9·K̄₂^{5/2}/K̄₀^{3/2} = 0` (redundant cross-form agreement, built directly).
- **INV5** literal anchors pinning `K̄₀,K̄₂,K̄₄,BT` to `54/5, 6/5, 8/15, 2/5` and the structural literals
  `27, 9, 5/2, 3/2` (`matchback.py:171-178`; `matchback.wl:150-157`) — independently-restated literals, so
  they defeat a coherent-rescale rig.

**Able-to-fail (the anti-rig backstop):** an 11-mutation probe feeds data-only corruptions through the SAME
computation path and asserts an expected caught-by matrix (`matchback.py:217-235`; `matchback.wl:186-208`).
Every mutation is caught; notably the **coherent-scale** `(K̄₀,K̄₂,K̄₄,BT)×2` — which passes INV1–INV4 — is
caught **only** by the INV5 anchors, and the **coupled** `K̄₀,K̄₂,K̄₄×2` (BT fixed) leaves INV4 zero while
INV2/INV3/anchors fire. Exact arithmetic with a no-float guard (`matchback.py:151-154`); honest provenance
block (moments CALIBRATED, `G = GENUINE_BLOCKED`, `27` upstream-earned, corpus form imported —
`matchback.py:199-205`); no runtime cross-consumption of any peer/report/source/note/`_scratch` file.

Verification trail: Codex design-review (caught a real coherent-rescale blind spot → INV5 anchors added) →
GLM-5.2 tertiary (`SOUND`; could not construct a passing rig) → Codex re-green `DIRECTIVE_CLEAN` → both
engines exit 0 → orchestrator arbiter re-run reproduced both → fidelity audit `FIDELITY_CLEAN`.

## 3. The `54/5 = 2 · 27/5` split (why the match is partly earned)

The magnitude that enters the normalization is `m̂₀² P₀ = 54 G c_s⁵/(5 a⁵ c⁵)`, i.e. the coefficient `54/5`.
Per the pathA_43 directive (`directive.md:67-68`):

> `54/5 = 2·27/5` — the **27 earned** from the outgoing ℓ=2 fingerprint (`1/9, 4/81, 1/27`, derived in
> pathA_43), the **2/5 · G calibrated**.

So the `27` (the density-Hankel tail that makes `Γ̄₅ = K̄₀ a⁵/(27 c_s⁵)` land on `2G/(5c⁵)`) is genuinely
derived on the density mode; the `2/5` Burke–Thorne factor and `G` are the calibrated / external inputs. (The
literal string `54/5 = 2·27/5` appears in the **directive**, not the reports. The explicit
`Γ̄₅ = K̄₀ a⁵/(27 c_s⁵)` formula (with the earned `27`) lives in the **tools** (`_sympy.py:709`; `.wl:359`);
the **reports** carry the resulting zero residual `Γ̄₅ − 2G/(5c⁵) = 0` (`results.yaml:193`).)

## 4. Honest scope — what this is and is NOT

**IS:** a dual-engine consistency check that pathA_43's density-mode ℓ=2 port, at its calibrated closure,
reproduces the 2.5PN Burke–Thorne normalization `research/4d_2_5pn` left open (`Γ̄₅ = 2G/(5c⁵)`) and the
moment-invariant `K̄₄ = 4K̄₂²/K̄₀`. This is the *cheapest decisive falsifier* for the density-mode port: had
the calibrated moments been mutually inconsistent, or inconsistent with Burke–Thorne, these residuals would
be nonzero and the gate would fail.

**IS NOT** a first-principles derivation of `Γ̄₅`. The moment coefficients `K̄₀, K̄₂, K̄₄` are **hardcoded
calibrated closure-overlay inputs** (`_sympy.py:705-707`), not solved from the port numerator `N0_den`. The
magnitude stays **CALIBRATED**: `results.yaml:171-187` (`calibrated_split`) tags
`CALIBRATED: [G, "2/5", "54/5"]` and `SIM_DEFERRED: ["literal Y2* l=2 moment magnitude", "Xi_Q exact branch
magnitude", eta_q, eta_phi, "lambda_c literal throat value"]`; the directive marks `G = GENUINE_BLOCKED`
(`directive.md:35`, `directive.md:226`). **The model delivers the FORM/branch, not Newton's `G`.**

**A full first-principles 1PN→4PN re-derivation from the throat interior stays SIM-DEFERRED (Gate 6)** — that
is out of reach (solver tractability, not hardware; `project-simulation-deferred-complete-pde-strategy`). A3
is the *reachable* consistency check, consistent with the whole program's calibrate-predict discipline and
the sim-deferred guardrail (completing the SPEC, not proving the theory).

## 5. Provenance / cross-refs

- **Supersedes nothing; duplicates nothing** — no prior 2.5PN-match-back note for pathA_43 existed under
  `notes/`. Related but distinct: `notes/25pn_notes.md` (pre-pathA_43 2.5PN working notes, Mar 28).
- **pathA_43 (Phase A1):**
  `software/stage1_solver/{directives,tools,reports}/pathA_43_density_quadrupole_port*`; verdict
  `DENSITY_PORT_HOSTED` (EARNED, tri-reviewed — caught+remediated a rig, re-tri-review + hardening CLEAN).
- **The PN corpus (separate, audited):** `research/4d_2_5pn` (conditional 2.5PN theorem); the broader ladder
  `research/4d_{1pn_full,2pn,3pn,4pn}` + `research/1pn_orbital_dynamics` (`project-pn-gravity-ladder` —
  calibrated / GR-matched controlled-reduction, `G` a genuine gap; **do not re-derive**).
- **Plan / front door:** `notes/ledger_v2_rebuild_plan.md` §4 (A3); `STATUS.md` ▶ RESUME HERE;
  `software/stage1_solver/decisions/13_emergent_constants_derivation.md` §0.
- **Next after A3:** Phase B — build `research/pde_ledger_v2/` (blueprint → machinery copy → assemble),
  where this note becomes a Part-II gravity stage; then Phase C (redteam).
