# G0 — Shared Minimal Closure Card: SCOPE

**Status: SCOPE DRAFT — 2026-07-20.** Not yet a build directive; this defines *what the card must fix and with what discipline*, for design-review (Codex → Grok → Codex) before the card's functional forms are authored. Requirement + acceptance only — it does **not** pre-design the functional forms (those are the card's postulation act).

## 0. Why this exists (the converged three-way conclusion)

We are computing the **electric force sign** by an analytic **four-channel far-field** "calibrate-and-exclude" pass, deferring simulation. A read-only harvest + independent advisories (Claude pre-registration, Codex xhigh, Grok 4.5 — all three converged) established:

- The existing analytic corpus is **silent on the throat surface rule `𝔅`**; every route to the sign runs through it (U2). ⇒ We must **postulate a closure** first — unavoidable, analytic or sim.
- The one existing "static electric" result (`U₊₊>0`, Coulomb-looking) is the **conservative `F_var` channel only**. The full force is `F = F_var + F_flux + F_𝔅 + F_rad`; the **flux (active drain)** and **`𝔅` (boundary)** channels — absent from every prior calc — can **flip** it. So a var-only kernel is the frozen-static artifact this whole pivot exists to escape.
- **The first *move* is ONE shared minimal closure card** where the parent-action pieces are *fixed* and **only `𝔅` differs** between branches. Two improvised closures would confound `𝔅` with parent-action choices and "relearn U2 by accident."

**This document scopes that shared card (G0).** The card is the prerequisite to grinding any branch; grinding is a separate downstream directive.

## 1. What the card FIXES (the shared parent pieces — postulate minimally, once)

The card supplies a **minimal** instantiation of the `§2.6` closure register sufficient to pose the **static one-throat and two-throat far-field** for the **electric (`h`-channel + drain) sector** in the two-sided `R_w` ambient. Every entry is tagged `[POSTULATE]` with its dimension and provenance; each is the *same* across all branches.

Essential (must be fixed):
1. **Wall-amplitude action** `r_B = |χ_B|`: kinetic (`Z_χ`), gradient stiffness / double-well potential (`κ_χ`, `V_χ`) — minimal regular form, non-singular where `r_B→0`.
2. **Sleeve geometry** — **frozen `Σ` for the MVP** (fixed sleeve shape, frozen `L/a=37/20`, `ε_r=1/20`), so the first pass is a *fixed-boundary* elliptic far-field, not a free-boundary solve. Frozen `Σ` means a **body-attached prescribed shape**, *not* a lab-pinned surface; the card must **book where the shape-maintaining reaction goes** (else a hidden support force conflicts with `F_𝔅=0`), and collar/interface conditions remain required even with bending/free-shape deferred. (Free-`Σ` is a later escalation, per the spec §6.3 R1-vs-R1-MVP split.)
3. **`χ↔h` coupling** (`g_χh`): how the sleeve amplitude sources the distinct `h` scalar (the electric mediator).
4. **Scalar-wave coefficients** (`M_h`, `K_h=M_h c_E²`, `C_hu`) with the **stability inequality `C_hu² < B_eff·K_h`, `M_h>0`** as a *hard admissibility constraint on the postulated values* (a card violating it is inadmissible, not a result).
5. **Drain + return, kept ACTIVE** (`g_χdrain`, `S_drain`, `S_leakage`): the throat sink and the balancing global return. **Non-zero by construction** — this is the anti-frozen discipline made concrete: the flux channel must be alive. Fixed via the finite-control-volume balance (`∫(S_drain+S_leakage)=0` at steady state).
6. **IR / return scheme**: two-sided far-field decay / return closure for the `R_w` ambient.
7. **Nondimensionalization**: `a=1` (or equivalent), the independent dimensionless groups, and the wall/mouth/sleeve/separation/outer scale hierarchy.

Deferrable / minimal-simplification (fix to a stated simple value, flag as a knob for later sensitivity):
8. **Phase representation + `θ_B`/`C_hu`** — one **consistent** choice, no double-counting: `Q_s(u_L,h)` vs the `(u_L,θ_B)` block (with bare `ρ_br` vs Schur-renormalized `A_eff`), and the `θ_B` contract (slaved/eliminated/independent — a *constant independent* phase is singular at `r_B→0`, so an elimination/slaving contract is required). State whether E3's "phase texture" is the bulk `θ` or the brane `θ_B`. `C_hu` may be fixed minimally (a decoupling limit), **but that is an explicit channel ABLATION, not neutral** (pathA39's pole residues depend on `C_hu`): a decoupled-card sign may **not** be promoted beyond that submodel without power-counting or a preregistered nonzero-`C_hu` sensitivity.
9. **Quantum pieces** `(Π_Q, ε_Q, 𝐣_ε^Q)` — one convention package (only needed for the energy-balance/finite-energy gate; a stated standard Madelung convention suffices for the MVP).

## 1a. Card-completeness — what the card must *fully specify*, not merely name

Codex's scope design-review (2026-07-20) confirmed the strategy is sound but flagged that a coefficient *list* does not by itself earn "no `UNRESOLVED(closure)` gaps." The card must additionally fix, as **explicit functionals** — with every omitted symmetry-allowed cross-coupling set to zero **or** a value:
1. **Complete bulk + surface + collar + interface + jump functionals.** `g_χh`, `g_χdrain` as *names* do not determine the Euler–Lagrange equations, jump laws, surface stresses, or momentum flux — the domain and surface action functionals must be written out (and the branch-specific completions: E2's action-derived traction, E3's permeability/texture + drain-source matching).
2. **The full drain/return momentum closure.** Not just `∫(S_drain+S_leakage)=0` (that closes *mass* only), but the spatial source laws, the momentum/energy sources `(𝐟_drain, 𝐟_leakage, w_drain, w_leakage)`, the declared **control surface** and **partition convention**, and **whether the return responds to incident `h`/pressure fields** (it feeds the flux force — the sign-deciding channel).
3. **The exact mouth functional.** The Neumann/Robin functional, annular support, normalization, sign convention, and conjugate natural condition — for the fixed-source ensemble (and later the fixed-defect exchange). Defined *independently* (see §3 anti-smuggling).
4. **The phase representation + zero-mode quotient** (§1 item 8).

These are the card's job (the postulation); they are enumerated here so the scope's acceptance (§6) is *earned*, not asserted.

## 2. What VARIES (only `𝔅`)

The **only** thing that differs between branches is the **boundary-operator completion on `Σ`**:
- **Branch 1 (lead) = E3** (permeable / phase-texture) — the far end of the permeability axis; it exercises the load-bearing **active-drain flux immediately**.
- **Branch 2 = E2** (impermeable + free-slip) — the clean contrast: it flips the *key* property (permeable ↔ impermeable) while keeping the **free-slip** condition natural to an inviscid, shear-free potential-flow superfluid. Both E2 and E3 have **`F_𝔅 = 0` by the committed bookkeeping**, so the comparison changes only the flow-contact property, not the channel partition.
- **E1 (impermeable + no-slip) is NOT a branch.** No-slip imposes the *full* surface velocity on a **shear-free, inviscid** potential flow — with no boundary layer, vortices, or extra surface DOF to support it, that is **over-determined**, and its reaction splits ambiguously between `F_var` and `F_𝔅`. E1 is retained only as a **later able-to-fail fixture**, not a physical baseline. *(3-advisor call: Codex + Claude for E2 over E1 on the over-determination + `F_𝔅=0`-clean argument; Grok preferred E1 on simplicity, overridden.)*
- **E4 (roll-lock) and E5 stay DEFERRED** to a later reduced *static elliptic* solve — **not** retired (see §5).

Each branch is otherwise identical (same card §1). This is the boundary-sensitivity experiment.

## 3. The mouth-ensemble policy (a co-load-bearing axis)

- Fix **one** ensemble for the first pass: the **fixed-source** case (pathA38-compatible, `~J²/K`), used for pipeline regression.
- The **fixed-defect** exchange (`~K·q_def²`, opposite stiffness scaling) is **required before any sign is promoted** — run as the diagnostic exchange (spec §5.3(a)). Never leave the ensemble free after seeing a sign.
- The ensemble choice is stated on the card label.
- **Anti-smuggling (load-bearing):** the mouth boundary functional is defined **independently** (its Neumann/Robin form, support, normalization, sign convention, conjugate condition) — pathA_38 is then used only as a **regression check** on the reduced kernel, **never** by fitting the parent source to reproduce `U₊₊>0`. The pathA_38 propagator itself enters only after the card supplies its Tier-2 embedding, normalization, and parent-operator match.

## 4. The force observable + exclusion gates (what the card must make computable)

- **Four-channel far-field**: `F = F_var + F_flux + F_𝔅 + F_rad`, computed **channel-resolved**, over the four `{s1,s2}` orientation cells (two-sided `R_w`). For the E2/E3 branches, `F_𝔅 = 0` and (stationary) `F_rad = 0` — but these zeros must be **demonstrated** from the stationary outgoing solution, not silently assumed. So the operative passive-branch force is `F_var + F_flux`; **dropping `F_flux` is forbidden.**
- **Only the TOTAL force is invariant** — the channel split depends on the declared control surface and partition convention (both stated on the card).
- **HARD RULE: a var-only kernel is NEVER promoted to a physical electric sign.** A var-only pass is a **regression test only** — it must reproduce pathA_38's relative `±w` sign matrix and verify the Hadamard projection; it is not evidence for the physical sign. `F_flux` (drain alive) enters before any attract/repel sentence, unless power-counting *proves* its orientation-pair component subleading.
- **Pre-built exclusion gates** the card must feed: pathA_38 FAIL_* forms (delocalized→`1/r³`, ghost, no-monopole, pinned-branon, Yukawa) + the §1(4) stability inequality. A branch whose reduced kernel hits a FAIL is **excluded** immediately.
- **`R_w`-conditional**: every sign is conditional on the two-sided ambient postulate (stated on the label).
- **Held-out**: the `~10⁴²` hierarchy magnitude and the detailed magnetic/Darwin structure are **kept out** of calibration (the real predictive test on survivors).

## 4a. Computation strategy — response-first, with a go/no-go and a fallback

Do **not** assume the two-body far-field is already analytic. pathA_38 solved a far-zone Green kernel for a *prescribed* compact source; it did **not** derive that source from the sleeve/drain/mouth boundary — so the **leading** monopole amplitude and the flux cross-term are **not** in hand (E3 needs the permeability/texture law + drain-source matching; E2 needs the action-derived traction; these set the leading term, not just short-range corrections). Therefore the pass is **response-first**:
1. Compute a **one-throat static affine T-matrix** per branch — `permanent sources + incident long-wavelength (h, pressure-phase, brane) fields → induced multipoles, flux read-outs, and reactions` — using the full coupled propagator, source normalization, multiple-scattering/self-consistency, and *separate* force/flux read-out operators. (An active, possibly non-self-adjoint drain cannot assume reciprocity; a *homogeneous* response matrix is insufficient — the affine "permanent-source" term is required.)
2. **Assemble the two-body far-field** by connecting two T-matrices with the already-known far-field propagators.
3. **Go/no-go (refined):** if the **total Hadamard `(1,1)` force coefficient at `1/R²`** is not uniquely determined, finite, and control-surface invariant under the completed card — or the one-body factorization fails — the analytic sign is **no-go**. (A *proved zero* is a null/range **result**, not a no-go.) **Missing closure sends the work back to the card, not to numerics.** If the obstruction is an unknown *local one-body response profile*, escalate to the fixed-`Σ` 2D R1-MVP solve; if it is a *shared nonlocal return* or a *pair-irreducible flux term*, a one-body 2D solve cannot cure it → escalate to the R3 two-center solve.

This gives the pass a clean exit criterion and a defined fallback ladder, not an open-ended grind.

## 5. Discipline / guards (baked in, from the advisor corrections)

1. **No frozen fields.** The drain is always active (§1.5); solve the throat's *response*, never a frozen snapshot.
2. **V→0 cross-check.** Any quantity labeled "steady-state" must be verified as the `V→0` limit of a slightly-moving version before it is trusted (the 1PN lesson as a test).
3. **Agreement does NOT retire E4.** If E2 and E3 agree, that is robustness **within the passive/holonomic family only** — E4 injects a qualitatively new nonholonomic constraint reaction with no counterpart in the passive rules and no coefficient bound tying it to them, so it can flip the sign regardless. E4 stays scheduled; the passive pass hands it a calibrated baseline + form targets. (Disagreement quantifies sensitivity; it does **not** "re-confirm U2" — U2 already established underdetermination.)
4. **Preregister** the branch, the mouth ensemble, and the `R_w` label **before** any sign is computed. No card entry may be varied after seeing a sign.
5. **Shared card invariance.** Branches 1 and 2 must use byte-identical parent pieces; only the `§2` `𝔅` slot differs. A grind that alters a parent piece per-branch is invalid as a sensitivity experiment.

## 6. Acceptance (for the card, at design-review + build)

- Every entry tagged `[POSTULATE]`/`[CALIBRATED]`/`[CONVENTION]` with dimension + one-line provenance; the stability inequality satisfied by the postulated values.
- The card poses a **fixed-`Σ`, static, two-sided-`R_w`** far-field problem for E2 and E3 that is **instantiable** (no `UNRESOLVED(closure)` gaps for the electric sector) — i.e. the four-channel far-field is *computable in principle* from the card alone. This requires the full **§1a completeness checklist**, not only the §1 coefficient list.
- The drain/return is non-zero and balanced (flux channel live).
- The FAIL_*/stability gates are wired to the card's reduced kernel.
- A one-line statement of what the card does **not** fix (free-`Σ`, E4/E5, magnetism, the magnitude) and where each is deferred.

## 7. Process (what happens after this scope)

1. **Design-review this scope** (Codex → Grok → Codex) — is the minimal parent set correct/sufficient? are the deferrals safe? is anything smuggled?
2. **Author the card** — postulate the §1 functional forms (the conceptual step; user-involved on the form choices).
3. **Grind directive** — the response-first four-channel far-field pass (§4a), **E3 (lead) then E2** under the one card, with §4 gates and §5 guards. Separate directive. Full sequence: **closure/power-counting audit → E3 one-body response matrix → matched E2 response matrix → analytic pair assembly (go/no-go) → fixed-`Σ` 2D R1-MVP solve if a `1/R²` coefficient is undetermined → E4 reduced solve if still viable.**
4. **Grind-directive also fixes** (two-body machinery, not the closure card): the fixed-separation force convention + holding-reaction sign, the error/verdict thresholds, and the four-orientation `{s1,s2}` run manifest.
5. **Deferred escalations**: fixed-defect ensemble exchange (before any sign); E4/E5 reduced static solve; free-`Σ`; magnetism (R2+); E1 as an able-to-fail fixture.

**Bottom line:** G0 is the shared, minimal, postulated closure that makes the boundary-sensitivity grind non-vacuous. It fixes the parent physics once; varies only `𝔅` (**E3 lead, E2 contrast**; E1 dropped as over-determined for a shear-free superfluid; E4/E5 deferred, not retired); keeps the drain alive; runs **response-first with a go/no-go and a fixed-`Σ` 2D-solve fallback**; forbids var-only signs; and treats E2/E3 agreement as passive-family robustness only.
