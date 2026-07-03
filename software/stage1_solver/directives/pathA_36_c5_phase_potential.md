# pathA_36 — C5 decidability: can a brane-order-parameter phase supply the MacCullagh scalar potential?

**Status:** DRAFT for Codex design-review. **Type:** focused, flat-brane, linearized consistency test (single crux).
**Frozen baseline:** `T0_SHEAR_FROZEN(d9520d3819c3)` (pathA_35 G0). **Author:** orchestrator (scaffolding). **Date:** 2026-07-03.

---

## §0. Provenance and scope

The brane sector's crux is **C5** — the MacCullagh longitudinal zero mode. MacCullagh's curl-only potential `½μ_R(∇×u)²` is
gauge-invariant under `u→u+∇χ`, but the kinetic term `½ρ_br(∂_t u)²` is not, so the EOM forces `∂_t²(∇·u)=0` — the longitudinal
sector is a **constrained physical zero mode**. Maxwell escapes via a scalar potential (`φ→φ−∂_tχ`); the frozen light package has
**no φ**, so pathA_35 Gate L returned `FAIL_COUPLE_STRESS_NOGO` (and on the slaved-rigid `P` branch, which already clears the
hidden-mode and closure horns, **C5 is the sole surviving killer**).

The material-state pivot (`docs/conceptual_foundation.md` §2 v5) proposes the brane is an ordered phase of the one medium described
by a **complex** order parameter `χ_B = |χ_B| e^{iθ}`, and that its **phase θ** may supply the missing φ. A GLM reframe
(`notes/rung_W_reframe.md`) found **rung W passes** (a double-well `χ_B` wall is stable and does not starve in-plane shear) but
declared the θ-as-φ sub-hypothesis **`FATAL_FLAW`**, on the ground that an order-parameter phase couples to its conjugate *density*,
not to the displacement divergence `∇·u`.

**This directive tests that NEGATIVE verdict — as a hypothesis to BREAK, not accept** (`feedback-negative-verdict-short-circuit`:
a "can't-be-done" gets *harder* verification than a "did"). The motivating counter-observation to be adjudicated: if the
brane-order density fluctuation is identified with the longitudinal displacement, `δρ_B ≈ −ρ_B0 ∇·u` (number conservation for the
ordered-phase constituents), then the Josephson/Berry coupling `(∂_t θ) δρ_B` integrates by parts (space + time, boundary ledger
preserved) to the Maxwell cross-term `C_J·(∇θ)·(∂_t u)` with a **derived** sign/normalization `C_J = −J ρ_B0` (not automatic — see
§2.2). So the electric energy `½ε(∂_t u + ∇θ)²` *may* assemble from three provenanced pieces: **brane inertia** `(∂_t u)²` +
**Josephson coupling** `(∂_t u)·(∇θ)` + **phase stiffness** `(∇θ)²`. Whether it **actually**
assembles into the gauge-invariant perfect square — with genuine mechanical provenance, without disturbing the two transverse
photons, and without introducing a stray longitudinal mode — is a decidable calculation. **That calculation is this directive.**

**Scope:** flat-brane, linearized — C5 is a flat-brane light-sector question (the same flat-brane analysis pathA_35 Gate L ran); the
full `χ_B` wall (rung W) and the couple-stress details (a-i/a-ii/b) are **out of scope** except as non-disturbance checks. Outcome
routing: C5 resolves **with provenance** → the flat-brane C5 necessary condition is met (the full material-state route still needs
rung W + the Gate-L re-run) → proceed. C5 cannot be
resolved → shear-surface/material-state light is falsified at C5 (a first-class no-go) → pose the conceptual fork (non-Maxwell toy
light vs the previously-declined lattice/Wen route) to the user.

---

## §1. Methodology (BINDING)

1. **Analog, not derivation.** Postulate the minimal θ structure freely; test internal **consistency**; a genuine no-go (θ cannot
   supply φ without an unprovenanced term or a stray mode) is a **first-class, welcome** result — do not rescue or soften.
2. **Provenance rule** (inherited from pathA_35 §2.3). Any term that resolves C5 must be a **local variational structure with
   independent mechanical provenance** — a real term in the action arising from the complex brane-order parameter and its coupling
   to the medium. A raw `∇·u=0` projector, or an "electric" term postulated purely to force the perfect square via a **free-tuned**
   coefficient with no provenance, **FAILS or degrades** the verdict (labeled a postulate).
3. **Able-to-fail is mandatory.** The computation MUST be able to return each of {C5 zero mode, stray longitudinal mode, disturbed
   transverse sector, extra scalar DOF}. A control that cannot fail is rejected (`feedback-decisive-test-not-tautological`).
4. **Dual-engine REQUIRED.** SymPy (.py) AND Mathematica (.wl) independently assemble the quadratic action, run the gauge check, and
   compute the physical-DOF count; `ENGINE_AGREE` on the **HEADLINE** quantities (the longitudinal-sector Hessian/constraint
   structure and the mode counts), each engine assembling the headline itself — **no `x−x` vacuity** (`feedback-dual-engine-required`,
   and the pathA_32 dual-engine-gaming lesson).
5. **No result-emitters.** The mode count, gauge-closure result, and verdict must be **DERIVED** (a real constraint/Hessian/dispersion
   computation from the action), never typed literals keyed on a config string (the pathA_35 result-emitter lesson).
6. **Dimensional firewall** (`feedback-dimensional-consistency-check`): every term units-restored, with ≥2 able-to-fail dim ablations
   that MUST fire. We do not trust any engine's dimensional bookkeeping without an ablation that can catch a dropped scale.

---

## §2. Frozen inputs and the candidate θ structure

### §2.1 Frozen (reuse pathA_35 G0 — do not re-open)
- Flat brane, linearized; Fourier basis with in-plane wavevector `k` (choose `k` along one in-plane axis; decompose `u` into the
  **2 transverse** components (⊥`k`) and the **1 longitudinal** component (∥`k`, carrying `∇·u`)).
- In-plane displacement `u^a`; brane inertia `½ρ_br(∂_t u)²`.
- MacCullagh curl-only potential `½μ_R(∇×u)²`; `c_γ² = μ_R/ρ_br`.
- Out-of-plane `u_w` gapped (frozen `Ω_w > 0`).
- `P^i` arrows in the **slaved-rigid** form `P_∥ = ŵ×(∇×u)` (the branch on which a-ii/b pass), or held as a fixed background — the
  couple-stress sector is not the object under test, but its **non-disturbance** by θ must be checked.

### §2.2 The candidate θ sector — freeze the PRIMITIVE action, then only diagnose (B1/B3/B4/B5)

Introduce the phase `θ` of a complex brane-order parameter `χ_B = |χ_B| e^{iθ}`. The engines **MUST begin from an unreduced
primitive quadratic Lagrangian with named, frozen coefficients** — square-completion is a **diagnostic of the output, NEVER an
input** (a pre-completed `½ε(∂_t u+∇θ)²`, or a free `ε` inserted at the start, is a result-emitter → reject, per the pathA_35
lesson):

`L = ½ρ_br(∂_t u)² − ½μ_R(∇×u)² − ½B(∇·u)² + J(∂_t θ)δρ_B + ½K_θ(∇θ)²  [+ conjugate-density/amplitude terms per the branch below]`

with `c_γ² = μ_R/ρ_br` frozen, and `K_θ` the **SIGNED Lagrangian coefficient** of the phase-gradient term. **Its sign is a DERIVED,
load-bearing output, not an input:** a conventional stable phase stiffness is `K_θ = −κ_phase < 0` (the ordinary "potential" sign —
gradient energy `+½κ_phase(∇θ)²` in the Hamiltonian), whereas Maxwell's *electric* structure requires `K_θ > 0` (the opposite sign),
which generically signals a ghost/unbounded sector unless the order-parameter action derives it WITH a bounded-below Hamiltonian.
**Curl-only requires `B=0`** (a nonzero brane compression modulus `B` is itself a Cauchy `(∇·u)²` term) — and `B=0` must be
**justified from the order-parameter action, not assumed**.

**The `δρ_B` / amplitude decision (the crux fork — B4/B5).** `θ` is the phase of a complex order parameter, hence canonically
conjugate to the order-parameter **density/amplitude**, not directly to `u`. Freezing `|χ_B|=1` while using its density fluctuation
as `θ`'s partner is inconsistent unless the amplitude/density mode is explicitly handled. **Both branches are required; branch (b)
is the main single-medium test:**
- **Branch (b) — slaved (single-medium-preferred; the "carried by the medium" firewall):** `δρ_B = −ρ_B0 ∇·u` imposed as a genuine
  variational constraint = number conservation of the ordered constituents advected by the displacement `u` (the same
  carried-by-the-medium coupling that turns Frank into MacCullagh). No independent density mode. **This is the main test.**
- **Branch (a) — independent (the "second-medium?" alternative):** `δρ_B` is an independent field with its own compressibility
  `½(δρ_B)²/χ_c` (amplitude stiffness) + continuity `∂_t δρ_B + ρ_B0 ∇·(∂_t u) = 0`, then **integrated out with an explicit proof**
  (not by convenience); report what `θ`-kinetic structure the integration induces.
- **Compressibility is load-bearing (B5):** under `δρ_B=−ρ_B0∇·u` a finite conjugate-density stiffness becomes a `(∇·u)²` **Cauchy**
  term directly bearing on C5. Whether it is zero, finite, or infinite must be **derived from the order-parameter action**, not
  chosen for convenience.

**Josephson cross-term (B3).** `J(∂_t θ)δρ_B = −Jρ_B0(∂_t θ)(∇·u)` integrates by parts (space + time, boundary ledger preserved) to
`C_J(∂_t u)·(∇θ)` with the **derived** sign `C_J = −Jρ_B0` — report it, do not assume it.

**Gauge transformation (S2).** Test `u→u+∇χ`, `θ→θ+g(χ)` with `g` **derived** from the action, not imposed, and required **local and
off-shell** (e.g. `∝ ∂_tχ`; no inverse-`k`, inverse-`ω`, or dispersion-shell-only closure). Maxwell's is `θ→θ−∂_tχ` up to
normalization.

**Ledger.** Count every new field (`θ`; `δρ_B` in branch (a)) and constant (`ρ_B0`, `J`, `K_θ`, `χ_c`, `B`, any others) as drift,
pathA_35 `SECOND_MEDIUM_DRIFT` style. Report `DRIFT(n)`.

---

## §3. The question, the computation, and the verdict grammar

**THE QUESTION.** In the flat-brane linearized theory (frozen §2.1 + candidate §2.2), does the **longitudinal (`∇·u`) + θ sector**
reduce to **ZERO physical propagating degrees of freedom** — the longitudinal displacement removed as pure gauge, `θ` acting as a
Gauss-law constraint (Maxwell) — **while the two transverse photons remain exactly massless and undisturbed**?

**Compute (route is Codex's to design — this directive fixes the WHAT, not the HOW):**
- the full quadratic action for `(u_transverse, u_longitudinal, θ [, u_w, P slaved])` in `(ω,k)`, from the **primitive** §2.2
  Lagrangian (not a pre-reduced one);
- the **physical propagating DOF count** per sector via **full Dirac–Bergmann constraint analysis on the original first-order
  Lagrangian** (B8): canonical momenta; primary/secondary constraints; the first-class vs second-class split; the gauge generator;
  Hamiltonian boundedness; pole residues; and the number of independent initial-data functions per finite `k`. **Hessian
  null-space / propagator poles ALONE are insufficient** for the first-order Josephson term — they misclassify a legitimate
  constrained Maxwell-like sector (the false-FAIL risk);
- whether the gauge symmetry `u→u+∇χ, θ→θ+g(χ)` leaves the **full** action invariant **off-shell for all `k`** (not just at `k=0` or
  on the dispersion shell), and if so whether it renders the longitudinal sector pure-gauge;
- the **provenance status** against the **objective Maxwell locus (B2/B9).** After IBP the longitudinal/electric sector is
  `½ρ_br(∂_t u)² + C_J(∂_t u)·(∇θ) + ½K_θ(∇θ)² − ½B(∇·u)²` with `C_J=−Jρ_B0` and `K_θ` the signed gradient coefficient (a
  conventional stable phase → `K_θ<0`). A gauge-invariant perfect square `½ρ_br(∂_t u + (C_J/ρ_br)∇θ)²` requires ALL of:
  **(i) `C_J² = ρ_br·K_θ`** (square closes) — note this **REQUIRES `K_θ>0`, the anti-conventional "electric" sign; a conventional
  positive phase stiffness (`K_θ<0`) cannot close the square → automatic `FAIL_C5_LONGITUDINAL_ZERO_MODE`**; **(ii) `ε = ρ_br`
  exactly** — the coefficient of `(∂_t u)²` in the electric square IS the frozen brane inertia, not a new knob (else `c_γ` shifts →
  false PASS); **(iii) `B = 0`** (no residual Cauchy `(∇·u)²`); **(iv)** `K_θ>0` must come with a **bounded-below Hamiltonian** (the
  electric-sign gradient can be a ghost → `FAIL_GHOST_OR_NEGATIVE_NORM`). **`WITH_PROVENANCE` iff (i)–(iv) FOLLOW from the
  previously-frozen definitions** of `ρ_br, ρ_B0, J, K_θ, |χ_B0|, χ_c` (in particular the order-parameter action must **derive the
  Maxwell sign `K_θ>0` with boundedness**); **`BY_TUNING`** if any constant or sign must be adjusted to satisfy them (report the
  exact condition). This is the decisive discriminator.

**Verdict grammar.**
- **`C5_RESOLVED_MAXWELL_WITH_PROVENANCE`** — longitudinal sector → 0 physical DOF via a genuine (first-class, off-shell, all-`k`)
  gauge symmetry; transverse = 2 massless photons undisturbed with `c_γ²=μ_R/ρ_br`; the Maxwell locus (i)–(iv) holds from frozen
  definitions — **including the derived electric sign `K_θ>0` with a bounded-below Hamiltonian**. **The strong pass** → the flat C5
  leg is viable.
- **`C5_RESOLVED_MAXWELL_BY_TUNING`** — Maxwell holds only on a measure-zero locus reached by tuning a constant **without**
  provenance; labeled an additional postulate, conditionality degraded. **Report the exact tuning condition.**
- **`FAIL_C5_LONGITUDINAL_ZERO_MODE`** — θ fails to couple, OR the coefficients miss the locus, leaving `∇·u` a constrained physical
  zero mode (pathA_35 outcome persists).
- **`FAIL_CAUCHY_STRAY_LONGITUDINAL`** — the longitudinal+θ sector has **1** propagating physical mode (a stray longitudinal "third
  photon" / Cauchy-type wave; e.g. from `B≠0` or `C_J²≠ρ_br·K_θ`).
- **`FAIL_GAPPED_NOT_GAUGE_REMOVED`** — the longitudinal/scalar mode is **gapped (massive)** rather than removed by gauge (not
  Maxwell).
- **`FAIL_PARTIAL_GAUGE_ONLY` / `FAIL_GAUGE_ON_SHELL_ONLY`** — the gauge symmetry holds only at `k=0`, or only on the dispersion
  shell, not off-shell for all `k` (not a true gauge symmetry).
- **`FAIL_GHOST_OR_NEGATIVE_NORM`** — the first-order Josephson structure yields a ghost / negative-norm / unbounded-below
  longitudinal-scalar sector.
- **`FAIL_SECOND_CLASS_NOT_MAXWELL`** — the constraints are **second-class** (a pole is removed, but the initial-data count / gauge
  structure is not Maxwell's first-class one).
- **`FAIL_EXTRA_SCALAR_DOF`** — θ is an independent **propagating** scalar that neither removes the longitudinal mode nor pairs with
  it (a hidden mode / fifth-force-like scalar).
- **`FAIL_TRANSVERSE_DISTURBED`** — adding θ changes the transverse sector (mode count ≠ 2, `c_γ` altered, or a mass generated).
- Sub-tags as applicable: `DRIFT(n)`; `AXIS_RE_ADMITTED` (θ's coupling drags in `ŵ`/a chiral term — the pathA_35 parity trap);
  `U_W_COLLISION` (θ must be identified with, or gap-collides with, `u_w`).

---

## §4. Able-to-fail controls (each a real computation; each must fire correctly)

1. **No-θ control (B6)** — remove the **entire θ sector** (field, Josephson coupling, AND phase stiffness) → must return
   `FAIL_C5_LONGITUDINAL_ZERO_MODE` (recovers pathA_35). *Proves the test can see the zero mode.*
2. **Cauchy control** — replace the curl-only potential with a Cauchy bulk modulus `½λ(∇·u)²` (no θ) → must return a propagating
   longitudinal mode. *Proves the test can see a stray longitudinal mode.*
3. **Mismatched-coefficient control** — detune the coefficients **off** the locus (`C_J²≠ρ_br·K_θ` — **including `K_θ≤0`, the
   conventional/wrong (potential) sign** — or `ε≠ρ_br`, or `B≠0`) → must return `FAIL_CAUCHY_STRAY_LONGITUDINAL` /
   `FAIL_GAPPED_NOT_GAUGE_REMOVED` / `FAIL_C5_LONGITUDINAL_ZERO_MODE`, **not** a spurious pass. *Proves the pass is not generic;
   measures the tuning width.*
4. **Decoupled-θ control (B6)** — θ present with its **complete physically-justified sector** (including whatever kinetic term the
   §2.2 compressibility/amplitude branch induces) but `C_J=0` → must leave the C5 zero mode AND, if θ is dynamical, return
   `FAIL_EXTRA_SCALAR_DOF`; if θ is non-dynamical it simply decouples (report which). *Proves the coupling is load-bearing.*
5. **Transverse control** — confirm the 2 transverse photons are exactly massless with `c_γ²=μ_R/ρ_br` in ALL cases. *Guards
   `FAIL_TRANSVERSE_DISTURBED`.*
6. **Provenance ablation (the GLM crux)** — evaluate the Maxwell locus (i)–(iii) with the coefficients **fixed by the order-parameter
   structure** vs **free**; report whether provenance **forces** the locus or whether it needs **tuning**. *Distinguishes
   `WITH_PROVENANCE` from `BY_TUNING` — the whole verdict turns on this.*
7. **Compressibility control (B5)** — run branch (b) and branch (a) with the conjugate-density/amplitude stiffness **included vs
   absent**; a finite stiffness must surface as a `(∇·u)²` Cauchy term (→ `FAIL_CAUCHY_STRAY_LONGITUDINAL` unless the θ structure
   gauge-removes it). *Guards the omitted-compressibility false result.*
8. **θ-mass control (S3)** — add a phase-locking mass `½m_θ²θ²`; it must **break** the gauge structure (→ a FAIL) and be counted. If a
   symmetry forbids `m_θ`, state the symmetry. *Confirms the gauge closure genuinely needs θ massless.*

---

## §5. Dual-engine + dimensional firewall (REQUIRED)

- **SymPy (.py) and Mathematica (.wl)** independently: assemble the quadratic action; perform the gauge-invariance check; compute the
  physical-DOF count / dispersion; run all §4 controls. `ENGINE_AGREE` on the longitudinal-sector Hessian/constraint structure, the
  mode counts, and the verdict — **each engine assembling the headline itself** (no `x−x`; both engines must independently reach the
  DOF count, not transcribe one literal).
- **Comprehensive inline dimensional firewall** — every term (brane inertia, MacCullagh potential, Josephson coupling, phase
  stiffness, the electric energy `½ε(...)²`, `c_γ²`, any gap/mass) units-restored; **≥2 able-to-fail dim ablations that MUST fire**.
- **Timeouts** — every script runner `timeout 600`; a timeout (124) is a failure → reformulate the math, never raise the cap.

---

## §6. Discipline / scope

- Flat brane only; full `χ_B` wall (rung W) and couple-stress details are out of scope except as non-disturbance checks.
- **Do not rescue.** If C5 cannot be resolved with provenance, report the FAIL plainly and **identify which term won't line up** (the
  physical reason) — that reason is the deliverable's most valuable content.
- Orchestrator (Claude) reviews + arbiter-re-runs existing scripts unchanged; Codex codes/runs/iterates; no orchestrator code
  mutation.
- **Deliverables:** `tools/pathA_36_c5_{sympy.py,.wl}`; `reports/pathA_36_c5_phase_potential.md` + `_results.yaml` (YAML for the
  machine-readable result); the verdict + the DERIVED mode counts + the provenance status + all §4 controls + the dim firewall.

---

## §7. Acceptance (what "verified" means, per the user)

**Nothing is accepted until:** (a) it is **derived via BOTH SymPy and Mathematica** with `ENGINE_AGREE` on the headline quantities;
AND (b) it is **double-reviewed by Claude AND Codex** — this directive design-reviewed to green (Codex + Claude) before execution;
then, post-execution, the standing tri-review: orchestrator arbiter re-run + a transliteration-fidelity audit on a **fresh clean
agent** + an adversarial-with-ablation review on a **fresh clean agent** + Claude's own substantive read. **A negative verdict (any
`FAIL`) gets the harder scrutiny** (the pathA_35 result-emitter lesson): confirm the FAIL is genuinely DERIVED and **able-to-pass**
(the §4 controls demonstrate the gate can flip), not a hardcoded / config-keyed emitter.

---

## §8. Changelog
- v1 (2026-07-03) — initial draft for Codex design-review. Tests the GLM `FATAL_FLAW` on θ-as-C5-φ as a hypothesis to break;
  flat-brane, dual-engine, able-to-fail; verdict distinguishes provenanced Maxwell from tuned Maxwell from the three failure modes.
- v2 (2026-07-03) — folded Codex design-review (NOT_SOUND → 9 blockers + 3 suggestions), Claude-reviewed: freeze the **primitive**
  Lagrangian (square-completion diagnostic-only, B1); soften the IBP claim + report derived `C_J=−Jρ_B0` (B3); the **δρ_B/amplitude
  fork** — branch (b) slaved single-medium main test vs branch (a) independent-integrated-out (B4) + compressibility as a load-bearing
  `(∇·u)²` Cauchy term (B5); the **objective Maxwell locus** `C_J²=ρ_br κ_θ`, `ε=ρ_br`, `B=0` with provenance-from-frozen-defs the
  WITH/BY discriminator (B2/B9); **full Dirac–Bergmann** DOF counting for the first-order Josephson term (B8); fixed controls 1/4 (B6)
  + added compressibility & θ-mass controls (B5/S3); expanded verdict grammar with gapped / partial-gauge / on-shell-only / ghost /
  second-class failure labels (B7); off-shell-all-`k` gauge requirement + local `g(χ)` (S2); softened §0 routing to "flat C5 leg
  viable, not the whole route" (S1). Next: Codex confirm-pass → SOUND before execution.
- v3 (2026-07-03) — folded the Codex confirm-pass (NOT_SOUND, 11/12 ADDRESSED; B2 PARTIAL), Claude-verified: fixed the sign
  contradiction between the primitive `(∇θ)²` term and the Maxwell locus. Now the phase-gradient term is a **signed** coefficient
  `+½K_θ(∇θ)²`; a conventional stable phase is `K_θ<0` (potential sign), Maxwell needs `K_θ>0` (electric sign) → the locus is
  `C_J²=ρ_br·K_θ` (auto-FAIL for conventional stiffness) with the added condition (iv) that `K_θ>0` carry a bounded Hamiltonian (else
  ghost). This makes the electric-vs-potential SIGN of the phase-gradient the crux discriminant. Next: Codex re-confirm → SOUND.
