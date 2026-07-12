# Build Spec — Analog Four-Force Phenomenology Visualizer

**Status:** DRAFT for user gate → then hand to Codex (xhigh).
**Author:** Claude (coordinator; owns this spec, does not write the code).
**Builder:** Codex (xhigh) builds ALL code. **Reviewer:** Grok (against the cited docs) → Codex folds → iterate to clean. **Final:** a fresh review agent, once the whole thing is done.
**Location:** `software/force_visualizer/`.

---

## 0. What this is (and is NOT)

A **2D phenomenology visualizer** that shows the four force sectors of the toy superfluid-analog model *in action*, using the **reduced / effective, CALIBRATED force laws** — NOT the full 4D PDE.

- It **visualizes calibrated predicted phenomenology.** It does **NOT** simulate emergence-from-the-medium (that is the deferred, intractable full 4D nonlinear PDE). Do not present it, name it, or comment it as "deriving" or "emerging" the forces. Magnitudes and cross-sector speed ratios are **inputs**, not outputs.
- **2D only**, matplotlib. No web, no 3D. (If 3D is ever wanted, the render layer will be swapped — so keep physics render-agnostic; see §2.)

## 1. Source of truth (hard rule)

The **ONLY** authoritative source for every physics law is the cited research doc below. **Do NOT use, copy, import, or depend on `superfluid_lib/`** — its correctness is not trusted (it predates much of the current research). You may glance at it as a *non-authoritative* reference, but anything you take from it must be independently re-derived from and verified against the cited docs. Same for the existing `research/1pn_*/scripts/*.py` sims: reference-only, verify against the papers.

**Read-only** the ledger material under `research/pde_ledger_v2/` (a separate session is actively writing there — do not write anything under it).

## 2. Architecture (requirements; you own the implementation details within these)

Hard separation of concerns:

- **`physics/`** — pure, deterministic force-law + integrator functions. **No matplotlib / no rendering imports whatsoever** (enforced by a test, §6). Each module/function docstring **cites its source doc** (path + section/result key). This is a small, unit-testable library.
- **`scenes/`** — one matplotlib scene per force that *calls* the physics core and animates it (FuncAnimation). Rendering only.
- **`params.py`** — the single shared parameter block (§4).
- **`tests/`** — the analytic golden checks (§6).
- **`notes/`** — this spec; you may add a short README.

Rendering = matplotlib is the default and expected choice; you own the details (integrator scheme, API shape, file layout). Do not pre-optimize (problem size is tiny — a handful of bodies, small grids; plain numpy is real-time-or-faster).

## 3. The four sectors — effective forms to implement

For each: the effective form below is the **acceptance target / pointer**; the **cited doc is authoritative** — verify against it, and if the doc disagrees with this summary, follow the doc and flag the discrepancy. Implement FORM faithfully; magnitudes may be chosen (labeled CALIBRATED).

### 3.1 Gravity
- **Form:** Newtonian `a = -GM r̂ / r²` (1/r² inflow/drain, attractive, p=2) **plus** the post-Newtonian ladder for a visible orbit precession (1PN) and optionally a slow inspiral (2.5PN radiation reaction, Burke–Thorne).
- **Form derived, magnitude CALIBRATED** (G is a free effective Newton constant; the 2.5PN normalization is a benchmark target — `G = GENUINE_BLOCKED`).
- **Departure to expose (§5):** the bounded monopole/dipole `c_s`-radiation residual ∝ drain strength `ε0/(1+ε0)` — GR forbids ℓ=0/1 GW; this is the genuine falsifiable prediction of the drain picture.
- **Docs:** `software/stage1_solver/reports/pathA_29_brane_bulk_return.md` (+ `pathA_29_results.yaml`) for the 1/r² localization + the radiation residual; PN two-body forms in `research/4d_1pn_full/` and `research/4d_2_5pn/` (papers). `research/1pn_orbital_dynamics/scripts/` = reference-only.

### 3.2 Light
- **Form:** transverse (MacCullagh rotational-shear) wave on the brane, `ρ_br ∂²ₜ u_T = μ_R ∇² u_T`, `c_γ² = μ_R/ρ_br`, **2 transverse polarizations, no longitudinal**, dispersion `ω = c_γ k`. Plus a **lensing** scene: photon ray bending through a refractive profile `n(r)` around a mass.
- **Form derived, magnitude CALIBRATED** (μ_R, ρ_br → c_γ chosen; `λγ = c_γ/c_s` is a free calibration input).
- **Departure to expose (§5):** the stray longitudinal mode (`FAIL_CAUCHY_STRAY_LONGITUDINAL`) — clean 2-DOF Maxwell is reachable only BY TUNING (`B_eff → 0`).
- **Docs:** `software/stage1_solver/reports/pathA_36_c5_phase_potential.md` (+ results yaml); `research/pde_ledger_v2/notes/stages/ledger_stage003_transverse_photons_stray_longitudinal.md`, `ledger_stage005_sound_speed_light_ratio.md`; lensing ref `research/1pn_optics/` (verify).

### 3.3 Electric charge
- **Form:** brane-localized `1/r²` Coulomb, `F ∝ q₁q₂ / r²`; **sign = ±w puncture orientation**; **like-repel / unlike-attract** (from `G₀ > 0`). Parity ODD (charge ≠ mass).
- **Form derived, magnitude CALIBRATED** (coupling `Q_E` chosen; `N0` itself is derived).
- **Departure to expose (§5):** exponentially-suppressed **Yukawa** short-range corrections (`~ exp(-√3 R/ℓ)/R`) from gapped partner modes, on top of the long-range Coulomb.
- **Docs:** `software/stage1_solver/reports/pathA_38_throat_body_electric_localization.md` (+ `pathA_38_results.yaml`).

### 3.4 Magnetism
- **Form:** velocity-dependent (`O(V₁·V₂)`) current–current force between moving charges (in the throat, not bulk vorticity): `1/R` potential → **`1/R²` point force**; transverse channel → **like currents attract** (correct EM sign; flips with `μ_R → -μ_R`).
- **Form derived, magnitudes CALIBRATED / sim-deferred** (`aT, aL`, and `c_E = c_γ` are free).
- **Departure to expose (§5):** the **unavoidable attractive longitudinal scalar-current admixture** (same sign, uncancelable) — full field content is transverse-vector + charge-coupled scalar (`h`-branon), NOT exact Maxwell (`FIELD_SCALAR_VECTOR_DEPARTURE`); preferred-frame unless `c_E = c_γ`.
- **Docs:** `software/stage1_solver/reports/pathA_39_magnetic_force.md` (+ `pathA_39_magnetic_force_results.yaml`), `pathA_39_stage4_field_classification.md`, `pathA_39_scalar_admixture_screen.md`.

## 4. Parameter block (`params.py`)

One **internally-consistent, hand-picked** set, each value **labeled DERIVED vs CALIBRATED** with a one-line provenance + doc ref:

- Speeds: `c_s`, `c_γ` (→ `λγ = c_γ/c_s`), `c_E`. **`λγ` and `c_E` are CALIBRATION inputs** (pathA_40, the `Δr=2` result — each an independent calibration, not derived). Default `λγ = 1`, `c_E = c_γ`.
- Couplings/magnitudes: `G`, `Q_E`, `aT`, `aL`.
- EM-family shared block (reused across light/charge/magnetism, so they are mutually consistent by construction): `μ_R`, `ρ_br`, `B_eff = ρ_B0²/χ_c`.
- Drift trio `{ρ_B0, χ_c, C_hu}` = FREE-UNREDUCED (pathA_41) — pick values; label as un-reduced.

Docs for the cross-sector calibration status: `pathA_40_cone_lock.md`, `pathA_41_ng5_second_medium_drift.md`.

## 5. Honest-scope requirement — expose the departures

Each scene must be able to show the model's **characterized departure** from textbook GR/Maxwell as a **toggle or annotation** (per §3.x): the monopole/dipole gravity radiation residual, light's stray longitudinal mode, charge's Yukawa short-range term, magnetism's scalar-current admixture. **Do NOT hide these to make it look like clean GR + Maxwell** — surfacing them is the point (and the project's falsification ethos requires it). Label every value as DERIVED-form vs CALIBRATED-magnitude in the scene text.

## 6. Tests (HARD acceptance) — analytic golden checks

- **Charge:** two like charges repel, unlike attract; measured force-vs-distance slope ⇒ exponent p = 2 (within tol).
- **Gravity:** a 0-PN two-body orbit is closed (≈ zero precession over N periods); with 1PN on, the measured perihelion precession sign **and** magnitude match the PN formula to tolerance.
- **Light:** measured dispersion ⇒ `ω = c_γ k`; exactly 2 transverse polarizations; ray-bending deflection has the correct sign (toward the mass).
- **Magnetism:** two parallel currents attract; measured force-vs-distance ⇒ p = 2.
- **Purity:** a test asserts the `physics/` package imports **no** matplotlib/rendering module.
- **Determinism:** same params/seed ⇒ identical trajectories (bitwise or tight tol).

## 6.1 Numeric verification report (headless, text) — REQUIRED

The **primary human-verification surface** for this project. A headless mode (e.g. `python -m force_visualizer.report` or a `report.py`) that runs each sector's **physics core only (NO rendering)** and prints a structured, **deterministic, plain-text** transcript of the physical quantities — so correctness can be checked **from text alone, no screenshots**. Also write it to `output/verification_report.txt`.

Per sector, print **measured value + expected value + a ✓/✗**:
- **Gravity:** sampled orbit state + total-energy drift over N periods (0-PN should conserve); measured perihelion precession per orbit vs the 1PN prediction; (2.5PN on) measured inspiral rate vs Burke–Thorne.
- **Light:** measured dispersion `ω(k)` at several k vs `c_γ k`; count of transverse polarizations (expect 2, no longitudinal); lensing deflection angle (sign + scale).
- **Charge:** force magnitude at several separations → **fitted exponent p (expect 2)**; sign of force for like vs unlike pairs.
- **Magnetism:** force between parallel vs antiparallel currents at several separations → fitted exponent (expect 2) + sign (parallel attract); scalar-admixture contribution magnitude.

Also print, as labeled numbers, the **characterized-departure** quantities (gravity monopole/dipole residual amplitude; light stray-longitudinal presence; charge Yukawa short-range correction; magnetism scalar-current admixture), the **parameter block in use** (DERIVED vs CALIBRATED labels, values with units), and a final summary line (checks passed / total). Values with units, deterministic run-to-run.

## 7. Acceptance criteria (summary — all required)

1. All §6 tests pass.
2. Runs end-to-end and produces the four scenes' animations.
3. Each force law is faithful to its cited doc (Grok-verifiable against §3 docs).
4. Departures are toggle-able (§5); no scene hides them.
5. Every calibrated value is labeled with provenance (§4).
6. `physics/` has zero rendering dependencies; core is deterministic + unit-tested.
7. No dependency on `superfluid_lib`; no writes anywhere under `research/pde_ledger_v2/`.
8. The headless numeric report (§6.1) runs, prints per-sector measured-vs-expected values as deterministic plain text, and writes `output/verification_report.txt`.
9. The live interactive app (§9) builds and its headless smoke test passes; the GIF export, numeric report, and all prior tests still pass (the app is an additional front-end, it must not regress them).
10. Per-sector field visualizations (§10) render, each faithfully labeled for what it physically is; charge recolors by sign live and defaults to opposite signs; throats render as breathing; the field-direction golden checks pass; GIFs + report + live-app fix + all prior tests still pass.

## 8. Out of scope / do-not

- No full PDE; no emergence/derivation claim; no 3D; no web.
- Don't invent physics beyond the cited docs. If a doc is ambiguous or silent, **flag it** — do not guess or fill from `superfluid_lib` or general textbook EM/GR.
- Don't touch `research/pde_ledger_v2/` (read-only).

## 9. Live interactive app (deliverable alongside the GIFs)

A live, interactive matplotlib app (e.g. `app.py`, runnable as `python -m force_visualizer.app`) that runs the sectors in **real time** and lets the user manipulate them. It is a **second front-end over the SAME verified `physics/` core** — no re-implementation or re-derivation of any force law; the GIF export and the live app share one engine.

Requirements:
- **Reuse `physics/` + the scene logic.** Physics stays render-agnostic (still zero rendering imports in `physics/`, per §2/§6).
- **Live real-time animation** of each sector, with a way to switch sectors (tabs / selector / menu / small dashboard — your choice).
- **Interactive controls (matplotlib widgets):** sliders for each sector's key parameters (e.g. masses & eccentricity for gravity; k and c_γ for light; charge signs & Q_E for charge; `aL` and current directions for magnetism) **and a departures on/off toggle per sector** — the animation updates live as controls change.
- **Optional (nice-to-have, not required):** mouse interaction to place/drag particles or set initial velocities.
- **Must NOT regress** the GIF export (`render_all.py`), the numeric report (§6.1), or any existing test.
- **Backend:** use an interactive backend when a display is present; the app is meant to be **run locally by the user** (document the run command in the README). Fall back gracefully / exit with a clear message if run headless.
- **Headless acceptance test (since the build environment has no display):** add a smoke test (Agg backend) that **constructs** the app, wires the widgets, and **exercises the update callbacks** — a slider change and a departure toggle — against the physics core **without requiring a real display**, asserting the callbacks run and the underlying quantities change as expected. The user runs the actual interactive window locally to see/poke it.

## 10. Field visualizations + throat dynamics (per-sector medium fields)

Add a per-sector spatial field overlay so each scene shows its associated medium field, computed from the SAME verified `physics/` core (evaluate the force / potential / drain-velocity on a grid — **no new force laws**). Each field MUST be labeled truthfully for what it physically IS — do NOT render one sector's mechanism as another's:

- **Gravity — in-brane inflow streamlines (a real medium flow).** Streamlines/tracers of the drain velocity field: the medium flowing INTO the masses (each mass is a sink). This is the model's literal gravity picture ("the flow between the drains"), a **one-way** inflow — streamlines terminate at the masses; do NOT draw them circulating back out (the return is the separate global leakage, not a local loop). Animated inward tracers are ideal. Label: *medium inflow (gravity = the drain)*.
- **Charge — Coulomb field lines (a FIELD, not a flow).** Electric field lines / field-vector overlay of the two charges' Coulomb field, with correct like-repel / unlike-attract geometry. Label explicitly: *electric field (throat-body interaction) — a field, not a medium flow*. Do NOT render flowing/streaming tracers for charge (that would falsely imply the mass-flow picture the model rejects — a flow carrying energy is a mass current = gravity). Optionally also indicate the brane's displacement into `w` (the `h`-bulge) as a scalar height/color field around each throat.
- **Magnetism — circulating magnetic field lines on the brane (END-ON currents).** Draw the two currents **end-on** (pointing into/out of the viewing plane — e.g. dot = out-of-page, cross = into-page), so the **in-plane circulation** field around each is the *true* magnetic field of an end-on current (no out-of-plane hand-wave). The field grid must **follow the particles** (never decouple). **Animate the transverse force** — parallel currents attract (drift together), antiparallel repel (drift apart), driven by the verified core force — **not** in-plane streaming along a co-moving axis (which makes antiparallel charges diverge off the field grid, the bug being fixed). Label: *magnetic field on the brane (the moving throat's swirl, felt via localization)*; the literal swirl lives in the throat's 4D body.
- **Light — already the field** (the transverse shear wave); optionally annotate as the brane shear field. No force-law change.

Also:
- **Recolor particles by sign/direction, live.** Color each charge by its actual sign (red = positive, blue = negative), updating as the sign sliders move; **default the charge scene to OPPOSITE signs** (one positive, one negative) so it opens on the attract case. For magnetism, color by current direction.
- **Throat = a breathing object, not a fixed dot.** Render each particle/throat with a subtle breathing animation (a gently oscillating radius/ring) to reflect that the throat is a dynamical, breathing balance — nothing here is static. Keep it subtle; a visual cue, not a quantitative claim.

**Faithfulness guards (headless-checkable — REQUIRED):** add field-direction golden checks to the numeric report and/or tests, computed from the core and able-to-fail: gravity field at sample points points TOWARD the masses (inflow); the charge E-field points AWAY from a positive charge and TOWARD a negative one; the magnetism field CIRCULATES around a current with the sign set by the current direction. Must NOT regress the GIF export, the numeric report, the live-app retention fix, or any existing test; `physics/` stays free of rendering imports.

