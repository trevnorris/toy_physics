## 1) The core modeling goal

The overall objective of this sequence is:

1. **Specify a “throat”** as a localized, spherically symmetric *throughput* region that removes mass from the 3D medium at a rate ( \dot M ) (and possibly removes momentum/energy depending on model choices).

2. **Add a compensating refill** (uniform in space in the simplest tests) so that in a periodic box the net source integrates to zero. This is the operational stand‑in for “bulk leak‑in / dark energy” in the toy model: it keeps the 3D background from draining away while still allowing a persistent monopole‑like through‑flow. This framing is explicitly used in the toy model notes and the early throat scripts. 

3. **Show that the far‑field response is Newtonian‑like**, i.e.

* a potential that scales like ( \psi(r) \sim 1/r ) (so slope (\approx -1) in log‑log),
* and a field/acceleration ( g(r) \sim 1/r^2 ) (so slope (\approx -2) in log‑log),

and that this behavior **emerges** from the fluid equations in a controlled limit rather than being asserted by fiat.

A key refinement we learned along the way: in a driven, compressible medium with a steady sink+refill, the *most robust* “gravity analog” observable is **not** the raw pressure gradient, but the **longitudinal (potential) component** of the flow, represented by a scalar (\psi) where the longitudinal velocity is ( \mathbf v_L = \nabla \psi ). StepE makes that design choice explicit. 

---

## 2) What “Poisson emerges” means in this toy model

There are two distinct “Poisson” roles that showed up in your progression:

### A) “Asserting Poisson” (a regression / calibration test)

This is StepB: you *define* a Poisson sector (\Phi_P) from a prescribed source and then check that:

* it produces the expected radial slopes,
* your radial binning + fitting pipeline behaves,
* and that the “lag”/wave sector can be damped away cleanly.

This corresponds directly to the “instantaneous Newtonian/Poisson sector + optional lag sector” decomposition summarized in your toy model notes. 

### B) “Poisson as a static limit of dynamics” (the emergence bridge)

This is Steps C → D → E:

* **C:** solve Poisson for a potential (\psi) from the *sink source term*.
* **D:** evolve a (damped) acoustic/potential equation and show it relaxes to the Poisson solution.
* **E:** run full Euler and show the longitudinal component of the velocity field matches the Poisson predictor, i.e. that the Poisson behavior is the quasi‑static limit of the hydrodynamics.

The clean bridge derivation is the one you recorded in `stepC_bridge_derivation.md`: linearize continuity+momentum about ((\rho_0, \mathbf v=0)), assume potential flow, eliminate (\rho'), and you get a forced wave equation for (\psi). In the static limit, it becomes Poisson:
[
\nabla^2 \psi = \frac{S_\rho}{\rho_0}\quad (\text{with mean subtraction in a periodic box})
]


That derivation is exactly the non‑hand‑wavy story a referee will want: it’s a standard controlled limit (low Mach / quasi‑static / potential‑flow dominated), and it tells you precisely which assumptions you’re using and what breaks if those assumptions fail.

---

## 3) Step A: `single_throat_monopole.py` — why the “pressure → g_eff” route was fragile

This script was trying to do something physically intuitive but numerically unforgiving:

* Evolve a compressible fluid with a localized sink (the “throat”) and a compensating refill.
* Read off an “effective gravity” via a pressure‑gradient relation (hydrostatic‑looking),
  something of the form ( g_{\rm eff} \sim -\frac{1}{\rho}\frac{dP}{dr} ),
  and target the slopes:

  * (\Delta P \sim r^{-1}),
  * (g_{\rm eff}\sim r^{-2}).
    Those targets are stated right in the script’s intent. 

What we learned from the behavior you observed (NaNs, window sensitivity, “settles but wrong”, resolution dependence, slope drift):

### Why it was hard

1. **Pressure is not the clean far‑field variable** in a driven, compressible flow.
   Even if a far‑field (1/r) structure exists, (P(r)) is small in amplitude compared to the background, and its gradient is even smaller—so it’s easily swamped by:

   * standing acoustic modes,
   * discrete differentiation noise,
   * refill transients,
   * and the fact the system is not truly hydrostatic (it’s a *through‑flow* problem).
     StepE later explicitly documents this as the motivation for switching to the (\psi) diagnostic instead of pressure gradients. 

2. **The “fit NaN” failure mode** (like the early `dP_npts: 0`) is a diagnostic issue, not a physics issue:
   your fit is conditioned on thresholds / cutoffs / windows; if the data in that window falls below the cutoff or becomes sign‑mixed/noisy, the regression set can become empty, yielding NaNs. That showed up when the signal got too small or the window clipped into noise.

3. **Resolution sensitivity was real here** because you were differentiating noisy profiles and trying to infer slopes from tiny residuals. When you later compared (N=256) vs (N=512), the change in behavior is consistent with “the observable you’re extracting is near the numerical noise floor at low resolution,” not “the model breaks.” (The later (\psi)-based method is much more robust under the same grids.)

### Bottom line from Step A

Step A didn’t prove the model “wrong”; it told us that **pressure‑gradient‑as‑gravity** is the wrong *first diagnostic* for a throat‑driven through‑flow in a periodic, compressible box. The right variable to test first is the *potential‑flow* structure ((\psi)) that continuity directly ties to the sink.

(And that’s exactly what Steps C–E did.)

---

## 4) Step B: `stepB_scalar_poisson_lag.py` — “assert Poisson, measure lag, validate diagnostics”

This script is the cleanest “unit test”:

* Build a source distribution (a defect / density perturbation) and solve Poisson for (\Phi_P).
* Optionally evolve a damped wave equation for the total (\Phi=\Phi_P+\Phi_L), where (\Phi_L) is the “lag” piece.
* Measure:

  * the slope of (|\Delta \Phi|) (expect (-1)),
  * slope of (g=|\nabla \Phi|) (expect (-2)),
  * and the ratio (\mathrm{rms}(\Phi_L)/\mathrm{rms}(\Phi_P)) to quantify how much “non‑Poisson” content remains.

The script explicitly sets this up as:

* a Poisson sector plus a damped wave/lag sector, with damping (\gamma) acting as an absorber. 
  And it enforces periodic solvability by subtracting the mean / setting the (k=0) mode. 

**What you observed (and what it means):**

* Your slopes were basically stable from the start ((g) slope around (-2.035) in your logs), while the lag ratio dropped by orders of magnitude.
* Increasing (\gamma) changed *how fast* the lag decayed but not the final slopes.

Interpretation:

* **Good:** the Poisson sector is doing what it should; your slope fitting machinery is sane.
* **Expected:** (\gamma) is *not* supposed to change the Poisson solution; it’s there to damp the “lag” dynamics. Your observation that slopes didn’t change much with (\gamma) is the correct behavior.

This is exactly the “Poisson + lag sector” decomposition described in the toy model summary. 
(And per your earlier toolchain, I’m also referencing the marker you asked me to keep using.) 

---

## 5) Step C: throat Poisson / potential flow — why `dpsi_clean_slope` mattered

Step C moved from “defect mass” to “throat throughput”:

* You define a **sink source field** (S_\rho(\mathbf x)) representing the throat removing mass locally.
* In a periodic domain, you subtract the mean so the net source integrates to zero:
  [
  \nabla^2 \psi = \frac{S_\rho}{\rho_0} - \left\langle \frac{S_\rho}{\rho_0}\right\rangle
  ]
* Then (g = |\nabla \psi|) should scale like (r^{-2}) far from the throat.

Your first Step C output had:

* (g) slope (\approx -2.0367) (already good),
* but raw (d\psi) slope (\approx -1.318) (too steep).

The v2 script added a “clean” potential metric:

* Fit (\psi(r)) to a **monopole + constant**: (\psi \approx A/r + C),
* then define (d\psi_{\rm clean} = \psi - C),
* and measure the slope of (|d\psi_{\rm clean}|).

That’s why you got:

* `dpsi_clean_slope ≈ -0.997` while the raw slope looked off.

**Lesson:** In a periodic box (and with finite‑width sinks), constants and reference choices easily pollute raw (\Delta \psi) slopes. Removing the fitted constant is the right “physics‑invariant” way to measure the (1/r) falloff.

This is the first place the “throat → (1/r) potential” mapping becomes explicit and clean.

---

## 6) Step D: `stepD_superfluid_acoustic_throat.py` — Poisson as the static limit of a damped wave equation

Step D implemented the dynamical bridge from the derivation note:

Start with linearized continuity and momentum (isothermal, (P=c_s^2\rho)):
[
\partial_t \rho' + \rho_0 \nabla^2\psi = S_\rho,\qquad
\partial_t \psi + \frac{c_s^2}{\rho_0}\rho' = 0
]
Eliminate (\rho') to get the forced wave equation:
[
\partial_t^2\psi - c_s^2\nabla^2\psi = -\frac{c_s^2}{\rho_0}S_\rho
]
Add damping (\gamma \partial_t\psi) to absorb transients; the static limit gives:
[
\nabla^2\psi = \frac{S_\rho}{\rho_0}
]
with mean subtraction in a periodic box. 

**What your Step D run showed:**

* The lag ratio (\mathrm{rms}(\psi-\psi_P)/\mathrm{rms}(\psi_P)) dropped to (\sim 8\times 10^{-3}).
* `dpsi_clean_slope` stayed ~(-1) and `g_slope` ~(-2).

**Lesson:** even with a time‑dependent acoustic dynamics, the system is explicitly relaxing toward the Poisson predictor. And (just like in Step B) changing (\gamma) mostly affects the transient decay, not the target Poisson scaling.

This is the “referee‑friendly” bridge between:

* **Poisson as a static constraint**, and
* **Poisson as a limit of causal dynamics**.

---

## 7) Step E: `stepE_euler_throat_longitudinal.py` — the real “emergence” test (and why it bounced before settling)

Step E is the strongest step because it uses the **full nonlinear isothermal Euler system** with sink+refill, then checks whether the *longitudinal* component of the flow matches the Poisson predictor.

The key methodological upgrade is:

### Extract (\psi) from the flow via Helmholtz/FFT

In Fourier space, take the velocity field (\mathbf v(\mathbf k)) and compute the longitudinal potential:
[
\psi(\mathbf k) = -,\frac{i,\mathbf k\cdot \mathbf v(\mathbf k)}{|\mathbf k|^2}
]
so that ( \mathbf v_L = \nabla \psi ). StepE tracks `frac_long` as “how potential/irrotational is the flow.” 

### Compute a Poisson predictor (\psi_P) from the throat source

For periodic solvability, StepE uses mean‑zero handling (k=0 set to zero; refill handled as uniform), and builds:
[
\nabla^2 \psi_P = \frac{S_\rho}{\rho_0} - \left\langle \frac{S_\rho}{\rho_0}\right\rangle
]
Then it monitors:
[
\frac{\mathrm{rms}(\psi-\psi_P)}{\mathrm{rms}(\psi_P)}
]
as `lag_rms_over_psiP_rms`. 

### Why it bounced in your first StepE run

When you saw the slopes “look good” at one diag, then drift badly at the next, that’s exactly what you expect if the system is still carrying significant **acoustic ringing** (non‑steady compressible modes). In that case:

* (\psi) (from (k\cdot v)) includes time‑dependent wave content,
* the Poisson predictor (\psi_P) is quasi‑static,
* so the lag ratio stays high and the fitted slopes jump around depending on what phase of the ringing you sample.

StepE’s docstring calls this out: sustained sink through‑flow is not hydrostatic; you should diagnose the longitudinal potential directly. 

### What changed in your second StepE run (the “success” run)

With:

* `gamma_drag = 0.03` (momentum drag / damping),
* ramped throat (`ramp_time = 2.0`),
* refill at rest,
* and mean‑zero periodic solvability handling,

the system finally entered a stable quasi‑static regime.

Your final stable diagnostic at step 10200 showed:

* `frac_long ≈ 0.9999953` → the flow is overwhelmingly longitudinal/potential (very “superfluid‑like” in the sense relevant here),
* `lag_rms_over_psiP_rms ≈ 0.0316` → only a few percent residual “non‑Poisson / wave” content,
* `dpsi_clean_slope ≈ -0.995` and `g_slope ≈ -2.045` → essentially the target far‑field scalings,
* Monopole amplitude close to the expectation: StepE explicitly uses (A_{\rm expected}=\dot M/(4\pi\rho_0)) for (\psi\approx -A/r + C). 

**Lesson:** Once you (i) diagnose the *right variable* ((\psi) from (v_L)), and (ii) add enough damping to remove ringing, the throat‑driven Euler flow self‑organizes into the expected Poisson far‑field.

This is the cleanest evidence in the chain that “Poisson can emerge from superfluid‑like hydrodynamics” in the exact sense you need: it’s not asserted; it’s the quasi‑static limit of the driven Euler dynamics with a sink source.

---

## 8) Cross‑cutting lessons (physics + numerics)

### (i) Mean‑zero solvability is not optional in periodic gravity analogs

Every Poisson solve in a periodic box requires zero mean source; you handled that explicitly by subtracting means / killing (k=0). That is baked into StepB and StepE.

### (ii) “The slope is wrong” can be a *measurement* issue, not a field issue

Your StepC experience (raw (d\psi) slope off, clean slope right) is the canonical example. Constants and reference choices contaminate log‑log slope fits unless you remove them (monopole+constant fit). StepE includes that cleaning as a first‑class diagnostic (`dpsi_clean_slope`). 

### (iii) Pressure gradients were a poor early proxy for “gravity” in this setup

Not because “gravity can’t work,” but because:

* this is a **through‑flow** problem (sink + refill),
* compressibility creates ringing,
* and small pressure residuals are noisy to differentiate.
  StepE’s choice to use (\psi) instead is exactly the fix. 

### (iv) Resolution matters, but for understandable reasons

Your observation that 256 vs 512 gave materially different fitted slopes is consistent with:

* needing separation between the throat core ((\sigma)) and the fit window,
* needing enough bins / points in the far‑field region,
* and needing gradients that are well‑resolved.

The fact that the *(\psi)-based* diagnostics became stable at 512 is a sign you were moving out of the “numerical noise dominated” regime.

### (v) GPU implementation gotchas matter for diagnostics

Your CuPy `scatter_add` fix for spherical averaging is exactly the kind of thing that can silently corrupt diagnostics if mishandled (wrong return semantics / empty bins). This is a “keep it robust across backends” lesson, not a physics lesson—but it matters because the physics conclusions depend on the bin‑averaged radial profiles being correct.

---

## 9) Where the 4D throat / EM / “dark energy” pieces fit (and what we did *not* yet prove)

Nothing in Steps B–E required the *full* 4D implementation to show Poisson‑like scaling in the 3D effective theory.

What we effectively did is:

* treat the throat as a **parameterized sink** with parameters like (\dot M) and (\sigma),
* treat the “bulk” as an external reservoir providing a compensating refill (your “dark energy leak‑in” proxy),
* and test whether the 3D dynamics produce the Newtonian far‑field you want.

The 4D completion becomes necessary when you want to stop treating (\dot M), (\sigma), refill law, and coupling constants as knobs—and instead **derive** them from:

* a dynamical throat geometry,
* cavity/backreaction,
* and (in your original plan) EM as a projection of 4D interactions.

That roadmap (4D degrees of freedom + Maxwell + coupling structure) is what you’ve outlined in your next‑steps planning docs. 

So the right conclusion from this script chain is:

* We are **not** proving the full 4D + EM + throat‑geometry model end‑to‑end yet.
* We **are** demonstrating that the **Newtonian/Poisson sector you want** is internally consistent and can arise as a controlled limit of the driven hydrodynamics in 3D, provided you diagnose the correct field ((\psi)) and handle periodic solvability correctly.

That should reduce—not increase—the worry that “the papers must have fudged something.” The thing that was “wrong” was mainly the *early observable choice* (pressure‑gradient proxy) for this particular driven setup, plus some diagnostic fragility.

---

## 10) What you can credibly claim now (based on StepE’s stable run)

With the parameters that gave you the stable StepE output, you can now write (in referee‑language):

1. **Derivation:** In the linearized, low‑Mach, potential‑flow regime of isothermal Euler with a source term (S_\rho), the velocity potential obeys a forced wave equation whose static limit is Poisson:
   (\nabla^2\psi = S_\rho/\rho_0) (mean‑subtracted for periodic solvability). 

2. **Numerical verification:** In a periodic 3D domain with a localized sink (throat) and uniform compensating refill, the measured longitudinal component of the velocity field is overwhelmingly potential (`frac_long→1`) and the extracted (\psi) matches the Poisson predictor with small residual lag (`lag_rms_over_psiP_rms ≪ 1`). 

3. **Far‑field scaling:** After removing an additive constant (monopole+constant fit), the potential falls as (1/r) and its gradient as (1/r^2) over the fit window, consistent with the Newtonian monopole. 

That is, in my view, a coherent “we learned something real” endpoint for this chain.

===

## Where we are now

### What the scripts established

**StepB: scalar-field Newtonian sector sanity check**

* We implemented a split “Poisson + finite-speed wave” scalar structure:

  * Instantaneous Poisson reference: (\nabla^2\Phi_P = 4\pi G\rho)
  * Total potential evolves with a sourced wave equation, with optional absorber:
    [
    \partial_t^2 \Phi + 2\gamma,\partial_t\Phi = S^2\left(\nabla^2\Phi - 4\pi G\rho\right)
    ]
  * Define lag field (\Phi_L\equiv \Phi-\Phi_P); on a periodic box, (\gamma>0) is the “radiate away” surrogate so the Poisson state is an attractor.

This is **not** the “superfluid derivation” yet — it’s just a controlled regression test: “if we build the scalar sector this way, do we measure (1/r) and (1/r^2) with our diagnostics?”

---

**StepC: throat-as-sink ⇒ Poisson for the velocity potential (\psi)**

Here we directly test the core identity that later becomes “gravity emerges”:

1. Start from continuity with a mass source/sink (S_\rho):
   [
   \partial_t\rho + \nabla\cdot(\rho\mathbf v)=S_\rho
   ]

2. In the **low-Mach, quasi-steady, nearly-uniform** regime, take (\rho\approx\rho_0) and (\partial_t\rho\approx 0), giving:
   [
   \nabla\cdot \mathbf v \approx \frac{S_\rho}{\rho_0}
   ]

3. If the flow is **irrotational** ((\mathbf v=\nabla\psi)), then:
   [
   \nabla^2\psi \approx \frac{S_\rho}{\rho_0}
   ]
   On a periodic FFT box we enforce solvability by subtracting the mean:
   [
   \nabla^2\psi = \left(\frac{S_\rho}{\rho_0}\right) - \left\langle\frac{S_\rho}{\rho_0}\right\rangle
   ]
   (exactly what your init printed). 

StepC_v2 additionally cleaned up the “constant offset” issue by reporting `dpsi_clean_slope` (reference-subtracted), which is the physically meaningful (1/r) metric.

---

**StepD: dynamic bridge (damped acoustic / wave relaxation) ⇒ Poisson static limit**

StepD is the “referee bridge” between StepB’s “assume Poisson+wave” and StepE’s “derive Poisson behavior from hydro”:

* It evolves a **damped acoustic-potential** equation in Fourier space (a wave-like evolution for (\psi)) and continuously compares it to the Poisson predictor (\psi_P), monitoring (|\psi-\psi_P|/|\psi_P|). 

This shows: if you start with a propagating scalar (sound/phase mode) and add mild damping, the system relaxes onto the Poisson solution sourced by the sink.

---

**StepE: full isothermal Euler + throat sink/refill ⇒ Poisson emerges in the longitudinal sector**

This is the big one.

* You evolve **isothermal Euler** (so (P=c_s^2\rho)) with:

  * A localized sink implementing the throat (Gaussian kernel width (\sigma))
  * A uniform refill chosen so the net source is mean-zero (periodic mass conservation + Poisson solvability)
  * A small drag (`gamma_drag`) to damp long-lived acoustic ringing.

* Diagnostics:

  1. **Helmholtz projection** in k-space to measure how longitudinal the velocity field is:

     * Build ( \mathbf v_k)
     * Longitudinal component via projection onto (k):
       [
       \mathbf v_{L,k} = \frac{\mathbf k(\mathbf k\cdot \mathbf v_k)}{k^2}
       ]
     * `frac_long` is the fraction of kinetic energy in (\mathbf v_L).
  2. Extract the scalar potential (\psi) from the longitudinal component:
     [
     \psi_k = -,i,\frac{\mathbf k\cdot \mathbf v_k}{k^2}
     ]
     (with (k=0) handled safely).
  3. Compute a **Poisson predictor** (\psi_P) from the source (S_\rho):
     [
     \psi_{P,k} = -,\frac{(S_\rho/\rho_0)_k}{k^2}
     ]
     (again, (k=0) set to 0).
  4. Compare:

     * Far-field slope of (\Delta\psi) (reference-subtracted): target (-1)
     * Far-field slope of (g=|\partial_r\psi|): target (-2)
     * Lag ratio (|\psi-\psi_P|/|\psi_P|): target small
     * Monopole amplitude: (A_{\rm expected} = \dot M/(4\pi\rho_0)) and compare to fit (A_{\rm fit}).

Your final “settled” run hits all of these in the expected ballpark, so StepE is a **successful emergence demonstration**: the sim does **not** enforce Poisson as an evolution law; Poisson shows up as a diagnostic identity of the longitudinal, quasi-steady limit of continuity + irrotationality.

---

## So… do we need StepF?

If “StepF” was intended as “do *even more* variants of single-throat steady state,” then I agree with your instinct: you’ve basically proven the single-source Newtonian sector already.

But there are still **three** next steps that are worth doing **before** 4D, because they turn this from “nice internal check” into “referee-resistant evidence” and they become regression tests for the 4D program.

---

## Next steps before full 4D

### Next Step 1 — Convergence & robustness sweep (make it undeniable)

**Goal:** Show the emergent Poisson behavior is not a lucky parameter choice.

**What to do (Python sims):**
Run StepE (or a copy “StepE_sweep”) across a small grid of:

* Resolution: (N\in{256, 384, 512, 768})
* Box size (L): e.g. (L\in{150, 200, 300})
* Throat width (\sigma): e.g. ({1,2,3,4})
* Drag (\gamma_{\rm drag}): e.g. ({0, 0.01, 0.03, 0.1})
* Ramp time: ({0, 1, 2, 5})

**What to record (per run):**

* `dpsi_clean_slope`, `g_slope`
* `A_fit/A_expected` and `Mdot_from_flux/Mdot_target`
* `lag_rms_over_psiP_rms`
* `frac_long`
* mass drift and mean momentum (should be ~0)

**Acceptance targets (suggested):**

* (|dpsi_clean_slope+1| \lesssim 0.05)
* (|g_slope+2| \lesssim 0.05)
* (|A_{\rm fit}/A_{\rm expected}-1| \lesssim 0.05)
* (|\dot M_{\rm flux}/\dot M_{\rm target}-1| \lesssim 0.05)
* `frac_long > 0.999` (or whatever you decide is “clean enough”)
* lag ratio (\lesssim 0.05) at late time (you got ~0.03 in the stable run)

**Deliverable:** one summary JSON/CSV + a few plots (slope vs dx, amplitude error vs dx).
This is the single best “pre-4D” thing you can do to make the Poisson-emergence claim hard to dismiss.

---

### Next Step 2 — Superposition and multipoles (two throats)

Right now you’ve proven: “one throat produces a monopole-like (1/r) potential in the longitudinal sector.”

A referee will immediately ask: “Does it superpose linearly?” Because Newtonian gravity is linear in the source.

**Goal:** Demonstrate that multiple throats produce a potential that is (approximately) the sum of single-throat potentials in the far field.

**What to do (Python sims):**
Create **two sink kernels** (and corresponding refills) separated by distance (d) (e.g. (d\sim 30)–(60) in your units). Options:

* Equal sinks: (\dot M_1=\dot M_2)
* Unequal sinks: (\dot M_2 = 2\dot M_1)

**Diagnostics to add:**

* Compute (\psi) extracted from velocity exactly as StepE does.
* Compute Poisson predictor (\psi_P) from the combined source (already the same formula).
* Instead of a single spherical average about the box center, do:

  * Spherical average around throat #1
  * Spherical average around throat #2
  * And a **line probe** along the axis connecting them (this is often more informative than radial averages in multi-source problems)

**Acceptance tests:**

* Far from both throats, (\psi \approx \psi_1+\psi_2) within a few percent
* Along the connecting axis, the gradient matches the expected “two monopoles” profile
* `frac_long` remains ~1 if you keep the system in the same irrotational regime

**Deliverable:** a “superposition plot” showing (\psi_{\rm measured}) vs (\psi_{P}) vs (\psi_1+\psi_2).

This is arguably *the* last “Newtonian-sector” check you want before 4D.

---

### Next Step 3 — Turn the field into an actual force law (test particles or mobile throats)

StepE currently demonstrates the *field* exists. The next credibility step is: does it *act like gravity* on something?

There are two levels:

#### 3A. Passive tracer test (quick and clean)

* Add “test particles” that experience acceleration:
  [
  \ddot{\mathbf x} = -\alpha,\nabla\psi(\mathbf x)
  ]
* Here (\alpha) is your “dictionary constant” converting (\psi) units into whatever you want to call “gravitational potential.” The exact value is not important yet; what matters is qualitative behavior (e.g., inverse-square attraction, stable circular orbit for tuned velocity).

This can be done without changing the fluid; you just interpolate (\nabla\psi) from the grid.

#### 3B. Two mobile throats (harder, but directly on-model)

* Make each throat center a dynamical DOF.
* Let them move under the force inferred from the local field (or from momentum-flux / impulse accounting).

This is more faithful, but you can postpone it until after 4D if you want. The tracer test already shows the field is usable as “gravity.”

**Deliverable:** a simple orbit or scattering demo that is clearly (1/r^2)-like.

---

## Optional (but important if you want to preserve the older “pressure = gravity” intuition)

Your earlier “single_throat_monopole” experiments were fitting slopes of (\Delta P) and (g_{\rm eff}) derived from pressure/enthalpy gradients, and you got confusing slopes (often trending toward 0, or weird steep values).

From the current sequence, the clean statement is:

* **The thing that naturally obeys Poisson in this sink-driven, irrotational regime is the velocity potential (\psi)** (extracted from the longitudinal velocity), not the pressure field directly.

If the papers (or your preferred ontology) want “gravity = pressure/enthalpy gradient,” then you should explicitly re-derive the mapping you want (e.g., via Bernoulli / enthalpy relations in the appropriate limit) and confirm it numerically.

I’d treat this as **optional pre-4D** because:

* it’s more about matching legacy interpretation,
* whereas the (\psi)-Poisson emergence is already a clean hydro identity.

But if you don’t reconcile it, you risk confusion later when you try to “connect” this to the n=5 GNLS/hydro in 4D (where enthalpy matters a lot). 

---

## Implementation notes to carry forward

### 1) The “constant offset” / reference subtraction is not cosmetic

Any Newtonian-like potential is defined only up to an additive constant, so you **must** fit slopes using a reference-subtracted quantity (`dpsi_clean`, `dPhi`, etc.). StepB explicitly does that by using an outer-shell reference.
That was the main reason early fits could look “wrong” even when the physics was fine.

### 2) Periodic Poisson solvability requires mean-zero source

Every time you compute a Poisson predictor on a periodic grid, you need
(\langle \text{rhs}\rangle=0). StepB does this for (\rho) by subtracting (\langle\rho\rangle).
Your throat scripts do the same conceptually for (S_\rho/\rho_0). 

Interpretation-wise, the uniform refill is “dark-energy-like background bookkeeping,” not a local field source.

### 3) CuPy radial binning: `cupyx.scatter_add` return value is version-dependent

You found a real interoperability hazard: some CuPy versions return `None` from `cupyx.scatter_add` because it’s in-place. The StepC code as written assumes it returns the array, which can break. 
Your patch (allocate `sums`, then `cupyx.scatter_add(sums, ...)`, then divide) is the robust pattern — keep that as a utility function and reuse it everywhere.

---

## When it’s reasonable to move on to full 4D

Once you have:

1. **Single throat** steady-state Poisson emergence (done; StepE success),
2. **Convergence sweep** showing it doesn’t depend on one lucky (N)/(L)/(\gamma_{\rm drag}),
3. **Two-throat superposition** showing linearity in the far field,

…then you’ve basically maxed out what 3D periodic-box hydro can do for the *Newtonian scalar sector*.

At that point, the remaining open questions are genuinely “4D questions”:

* How the throat geometry ((a,L)) is defined and evolves self-consistently.
* How “brane observables” are projected from a unified 4D field with confinement.
* How the effective mouth operator (Z^{\rm eff}(\omega)) emerges (instead of being imposed).
* How the fixed microphysics constraint (the stiff (n=5) EOS) changes the story compared with the isothermal placeholder used in StepE. 
* How EM arises from the transverse sector / cavity selector, and how geometry picks (L/a\approx 1.85) in the relevant regime. 

Those are exactly the “hard-mode unified 4D” deliverables already laid out: single 4D PDE, brane emerges via confinement in (w), geometry as a potential (V_{\rm conf}(\mathbf X; a,L)), and extracting an effective response operator operationally rather than by boundary matching.

---

## Concrete “new session” task list

If we start a new chat, here’s the exact checklist I’d want to pick up with:

### A. Package the current 3D result as a regression test

* Freeze one “golden” StepE parameter set (like the one that settled at step ~10200).
* Add a small script that runs it and asserts thresholds on:

  * `dpsi_clean_slope`, `g_slope`, `lag_rms_over_psiP_rms`, `frac_long`, `A_fit/A_expected`, `Mdot_from_flux`
* Purpose: you’ll use this as a regression when you later modify anything.

### B. Write `stepE_sweep_convergence.py`

* Runs the small sweep over (N,L,\sigma,\gamma_{\rm drag},) ramp time.
* Outputs one JSON line per run.
* Postprocess to plots.

### C. Write `stepE_two_throats_superposition.py`

* Two sinks at known separation.
* New diagnostics: line probes + per-throat radial averages.
* Show superposition and far-field slopes.

### D. Write a tiny “dictionary” note (1–2 pages)

* What is the gravitational potential in this hydro limit?

  * (\psi) (velocity potential) is the clean Poisson object here.
* What is the gravitational acceleration?

  * (\mathbf g \propto -\nabla\psi).
* How does “mass” relate to throat parameters?

  * Monopole coefficient (A \sim \dot M/(4\pi\rho_0)) (already in StepE’s monopole fit logic).
* What assumptions are required?

  * subsonic, nearly uniform (\rho), irrotational flow, quasi-steady, small drag to remove standing waves.

This “dictionary page” becomes the bridge to the 4D derivations (where (\dot M), (a), (L), and confinement determine the effective source seen on the brane).

### E. Only then: begin the 4D derivations/sims

Use the “hard-mode” plan:

* Unified 4D model in ((x,y,z,w)); brane emerges via confinement in (w); throat is a localized relaxation of confinement with geometry DOFs ((a,L)).
* Keep the (n=5) EOS as non-negotiable microphysics when you switch from isothermal placeholders. 
* Define and measure an operational (Z^{\rm eff}(\omega)) on/near the brane rather than imposing a boundary impedance.

---

If you want, in the next session we can start directly with **B + C** (convergence sweep + two-throat superposition), because those are the last “prove it” steps that will pay off massively once you jump into the much more coupled 4D geometry work.

