## 0. Executive Summary

This session took the “Hutchison Effect” idea and forced it into a **tight, falsifiable mathematical box** using the superfluid-defect toy model you’ve been developing. The goal was not to argue that any controversial claim is true, but to ask a sharper question:

**If** a Hutchison-like phenomenology could arise from the toy model’s “vacuum acoustic / enthalpy” sector, **what would it have to look like in a reproducible lab test?**

We made one key modeling choice to keep things disciplined: **Model A (collective/effective mode)**. In Model A, the object’s response is dominated by a **single collective “breathing/compressional” mode** set primarily by its size (L), with a quality factor (Q), and driven by a scalar/enthalpy standing wave characterized by an effective pressure amplitude (p_A) and wave number (k). This gave us two channels with different signatures:

* **Regime I (ponderomotive / radiation-pressure channel):** a **time-averaged force** from gradients in (\langle p'^2\rangle), leading to a “weight-change / lift capability” threshold that is *quadratic* in drive amplitude and largely insensitive to fine frequency tuning.
* **Regime II (resonant softening channel):** a **narrowband internal resonance** where the object’s effective rigidity collapses when a strain proxy (\epsilon_{\mathrm{rms}}) crosses a critical value (\epsilon_{\mathrm{crit}}). This channel is *linear* in drive amplitude but amplified by (Q), and is sharply frequency-selective.

From these, the session produced the main deliverables:

1. **Closed-form thresholds and scalings (the “replicability spine”)**

* A size-set resonance:
  [
  f_0(L)=\frac{\chi c_*}{2L}\quad\Rightarrow\quad f_0\propto \frac{1}{L}
  ]
* A resonant softening threshold (on resonance):
  [
  p_{\mathrm{jelly,on-res}}\propto \frac{1}{Q}
  ]
* A lift/weight-change capability threshold (under the geometry lock (k\simeq \pi/L)):
  [
  p_{\mathrm{lev}}\propto \sqrt{L}
  ]
  These three relations became explicit **falsifiers**: if a dataset doesn’t follow them (within geometry factors), Model A is wrong.

2. **Phase maps in ((f,p_A))**
   We converted the model into **2D phase diagrams** separating regions where:

* nothing happens,
* only softening happens,
* only lift capability exists,
* both occur.

This directly answers the replication critique: the model predicts that the “effect” should live in **a narrow wedge around (f_0)** (Regime II), while the ponderomotive channel requires crossing a distinct amplitude threshold (Regime I). In other words, the toy model naturally supports “softening without levitation” as a clean regime.

3. **A coupled dynamical simulation (COM + internal mode)**
   We upgraded from purely analytic thresholds to a minimal ODE model coupling:

* center-of-mass motion in a ponderomotive potential, and
* an internal driven mode represented by slow quadratures (envelope dynamics).

We verified the analytic resonance rolloff against the dynamical envelope model, showing the narrowband structure is not an artifact of a single derivation style.

4. **Nonlinear “endgame” upgrade: Duffing + saturation + hysteresis**
   To model the “unreplicable / sometimes it works” complaint honestly, we introduced standard nonlinear resonance physics:

* **Duffing nonlinearity** (amplitude-dependent resonance shift),
* **saturating coupling** (prevents runaway),
* **hysteresis** via branch selection and continuation (up-sweep vs down-sweep differ).

This produced a new, extremely specific, testable signature:

> In certain amplitude bands, the outcome at the same ((f,p_A)) can depend on sweep history (up vs down), yielding a measurable hysteresis window.

In short: we ended the session with a compact mathematical framework that can generate **hard predictions** (threshold curves, scalings, hysteresis bands) that—if ever compared to real measurements—would either validate the model’s structure or falsify it quickly.

---

## 1. Scope, Framing, and Guardrails

### 1.1 Scope: what we are doing

We are treating the Hutchison Effect as a **thought experiment constraint problem**: take the toy model seriously on its own terms and ask whether it can produce Hutchison-like phenomenology (levitation/weight anomalies, rigidity loss, apparent “merging” behavior) in a way that is **quantitatively testable**.

The priority is *not* to preserve the claim. The priority is to produce:

* explicit parameter dependencies,
* predicted thresholds,
* narrowband/geometry conditions,
* “if–then” falsifiers.

This is the opposite of story-first reasoning: it is *boundary-first* reasoning.

### 1.2 Guardrails: what we are not doing

* We are **not** asserting that Hutchison reports are accurate.
* We are **not** asserting a working engineering pathway or giving instructions to create any apparatus.
* We are **not** claiming the toy model is “true physics.” It is a sandbox model whose value comes from the clarity of its predictions.

The output of this effort is best interpreted as a **hypothesis generator with built-in failure modes**: it tells you what would have to be observed if the model is relevant, and it tells you how to refute it.

### 1.3 Why controversial claims demand stricter testability

The major complaint about Hutchison-like reports is that they are “not replicable.” Within this project, we treat “not replicable” as a design requirement for the *model*:

If the toy model can explain anything here, it must explain it in a way that makes the phenomenon **more replicable**, not less—meaning:

* clear resonance conditions,
* clear geometry dependence,
* clear scaling with size and damping,
* clear sweep-direction effects if nonlinearity is essential.

So our standard for “progress” is: *the model compresses vague anecdotes into sharp boundaries.*

### 1.4 Modeling strategy: why we chose Model A

We chose **Model A (collective/effective mode)** because it is the smallest nontrivial model that can produce hard experimental-style predictions without needing detailed microphysics.

Model A assumes:

* The object responds predominantly through a **single collective mode** set by size (L).
* The drive can be represented as an effective scalar/enthalpy pressure amplitude (p_A) with a standing-wave geometry.
* The object’s “rigidity loss” can be represented by a strain proxy (\epsilon_{\mathrm{rms}}) crossing a threshold (\epsilon_{\mathrm{crit}}).

This is deliberately minimal. It buys three advantages:

1. **Closed-form scalings**: (f_0\propto 1/L), (p_{\mathrm{jelly}}\propto 1/Q), (p_{\mathrm{lev}}\propto \sqrt{L}) (under a specific geometry lock).
2. **Clean separation of channels**: ponderomotive vs resonant softening.
3. **Clear extension path**: once the minimal version works, you add realism (spatial mode structure, multiple modes, noise/drift, nonlinearities). We did exactly that by adding Duffing/saturation/hysteresis.

### 1.5 What counts as “derived” vs “assumed” in this session

**Derived within the model’s rules**

* The two-channel structure (ponderomotive vs resonance) and their different amplitude/frequency signatures.
* The phase-map construction in ((f,p_A)).
* The scaling relations as falsifiers (given Model A assumptions and geometry lock).

**Assumed / phenomenological inputs**

* The existence of an effective coupled density (\rho_*) and wave speed (c_*).
* The reduction of “rigidity loss” to a single threshold (\epsilon_{\mathrm{crit}}).
* The geometry lock (k\simeq \pi/L) as a simplifying, testable choice.
* The nonlinear form (Duffing + saturation) as a representative minimal nonlinearity.

These assumptions are not weaknesses if they’re explicit: they are knobs that an eventual dataset would either constrain or force us to replace (e.g., by moving to Model B/C).

---

## 2. Minimal Definitions and Dictionary

This section fixes notation and the “translation layer” between (i) the superfluid-defect ontology and (ii) the reduced Model A variables we actually simulate and map.

### 2.1 Effective field variables

We represent the driven “vacuum acoustic / enthalpy” sector using an effective scalar pressure-like fluctuation (p'(\mathbf r,t)) (units Pa). In Model A we do **not** resolve microscopic electrodynamics; we treat the drive as an effective compressive mode that can couple to a macroscopic defect aggregate (the sample).

We also use an enthalpy-like potential (h(\mathbf r,t)) (units m(^2)/s(^2)) related by an effective constitutive map
[
p'=\rho_*,h,
]
where:

* (\rho_*) is an **effective coupled density** for the mode that the sample “feels.” It is a phenomenological parameter to be inferred or constrained, not assumed from air or material density.
* (c_*) is an **effective propagation speed** for this coupled mode (the “sound speed” in the toy sector). Again, it is treated as a fit parameter in the thought model.

### 2.2 Sample variables

We model the test object as a single effective “lumped” oscillator characterized by:

* (L): characteristic size (for cubes, the side length; for other shapes, an equivalent dimension).
* (\rho_m): material mass density (e.g., aluminum (\rho_m\approx 2700) kg/m(^3)).
* (M=\rho_m V): total mass.

We intentionally collapse the object’s internal mechanics into a single collective coordinate (the “breathing”/compressional mode). Shape dependence is collected into a dimensionless factor (\chi=\mathcal O(1)).

### 2.3 Drive variables and standing-wave ansatz

We assume the drive produces an approximately standing field in one dimension (vertical (z) axis for convenience):
[
p'(z,t)=p_A\cos(kz)\cos(\omega t),
]
where:

* (p_A) is the **pressure amplitude** (Pa).
* (k=2\pi/\lambda) is the wavenumber.
* (\omega=2\pi f) is the drive angular frequency.

A crucial modeling choice in this session was the **geometry lock**:
[
\boxed{k \simeq \pi/L \quad (\lambda \simeq 2L)}
]
This is not claimed as universally correct; it is a *testable simplifying hypothesis* that ties the standing-wave geometry to the sample size so the scaling predictions become sharp.

### 2.4 Dissipation and resonance sharpness

We describe resonant selectivity by a quality factor (Q) and linewidth (\Delta f):
[
Q \equiv \frac{f_0}{\Delta f},
\qquad
\Delta f \simeq \frac{f_0}{Q}.
]
We use (Q) as an input knob because it is the most direct mathematical handle on “replicability difficulty.” If a claimed effect is real and resonant, its onset should vary strongly with (Q).

### 2.5 Internal response variable: strain proxy and threshold

Rather than modeling a full elastic tensor, Model A uses a single strain-like amplitude:
[
\epsilon_{\mathrm{rms}} ;;\text{(dimensionless)}.
]
In the envelope oscillator formulation, the internal displacement amplitude is (a) (peak), and we map it to a strain proxy by
[
\epsilon_{\mathrm{rms}} \approx \frac{a}{\sqrt{2},L}.
]

We then define a phenomenological “softening onset” threshold:
[
\boxed{\epsilon_{\mathrm{rms}} \gtrsim \epsilon_{\mathrm{crit}}.}
]
This (\epsilon_{\mathrm{crit}}) is **not** assumed to be a universal constant; it’s an effective parameter encoding “when rigidity meaningfully degrades” in this toy mapping.

### 2.6 Measurement proxy for “rigidity”

To connect the toy model to something that looks like data, we introduced a monotone “rigidity signal”:
[
\boxed{R(\epsilon)=\frac{1}{1+(\epsilon/\epsilon_{\mathrm{crit}})^2}.}
]
This gives:

* (R\approx 1): rigid/normal
* (R\ll 1): strongly softened
* The contour (R=0.5) corresponds to (\epsilon=\epsilon_{\mathrm{crit}}) and becomes our operational “threshold curve” in ((f,p_A)) maps.

### 2.7 Regime I vs Regime II: two channels

We explicitly separate two mechanisms:

**Regime I (ponderomotive / radiation-pressure channel):**
A time-averaged force on the center of mass from gradients in the intensity (\langle p'^2\rangle).

**Regime II (resonant softening channel):**
A narrowband internal resonance where (\epsilon_{\mathrm{rms}}) is amplified near (f_0), producing softening when (\epsilon_{\mathrm{rms}}) crosses (\epsilon_{\mathrm{crit}}).

This two-channel separation is central: it makes “softening without levitation” a natural and testable outcome rather than a contradiction.

---

## 3. Model A Core: Collective Mode Frequency and Geometry Locking

This section states the minimal Model A structure and extracts the scalings that later become falsifiers.

### 3.1 Collective resonance frequency

We postulate that the dominant internal response of the object is a single collective compressional mode with frequency set by size:
[
\boxed{f_0(L)=\frac{\chi,c_*}{2L}.}
]
Interpretation:

* (c_*) is the effective wave speed of the coupled mode in the toy vacuum sector.
* (\chi) collects geometry/mode-shape factors (order unity).

**Immediate falsifier:**
[
\boxed{f_0 \propto \frac{1}{L}.}
]
If you change the sample size and the resonance does not shift like (1/L), Model A is structurally wrong.

### 3.2 Bandwidth and the role of (Q)

With quality factor (Q), the resonance bandwidth is
[
\Delta f \simeq \frac{f_0}{Q}.
]
So as (Q) increases, the resonance becomes narrower and peak response grows. This alone provides a clean mathematical basis for “hard to reproduce”: small detunings in frequency can move the system from “on” to “off.”

### 3.3 Regime II (linear) amplitude scaling near resonance

In the linear envelope picture, the internal strain proxy near resonance takes the generic form
[
\epsilon_{\mathrm{rms}}(f,p_A)
\simeq \epsilon_{\mathrm{res}}(p_A),
\frac{1}{\sqrt{1+\left(2Q\delta\right)^2}},
\qquad
\delta\equiv \frac{f-f_0}{f_0}.
]
At exact resonance ((f=f_0)), the strain scales approximately as
[
\boxed{\epsilon_{\mathrm{res}}(p_A)\propto p_A,Q.}
]
Imposing the onset condition (\epsilon_{\mathrm{res}}=\epsilon_{\mathrm{crit}}) gives the on-resonance threshold
[
\boxed{p_{\mathrm{jelly,on-res}} \propto \frac{1}{Q}.}
]
This is our **Regime II falsifier**: the onset amplitude must fall like (1/Q) if the phenomenon is truly resonant.

### 3.4 Regime I (ponderomotive) force scale and geometry lock

From the standing-wave intensity,
[
\langle p'^2\rangle=\frac{p_A^2}{2}\cos^2(kz),
]
the gradient produces a time-averaged force density whose maximum acceleration scale (for the COM) is
[
a_0 \sim \frac{p_A^2,k}{4\rho_m\rho_*c_*^2}.
]
Define “lift capability” as (a_0/g). The “weight-change / levitation capability” threshold is the condition (a_0 \gtrsim g), giving
[
p_{\mathrm{lev}} \sim 2c_*\sqrt{\frac{\rho_m\rho_*g}{k}}.
]

Under the **geometry lock** (k\simeq \pi/L),
[
\boxed{
p_{\mathrm{lev}}(L)\propto \sqrt{L}.
}
]
This is our **Regime I falsifier** under the lock assumption. If the geometry lock is relaxed, the dependence shifts to (k) explicitly, which is why we keep this assumption visible: it’s a choice you can test.

### 3.5 Why Model A produces “strict testable boundaries”

With only these ingredients, Model A produces two different “activation rules”:

* **Resonant softening** is narrowband in frequency:
  [
  \epsilon_{\mathrm{rms}}(f,p_A)\ \text{peaks sharply at } f_0,\ \text{width }\Delta f\sim f_0/Q.
  ]
* **Ponderomotive lift capability** is mostly amplitude- and geometry-controlled:
  [
  a_0/g \propto p_A^2 k.
  ]

This yields a clean phase structure in ((f,p_A)):

* a narrow “wedge” of softening around (f_0),
* plus a separate high-amplitude threshold for lift capability,
* and a large region where nothing happens.

That is exactly the kind of structure that can turn a controversial anecdote into a falsifiable program: **change (L), change (Q), sweep (f) up/down**, and the model makes sharp, quantitative claims about what should move and how.

---

## 4. Regime I: Ponderomotive / Radiation-Pressure Channel

This regime models “weight anomaly / lift capability” as a **time-averaged** (cycle-averaged) force produced by spatial gradients in the **drive intensity**. The key point is that this channel is primarily **amplitude–geometry controlled** (quadratic in the drive), and only weakly dependent on being exactly on resonance.

### 4.1 Standing-wave intensity and the cycle average

We model the driven vacuum/enthalpy disturbance as an effective scalar pressure field
[
p'(z,t)=p_A\cos(kz)\cos(\omega t).
]
The cycle-averaged intensity is
[
\langle p'^2\rangle
===================

# \left\langle p_A^2\cos^2(kz)\cos^2(\omega t)\right\rangle

\frac{p_A^2}{2}\cos^2(kz),
]
since (\langle \cos^2(\omega t)\rangle=1/2).

### 4.2 Ponderomotive acceleration scale

In the toy superfluid mapping, gradients in the squared field produce an effective time-averaged “radiation pressure” push. In Model A we parameterize the resulting center-of-mass (COM) acceleration scale as
[
a_{\rm rad}(z);\sim;
+\frac{k}{4\rho_m\rho_*c_*^2},p_A^2,\sin(2kz),
]
i.e. the force is proportional to the gradient of (\cos^2(kz)). The peak acceleration scale is therefore
[
\boxed{
a_0 \equiv \max_z |a_{\rm rad}(z)|
==================================

\frac{p_A^2,k}{4\rho_m\rho_*c_*^2}.
}
]
This expression captures the two core features of the ponderomotive channel:

* **Quadratic drive dependence:** (a_0\propto p_A^2).
* **Geometry dependence:** (a_0\propto k), i.e. shorter spatial scale (larger (k)) increases the gradient and thus the force.

### 4.3 “Lift capability” and the operational threshold

We define the dimensionless lift capability as
[
\boxed{
\Lambda \equiv \frac{a_0}{g}.
}
]
Interpretation:

* (\Lambda<1): the available upward acceleration cannot exceed gravity, so sustained lift/weight-cancellation is not available in this simplified channel.
* (\Lambda\gtrsim 1): the ponderomotive channel is strong enough that, in principle, a configuration exists where upward radiation-pressure acceleration can balance (g).

The **Regime I threshold** is therefore the condition (\Lambda=1), i.e.
[
\frac{p_A^2,k}{4\rho_m\rho_*c_*^2} = g.
]
Solving for the amplitude gives
[
\boxed{
p_{\rm lev}(k)=2c_*\sqrt{\frac{\rho_m\rho_*g}{k}}.
}
]
This is the “weight-change / levitation capability threshold” used throughout our phase maps.

### 4.4 Geometry lock and the (\sqrt{L}) scaling

To convert (k) into a clean size dependence, we adopted the session’s simplifying, testable hypothesis:
[
k\simeq \frac{\pi}{L}.
]
Substituting into the threshold yields
[
\boxed{
p_{\rm lev}(L)
;=;
2c_*\sqrt{\frac{\rho_m\rho_*g}{\pi/L}}
;=;
2c_*\sqrt{\frac{\rho_m\rho_*g,L}{\pi}}
;\propto;
\sqrt{L}.
}
]
This is **Falsifier C** (under geometry lock): if an observed “weight anomaly threshold” does not scale like (\sqrt{L}), then either the geometry lock is wrong, or the ponderomotive mechanism as modeled here is not the right channel.

### 4.5 Why Regime I is not automatically “levitation in practice”

Even if (\Lambda\gtrsim 1), the actual COM motion depends on where the object sits relative to the standing-wave nodes, damping, and capture dynamics. Model A treats (p_{\rm lev}) as a **capability threshold**, not a guarantee that arbitrary placement yields levitation. This distinction matters:

* Regime I provides a *maximum available acceleration scale*.
* Whether the object is actually stabilized requires separate capture/trapping considerations (addressed later in the dynamical upgrade).

This is one reason the model can simultaneously predict “hard thresholds” and still allow for “finicky” behavior without abandoning physics: you can have a large capability but poor capture robustness if the node structure drifts or the damping is unfavorable.

---

## 5. Regime II: Resonant Softening / “Jellification” Channel (Linear Response)

This regime models the claimed “loss of rigidity” as a **narrowband resonant amplification** of an internal collective mode. The key point is that this channel is **frequency-selective** and strongly sensitive to (Q), creating strict test boundaries in ((f,p_A)).

### 5.1 Internal mode as a driven damped oscillator

Model A represents the object’s dominant internal response by a single coordinate (x(t)) (a collective compressional/breathing displacement). In the linear regime:
[
x'' + 2\gamma x' + \omega_0^2 x = F_{\rm drv}\cos(\omega t),
]
where:

* (\omega_0 = 2\pi f_0) is the collective mode frequency,
* (\gamma=\omega_0/(2Q)) is the damping rate (equivalently (Q) sets linewidth),
* (F_{\rm drv}) is the effective drive per unit mass (or per effective inertia), proportional to (p_A).

In our session’s reduced mapping, the drive amplitude scales as
[
F_{\rm drv}\propto \frac{1}{\alpha/\beta}\frac{p_A}{\rho_m L},
]
where (\alpha/\beta) is a dimensionless coupling/overlap ratio (kept explicit so it can be inferred rather than assumed).

### 5.2 Strain proxy and resonance enhancement

We map the internal amplitude into a dimensionless strain proxy
[
\epsilon_{\rm rms}\approx \frac{a}{\sqrt{2},L},
]
where (a) is the peak displacement amplitude of the internal mode.

For a driven damped oscillator, the steady-state amplitude has the standard Lorentzian rolloff. In a compact form we use
[
\boxed{
\epsilon_{\rm rms}(f,p_A)
=========================

\epsilon_{\rm res}(p_A),
\frac{1}{\sqrt{1+\left(2Q\delta\right)^2}},
\qquad
\delta\equiv\frac{f-f_0}{f_0}.
}
]
At exact resonance ((\delta=0)), the response scales like
[
\boxed{
\epsilon_{\rm res}(p_A)\propto p_A,Q,
}
]
which is the “gain” statement: increasing (Q) increases the peak internal strain linearly in this linear-response approximation.

### 5.3 Jellification onset condition and on-resonance threshold

We define “softening onset” by the criterion
[
\epsilon_{\rm rms}\gtrsim \epsilon_{\rm crit}.
]
At resonance this becomes
[
\epsilon_{\rm res}(p_A)\gtrsim \epsilon_{\rm crit}.
]
Using the proportionality derived and keeping the coupling ratio explicit, the session used the concrete mapping
[
\epsilon_{\rm res}(p_A)
\approx
\frac{p_A,Q}{(\alpha/\beta),\pi^2,\rho_m,c_*^2}.
]
Setting (\epsilon_{\rm res}=\epsilon_{\rm crit}) yields the on-resonance threshold:
[
\boxed{
p_{\rm jelly,on-res}
====================

(\alpha/\beta),\pi^2,\rho_m,c_*^2,\frac{\epsilon_{\rm crit}}{Q}.
}
]
This is **Falsifier B**:
[
\boxed{
p_{\rm jelly,on-res}\propto \frac{1}{Q}.
}
]
If the softening onset does not move like (1/Q) as damping changes, the “narrowband resonance-driven” hypothesis fails.

### 5.4 Off-resonance threshold curve (p_{\rm jelly}(f))

Using the Lorentzian rolloff, the amplitude required to reach (\epsilon_{\rm crit}) increases away from resonance. Solving
[
\epsilon_{\rm rms}(f,p_A)=\epsilon_{\rm crit}
]
gives the frequency-dependent threshold curve
[
\boxed{
p_{\rm jelly}(f)
================

# \frac{p_{\rm jelly,on-res}}{\sqrt{1+\left(2Q\delta\right)^2}}^{-1}

p_{\rm jelly,on-res},\sqrt{1+\left(2Q\delta\right)^2}.
}
]
Equivalently: the “jellification region” in ((f,p_A)) is a narrow wedge centered at (f_0), with width set by (\Delta f\sim f_0/Q).

This is the simplest mathematical encapsulation of “hard to reproduce” in a resonant theory: the effect is not merely high-threshold; it is **highly frequency-selective**.

### 5.5 Separation of regimes: softening vs lift

Because:

* Regime II is **linear in (p_A)** but enhanced by (Q),
* Regime I is **quadratic in (p_A)** and depends on geometry via (k),

Model A naturally predicts regions where:

* **softening can occur without lift capability** (on resonance but below (p_{\rm lev})),
* **lift capability can exist without strong softening** (very high (p_A) but detuned from (f_0)),
* and a narrow band where both occur.

This regime separation is a key output of the session: it produces a structured “phase space” with strict boundaries instead of a vague “anomaly.”

---

## 6. Analytic Phase Boundaries and “Strict Testable Boxes”

This section describes how we converted the two-channel Model A physics into explicit **phase diagrams** in the space of drive frequency and drive amplitude, ((f,;p_A)). The key outcome is that the model does not predict a vague “sometimes anomaly.” It predicts **well-defined regions** separated by hard boundaries, with simple scaling laws that shift those regions in controlled ways.

### 6.1 Why a phase plane is the right object

The replication complaint about Hutchison-style reports is essentially a statement about *parameter space*: people don’t know which combination of frequency content, geometry, power density, damping, and material state matters.

Model A addresses this by making the minimal parameter space explicit:

* **frequency** (f) (how close you are to a collective resonance (f_0)),
* **amplitude** (p_A) (how strong the effective scalar/enthalpy disturbance is),
* plus the sample parameters (L, Q) that determine where resonance lives and how narrow it is.

So the natural diagnostic object is a **2D map** in ((f,p_A)), with boundaries marking where each channel turns on.

### 6.2 Two boundaries define four regions

Model A defines two threshold surfaces:

**(i) Resonant softening threshold**
[
p_A \ge p_{\rm jelly}(f)
\quad\text{where}\quad
p_{\rm jelly}(f)=p_{\rm jelly,on-res}\sqrt{1+(2Q\delta)^2},
;;\delta=\frac{f-f_0}{f_0}.
]

**(ii) Lift/weight-change capability threshold**
[
p_A \ge p_{\rm lev}(k)
\quad\text{where}\quad
p_{\rm lev}(k)=2c_*\sqrt{\frac{\rho_m\rho_*g}{k}}
\quad (\text{or }p_{\rm lev}(L)\propto \sqrt{L}\text{ if }k\simeq\pi/L).
]

These two inequalities partition the plane into four distinct regimes:

* **Region 0 (none):** (p_A < p_{\rm jelly}(f)) and (p_A < p_{\rm lev})
  No softening, no lift capability.

* **Region 1 (softening only):** (p_A \ge p_{\rm jelly}(f)) but (p_A < p_{\rm lev})
  Narrowband rigidity loss possible without weight anomaly.

* **Region 2 (lift only):** (p_A < p_{\rm jelly}(f)) but (p_A \ge p_{\rm lev})
  Weight-change capability exists without strong resonant softening (typically detuned).

* **Region 3 (both):** (p_A \ge p_{\rm jelly}(f)) and (p_A \ge p_{\rm lev})
  Parameter overlap where both channels may appear.

This structure is one of the most important outcomes of the session because it gives you a clean way to interpret seemingly inconsistent reports. In this model, “I saw softening but no levitation” is not paradoxical—it corresponds to a well-defined region.

### 6.3 Shape of the softening boundary: the “narrow wedge”

Because (p_{\rm jelly}(f)) grows like (\sqrt{1+(2Q\delta)^2}), the softening region is always a narrow band around (f_0), with width on the order of
[
\Delta f \sim \frac{f_0}{Q}.
]

Two immediate, testable consequences follow:

1. **Narrowband requirement:** If an effect exists in this channel, a wideband source should produce it only fleetingly or intermittently (only when spectral content lands near (f_0)).
2. **Damping dependence:** Increasing (Q) should both narrow the band and lower the required on-resonance amplitude (p_{\rm jelly,on-res}\propto 1/Q).

This is the model’s built-in “replicability explanation” that does not rely on mystique: it’s just resonance.

### 6.4 How size (L) moves the map (three-size comparison)

We explicitly ran the analytic phase boundaries for three representative sizes:
[
L={5,;10,;20}\text{ mm}.
]

The map shifts in two distinct ways:

* The resonance center shifts as
  [
  f_0(L)\propto \frac{1}{L},
  ]
  so doubling (L) halves the predicted resonance frequency.

* Under the geometry lock (k\simeq \pi/L), the lift threshold shifts as
  [
  p_{\rm lev}(L)\propto \sqrt{L},
  ]
  so larger samples require higher (p_A) to reach the same lift capability.

This combination is particularly constraining: if you rescale the sample size and the “anomalous window” does not move in frequency like (1/L), Model A is wrong.

### 6.5 What makes these boundaries “strict and useful”

The strictness comes from the fact that each boundary is tied to a **different functional dependence**:

* softening boundary: narrowband, set by ((f-f_0)) and (Q),
* lift boundary: mostly amplitude + geometry, set by (k) (or (\sqrt{L}) under lock) and (\rho_*).

So the model gives you immediate “if–then” tests:

* If an effect is reported but does not sharpen with increasing (Q), the resonant channel is suspect.
* If the onset frequency does not move with sample size like (1/L), the collective-mode hypothesis is suspect.
* If “weight change” occurs but has no geometry dependence (via (k) or (\sqrt{L})), the ponderomotive mapping is suspect.

That’s the core deliverable: a controversial claim transformed into falsifiable constraints.

---

## 7. Upgraded Dynamics: Coupled COM + Internal Mode ODE Model

The analytic phase boundaries are powerful, but they do not address two practical issues that matter for interpreting “hard to reproduce” reports:

1. **capture/trapping:** even if a potential exists, does the object settle into it?
2. **coupling consistency:** does the narrowband envelope picture match a time-domain evolution?

So we upgraded to a minimal dynamical model that couples the object’s **center-of-mass motion** to its **internal mode amplitude**.

### 7.1 Why the upgrade was necessary

Analytic thresholds describe “capability” (e.g., (a_0/g)) and “onset” (e.g., (\epsilon_{\rm rms}=\epsilon_{\rm crit})). But in a standing-wave landscape, an object can still “fall through” nodes or fail to settle into a stable point depending on damping and initial conditions.

Because “replication” is fundamentally about *dynamics in parameter space*, we needed an ODE model that answers: *given (f, p_A, L, Q), what does the system actually do over time?*

### 7.2 State variables and decomposition

We modeled two coupled subsystems:

**(i) Center-of-mass vertical motion**

* position (z(t))
* velocity (\dot z(t))

**(ii) Internal collective mode envelope**

* quadratures (u(t)) and (v(t))

The quadratures represent the slow amplitude/phase evolution of the internal oscillator near resonance (rotating-wave / envelope approximation). The internal peak amplitude is
[
a(t)=\sqrt{u(t)^2+v(t)^2},
\qquad
\epsilon_{\rm rms}(t)\approx \frac{a(t)}{\sqrt{2},L}.
]

### 7.3 COM dynamics under the averaged ponderomotive force

Using the same standing-wave geometry, the COM experiences gravity plus an averaged radiation-pressure acceleration:
[
\ddot z = -g + a_0\sin(2kz) - \gamma_z \dot z,
\quad
a_0=\frac{p_A^2 k}{4\rho_m\rho_*c_*^2}.
]

This captures the key feature that the COM force is **spatially periodic** and depends on where the object sits relative to nodes/antinodes.

### 7.4 Internal mode drive depends on position

Because the standing wave has spatial structure, the internal mode is driven by the *local* pressure amplitude:
[
p_{\rm local}(z)\approx p_A\cos(kz).
]

So the internal quadratures evolve under a driven envelope system (schematically):
[
\dot u = -\gamma_{\rm int}u + \Delta v + (\text{drive})\cdot p_A\cos(kz),
]
[
\dot v = -\gamma_{\rm int}v - \Delta u,
]
with (\gamma_{\rm int}=\omega_0/(2Q)) and detuning (\Delta) set by (f-f_0).

This coupling is important: it creates a feedback loop where COM motion changes local drive strength, which changes internal excitation, which (in nonlinear extensions) can change effective response.

### 7.5 Cross-check: envelope ODE vs analytic rolloff

We explicitly compared the frequency response from:

* the analytic Lorentzian rolloff formula for (\epsilon_{\rm rms}(f)), and
* the numerically integrated envelope ODEs for the internal quadratures.

They matched in shape and bandwidth across representative amplitudes, confirming that the narrowband prediction (width (\sim f_0/Q)) is not an artifact of a single derivation route.

### 7.6 Capture and “punch-through” interpretation

The dynamical model clarified an important conceptual distinction:

* The analytic (p_{\rm lev}) is a **capability threshold**: it says the field can supply acceleration comparable to (g).
* But “levitation observed” depends on whether the object is **captured** into a stable node rather than passing through (“punch-through”).

Even with strong ponderomotive gradients, capture can fail if:

* damping is too small (object oscillates through nodes),
* the node pattern drifts,
* initial conditions inject too much kinetic energy.

This provides a clean toy explanation for erratic outcomes without adding mysticism: the phenomenon can be thresholded but still dynamically sensitive.

### 7.7 What the coupled model adds to the program

The coupled COM+internal model gives you three new, testable handles:

1. **time-to-settle:** how quickly the object approaches a stable location (a measurable timescale)
2. **position dependence:** outcomes should depend on node placement and spatial phase (\cos(kz))
3. **dynamic vs static signatures:** you can distinguish a static threshold effect from a transient crossing of a resonance window

This set the stage for the nonlinear “endgame” upgrades (Duffing shift, saturation, hysteresis), where the same ((f,p_A)) point can have multiple stable responses depending on history.

---

## 8. Nonlinear Endgame Upgrade I: Duffing + Saturating Coupling

The linear Regime II picture (a clean Lorentzian resonance) is already enough to produce strict, testable boundaries. But it does *not* capture a key empirical feature often associated with “hard to reproduce” resonance phenomena: **the response can distort, shift, or saturate at higher drive**, and the “onset curve” can become asymmetric.

To model that, we introduced two minimal nonlinear mechanisms that are standard in driven resonators:

1. **Duffing nonlinearity** (amplitude-dependent stiffness / resonance shift)
2. **Saturating coupling** (the effective drive stops increasing linearly with amplitude)

The purpose of this upgrade was not to chase complexity, but to add the *smallest* nonlinear structure capable of producing:

* shifted/tilted onset contours,
* non-symmetric response maps,
* and the conditions required for hysteresis (addressed in Section 9).

### 8.1 Duffing model for internal mode

We replace the linear internal oscillator with a Duffing oscillator:
[
x'' + 2\gamma x' + \omega_0^2 x + \beta x^3 = f_{\rm eff}\cos(\omega t),
]
where:

* (\gamma=\omega_0/(2Q)) as before,
* (\beta) is the Duffing coefficient (sign matters):

  * (\beta>0): **hardening** (resonance shifts upward with amplitude),
  * (\beta<0): **softening** (resonance shifts downward with amplitude).

In steady state, the Duffing response implies that the amplitude is no longer a simple Lorentzian in frequency: the resonance peak bends and can become multi-valued in certain parameter ranges.

### 8.2 Saturating coupling (drive soft-clipping)

We also introduced a saturating form for the effective drive:
[
f_{\rm eff}=\frac{f_0}{1+(a/a_{\rm sat})^2},
]
where (a) is the internal displacement amplitude and (a_{\rm sat}) is a saturation scale (chosen relative to (\epsilon_{\rm crit}) so it remains interpretable).

This does two things:

* prevents runaway response as (p_A) grows,
* produces a more realistic onset behavior where increasing drive eventually yields diminishing returns.

### 8.3 Steady-state amplitude equation and solution strategy

The Duffing steady-state response can be written as an algebraic equation for the amplitude (a). In the formulation we used, it becomes a cubic equation in (y=a^2) of the schematic form:
[
A^2 y^3 + 2AB y^2 + (B^2+C^2)y - f_{\rm eff}^2 = 0,
]
with (A\propto \beta), (B=\omega_0^2-\omega^2), and (C\propto \gamma\omega).

Because (f_{\rm eff}) itself depends on (a) via saturation, we solved this with a simple fixed-point outer loop:

1. assume an amplitude (a),
2. compute (f_{\rm eff}(a)),
3. solve the cubic for (y=a^2),
4. update (a) and iterate to consistency.

This gave us a fast way to populate full ((f,p_A)) response maps without time-domain integration.

### 8.4 Rigidity proxy response maps (R(f,p_A))

To make the nonlinear results “look like data,” we used the rigidity proxy introduced earlier:
[
R(\epsilon)=\frac{1}{1+(\epsilon/\epsilon_{\rm crit})^2}, \qquad
\epsilon_{\rm rms}\approx \frac{a}{\sqrt{2}L}.
]

We then generated response maps (R(f,p_A)) for three sizes:
[
L={5,;10,;20}\text{ mm},
]
plotting:

* the full heatmap of (R),
* the onset contour (R=0.5) (equivalently (\epsilon=\epsilon_{\rm crit})),
* and the lift threshold line (p_{\rm lev}).

**Key qualitative outcome:** the onset contour is no longer necessarily symmetric around (f_0). With (\beta>0) (hardening), the response can “lean” to higher frequency; with (\beta<0) it would lean to lower frequency. This produces a specific, falsifiable signature: **the onset band may shift with amplitude even at fixed (L)**.

### 8.5 What this upgrade adds to falsifiability

The nonlinear map adds two new test hooks beyond the linear model:

1. **Amplitude-dependent detuning:** the effective resonance frequency can shift as drive increases.
   If a real system shows “it only works when I tune a little above/below the naive resonance,” this is an exact mechanism for that.

2. **Onset contour distortion:** instead of a clean Lorentzian wedge, the region can form a bent “hook.”
   That is not handwaving; it is a geometrical feature of nonlinear resonance.

This creates a tighter prediction set than the linear model, not a looser one: once you include Duffing structure, the model predicts not just where onset occurs, but how that onset curve should deform as you vary (p_A).

---

## 9. Nonlinear Endgame Upgrade II: Hysteresis and Memory

After introducing Duffing nonlinearity, the model can enter a regime where, for the *same* ((f,p_A)), **multiple steady-state amplitudes are mathematically allowed**. This is the technical basis for hysteresis and “history dependence.”

We implemented hysteresis explicitly in the model and extracted a new high-value falsifiable signature: **sweep direction matters**.

### 9.1 Multistability and multiple amplitude branches

In Duffing systems, the steady-state amplitude equation can admit multiple real positive solutions for (a) at a given ((f,p_A)). Physically (in standard nonlinear oscillator theory), one typically finds:

* a “low-amplitude” branch,
* a “high-amplitude” branch,
* and an unstable intermediate branch.

In our implementation we do not assume the full stability analysis; instead, we treat the appearance of multiple roots as the mathematical marker that **branch selection** becomes relevant.

### 9.2 Branch selection as “up-sweep vs down-sweep”

To convert multistability into a prediction you can test, we implemented a continuation rule:

* **Up-sweep (in frequency):** start on the low branch at low (f) and track the solution closest to the previous amplitude as (f) increases.
* **Down-sweep:** start on the high branch at high (f) and track the closest solution as (f) decreases.

This is the standard numerical analog of how real nonlinear resonators exhibit hysteresis: the observed state tends to persist until the branch disappears (a jump).

### 9.3 Direction-dependent rigidity response (R(f))

Using the rigidity proxy (R), the model produces two distinct curves at the same (p_A):

* (R_{\uparrow}(f)): measured while sweeping up
* (R_{\downarrow}(f)): measured while sweeping down

In certain amplitude windows, these differ strongly over a finite interval in frequency, producing a **hysteresis window** where:
[
R_{\uparrow}(f)\neq R_{\downarrow}(f).
]

We quantified this with a simple scalar metric:
[
\mathcal{M}*{\rm hyst} \equiv \int |R*{\downarrow}(f)-R_{\uparrow}(f)|,df,
]
which becomes nonzero when “memory” is active.

### 9.4 Hysteresis bands in ((f,p_A)): direction-dependent thresholds

To make hysteresis even more “lab-test-ready,” we constructed **direction-dependent threshold curves** by defining onset as (R=0.5) and solving for the required (p_A) at each frequency:

* (p_{\rm th}^{\uparrow}(f)): threshold on up-sweep (low branch criterion)
* (p_{\rm th}^{\downarrow}(f)): threshold on down-sweep (high branch criterion)

Where these two curves separate, the model predicts a **hysteresis band**: a region in which the apparent onset depends on history.

This is a particularly strong falsifier because it is not a vague “it depends.” It is a concrete prediction:

> In the nonlinear regime, there exists a window in ((f,p_A)) where sweeping frequency upward and downward yields systematically different onset behavior, even if the final ((f,p_A)) point is the same.

### 9.5 Why hysteresis is an “endgame” signature for replication

Hysteresis is exactly the kind of feature that can make an effect look “unreplicable” to casual observers while still being fully deterministic:

* If one experiment “walks into” the region from one side of parameter space, it lands on one branch.
* Another experiment approaching from the other side lands on a different branch.
* Both think they used “the same settings” because they report the same final (f) and (p_A), but the system’s state differs because the path differed.

So the model turns “unreplicable” into a testable hypothesis:

* **If** the phenomenon is governed by nonlinear resonance in this way,
* **then** sweep direction and path history should predictably change outcomes,
* **and** the size/damping scalings from the linear backbone should still hold outside the nonlinear window.

This is the tightest “replication story” we produced in the session because it doesn’t relax falsifiability—it increases it by adding a new crisp experimental signature: **direction-dependent thresholds.**

---

## 10. Consolidated Predictions and Falsifiers

This section collects the session’s outputs into a small set of **hard predictions**. The point is to make Model A easy to kill (or constrain) with clean measurements, rather than allowing it to “explain anything.”

### 10.1 Falsifier A: size scaling of the internal resonance

Model A asserts a dominant collective mode with
[
\boxed{f_0(L)=\frac{\chi c_*}{2L};;\Rightarrow;; f_0\propto \frac{1}{L}.}
]
**Test implication:** change sample size (L) while holding geometry as similar as possible; the onset window in frequency must shift like (1/L). If it doesn’t, the collective-mode hypothesis (Model A) fails.

### 10.2 Falsifier B: damping (Q) scaling of resonant softening threshold

In the linear regime,
[
\epsilon_{\mathrm{res}}(p_A)\propto p_A Q,
]
so the on-resonance threshold for softening is
[
\boxed{p_{\mathrm{jelly,on-res}}\propto \frac{1}{Q}.}
]
**Test implication:** if (Q) increases (narrower resonance), the required on-resonance drive amplitude should fall linearly as (1/Q). If “softening” onset is insensitive to (Q), then the resonance-driven explanation is wrong (or (\epsilon_{\rm crit}) is not the right reduction).

### 10.3 Falsifier C: lift/weight-change capability scaling under geometry lock

For the ponderomotive channel, the peak acceleration scale is
[
a_0=\frac{p_A^2k}{4\rho_m\rho_*c_*^2}.
]
The capability threshold (a_0=g) gives
[
p_{\mathrm{lev}}(k)=2c_*\sqrt{\frac{\rho_m\rho_*g}{k}}.
]
Under the session’s geometry lock (k\simeq \pi/L),
[
\boxed{p_{\mathrm{lev}}(L)\propto \sqrt{L}.}
]
**Test implication:** if a “weight anomaly” threshold exists in this channel, it should show a strong, predictable dependence on geometry/length scale via (k) (or (\sqrt{L}) under lock). If no such dependence appears, this channel is not the right one.

### 10.4 Falsifier D: sweep-direction dependence (hysteresis window)

With Duffing nonlinearity + saturation, the steady-state amplitude can be multivalued and the observed response depends on history. That produces a direct experimental-style signature:

[
\boxed{\text{There exists a region where }R_{\uparrow}(f)\neq R_{\downarrow}(f)\text{ at the same }(f,p_A).}
]

Equivalently, the onset threshold in amplitude becomes direction-dependent:

* (p_{\rm th}^{\uparrow}(f)) for up-sweeps,
* (p_{\rm th}^{\downarrow}(f)) for down-sweeps,
  and these define a **hysteresis band**.

**Test implication:** if the model’s nonlinear regime is relevant, experiments should show systematic, repeatable differences between up-sweep and down-sweep scans—not random inconsistency.

### 10.5 “If–then” decision tree for Model A

A practical way to use these falsifiers:

1. **First check (f_0\propto 1/L).**
   If false → Model A fails at the foundation.

2. **Then check (\Delta f\sim f_0/Q)** (linewidth), and **(p_{\rm jelly}\propto 1/Q)** near resonance.
   If false → resonance-driven softening channel fails.

3. **Then check lift capability scaling with geometry** (via (k) or (\sqrt{L})).
   If false → ponderomotive channel fails.

4. **Only then check hysteresis.**
   If hysteresis is absent but (1–3) hold, the system may simply be in a linear regime.
   If hysteresis is present, it becomes a strong confirming signature of nonlinear resonance structure.

---

## 11. Parameter Inference Plan: How the Model Becomes “Hard Values”

Model A has a small number of effective parameters. The goal is not to guess them; it’s to show how they would be inferred from clean measurements if one ever had data consistent with the model’s structure.

### 11.1 Minimal observables to measure (conceptually)

To constrain Model A, you’d want measurements of:

* **Resonant frequency** (f_0) as a function of size (L).
* **Linewidth** (\Delta f) (or (Q=f_0/\Delta f)).
* An internal-response proxy (anything monotone in internal vibration / compliance) to define an onset contour analogous to our (R=0.5) line.
* Any COM-response proxy (position stabilization, apparent weight change, etc.) to compare to the (a_0/g) channel.

### 11.2 Inferring (c_*) and (\chi) from size scaling

From
[
f_0(L)=\frac{\chi c_*}{2L},
]
a fit of (f_0) vs (1/L) yields slope
[
m \approx \frac{\chi c_*}{2}.
]
If (\chi) is known/assumed from geometry, (c_*) follows. If not, the product (\chi c_*) is still tightly constrained.

### 11.3 Inferring the effective “softening parameter” ((\alpha/\beta)\epsilon_{\rm crit})

The on-resonance softening threshold is
[
p_{\rm jelly,on-res}=(\alpha/\beta),\pi^2\rho_m c_*^2\frac{\epsilon_{\rm crit}}{Q}.
]
So at fixed size (fixed (f_0)) and varying (Q), a plot of (p_{\rm jelly,on-res}) vs (1/Q) gives slope
[
m_J \approx (\alpha/\beta),\pi^2\rho_m c_*^2,\epsilon_{\rm crit}.
]
This means the model does not require knowing (\epsilon_{\rm crit}) and (\alpha/\beta) separately at first—only their product matters for the onset boundary.

### 11.4 Inferring (\rho_*) from the lift threshold scaling

From the lift capability threshold
[
p_{\rm lev}(k)=2c_*\sqrt{\frac{\rho_m\rho_*g}{k}},
]
if (c_*) is already constrained from (f_0(L)) and (k) is known (or locked), then (p_{\rm lev}) measurements can constrain (\rho_*). Under geometry lock (k=\pi/L),
[
p_{\rm lev}(L)=2c_*\sqrt{\frac{\rho_m\rho_*g,L}{\pi}}.
]
So (p_{\rm lev}^2) vs (L) is linear, and the slope infers (\rho_*).

### 11.5 Minimal dataset requirements

Conceptually, Model A can be meaningfully constrained with surprisingly small datasets:

* 3–5 sizes (L) to test (f_0\propto 1/L),
* a few controlled (Q) values (or effective damping states) to test (p_{\rm jelly}\propto 1/Q),
* a geometry scan (or at least a (k) estimate) to check the ponderomotive scaling,
* up/down sweeps to test hysteresis only if nonlinear behavior is suspected.

The point is that Model A is not flexible in how these dependencies should look—so even sparse data can be decisive.

---

## 12. Limitations, Non-uniqueness, and Next Modeling Extensions

This section makes explicit what was assumed, where the toy model is deliberately crude, and what the next modeling steps would be if you wanted to tighten the framework while keeping falsifiability.

### 12.1 What is phenomenological vs derived

**Derived (given Model A assumptions):**

* Two-channel structure (ponderomotive vs resonant softening).
* Phase plane partitioning into four regions.
* Primary scalings (f_0\propto 1/L), (p_{\rm jelly}\propto 1/Q), (p_{\rm lev}\propto 1/\sqrt{k}).
* Hysteresis as a consequence of Duffing multistability.

**Phenomenological inputs:**

* Effective coupling density (\rho_*) and propagation speed (c_*).
* Reduction of complex internal mechanics to a single collective coordinate.
* The single-parameter onset rule (\epsilon_{\rm crit}).
* Geometry lock (k\simeq \pi/L) as a simplifying test case.

These phenomenological elements are not hidden—each one is a knob that data would either constrain or force us to replace.

### 12.2 Sensitivity to geometry and standing-wave idealization

We assumed a clean standing wave with a single wavenumber (k) and ideal spatial dependence (\cos(kz)). Real systems would have:

* mode mixing,
* spatial inhomogeneity,
* multiple standing-wave components,
* drift in node location.

In Model A, that complexity compresses into “effective (k)” and reduced coupling ((\alpha/\beta)). That is useful for prediction, but it means the model should be tested first on the simplest geometries where (k) is plausibly well-defined.

### 12.3 What “softening” means here vs real materials

Our “jellification” is an **effective rigidity loss proxy**, not a microscopic theory of solids. The toy model replaces complicated material response with a single threshold in (\epsilon), and maps it to a signal (R(\epsilon)).

This is appropriate for a first-pass falsifiable framework, but it does not distinguish:

* shear vs bulk modulus changes,
* plastic deformation vs reversible compliance,
* thermal vs nonthermal mechanisms.

Those distinctions would matter in any real interpretation, but they are orthogonal to the model’s main job here: predicting strict resonance geometry/damping structure.

### 12.4 Candidate extensions that preserve falsifiability

If you want to upgrade without losing the “strict boundary” character:

1. **Multi-mode internal structure:**
   Replace the single oscillator with 2–5 modes. Prediction: multiple narrow onset bands at size-dependent frequencies.

2. **Stochastic drift / noise model:**
   Add slow drift in (f_0), (k), or phase to model day-to-day variability. Prediction: broadened onset probability but preserved scalings.

3. **Spatial field model:**
   Solve for a realistic standing-wave pattern (or multiple axes) and compute (\langle p'^2\rangle) gradients more faithfully. Prediction: lift capability becomes geometry-sensitive in a more detailed, but still testable way.

4. **Nonlinear mode coupling:**
   Add coupling between internal mode amplitude and effective detuning (beyond Duffing), producing richer hysteresis structure. Prediction: sweep-rate dependence emerges as an additional signature.

### 12.5 When to abandon Model A

Model A should be retired in favor of a different modeling choice if any of these happen:

* the dominant response is not size-set ((f_0) does not scale like (1/L)),
* onset is not narrowband (no (Q) dependence),
* thresholds do not show consistent geometry dependence,
* observed behavior is strongly broadband or dominated by non-resonant dissipative effects.

In that case, the correct response is not to “patch” Model A indefinitely, but to move to a different reduced model (e.g., explicit spatial modes, multiple internal degrees of freedom, or a different coupling mechanism).

---


## Appendix A. Key Equations (One-Page Reference)

### A.1 Standing-wave driver (Model A ansatz)

Effective pressure-like field:
[
p'(z,t)=p_A\cos(kz)\cos(\omega t),
\qquad \omega=2\pi f .
]
Cycle average:
[
\langle p'^2\rangle=\frac{p_A^2}{2}\cos^2(kz).
]

### A.2 Collective resonance frequency and linewidth

Dominant collective mode:
[
\boxed{f_0(L)=\frac{\chi c_*}{2L}}
\qquad\Rightarrow\qquad
\boxed{f_0\propto \frac{1}{L}}.
]
Quality factor and linewidth:
[
Q=\frac{f_0}{\Delta f},
\qquad
\boxed{\Delta f\simeq \frac{f_0}{Q}}.
]
Detuning:
[
\delta \equiv \frac{f-f_0}{f_0}.
]

### A.3 Regime I: ponderomotive / radiation-pressure acceleration scale

Peak acceleration scale:
[
\boxed{
a_0=\frac{p_A^2,k}{4\rho_m\rho_*c_*^2}
}
\qquad\Rightarrow\qquad
\Lambda\equiv \frac{a_0}{g}.
]
Capability threshold ((a_0=g)):
[
\boxed{
p_{\rm lev}(k)=2c_*\sqrt{\frac{\rho_m\rho_*g}{k}}
}.
]
Geometry lock used in session:
[
\boxed{k\simeq \frac{\pi}{L}}
\quad\Rightarrow\quad
\boxed{
p_{\rm lev}(L)=2c_*\sqrt{\frac{\rho_m\rho_*g,L}{\pi}}
\propto \sqrt{L}
}.
]

### A.4 Regime II (linear): resonant strain proxy and threshold

Strain proxy from peak internal displacement amplitude (a):
[
\boxed{\epsilon_{\rm rms}\approx \frac{a}{\sqrt{2},L}}.
]
Lorentzian rolloff (linear model):
[
\boxed{
\epsilon_{\rm rms}(f,p_A)
=========================

\epsilon_{\rm res}(p_A),
\frac{1}{\sqrt{1+(2Q\delta)^2}}
}.
]
On-resonance scaling:
[
\epsilon_{\rm res}(p_A)\propto p_A Q.
]
Session’s explicit mapping:
[
\boxed{
\epsilon_{\rm res}(p_A)
\approx
\frac{p_A,Q}{(\alpha/\beta),\pi^2\rho_m c_*^2}
}.
]
Jellification/softening onset:
[
\boxed{\epsilon_{\rm rms}\gtrsim \epsilon_{\rm crit}}.
]
On-resonance threshold:
[
\boxed{
p_{\rm jelly,on!-!res}
======================

(\alpha/\beta),\pi^2\rho_m c_*^2,\frac{\epsilon_{\rm crit}}{Q}
\propto \frac{1}{Q}
}.
]
Off-resonance threshold curve:
[
\boxed{
p_{\rm jelly}(f)=p_{\rm jelly,on!-!res},\sqrt{1+(2Q\delta)^2}
}.
]

### A.5 Rigidity measurement proxy used for response maps

Define a monotone “rigidity signal”:
[
\boxed{
R(\epsilon)=\frac{1}{1+(\epsilon/\epsilon_{\rm crit})^2}
}.
]
Operational onset contour:
[
\boxed{R=0.5 \iff \epsilon=\epsilon_{\rm crit}}.
]

### A.6 Nonlinear internal response: Duffing + saturation (steady state)

Duffing oscillator:
[
x'' + 2\gamma x' + \omega_0^2 x + \beta x^3 = f_{\rm eff}\cos(\omega t),
\qquad
\gamma=\frac{\omega_0}{2Q}.
]
Saturating coupling:
[
\boxed{
f_{\rm eff}=\frac{f_0}{1+(a/a_{\rm sat})^2}
}.
]
Steady-state amplitude equation becomes a cubic in (y=a^2):
[
A^2 y^3 + 2AB y^2 + (B^2+C^2)y - f_{\rm eff}^2 = 0,
]
with (A=(3/4)\beta), (B=\omega_0^2-\omega^2), (C=2\gamma\omega).

### A.7 Hysteresis implementation (branch selection / continuation)

When multiple real positive roots exist, define:

* **Up-sweep:** track the root closest to the previous amplitude starting from the low branch.
* **Down-sweep:** track the root closest to the previous amplitude starting from the high branch.

Direction-dependent thresholds from (R=0.5):

* (p_{\rm th}^{\uparrow}(f)) (up-sweep),
* (p_{\rm th}^{\downarrow}(f)) (down-sweep),
  defining a hysteresis band where they differ.

---

## Appendix B. Numerical Methods and Scripts

### B.1 Script inventory (this session)

We produced three Python artifacts, each with a specific role:

1. **Coupled dynamics (COM + internal envelope ODE):**
   `hutchison_modelA_upgrade.py`

* Implements a minimal time-domain model for:

  * COM vertical motion in the averaged ponderomotive landscape
  * internal collective mode via envelope quadratures ((u,v)) near resonance
* Used to:

  * cross-check analytic resonance rolloff
  * explore “capture/punch-through” behavior in a controlled toy setting
    **File:** [Download hutchison_modelA_upgrade.py](sandbox:/mnt/data/hutchison_modelA_upgrade.py)

2. **Fast nonlinear steady-state maps (Duffing + saturation):**
   `hutchison_modelA_nonlinear_steady.py`

* Implements:

  * Duffing nonlinear response as a cubic in (a^2)
  * saturating coupling fixed-point
  * rigidity proxy (R(\epsilon))
* Used to:

  * generate (R(f,p_A)) response maps quickly
  * extract onset contours (e.g., (R=0.5))
    **File:** [Download hutchison_modelA_nonlinear_steady.py](sandbox:/mnt/data/hutchison_modelA_nonlinear_steady.py)

3. **Endgame hysteresis (branch selection + continuation):**
   `hutchison_modelA_endgame_hysteresis.py`

* Implements:

  * multi-root Duffing solver
  * continuation sweeps in frequency
  * direction-dependent response (R_{\uparrow}(f)), (R_{\downarrow}(f))
  * direction-dependent threshold curves (p_{\rm th}^{\uparrow}(f)), (p_{\rm th}^{\downarrow}(f))
    **File:** [Download hutchison_modelA_endgame_hysteresis.py](sandbox:/mnt/data/hutchison_modelA_endgame_hysteresis.py)

### B.2 Numerical choices (what we actually did)

**Analytic phase boundaries**

* Computed (p_{\rm jelly}(f)) from the Lorentzian rolloff using (Q) and (f_0(L))
* Computed (p_{\rm lev}) from (a_0=g) using (k\simeq \pi/L)
* Classified each ((f,p_A)) point into 0/1/2/3 regions

**Envelope ODE cross-check**

* Integrated slow quadrature equations ((u,v)) (rotating-wave/envelope)
* Compared frequency response curves to analytic Lorentzian rolloff

**Nonlinear steady-state maps**

* Solved Duffing amplitude equation as a cubic in (y=a^2)
* Wrapped a fixed-point iteration to enforce saturation (f_{\rm eff}(a))
* Converted amplitude (a\to \epsilon_{\rm rms}\to R)

**Hysteresis**

* Per frequency point, solved for multiple amplitude branches (when present)
* Per sweep direction, used continuation to follow the nearest branch
* Compared (R_{\uparrow}(f)) and (R_{\downarrow}(f))
* Extracted direction-dependent threshold curves where (R=0.5)

### B.3 Parameter defaults used (session baseline)

These were the working defaults in the session (explicitly “toy”):

* (L \in {5,10,20}) mm (size sweep)
* (Q=300) (baseline)
* (\epsilon_{\rm crit}=10^{-4}) (toy softening threshold)
* (k\simeq \pi/L) (geometry lock)
* (c_*) and (\rho_*) treated as effective coupling parameters
* Duffing (\beta) chosen so nonlinear shift at (\epsilon_{\rm crit}) is comparable to half-linewidth (to make hysteresis demonstrable in the thought model)
* Saturation scale (a_{\rm sat}) set relative to (5,\epsilon_{\rm crit},L) (so saturation appears after onset)

(These are not claims about reality; they’re a consistent internal baseline to produce plots and scalings.)

### B.4 “Standard run” checklist (recommended plot set)

If you want a repeatable report pack for any new parameter set:

1. **Size scaling panel**

* plot (f_0) vs (1/L)
* plot (\Delta f) vs (f_0/Q)

2. **Linear phase boundaries**

* for each (L), plot (p_{\rm jelly}(f)) and (p_{\rm lev})

3. **Nonlinear response maps**

* (R(f,p_A)) heatmap + (R=0.5) contour + (p_{\rm lev})

4. **Hysteresis panel**

* for selected (p_A): overlay (R_{\uparrow}(f)) and (R_{\downarrow}(f))
* plot (p_{\rm th}^{\uparrow}(f)) and (p_{\rm th}^{\downarrow}(f)) (hysteresis band)

### B.5 Notes on reproducibility within the toy model

* The analytic boundary model should be identical run-to-run.
* ODE results depend on integration step and initial conditions, but qualitative structure (narrowband peak, (Q)-width scaling) should persist.
* Hysteresis depends on sweep direction by construction; it is itself a reproducibility test (run the same sweep twice and check the same branch is followed).
