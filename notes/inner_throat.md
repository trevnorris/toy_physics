Level A is a really good target: **solve the 4D “inner throat” problem locally**, and **match it onto the 3D “outer brane” effective fields**. That *does* give you a genuine 4D model of what’s happening in/near the throat, without committing to a full global 4D simulation.

Below is a concrete roadmap for the next paper, built to “close the missing piece” you’re pointing at: **a clean near–far matching layer** that quantifies *how much* of the 4D interior is visible in 3D, including near black-hole-like (acoustic horizon) regimes.

---

## What the next paper should accomplish (numbered, pass/fail style)

1. **Define a single throat geometry family** that resolves the “sphere–cylinder tension”

   * Outer mouth looks spherical at large (r) (needed for the 1PN gravity sector), but transitions to a cylindrical interior (needed for the EM cavity modes and (L/a\simeq 1.85) selection).
   * Deliverable: a parametric “rounded funnel” (you already used this idea in the ontology) with a small set of shape parameters.

2. **Write the inner 4D PDE problem explicitly (the throat-local model)**

   * Minimal inner model: 4D acoustic/enthalpy wave equation in the bulk throat region, with clear boundary conditions (walls, mouth, bottom).
   * Deliverable: a well-posed BVP/IBVP on the throat manifold (\mathcal{T}).

3. **Compute the throat mode spectrum and identify the EM-relevant branch**

   * Show the cylindrical/Bessel separation inside the throat (your ontology already points to this).
   * Deliverable: eigenfrequencies (\omega_n(L/a,\text{shape})), mode shapes, and the “selected” (L/a\approx 1.85) mechanism as an actual minimization/variational result (not just narrative).

4. **Derive a Dirichlet-to-Neumann (DtN) / impedance map at the mouth**

   * This is the key “projection” math: the 4D interior becomes a boundary operator that the 3D brane solver can use.
   * Deliverable: a frequency-dependent boundary condition at the mouth, e.g.
     [
     \partial_n h\big|*{\text{mouth}} = \mathcal{Z}(\omega), h\big|*{\text{mouth}}
     ]
     (or the appropriate pair of conjugate variables in your chosen formulation).

5. **Re-derive the outer 3D far-field multipole expansion from the matched geometry**

   * Show (again, but now as part of a systematic matching expansion) that the brane far-field is monopolar with corrections suppressed by ((a/r)^2).
   * Deliverable: explicit coefficients (\alpha_\ell) from the *same* throat geometry family; reproduce the small quadrupole like the ontology’s (\alpha_2\sim 10^{-2}) example.

6. **Define a clean “near–far mixing diagnostic” (“the gradient”)**
   One simple, brutally clear diagnostic (gravity sector) is:
   [
   \mathcal{R}*\text{geom}(r) \equiv
   \frac{|\Phi*{\text{non-spherical}}|}{|\Phi_{\text{monopole}}|}
   \approx |\alpha_2|\left(\frac{a}{r}\right)^2
   \quad\Rightarrow\quad
   \frac{d\ln\mathcal{R}*\text{geom}}{dr}=-\frac{2}{r}.
   ]
   For EM, you do the analogous thing with the DtN map:
   [
   \mathcal{R}*\text{EM}(r,\omega)\equiv
   \frac{\text{power/energy in throat modes}}{\text{power/energy in brane field}}
   ]
   evaluated in the overlap region and tracked toward strong-field radii.

7. **Make the wake-mixing (\alpha^2=3/4) a geometric throat statement (not a fitting constant)**
   Paper III’s translational wake ansatz is explicitly built from longitudinal/transverse projectors and yields (\alpha^2=3/4). This paper should show how that emerges from “two-channel” response: brane-parallel (transverse) + bulk-directed (longitudinal) components induced by throat motion, with the DtN map providing the coupling.

8. **Build the 2PN effective action for a throat as an object with internal DOF**

   * Your ontology already frames 2PN as the regime where finite-size ((a/r)^2) effects become physical.
   * Deliverable: an EFT-style object action with:

     * pointlike terms (recovering the calibrated 1PN story),
     * finite-size operators with coefficients fixed by matching,
     * *and* a small set of internal mode coordinates (a_n(t)) (optional but extremely valuable).

9. **Strong-field application: evaluate the near–far mixing at horizon / photon sphere / ISCO radii**
   `1pn_hybrid.tex` already defines an acoustic-horizon picture and strong-field observables (photon sphere, shadow size, etc.).
   Deliverable: explicit plots/estimates of (\mathcal{R}*\text{geom}(r)) and (\mathcal{R}*\text{EM}(r,\omega)) evaluated at (r=r_H, r_{\text{ph}}, r_{\text{ISCO}}) as functions of (M) (or your defect parameterization).

10. **Stress-test/failure scans** (the “make it fail” section)
    Deliverables:

* find parameter regions where ((a/r)^2) suppression breaks early (large (\alpha_2)),
* identify resonant frequencies where brane behavior becomes strongly nonlocal in time (big phase lag, hysteresis),
* locate where plasma micro-scales (Debye/skin depth, if you’re running kinetic) approach (a) and the “point throat” approximation collapses.

11. **Energy accounting closure including throat leakage**
    Explicitly show where energy goes: brane field ↔ particle motion ↔ throat modes ↔ (optional) deep-bulk leakage.
    This is the cleanest way to *diagnose* “hidden 4D interactions” without mysticism.

12. **A minimal inner/outer numerical coupling scheme that can be implemented**

* Inner solver: 4D throat PDE (or mode-truncated spectral surrogate).
* Outer solver: 3D brane fields + particle/throat motion.
* Coupler: the DtN/impedance map at the mouth + multipole matching for far-field consistency.

---

## Proposed paper structure (roadmap)

### 0. Executive summary and what is new

* “We provide the missing near–far matching layer: a local 4D throat PDE matched to 3D brane fields, yielding explicit 2PN finite-size coefficients and a measurable near–far mixing diagnostic.”

### 1. Series recap and constraints already fixed

Pull in the *non-negotiables* already established:

* 1PN orbital/optical calibration (and strong-field consequences from `1pn_hybrid.tex`),
* translational wake projector structure + (\alpha^2=3/4) from `1pn_spin_and_nbody.tex`,
* EM charge mapping (q\propto \rho_0 \pi a^2\Gamma) and field-theory mode option from `em_fields.tex`,
* “sphere–cylinder tension” and ((a/r)^2) suppression logic from `brane_bulk_ontology.tex`.

### 2. Geometry model of a brane–bulk throat (rounded funnel family)

* define the geometry with a few parameters
* define the matching interface (mouth)
* define what “near field” and “far field” mean operationally

### 3. Inner 4D throat dynamics

* write the governing PDE and boundary conditions
* classify the bottom condition (closed cavity vs leaky vs connected)
* solve for modes, show Bessel/cylindrical structure

### 4. Mouth matching: DtN / impedance map

* derive the boundary operator
* show low-frequency and near-resonant limits
* interpret as “hidden 4D degrees of freedom” that appear as dispersion/memory in 3D

### 5. Outer 3D brane fields and multipole expansion

* recover monopole + ((a/r)^2) corrections
* compute (\alpha_2) etc as functions of geometry parameters

### 6. Near–far mixing diagnostics (“the gradient”)

* define (\mathcal{R}*\text{geom}(r)) and (\mathcal{R}*\text{EM}(r,\omega))
* show how these scale and where they become (\mathcal{O}(1))

### 7. 2PN effective theory for throats with internal modes

* write (L_\text{eff}=L_N+\epsilon L_{1PN}+\epsilon^2 L_{2PN}+\delta^2\tilde L_{2PN}+\cdots)
* show which coefficients are fixed by matching (not fitted)
* show how this reduces to prior 1PN results in the appropriate limits

### 8. Strong-field application using the hybrid paper’s horizon/optics machinery

* evaluate the mixing ratios at (r_H), (r_\text{ph}), (r_\text{ISCO})
* identify observable consequences and “failure radii”

### 9. Failure scans

* where the multipole suppression fails
* where throat resonances dominate
* where nonlocal time response cannot be approximated by local closures
* where micro-scales cross (a)

Appendices: full DtN derivation, numerical implementation notes, parameter tables.

---

## Should this paper connect to 2PN GR now, or later?

**Recommendation:** do a *structural* connection now, and save *coefficient-level matching* to 2PN GR for a later paper.

What “structural connection” means (safe + useful now):

* Write your 2PN throat effective action in the same *operator language* as standard PN/EFT approaches: point-particle terms + finite-size multipoles/tidal operators.
* Emphasize: “our coefficients are determined by throat geometry and inner PDE matching.”

What to postpone (unless you want a much heavier paper):

* a full, coefficient-by-coefficient identification with the GR 2PN EIH/NRGR coefficients across all channels (and especially anything involving radiation reaction).

This keeps the next paper laser-focused on the missing toy-model piece (4D throat inner problem + matching), while making it *easy* to compare to 2PN GR later because your result will already be expressed in the right basis.

===

# Paper 7 — Cylinder throat DtN: Step 2 writeup (what we learned)

## What Step 2 was trying to establish

Step 2 is **not** the final “find (L/a)” mechanism. It’s a **sanity + feasibility step** to confirm we can treat the inner throat as:

1. **Local / PN-expandable** (at least for the (\lambda=0) “monopole” channel), and
2. A **well-defined boundary operator** (Z(\omega)) that the outer 3D solver/EFT can consume, and
3. A **diagnostic scan** that doesn’t accidentally “select (L/a)” just because we picked a damping model or a bad scan band.

So Step 2 is mostly about **making sure the math behaves** (local expansion exists, matching criteria aren’t biased), and about **learning which modes are numerically/physically fragile** (resonances/cutoffs).

---

## Model recap: what the script is computing

### Cylinder throat DtN eigenvalues

We use a separable throat wave equation in a cylinder with:

* Radius (a), length (L)
* Wall boundary at (r=a): Neumann or Dirichlet → fixes Bessel roots
* Bottom boundary at (z=-L): Neumann or Dirichlet → fixes standing-wave condition

For each azimuthal/radial pair ((m,n)), define the radial separation constant
[
\lambda_{mn}=\frac{x_{mn}}{a},
]
where (x_{mn}) is:

* Dirichlet: a root of (J_m(x)=0)
* Neumann: a root of (J_m'(x)=0) (with the special Neumann (m=0,n=0\Rightarrow \lambda=0) mode)

Define the axial wavenumber
[
\kappa_{mn}(\omega)=\sqrt{\Big(\frac{\omega}{c}\Big)^2-\lambda_{mn}^2},
]
(using a consistent branch choice).

Then the **DtN eigenvalue** at the mouth has the standard waveguide/cavity form:

* bottom Dirichlet: (Z_{mn}=\kappa\cot(\kappa L))
* bottom Neumann: (Z_{mn}=-\kappa\tan(\kappa L))

For the specific case you ran most often (Neumann wall + Neumann bottom):

* ((0,0)): (\lambda=0\Rightarrow \kappa=\omega/c\Rightarrow Z_{00}(\omega)=-\frac{\omega}{c}\tan!\big(\frac{\omega}{c}L\big))
* ((0,1)): (\lambda_{01}\approx 3.8317/a) and (\kappa) becomes small near cutoff.

---

## Step 1 result we carried into Step 2 (quick)

The Step 1 script already validated:

* **Mode table sanity**: for Neumann wall, (m=0,n=0\Rightarrow \lambda=0); (m=0,n=1\Rightarrow \lambda\approx 3.8317).
* **Resonance prediction**: for ((0,0)) Neumann bottom, resonances at
  [
  \omega_q \approx \frac{(q+\tfrac12)\pi c}{L},
  ]
  and your output matched the predicted (\omega_0=\frac{\pi c}{2L}) and spacing (\Delta\omega=\frac{\pi c}{L}) exactly.
* **Power proxy sign check**: (\omega,\mathrm{Im}(Z)) came out negative for the normal you used, and (-\omega,\mathrm{Im}(Z)) positive. That’s not “wrong”; it’s just the sign convention for the outward normal / conjugate variables. The important point is: we can choose a convention where the throat is **passive** (dissipative sign definite).

---

## Step 2A: Locality / PN-expandability of the monopole channel ((0,0))

### What we tested

We tested whether the ((0,0)) DtN map is **analytic at (\omega=0)** and admits a **small-(\omega)** expansion accurate on a “PN-ish” band.

### Final locality result (v4–v6 consistent)

For ((m,n)=(0,0)), (L/a=1.85), you got:

* Test range: (\omega \in [10^{-3}\omega_0,;0.25,\omega_0])
* Series order: 8
* Tolerance: 0.1
* **Max relative error over band:** (\sim 1.25\times 10^{-5})
* **Largest ω meeting tol:** the full tested max (0.25 ω0)

The printed series:
[
Z_{00}(\omega)\approx
-1.85,\omega^2
-2.1105416667,\omega^4
-2.8893315417,\omega^6
-4.0025841053,\omega^8
]
(consistent with the Taylor series of (-\omega\tan(\omega L)) when (c=1)).

### What this means (paper-relevant)

This is the key Step-2 deliverable:

* The ((0,0)) channel is **strongly local** in the low-frequency regime.
* In EFT language, the inner throat contributes a **derivative expansion** at the mouth (a local series of operators), rather than a mandatory long-memory kernel.
* This is exactly the “sanity check” needed before trying to write a **2PN effective action with finite-size terms**, because it shows at least one physically important channel (monopole) behaves like a controlled PN expansion.

---

## Step 2B: Criteria scans over (L/a) (what they were *actually* doing)

### What we were scanning

We scanned (L/a) over a range (typically 1.2 → 2.6) and computed:

1. A **DtN mismatch** criterion: compare inner (Z_{mn}(\omega;L)) to an outer “reference” impedance (Z_{\text{out}}(\omega)) on a chosen frequency band, using a normalized mean-square mismatch.
2. An **energyRed** proxy: a reduced stored-energy metric (normalized to remove trivial scaling).

These are *diagnostics*, not a guaranteed selector.

### The biggest pitfall we found (and fixed)

**Absolute damping (\gamma) breaks scale invariance.**

* For the ((0,0)) band, (\omega\sim\omega_0(L)\propto 1/L).
* If you keep (\gamma) fixed in absolute units, then (\gamma/\omega) changes with (L), and the “best” (L/a) will often run to the scan boundary as a pure artifact.

**Fix:** switch to
[
\omega\to \omega,(1+i\delta) + i\gamma_{\rm floor}
]
so (\gamma_{\rm eff}/\omega\approx\delta) is constant. In your v6/v7 outputs, the script reported
[
\mathrm{median}(\gamma_{\rm eff}/\omega)\approx 0.02
]
and it stayed essentially constant across the scan. That was the right move.

### The second pitfall (partially fixed)

**Resonance exclusion can erase your ω-grid, especially near cutoff.**

We excluded ω values too close to poles (to avoid huge spikes from (\tan(\kappa L)) / (\cot(\kappa L))). For (\lambda>0) modes, resonance spacing in ω is not (\pi c/L); it’s smaller near cutoff because
[
\Delta\omega \approx \frac{\pi}{L}\frac{c^2|\kappa|}{|\omega|}.
]
When (|\kappa|) is small (near cutoff), resonances get dense in ω, and an exclusion rule can remove most of the band.

We improved this (v7) by:

* using a **mode-correct spacing estimate**, and
* capping exclusion width to a fraction of band width.

---

## What the final v7 scan showed (from your `criteria-scan5.out`)

### Mode (0,0): “no size selection” (expected, good)

* `mismatch(min)` ~ 8.4424 and basically **flat** vs (L/a)
* `energyRed(min)` ~ 7.9570 and also effectively **flat** (tiny drift only)

**Interpretation:** the monopole channel does not intrinsically pick (L/a) under these scale-free criteria — and it **shouldn’t**. Once we removed the damping bias, the scan behaved properly.

This is a positive result: it means Step 2 is no longer “secretly selecting geometry” due to numerics.

### Mode (0,1): a real-looking optimum *before* NaNs

The scan reported:

* `best_L/a ≈ 1.92` for **both** mismatch(min) and energyRed(min)
* Past ~2.34, values become **NaN**.

**Interpretation of the optimum:**

* This looks like a genuine “best matching” region for that band + outer model choice.
* It’s not yet a universal throat-size selector; it is a **mode-dependent matching optimum**.

**Interpretation of the NaNs:**

* Not “physics blowing up.” It means: after resonance exclusion, **too few ω points survived** (or too many evaluations returned Indeterminate near poles), so the script correctly refused to compute a mean and printed NaN.
* In other words, for that mode and chosen band, large (L/a) lives in a regime where the response is **too resonance-dominated** to summarize with a “smooth band average” unless you change the band or the exclusion strategy.

This is still useful: it marks a **nonlocal / resonance-dense regime** where an outer local closure is fragile.

---

## What Step 2 accomplishes for the paper (clean “conclusions”)

This is the “what we now know” list:

1. **The DtN map is explicit and mode-resolved.**
   We now have a closed-form cylinder throat DtN eigenvalue (Z_{mn}(\omega)) with validated resonance structure.

2. **Monopole channel ((0,0)) is PN-local and operator-expandable.**
   The small-(\omega) series matches the exact impedance extremely well up to at least (0.25\omega_0) (and we can push that if needed). This is exactly what we need to justify a PN/EFT-style local expansion at the mouth.

3. **Scale invariance matters: damping must scale with ω.**
   Using absolute (\gamma) makes scans “select (L/a)” by accident. Using (\gamma_{\rm eff}\propto \omega) fixes this. This is an important methodological lesson to include in the “make it fail / stress test” section.

4. **Higher modes ((\lambda>0)) expose resonance-density and cutoff fragility.**
   The ((0,1)) mode shows a meaningful optimum (~1.92) in the “safe” region, but also shows that at larger (L/a) the response becomes resonance-dense enough that band-averaged criteria can fail (NaNs) unless you redesign the diagnostic.

5. **Therefore Step 2 is a feasibility filter, not the final size selector.**
   ((0,0)) proves “we can write local operators.” ((0,1)) shows “watch out: near cutoff you can get nonlocal/resonant behavior.”
   The **actual selection of (L/a)** still needs Step 3 (frequency anchor, variational selection, or external scale).

---

## What we should do next (Step 3 direction)

Step 3 should introduce an **external anchor** so (L/a) becomes physically selectable, e.g.:

* a fixed drive frequency (\omega_*) (from EM/plasma or experiment), and show which (L/a) aligns (\omega_*) with a mode family (\omega_{mnq}(L/a)), or
* a variational/coupling functional (stored energy or mouth admittance under fixed forcing) that yields a robust optimum.

Step 2 already gave us the tools we need to do that credibly (explicit (Z_{mn}), resonance structure, locality regime, and failure regimes).

---

## Do we “have what we need” despite NaNs?

**Yes.** For Step 2, NaNs are not a blocker; they’re part of the story:

* Step 2’s core result is the ((0,0)) locality + the scale-invariant scan behavior.
* The ((0,1)) NaNs are a useful red flag: “resonance-dense / nonlocal regime where coarse band diagnostics break.” That’s *exactly* the type of failure-mode we wanted to identify.

If you want, the next incremental improvement (not required for Step 2) would be: instead of averaging over a fixed ω band for ((0,1)), compute **resonance-resolved metrics** (evaluate at controlled offsets from each pole), so exclusion never wipes out the dataset.

===

## Paper 7 — Cylinder throat DtN: Step 3 writeup (what we learned)

### What Step 3 was trying to do (the “why”)

Step 2 showed two key facts:

1. **Monopole channel (m,n)=(0,0)** is *PN-local* at low frequency (clean Taylor/derivative expansion exists).
2. **Higher modes (λ>0)** can become *resonance-dense / fragile* (band-averaged criteria can collapse or bias).

So Step 3’s job was to introduce an **external anchor** (a dimensionless band or resonance family) so that **(L/a)** becomes *selectable*, and to do it in a way that doesn’t accidentally “pick L/a” due to numerics.

In other words: Step 3 is our first serious attempt at a **geometry-selection diagnostic** that can be fed into the paper as a credible procedure (and that we can later swap in a better outer model for).

---

## Model objects used in Step 3 (recap)

We keep the same cylinder throat DtN eigenvalues:

* Radial separation constant:
  [
  \lambda_{mn}=\frac{x_{mn}}{a}
  ]
* Axial wavenumber:
  [
  \kappa_{mn}(\omega)=\sqrt{\Big(\frac{\omega}{c}\Big)^2-\lambda_{mn}^2}
  ]
* DtN eigenvalue at the mouth (Neumann bottom in your runs):
  [
  Z_{mn}(\omega)=-\kappa\tan(\kappa L)
  ]

And Step 3 introduces:

* A **dimensionless frequency anchor** (for λ=0): a fixed **ka band**.
* A **near-cutoff / propagating anchor** (for λ>0): a band in
  [
  \eta(\omega)\equiv \frac{k}{\lambda}=\frac{\omega/c}{\lambda}
  ]
  (your `etaCutBand={...}`), which effectively chooses a **controlled region relative to cutoff**.

We also keep the Step-2 lesson: damping must be scale-aware, so we use **relative-ω damping** (`gammaMode=RelativeOmega`, `deltaRel=0.02`, with a tiny floor).

---

## Step 3.0 (“grid band average”) — what happened

### What it did

For each candidate (L/a) on a scan grid:

* pick a frequency set (\omega) inside the chosen anchor band (ka-band for monopole; eta-band for λ>0),
* compute (Z_{\text{in}}(\omega;L)),
* compare it to a simple **outer reference impedance** (Z_{\text{out}}(\omega)) (your script’s `outerModel=AutoScaled`),
* compute a mismatch metric averaged over the band,
* run robustness tests: band scaling by `s∈{0.95,1,1.05}` and a small band jitter.

### Result: (0,0) behaves “correctly” (but doesn’t select L/a)

For (0,0), mismatch was smooth and monotonic with (L/a), and the “best” value always landed at the **lower scan boundary** (e.g. 1.2). That’s *not* physics selecting a special (L/a); it’s telling us:

* **Monopole is too scale-free under these criteria**, and/or
* the particular “outerModel” reference is not imposing a meaningful size anchor (so the metric prefers “smaller L” because the impedance magnitude trend is monotonic).

This is consistent with Step 2: **(0,0) is a good PN-local channel, but not a reliable size selector by itself**.

### Problem: (0,1) selection was not robust

For (0,1), Step 3.0 gave a plausible best (L/a) (earlier you saw values like ~1.6–1.7), but the **robustness test drifted badly**: small changes in the band endpoints could move the optimum by ~0.7–1.0 in (L/a).

That’s the red flag:

* **band-averaged metrics are not stable for λ>0** because the response is dominated by **where your ω-band lands relative to resonances** and cutoff structure.

---

## Step 3.1 and 3.1b (“resonance-aware sampling”) — what changed and what we learned

### Why we changed it

The λ>0 failure mode is: *a fixed ω-grid is not the right object near resonance structure.* A slight change in band endpoints can include/exclude poles (or near-pole spikes), which dominates the average.

So 3.1 tried to fix that by sampling **controlled offsets from resonances** (“detunes”), instead of brute averaging a dense ω-grid.

### What initially went wrong (3.1)

The first 3.1 output produced **all NaNs** for (0,1): it was filtering out *everything*. That told us our resonance selection logic was too strict (or the anchor filtering was incompatible with how we were choosing sample points).

### 3.1b produced finite tables, but still drifted

In 3.1b you got a clean scan table for (0,1), with:

* (N_\omega\sim 20) up to (\sim 47),
* (N_{\text{res}}\sim 2) or 3,
* smooth mismatch vs (L/a).

But even then, the robustness test still drifted:

* Example: best (L/a) might be ~1.74 at s=1.00, but shift to ~2.14 at s=0.95 and ~1.46 at s=1.05.

**Interpretation:** resonance-aware sampling is the right direction, but if the sampling is still “band-defined” in a way that changes *which resonances* are included, the optimum will still move.

That led directly to Step 3.2.

---

## Step 3.2 (“fixed-q resonance family”) — the key Step-3 outcome

### The core idea

For λ>0 modes, instead of:

* “sample whatever resonances happen to fall into this band,”

we do:

* “always sample the **same axial resonance indices** q (e.g. q=0,1,2…) for every L/a, and compare like-to-like.”

This is exactly what `resSamplingMode=FixedQ` and `qFixedList={...}` did.

### What it fixed

For the (0,1) mode, Step 3.2 produced **stable optima under band scaling**:

* In your Step 3.2 output with `qFixedList={0,1,2}`, the best (L/a) was **1.46**, and the robustness test gave:

  * s=0.95 best = 1.46
  * s=1.00 best = 1.46
  * s=1.05 best = 1.46
    → **Band-jitter drift Δ(L/a)=0**

That’s a big deal: it means the selection is no longer an artifact of “which resonances got included.”

### What it reveals (important)

Even though 3.2 makes the selection robust, it also shows something unavoidable:

> **The selected (L/a) depends strongly on which resonance family (which q list) you assume is physically relevant.**

You demonstrated this explicitly:

* `qFixedList={0}` → best (L/a = 2.6) (edge of scan)
* `qFixedList={0,1}` → best (L/a = 1.44)
* `qFixedList={1,2}` → best (L/a = 1.62)
* `qFixedList={0,1,2}` → best (L/a = 1.46)
* `qFixedList={0,1,2,3}` → best (L/a = 1.38)

And these were **internally robust** (no drift) once you fix the q-family.

**Interpretation:** Step 3.2 is telling us that geometry selection is *not unique* until we specify which axial resonances the outer system actually excites or couples to.

That is not a bug — it’s physics/identifiability:

* The throat has a **mode family**.
* The outer world decides **which member(s)** matter.

So Step 3.2 gives us a *clean interface*:
**“Given a chosen mode family (m,n,q set) + anchor definition + outer impedance model, we can stably compute best L/a.”**

---

## What we can safely conclude from Step 3 (paper-ready)

### 1) Monopole channel (0,0): good for PN locality, not for size selection

* It stays well-behaved under ka-band anchoring and robustness tests.
* But the mismatch is monotonic and “best L/a” hits the scan boundary (1.2) → **no intrinsic size-selection here** under this outer-model/metric.

This is consistent with Step 2 and actually strengthens the story: monopole gives you **local operators**, not a geometry pin.

### 2) Higher mode (0,1): band-averaged scans are unreliable unless resonance structure is controlled

* Step 3.0/3.1b drift under small band changes shows the failure mode clearly.
* This is the “make it fail” lesson we should explicitly include in the paper.

### 3) Fixed-q resonance sampling (Step 3.2) is the first robust selector we have

* When we compare like-to-like (same q indices across L/a), the optimum becomes **stable** under band scaling/jitter.
* This is the strongest Step-3 deliverable: it’s a **methodologically clean** way to select geometry from resonance families.

### 4) But Step 3.2 also proves an identifiability point: you must specify which resonance family is physically excited

* Different q lists give different best L/a.
* Therefore: **the final selection requires an input from the outer physics** (drive spectrum, coupling strength, which resonances dominate).

That tells us exactly what Step 4 must add.

---

## What Step 4 needs to do (so L/a becomes “real”)

Step 4 should supply the missing physical information that Step 3.2 makes explicit:

1. **A physically grounded outer model** (Z_{\text{out}}(\omega)) (from the outer 3D solver / EFT / boundary condition model), not the placeholder AutoScaled reference.
2. **A coupling/weight model** telling us which throat channels and which q’s matter:

   * e.g. weights (w_{mnq}) from mouth overlap integrals / outer field symmetry / driving configuration,
   * or directly: “the drive is narrowband near ω* so only q≈q* is relevant.”
3. Then run Step 3.2 as:
   [
   \text{choose } L/a \text{ to minimize } \sum_{mnq} w_{mnq},\mathcal{M}[Z_{mn}(\omega_{q}(L)), Z_{\text{out}}(\omega_{q}(L))]
   ]
   where (\omega_q(L)) is the predicted resonance location for that fixed-q family.

At that point, “best L/a” stops being an internal numerical artifact and becomes a real inference from outer physics.

---

## Short “what to paste after Step 2” summary

* Step 3 introduced dimensionless anchors (ka-band for λ=0; eta-band / resonance sampling for λ>0) to attempt geometry selection.
* (0,0) stayed stable and PN-friendly but did not uniquely select L/a (best ran to scan boundary).
* (0,1) showed that naive band averaging is unstable (optimum drifts under tiny band changes) because resonance structure dominates.
* Step 3.2 fixed this by switching to **fixed-q resonance-family sampling**, producing **robust** best L/a values under band scaling/jitter.
* However, Step 3.2 also shows geometry selection depends on which q-family you assume is relevant; therefore Step 4 must provide the outer-physics input that determines which resonances are actually excited/coupled, and should replace the placeholder outer impedance model with something derived from the outer solver/EFT.

===

## Paper 7 — Cylinder throat DtN: Step 4 writeup (what we learned)

### What Step 4 was trying to do (the “why”)

Steps 1–3 established a usable **inner throat DtN map** (Z_{mn}(\omega)) and showed:

* the **monopole** ((m,n)=(0,0)) channel is PN-local (admits a clean derivative expansion),
* higher modes ((\lambda>0)) can become **resonance-dense / fragile**, and
* **fixed-q** sampling is the first way to make “geometry selection” numerically robust (no band-endpoint artifacts).

Step 4’s job was to add the missing second half:

> **Build and test a concrete “outer” impedance/DtN map (Z_{\text{out}}(\omega))** and check whether *inner/outer matching* can (i) behave sanely, (ii) be stable under perturbations, and (iii) *actually select* (L/a) in a way that is compatible with a future **2PN / EFT-style** formulation (i.e., conservative-static limit is meaningful).

So Step 4 is the first “inner PDE ↔ outer field” coupling test that the eventual paper needs.

---

## Model objects used in Step 4 (recap)

### Inner: cylinder throat DtN eigenvalues (same as Steps 1–3)

For each ((m,n)):

* Radial separation constant:
  [
  \lambda_{mn}=\frac{x_{mn}}{a},
  ]
  where (x_{mn}) is a Bessel root determined by the wall BC:

* Dirichlet wall: (J_m(x_{mn})=0)

* Neumann wall: (J_m'(x_{mn})=0) (with special Neumann mode ((m,n)=(0,0)\Rightarrow \lambda=0))

* Axial wavenumber:
  [
  \kappa_{mn}(\omega)=\sqrt{\Big(\frac{\omega}{c}\Big)^2-\lambda_{mn}^2}
  ]
  (using a consistent branch choice).

* DtN eigenvalue at the mouth (Neumann bottom in your runs):
  [
  Z_{mn}(\omega)= -\kappa \tan(\kappa L),
  \qquad L=(L/a),a.
  ]

We keep **scale-aware damping**:
[
\omega\to \omega(1+i\delta) + i\gamma_{\rm floor},
]
with (\delta=0.02) and (\gamma_{\rm floor}=10^{-6}), to avoid the Step-2 pitfall where fixed absolute damping biases the scan.

### Outer: spherical DtN map at (r=a)

Step 4 introduced an explicit outer boundary operator that can be evaluated at the mouth radius (r=a), multipole by multipole:

* **ConservativeStatic** (Laplace/static limit):
  [
  Z_\ell^{\rm static}(a)= -\frac{\ell+1}{a}.
  ]

* **Radiating** (Helmholtz outgoing-wave DtN):
  [
  Z_\ell^{\rm rad}(\omega;a)= k\frac{h_\ell^{(1)\prime}(ka)}{h_\ell^{(1)}(ka)},
  \qquad k=\omega/c,
  ]
  using spherical Hankel functions.

Then define an effective outer impedance as a weighted multipole sum:
[
Z_{\text{out}}(\omega)=\sum_{\ell\in\mathcal{L}} w_\ell(\omega),Z_\ell(\omega).
]

### Key symmetry/coupling fact discovered (very important)

For the Neumann-wall cylinder, the ((m,n)=(0,1)) radial mode has **zero monopole (disk-average) overlap** at the mouth. Operationally: it does not couple to (\ell=0) in the simplest “uniform mouth” projection.

That means: if you try to match ((0,1)) using only (\ell=0), you are matching the wrong channel. The first relevant outer channel is (\ell=2) (or higher).

Step 4 operationalized this by running (\ell=2) outer matching:

* `lList=(2,)`
* weights `wOuter={2:1, 0:0}`

---

## Step 4A: Sanity checks that the outer operator is correct

### Static limits (passed)

Your script prints:

* (Z_0^{\rm static}=-1/a)
* (Z_2^{\rm static}=-3/a)

and the output matched expectations exactly:

* `Zout0_static = -1 (expected -1/a)`
* `Zout2_static = -3 (expected -3/a)`

This is a key correctness check: the outer operator reduces to the right conservative multipole DtN.

### Passivity sign proxy (consistent)

The script prints a “power proxy” sign check:

* (\omega,\mathrm{Im}(Z_{\rm in})) and (-\omega,\mathrm{Im}(Z_{\rm in}))

For your runs, (\omega,\mathrm{Im}(Z_{\rm in})<0) at sample points, and the flipped sign is (>0). This is not an error; it is a **convention** (outward normal / conjugate variables). It confirms we can choose a sign convention where the throat is passive.

---

## Step 4B: A real bug we hit and fixed (and why it matters)

When you switched `lList=(2,)` but left `wOuter={0:1, 2:0}`, the effective outer impedance became:
[
Z_{\text{out}}(\omega)\equiv 0,
]
so the **RelL2** mismatch (which divides by (\langle|Z_{\text{out}}|^2\rangle)) exploded to (\sim 10^{13}-10^{14}).

Fix: weights must be consistent with chosen multipoles:

* for `lList=(2,)`, set `wOuter={2:1, 0:0}`.

This is worth mentioning in the “make it fail / implementation pitfalls” section: *Rel-normalized mismatch will blow up if the reference channel is weighted out.*

---

## Step 4C: Main matching scans and what they show

### Setup used for the key scans

* Inner mode: ((m,n)=(0,1))
* Inner BCs: Neumann wall, Neumann bottom
* Outer multipole: (\ell=2) only
* Weights: (w_2=1), (w_0=0)
* Sampling: **fixed-q resonance-family sampling**, with detunes:

  * `qFixedList=(0,1,2)` or `qFixedList=(0,)`
  * `detuneRel=0.02`, detune signs (\pm)
* Damping: RelativeOmega with (\delta=0.02)

Mismatch metrics used:

* `RelL2`: normalized mean-square mismatch
* `AbsL2`: absolute mean-square mismatch (to rule out normalization artifacts)

---

## Step 4C.1: Multi-q matching does NOT select an interior (L/a) (robust boundary drift)

With `qFixedList=(0,1,2)` and (\ell=2) matching:

* In **Radiating** mode, (J(L/a)) decreased monotonically as (L/a) increased.
* Extending the scan to (L/a=4) showed the same trend.
* Best (L/a) always landed on the **upper scan boundary**, with robustness drift (\Delta(L/a)=0).

You repeated the same style of scan in **ConservativeStatic** and again:

* (J(L/a)) decreased monotonically, best at the largest (L/a) scanned.

**Interpretation (physics, not numerics):**

Fixed-q resonance sampling has an inherent “cutoff drift” property for (\lambda>0) modes:
[
\omega_{mnq}(L)=c\sqrt{\lambda_{mn}^2+\kappa_q^2},
\qquad \kappa_q\propto 1/L
\quad\Rightarrow\quad
\omega_{mnq}(L)\to c\lambda_{mn}
\text{ as } L\to\infty.
]

So as (L/a) grows, the sampled frequencies slide toward cutoff ((\omega\to c\lambda)), and the mismatch tends to improve monotonically. That yields a stable but non-identifying “best at boundary” outcome.

**Paper-safe conclusion:**

> Inner/outer DtN matching (with multi-q fixed-q sampling) is stable, but **does not uniquely determine (L/a)**; it is underdetermined without additional outer-physics anchoring.

This is exactly the kind of “make it fail” result we want: it prevents us from overclaiming that conservative PN structure will automatically pin the throat length.

---

## Step 4C.2: Single-q (q=0) gives a weak interior optimum near (L/a \sim 3.1–3.2)

When you assume a **q=0-dominant** scenario (`qFixedList=(0,)`), you observed a genuine interior minimum:

### Local refined scan (RelL2)

In a narrow window 3.04–3.14:

* best (L/a = 3.08) with (J\approx 1.4870388)
* the minimum was real but extremely shallow.

### Broad scan (AbsL2)

Across 1.2–4.0:

* best (L/a \approx 3.2), (J\approx 17.0966)
* this demonstrated the feature is **not purely a RelL2 normalization artifact**.

### Robustness for q=0 (both metrics)

Scaling detuneRel and deltaRel together by (s\in{0.95,1.00,1.05}):

* **RelL2:** best (L/a) moved (3.16 \to 3.08 \to 3.00)
* **AbsL2:** best (L/a) moved (3.26 \to 3.18 \to 3.10)
* Drift in both cases: (\Delta(L/a)=0.16)

**Interpretation:**

* The q=0 interior optimum exists, but it is **soft/phase-sensitive**: it shifts measurably under small damping/detune changes.
* Therefore it can be used as a *conditional inference* (“if q=0 dominates”), but it is not a universal geometry selector.

**Paper-safe phrasing:**

> Under a q=0-dominant hypothesis, radiating outer matching yields a shallow optimum near (L/a\approx 3.1–3.2), with sensitivity (\Delta(L/a)\approx 0.16) under ±5% damping/detune scaling.

---

## Step 4D: Conservative-static is the “2PN sanity” filter — and it fails to select (L/a)

Your ConservativeStatic scan for q=0 over the local window still drifted to the boundary (best at the maximum scanned (L/a)), with (J) still decreasing as (L/a) increased.

This is a crucial “sanity” point relative to future 2PN matching:

* In the conservative (static) limit, there is no strong mechanism here that pins a unique (L/a).
* That’s consistent with EFT/PN expectations: conservative outer structure organizes operators and coefficients but typically does **not** fix a cavity’s internal length without additional physics.

So Step 4’s conservative-static result is exactly what we wanted to know early: **we won’t be able to “derive (L/a)” from 2PN structure alone**.

---

## What Step 4 accomplishes for the paper (clean conclusions list)

1. **Outer DtN operator is explicit and correct.**
   Static limits (Z_\ell^{\rm static}=-(\ell+1)/a) pass; radiating DtN is well-defined.

2. **Channel selection matters (symmetry/overlap).**
   The ((0,1)) inner mode does not couple to (\ell=0) under a naive mouth projection; matching must include (\ell=2) (or higher). This is a major conceptual constraint.

3. **Matching is numerically stable once weights are consistent.**
   The earlier “huge J” blowups were caused by (w_\ell) accidentally zeroing the chosen (\ell)-channel while using a normalized (RelL2) metric.

4. **Multi-q fixed-q matching does not uniquely select (L/a).**
   With ((0,1)) and (\ell=2), the mismatch decreases monotonically with (L/a) and runs to scan boundaries. This persists in both radiating and conservative-static outer models.

5. **A weak interior optimum appears only under restrictive excitation assumptions (q=0).**
   For q=0 only, an interior optimum near (L/a\sim 3.1–3.2) exists, but it is shallow and shifts by (\Delta(L/a)\sim 0.16) under small perturbations.

6. **Therefore: geometry selection requires an external anchor.**
   Without specifying drive spectrum (\omega_*), coupling weights (w_{mnq}), or a fixed “distance above cutoff” anchor, inner/outer DtN matching is underdetermined in the conservative regime. This is the main Step-4 “make it fail” result.

---

## What Step 4 implies for the next step

Step 4 makes the identifiability gap explicit:

To turn matching into a credible (L/a) selector (especially in conservative/2PN-friendly settings), Step 5 must introduce at least one of:

* a **fixed drive frequency** (\omega_*) (or measured/narrowband spectrum),
* a **fixed (\eta=\omega/(c\lambda)) band** (distance above cutoff, to prevent cutoff drift),
* or **physics-derived weights** (w_{mnq}) (overlap integrals / mouth forcing symmetry).

This is not optional if we want geometry selection to survive the conservative-static filter.

---

## Implementation notes you can keep with Step 4

* **Do not use RelL2 unless (Z_{\text{out}}\neq 0)** on the sampled set. If the chosen multipole channel is weighted out, RelL2 will blow up.
* Always print:

  * static limits for (Z_\ell),
  * a passivity proxy,
  * the number of ω samples (N_\omega), and ωmin/ωmax,
  * robustness drift under small scaling of detune/damping.

That set of checks is exactly what prevented Step 4 from fooling us.

===

## Where (1.85) comes from in the existing series

### In `brane_bulk_ontology.tex`

It states explicitly that the EM cavity analysis picks
[
\frac{L}{a}=\frac{\sqrt{2},\pi}{x_{01}}\approx 1.85,
]
**where (x_{01}) is the first zero of (J_0)** (i.e. Dirichlet radial condition (J_0(x_{01})=0)). This is the classic (x_{01}\approx 2.40482556) root, so
[
\frac{L}{a}=\frac{\sqrt{2}\pi}{2.40482556}=1.84748658.
]

### In `em_fields.tex`

You actually have the more detailed derivation in the enthalpy model:
[
H(a,L)=A,L + B,\frac{a^2}{L}+C,a^2L,
]
with (A\propto x_{01}^2) and (B\propto \pi^2). Extremizing w.r.t. (a) and (L) gives
[
\frac{L}{a}=\sqrt{\frac{2B}{A}}=\frac{\sqrt{2},\pi}{x_{01}}\simeq 1.85,
]
and **the (C) term cancels out of the ratio**, so the aspect ratio is robust within that model structure.

**Key point:** this (1.85) result is a *variational/enthalpy selection* for a **specific cavity mode family** with **Dirichlet radial root (x_{01})** and a **standing-wave axial (k\sim \pi/L)**.

## What we did in Steps 1–4, and why it doesn’t automatically reproduce 1.85

In Steps 1–4 we mostly used:

* **Neumann wall** for the cylinder (so radial eigenvalues use (J_0'(x)=0\Rightarrow J_1(x)=0), with first root (3.83170597)),
* an **open mouth** (DtN/impedance port), not a closed “node at the mouth” cavity,
* and a **matching objective** (inner (Z_{\text{in}}) vs outer (Z_{\text{out}})), not “minimize enthalpy at fixed charge.”

If you take the *same* aspect-ratio formula but swap the radial root from Dirichlet to Neumann, you get:
[
\frac{L}{a}=\frac{\sqrt{2}\pi}{3.83170597}=1.15950518,
]
which is **nowhere near 1.85**. So unless we change the wall condition / variable definition, **we should not expect the EM-cavity-derived 1.85 to appear**.

And Step 4 added an additional reason matching won’t magically pin (L/a): for (\lambda>0) modes, fixed-(q) resonance sampling tends to drift toward cutoff as (L) grows, producing the monotone “best at boundary” behavior you observed. That’s a *different identifiability issue* than the EM enthalpy minimization.

## So do we have a “problem” if Step 4 doesn’t give 1.85?

Not necessarily — but **Paper 7 must explicitly reconcile the two regimes**. Right now, Step 4 is not the same model that produced 1.85 in the series.

To keep the toy model consistent, Paper 7 should do this:

### Required consistency check (high priority)

Show that **in the EM-cavity limit assumed in `em_fields.tex`**, your throat formalism reproduces:
[
L/a=\sqrt{2}\pi/x_{01}\approx 1.85.
]

That means, at minimum, Paper 7 needs a subsection like:

* **“Closed-cavity / high-impedance mouth limit”** (mouth behaves approximately like a node)
* **Dirichlet-type radial condition for the EM scalar** (so (x_{01}) is the first zero of (J_0))
* **Enthalpy functional (H(a,L))** → extremize → recover the same ratio

Once you do that, Steps 1–4 can be presented as the **open-mouth / coupling-to-outer-field generalization**, where geometry selection is underdetermined unless you supply an external anchor (drive spectrum, coupling weights, fixed (\eta), etc.).

## Concrete comparison you can paste into the paper

**Old (series EM selection):**

* Mode: (J_0(x_{01}r/a)\times \sin(\pi z/L))
* Radial root: (J_0(x_{01})=0 \Rightarrow x_{01}=2.4048)
* Selection principle: extremize enthalpy at fixed microphysics/charge
* Output: (L/a=\sqrt{2}\pi/x_{01}\approx 1.85)

**New (Paper 7 DtN matching work so far):**

* Mode basis: DtN eigenvalues (Z_{mn}(\omega)) with an **open mouth**
* Often used wall BC: Neumann ((J_0'(x)=0\Rightarrow x=3.8317) for first nontrivial)
* Selection principle: match (Z_{\text{in}}) to (Z_{\text{out}}) under sampling assumptions
* Output: stable matching exists, but conservative matching doesn’t uniquely select (L/a) without an anchor

## What I would do next (fastest path to resolve this cleanly)

1. **Re-run the cylinder throat scripts for the EM branch with `wallBC="Dirichlet"`** (so (x_{01}=2.4048)).
2. Add a short “EM cavity limit” derivation directly in Paper 7 reproducing (L/a=\sqrt{2}\pi/x_{01}). (This is essentially already in `em_fields.tex`; Paper 7 can reframe it as the *closed-mouth/high-impedance limit* of the throat model.)
3. Then explicitly state: **the open-mouth DtN generalization does not, by itself, guarantee the enthalpy-based selector**—it needs the same physical anchor (drive/coupling/termination) that made the cavity variational problem well-posed.

===


