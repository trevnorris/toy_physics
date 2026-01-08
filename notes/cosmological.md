## Slab vacuum mechanism: statement of the problem

We consider the scalar potential (\Phi) sourced by a localized mass/defect on the brane, embedded in a bulk that is **finite in one transverse direction**, with “walls” at (z=\pm H). The defining assumption of the *slab vacuum* is a **no-flux (Neumann) boundary condition** at the walls:
[
\partial_z \Phi\big|_{z=\pm H}=0,
]
which captures the idea that the bulk is bounded (or effectively bounded) in the transverse direction, so scalar flux cannot spread into an infinite 3D volume at arbitrarily large radius.

This is a controlled, global/topological modification: local physics near the source is unchanged, but the long-range Green’s function changes because of the global boundary condition.

---

## Two equivalent representations of the slab Green’s function

### 1) Method of images (direct sum)

On the brane ((z=0)), the Neumann slab can be implemented by an infinite stack of image sources at
[
z_n = 2nH,\qquad n\in\mathbb{Z}.
]
For a point source, the brane potential can be written as the image sum
[
\Phi_{\rm img}(r) \propto \sum_{n=-\infty}^{\infty}\frac{1}{\sqrt{r^2+(2nH)^2}},
]
and the radial force (magnitude) is
[
F_{\rm img}(r) \propto \sum_{n=-\infty}^{\infty}\frac{r}{\left(r^2+(2nH)^2\right)^{3/2}}.
]

This representation is conceptually transparent (flux confinement as “many copies”) but converges slowly for force at large (r) if truncated naïvely.

### 2) Mode / KK expansion (fast Bessel sum)

The same boundary value problem admits a Neumann mode expansion with transverse wave numbers
[
k_m=\frac{m\pi}{H},\qquad m=0,1,2,\dots
]
On the brane, the **force** can be written as a rapidly convergent Bessel-(K) series:
[
F_{\rm mode}(r)=\frac{GM}{H}\left[\frac{1}{r}+2\sum_{m=1}^{\infty}k_m,K_1(k_m r)\right].
]
This representation is numerically excellent at moderate/large (r) because the (m\ge1) terms decay exponentially (\sim e^{-k_m r}).

We also implemented a **finite potential difference** (to handle the 2D/logarithmic zero mode cleanly):
[
\Phi(r)-\Phi(r_0)=\frac{GM}{H}\left[-\ln!\left(\frac{r}{r_0}\right) +2\sum_{m=1}^{\infty}\big(K_0(k_m r)-K_0(k_m r_0)\big)\right],
]
so the log divergence is treated as gauge/offset (only differences matter).

**Key result:** The image sum and mode sum are two equivalent constructions of the same Green’s function for the Neumann slab.

---

## Asymptotic behavior: what is proved/verified

### Near-field: recovery of Newtonian scaling

For (r\ll H), the (n=0) term dominates the image sum and one recovers
[
F(r)\sim \frac{GM}{r^2}.
]
In the script we validate this numerically using the dimensionless diagnostic
[
r^2 F(r)\to GM\quad\text{(here (G=M=1\Rightarrow r^2F\to 1)).}
]

**Observed:** at (H=10), (r=0.5), we obtain
[
r^2F(r)=1.000037534
]
for both images and modes, confirming the Newtonian near-field limit to (\sim 4\times 10^{-5}) at a radius already “deep” but not asymptotically tiny.

### Far-field: dimensional crossover to (1/(Hr))

For (r\gg H), the slab enforces an effective dimensional reduction: flux cannot spread spherically forever, and the dominant “zero mode” behaves as a 2D Green’s function (logarithmic potential), yielding a force that scales as (1/r).

The force asymptote is
[
F(r)\sim \frac{GM}{H,r}.
]
The paper-grade normalization check is
[
H,r,F(r)\to GM\quad(\text{here }G=M=1\Rightarrow HrF\to1).
]

**Observed:** at (H=10), (r=1000), we obtain
[
H,r,F(r)=1.000000000
]
for both methods (with tail correction in the image sum), i.e. the far-field prefactor is confirmed at essentially machine precision.

---

## Rotation curve prediction from the slab force

For circular motion in the brane plane,
[
\frac{v^2(r)}{r}=F(r)\quad\Rightarrow\quad v(r)=\sqrt{r,F(r)}.
]

* In the Newtonian regime (F\sim GM/r^2), one gets Keplerian decay:
  [
  v(r)\sim \sqrt{\frac{GM}{r}}.
  ]
* In the slab regime (F\sim GM/(Hr)), one gets a **flat asymptote**:
  [
  v^2(r)\to \frac{GM}{H},\qquad v_\infty=\sqrt{\frac{GM}{H}}.
  ]

This is a strict consequence of the slab Green’s function and does not require any additional “dark matter” source term.

---

## Transition scale (“knee”) and log-slope diagnostic

To quantify where the crossover occurs, we computed the local log-slope of the force:
[
\frac{d\ln F}{d\ln r},
]
which should interpolate from (-2) (Newton) to (-1) (slab).

We defined a knee scale (r_{\rm knee}) by the slope criterion
[
\frac{d\ln F}{d\ln r}\Big|*{r=r*{\rm knee}}=-\frac{3}{2}.
]

**Observed:** for (H=10),
[
r_{\rm knee}\approx 10.496773,
]
i.e. the crossover occurs at (r\sim H) with an (\mathcal{O}(1)) coefficient. This gives a clean interpretation of (H): it sets the geometric transition scale.

---

## Numerical validation strategy and why it is robust

### Cross-check of two independent methods

A central reliability feature is that we computed the same physical force in two mathematically distinct ways:

1. image sum (physical intuition; slow convergence if truncated),
2. mode sum (spectral/Kaluza–Klein; fast convergence).

After adding a tail correction to the image force, the two agree at essentially floating-point limits:

Example cross-check table (relative error (|F_{\rm mode}-F_{\rm img}|/F_{\rm img})):

* (\sim 10^{-15}) at (r\sim 10)–30,
* rising to (\sim 10^{-13}) by (r\sim 1000),

which is consistent with machine precision effects (cancellation/roundoff) rather than physics differences.

### Tail correction for image force

Because the truncated image force has a slowly convergent tail at large (r), we added a controlled analytic tail estimate (integral approximation) to the images force. This collapses the large-(r) truncation error and makes the images–modes comparison clean across the entire range without extreme (n_{\max}).

**Paper point:** this is not “fitting” — it is the standard analytic remainder estimate for replacing a discrete sum tail by an integral.

### Convergence tables

We produced convergence tables vs truncation:

* image+tail error vs (n_{\max}),
* mode sum error vs (m_{\max}),

demonstrating that:

* the mode method converges extremely rapidly (exponential decay of higher modes),
* the images method becomes comparably accurate once the tail correction is included.

These tables justify the numerical settings used to generate the plots and any quoted values (e.g. knee scale).

### “Weird discrepancy plot” interpretation

After the tail fix, the images–modes discrepancy is so small that a relative-error plot can look noisy/jagged because it is dominated by numerical precision and subtraction of nearly equal numbers. This is expected and not a concern. For paper figures, an absolute error plot or a smoothed log-grid relative error plot is preferable.

---

## Extended source: far-field axisymmetric check (disk)

To show that the slab mechanism is not an artifact of a point source, we included an **extended axisymmetric disk** sanity check using the far-field slab law.

In the slab regime, for an axisymmetric mass distribution,
[
F(r)\approx \frac{G,M_{\rm enc}(r)}{H,r}
\quad\Rightarrow\quad
v^2(r)\approx \frac{G,M_{\rm enc}(r)}{H}.
]
Thus as (r\to\infty), (M_{\rm enc}(r)\to M_{\rm tot}) and
[
v_\infty^2\to \frac{G M_{\rm tot}}{H}.
]

We demonstrated this explicitly for an exponential disk using the closed-form enclosed mass:
[
M_{\rm enc}(r)=M_{\rm tot}\Big[1-e^{-r/R_d}\big(1+r/R_d\big)\Big].
]
The resulting (v(r)) approaches (\sqrt{GM_{\rm tot}/H}) as expected, validating that the asymptote is controlled by total mass and the slab scale (H).

(We also left an optional “full near-field disk via mode integrals” implementation scaffold, but the far-field check is sufficient for establishing the robustness of the scaling in the slab regime.)

---

## Wake-force toy model: kept separate as a distinct mechanism

We retained a separate test demonstrating a different flattening mechanism:
[
\frac{v^2}{r}=\frac{GM}{r^2}+C\frac{v}{r}
\quad\Rightarrow\quad
v^2-Cv-\frac{GM}{r}=0,
]
which yields (v\to C) at large (r). The script verifies this numerically.

**Paper positioning:** this is not part of the slab Green’s function claim; it is an optional additive/dynamical mechanism (potentially tied to other sectors of the model) and should be presented as such.

---

## What you can now claim in the paper (tight, defensible)

1. **Boundary-condition mechanism:** A finite transverse bulk extent with Neumann boundary conditions modifies the scalar Green’s function at large radius.

2. **Controlled asymptotes:** The resulting brane-plane force transitions from
   [
   F\sim GM/r^2 \quad (r\ll H)
   \qquad\text{to}\qquad
   F\sim GM/(Hr)\quad (r\gg H),
   ]
   with the far-field prefactor (\frac{1}{H}) verified numerically.

3. **Rotation curve consequence:** The force law implies a crossover from Keplerian (v\propto r^{-1/2}) to a flat asymptote
   [
   v_\infty=\sqrt{GM/H}.
   ]

4. **Transition scale:** The crossover occurs at (r\sim H) (quantified by the slope criterion), consistent with interpreting (H) as the geometric scale controlling the regime change.

5. **Independent computational confirmation:** The same force is computed via (i) image sums and (ii) mode/Bessel sums, and the two agree to floating precision once the image-tail remainder is included.

6. **Not point-mass-specific:** In the slab regime, the far-field law generalizes to axisymmetric extended sources via (M_{\rm enc}(r)), and the asymptote is controlled by (M_{\rm tot}) and (H).

---

## How to incorporate this into the paper structure

### Main text (core mechanism)

* Define the slab boundary value problem and the Neumann condition.
* Present the two representations (images and modes) at least briefly, emphasizing equivalence.
* Derive the two asymptotes (1/r^2) and (1/(Hr)).
* Derive the rotation curve consequence (v_\infty^2=GM/H).
* Define and interpret the transition scale (r\sim H).

### Appendix (validation)

* Include the critical normalization checks:

  * (r^2F\to 1) near field (in units (G=M=1)),
  * (HrF\to 1) far field.
* Include convergence tables vs truncations.
* Include cross-check table images vs modes.
* Include the log-slope plot and the knee estimate procedure.
* Include the far-field exponential disk check.
* Mention that relative-error discrepancy plots at (10^{-15}) reflect numerical noise floors.

### Separate “ideas/what-if” section

* Wake-force flattening as an alternative/additive mechanism (clearly labeled speculative or model-dependent).
* Any cosmological extrapolations built on the same slab Green’s-function logic.


---

# 4) Boundary sensitivity / leaky slab (paper-ready)

## 4.1 Why boundary conditions are the whole story (zero-mode lever)

The strict slab far-field law
\[
F(r)\sim \frac{GM}{Hr}
\qquad (r\gg H)
\]
exists because Neumann walls at \(z=\pm H\) admit a **true zero mode** \(k_0=0\) in the transverse spectrum. That zero mode reduces the far-field problem to an effectively 2D Green’s function (logarithmic potential), hence \(F\propto 1/r\).

Changing the wall physics changes the spectrum:

* **Neumann (sealed walls):** \(k_0=0\) → 2D/log potential → \(F\propto 1/r\).
* **Dirichlet (strong leakage / absorbing walls):** no zero mode → leading mode \(k\sim \pi/H\) → Yukawa-like decay \(\sim e^{-\pi r/H}\).
* **Robin (partially leaky):** zero mode lifted to small \(k_0>0\) → intermediate \(1/r\) window then screening at \(r\gtrsim 1/k_0\).

This is an honest robustness taxonomy: **flat rotation curves are a diagnostic of bulk confinement.**

## 4.2 Model A (fast): partial reflections via weighted images

Parameterize leakiness by damping the image stack:
\[
F_w(r) \propto \sum_{n=-\infty}^{\infty} w^{|n|}\,
\frac{r}{\left(r^2+(2nH)^2\right)^{3/2}},
\qquad 0<w\le 1.
\]

Interpretation:

* \(w=1\): perfect reflection → exact Neumann slab (your validated result).
* \(w<1\): partial leakage → effective number of images finite → Newton at small \(r\), approximate \(1/r\) regime over some range, and eventual departure from strict flattening at very large \(r\).

What to plot/quote:

* log-slope \(d\ln F/d\ln r\) vs \(r\) for \(w=1,0.99,0.95,0.9,\dots\),
* “flatness window” size: range of \(r\) where slope \(\approx -1\).

## 4.3 Model B (more formal): Robin BC and lifted zero mode (screening length)

Impose Robin walls:
\[
\partial_z\Phi+\alpha \Phi = 0 \quad \text{at}\quad z=\pm H.
\]

Even modes satisfy
\[
\tan(kH)=\frac{\alpha}{k}.
\]
For small \(\alpha\), the lowest eigenvalue scales as
\[
k_0 \simeq \sqrt{\alpha/H}.
\]

The brane zero-mode potential becomes the 2D massive Green’s function:
\[
\Phi_0(r)\propto K_0(k_0 r),
\qquad
F_0(r)\propto k_0 K_1(k_0 r)\sim \sqrt{\frac{\pi k_0}{2r}}e^{-k_0 r}.
\]

So the \(1/r\) force holds only for \(r\lesssim 1/k_0\), then screens.

## 4.4 Paper language (tight + falsifiable)

* The \(1/r\) far-field is the **Neumann zero mode** of a finite geometry.
* If confinement is imperfect, the model predicts a **finite** flat window and/or eventual screening.
* Measuring how “flat” galaxies remain at very large radii constrains the bulk leakage parameter (e.g. \(w\) or \(k_0\)).

---

# 5) Lensing / optics predictions of the slab regime (paper-ready)

## 5.1 Reuse the Paper II optics mapping

From `1pn_optics.tex`, weak-field deflection is
\[
\Delta\theta \simeq \int_{-\infty}^{+\infty}\nabla_\perp \ln N(r(z))\,dz,
\qquad r(z)=\sqrt{b^2+z^2}.
\]
In the weak-field mapping used there,
\[
\ln N(r)\simeq -\frac{2\Phi(r)}{c^2}.
\]

Only gradients matter, so the additive gauge constant in a log potential is irrelevant.

## 5.2 Deflection expressed directly in terms of the slab force

Let \(F(r)\) be the **inward** acceleration magnitude on the brane (the one computed by the slab code). Then
\[
\Delta\theta(b)
\simeq \frac{2}{c^2}\int_{-\infty}^{+\infty}F(r(z))\frac{b}{r(z)}\,dz.
\]

A numerically stable form uses \(z=b\tan u\) so \(r=b\sec u\):
\[
\Delta\theta(b)
= \frac{2b}{c^2}\int_{-\pi/2}^{+\pi/2}F(b\sec u)\,\sec u\,du.
\]

## 5.3 Two asymptotes (the headline)

### Near field: \(b\ll H\)

If \(F\simeq GM/r^2\), then
\[
\Delta\theta(b)\to \frac{4GM}{b c^2}.
\]

### Deep slab regime: \(b\gg H\)

If \(F\simeq GM/(Hr)\), then
\[
F(r)\frac{b}{r}\sim \frac{GM}{H}\frac{b}{r^2}
\quad\Rightarrow\quad
\Delta\theta(b)\to \frac{2\pi GM}{H c^2},
\]
independent of impact parameter.

Equivalently, using \(v_\infty^2=GM/H\),
\[
\boxed{\ \Delta\theta_\infty = \frac{2\pi v_\infty^2}{c^2}\ }.
\]

This gives a clean lensing–rotation consistency test.

(If the slab is leaky, the plateau becomes a finite window and then weakens/screened.)

## 5.4 What to plot for the appendix

Using your existing \(F_{\rm mode}(r)\) function:

1. \(\Delta\theta(b)\) vs \(b\), with overlays:
   * \(4GM/(b c^2)\) near-field,
   * \(2\pi GM/(H c^2)\) far-field plateau.
2. log-slope \(d\ln\Delta\theta/d\ln b\):
   * \(-1\) at small \(b\),
   * \(0\) at large \(b\) (plateau).
3. Boundary sensitivity lensing: repeat for several \(w\) (or \(k_0\)).

## 5.5 Minimal Mathematica implementation sketch

Assuming you already have `FMode[r_]`:

* `cVal = 1.0;`
* `Deflection[b_] := (2*b/cVal^2) * NIntegrate[FMode[b*Sec[u]]*Sec[u], {u, -Pi/2, Pi/2}]`

Then plot `Deflection[b]` across `b∈[0.1 H, 100 H]`.
