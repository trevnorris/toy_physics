# Moving-Throat PDE — Stage 250: Crossing-vs-Collapse / Goldilocks-Window Compiler from the Stage-248 Event Chain and the Relaxed Wall-Timescale Closure

## Status

**Exact within**

1. the Stage-248 dynamic event chain
   \[
   E=\frac{m_s}{2}\dot r^{\,2}+V_{\rm eff}^{(247)}(r),
   \]
2. the Stage-248 trigger-width packet
   \[
   \lambda_{\rm th}(E)=\left|\frac{E}{V'\bigl(r_+(E)\bigr)}\right|,
   \]
3. the Stage-245 non-rigid dressing-leg data \((g_{UV},\mu_\eta)\),
4. the Session-III characteristic-collapse closure for the active dressing leg,
5. and the Session-III proton-proxy / aligned-branch benchmark specialization.

This stage is a **timescale compiler**, not a new barrier law.
It does not modify the Stage-247 stationary lowered branch, and it does not add the Stage-251 damping/export kernel yet.  Its job is narrower: turn the relaxed branch into an explicit survival criterion of the form
\[
\text{crossing outruns collapse}.
\]

---

## Purpose

Stage 248 supplied the dynamic crossing branch:
\[
(r_{\rm peak},V_{\rm peak}),
\qquad
v_{\rm crit,new},
\qquad
r_\pm(E),
\qquad
\lambda_{\rm th}(E),
\qquad
\Xi_{\rm turn}(E).
\]

Stage 249 then attached the aligned-vs-anti-aligned hidden-channel export diagnostic to that same event chain.
What was still missing was the next reduced question actually used in Session III:

> on the already-lowered and dynamically traversable branch, when does the barrier-region transit finish before the relaxed dressing leg collapses?

That is exactly what Stage 250 compiles.

The output of this stage is the explicit **Goldilocks-window ledger**:

1. an exact event-chain transit-time integral,
2. the characteristic reduced crossing-time compiler,
3. the relaxed dressing-leg collapse-time compiler,
4. the stability ratio
   \[
   \mathcal S(E)=\frac{t_{\rm cross}(E)}{t_{\rm collapse}},
   \]
5. the lower safe edge in both energy and speed form,
6. the heavy-throat scaling theorem,
7. and the trigger-width sensitivity theorem that explains why the Session-III window was one-sided and geometry-sensitive.

---

## Provenance

This stage is the correct continuation of Stages 245–249.

- Stage 245 activated the non-rigid \(U/V\) lane and made the dressing leg into a live reduced degree of freedom.
- Stage 248 supplied the dynamic crossing branch, the peak energy \(V_{\rm peak}\), and the turning-point trigger width \(\lambda_{\rm th}\).
- Stage 249 fixed the aligned branch as the preferred hidden-channel export branch for the benchmark dynamic check, but deliberately left the scalar crossing law unchanged.
- Session III then compared barrier-region transit against dressing-leg collapse on that already-compiled relaxed branch.

So Stage 250 is not a fresh phenomenology layer.  It is the derivation-stack object corresponding to the Session-III “Goldilocks zone.”

---

## 0. Why this stage is needed

Before this step, the stack could already say:

- the lowered branch has a finite threshold speed,
- a classical crossing window exists relative to Coulomb,
- the turning-point branch carries the diagnostics \(\Xi_{\rm turn}\) and \(\lambda_{\rm th}\),
- and the aligned branch exports hidden helicity more efficiently than the anti-aligned branch in the chosen reduced closure.

But it still could **not** say:

- whether a successful crossing event survives its own relaxed dressing-leg drain,
- where the lower survival edge sits,
- why the safe set is one-sided on the chosen closure,
- or which quantities actually control that edge.

Stage 250 closes that gap.

---

## 1. Event-chain transit time from the Stage-248 branch

Let
\[
V(r):=V_{\rm eff}^{(247)}(r)
\]
be the Stage-247 lowered barrier front end and let the Stage-248 event chain satisfy
\[
E=\frac{m_s}{2}\dot r^{\,2}+V(r).
\]
Then along an inward same-charge trajectory,
\[
\dot r=-\sqrt{\frac{2}{m_s}\bigl(E-V(r)\bigr)}.
\]
So the exact path-time between any two radii on the already-compiled branch is
\[
\boxed{
T_{\rm traj}(E;r_a\to r_b)
=
\sqrt{\frac{m_s}{2}}
\int_{r_b}^{r_a}
\frac{dr}{\sqrt{E-V(r)}}.
}
\]

This is the exact event-chain transit object inherited from Stage 248.

### 1.1 Characteristic barrier-region transit compiler

Session III did not use the full path-time integral as its stress-test clock.
Instead, it compressed the active barrier region to one effective width
\[
\lambda_{\rm eff}
\]
and used the peak-over-barrier kinetic speed
\[
\boxed{
 v_{\rm bar}(E)=\sqrt{\frac{2(E-V_{\rm peak})}{m_s}}.
}
\]
Hence the characteristic crossing-time compiler is
\[
\boxed{
 t_{\rm cross}(E)
 =
 \frac{\lambda_{\rm eff}}{v_{\rm bar}(E)}
 =
 \lambda_{\rm eff}\sqrt{\frac{m_s}{2(E-V_{\rm peak})}}.
}
\]

This is exact **within the declared characteristic-width closure** once the active width \(\lambda_{\rm eff}\) has been chosen.

### 1.2 Monotonicity

For fixed \((m_s,\lambda_{\rm eff},V_{\rm peak})\),
\[
\frac{d t_{\rm cross}}{dE}
=
-\frac{\lambda_{\rm eff}\sqrt{m_s}}{2\sqrt{2}\,(E-V_{\rm peak})^{3/2}}<0.
\]
So higher incident energy shortens the barrier-region transit monotonically on this closure.

---

## 2. Relaxed dressing-leg collapse-time compiler

Stage 245 already turned the dressing leg into a live reduced degree of freedom.
Session III then adopted the local unstable normal form
\[
\mu_\eta\,\ddot V = g_{UV}\,\chi_{\rm peak}\,V,
\]
where

- \(\mu_\eta\) is the wall-inertia scale of the active dressing/geometry leg,
- \(g_{UV}\) is the non-rigid transfer-to-dressing coupling carried from Stage 245,
- and
  \[
  \chi_{\rm peak}:=\max\chi_\lambda(r)
  \]
  is the steepest logarithmic barrier gradient on the active branch.

The corresponding local growth rate is
\[
\Gamma_{\rm coll}=\sqrt{\frac{g_{UV}\chi_{\rm peak}}{\mu_\eta}},
\]
so the characteristic collapse time is
\[
\boxed{
 t_{\rm collapse}
 =
 \Gamma_{\rm coll}^{-1}
 =
 \sqrt{\frac{\mu_\eta}{g_{UV}\chi_{\rm peak}}}.
}
\]

This is the exact Stage-250 collapse compiler **within the declared local unstable-envelope closure**.

---

## 3. Stability ratio and lower safe edge

Define the characteristic survival ratio
\[
\boxed{
\mathcal S(E)=\frac{t_{\rm cross}(E)}{t_{\rm collapse}}.
}
\]
Substituting the previous two compilers gives
\[
\boxed{
\mathcal S(E)
=
\lambda_{\rm eff}
\sqrt{\frac{m_s g_{UV}\chi_{\rm peak}}{2\mu_\eta\,(E-V_{\rm peak})}}.
}
\]

### 3.1 Exact lower-edge theorem

The characteristic-safe branch is defined by
\[
\mathcal S(E)<1.
\]
Because \(\mathcal S(E)\) is monotone decreasing in \(E\), there is at most one lower survival edge and it is fixed by
\[
\mathcal S(E_{\rm edge})=1.
\]
Solving gives
\[
\boxed{
E_{\rm edge}
=
V_{\rm peak}
+
\frac{m_s}{\mu_\eta}
\frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2}.
}
\]
So the reduced Goldilocks window is
\[
\boxed{
E>E_{\rm edge}
}
\]
intersected with whatever energy band the scan or branch realization actually samples.

### 3.2 Why the window is one-sided on this closure

Under the present closure:

1. \(t_{\rm collapse}\) is energy-independent,
2. \(t_{\rm cross}(E)\) decreases monotonically with \(E\),
3. no additional high-energy collapse channel is turned on.

Therefore once \(E\) crosses the lower edge, the reduced timescale inequality stays satisfied for all larger energies on the same branch.
So the safe set is a **half-line**, not a closed island.
A finite upper edge can only appear if one adds extra energy dependence to the collapse law or an additional high-energy failure mechanism.  Those are outside Stage 250 and belong to later damping/export work.

---

## 4. Speed-space compiler

Let the reduced trajectory start at radius \(r_0\) with inward speed \(v_0\).
Then
\[
E=\frac{m_s}{2}v_0^2+V(r_0).
\]
Substitute this into the Stage-250 safe-edge law.

### 4.1 Reduced safe speed at fixed launch radius

The lower safe speed is
\[
\boxed{
 v_{\rm safe,min}
 =
 \sqrt{\frac{2(E_{\rm edge}-V(r_0))}{m_s}}.
}
\]
Using
\[
 v_{\rm crit,new}^2
 =
 \frac{2(V_{\rm peak}-V(r_0))}{m_s},
\]
this becomes the exact speed-space identity
\[
\boxed{
 v_{\rm safe,min}^2
 =
 v_{\rm crit,new}^2
 +
 \frac{\lambda_{\rm eff}^2 g_{UV}\chi_{\rm peak}}{\mu_\eta}.
}
\]
So Stage 250 turns the Goldilocks edge into a direct shift above the Stage-248 classical threshold speed.

### 4.2 Crossing-time compiler in speed form

At fixed launch radius,
\[
E-V_{\rm peak}=\frac{m_s}{2}\bigl(v_0^2-v_{\rm crit,new}^2\bigr),
\]
so the characteristic crossing time can be written exactly as
\[
\boxed{
 t_{\rm cross}(v_0)
 =
 \frac{\lambda_{\rm eff}}{\sqrt{v_0^2-v_{\rm crit,new}^2}}.
}
\]
The corresponding survival ratio is
\[
\boxed{
 \mathcal S(v_0)
 =
 \frac{\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}/\mu_\eta}}{\sqrt{v_0^2-v_{\rm crit,new}^2}}.
}
\]
Thus the reduced window is simply
\[
\boxed{
 v_0>v_{\rm safe,min}.
}
\]
Again, this is one-sided unless additional high-speed failure channels are added.

---

## 5. Heavy-throat scaling theorem

Session III emphasized the case
\[
\mu_\eta=\alpha m_s.
\]
Substituting into the stability ratio gives
\[
\boxed{
\mathcal S(E)
=
\frac{\sqrt{2}\,\lambda_{\rm eff}\sqrt{g_{UV}\chi_{\rm peak}}}{2\sqrt{\alpha}}
\,(E-V_{\rm peak})^{-1/2}.
}
\]
The lower edge becomes
\[
\boxed{
E_{\rm edge}
=
V_{\rm peak}
+
\frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2\alpha}.
}
\]
Therefore the simultaneous scaling of particle mass and wall inertia cancels out of the survival edge.
In particular,
\[
\frac{\partial E_{\rm edge}}{\partial m_s}=0
\qquad
\text{when}
\qquad
\mu_\eta=\alpha m_s.
\]
So heavy-throat scaling by itself does **not** move the Goldilocks edge if mass and wall inertia are scaled together.
What still matters is
\[
\lambda_{\rm eff},\qquad \chi_{\rm peak},\qquad g_{UV},\qquad \alpha.
\]

---

## 6. Trigger-width specialization from Stage 248

Session III used the Stage-248 turning-point trigger width as the active barrier width,
\[
\boxed{
\lambda_{\rm eff}=\lambda_{\rm th}\bigl(r_{\rm turn}\bigr),
\qquad
\lambda_{\rm th}(r_{\rm turn})
=
\left|\frac{V\bigl(r_{\rm turn}\bigr)}{V'\bigl(r_{\rm turn}\bigr)}\right|.
}
\]
This is the cleanest front-end specialization because it ties the Stage-250 timescale ledger directly to the Stage-248 event-chain diagnostic rather than to an unrelated geometric guess.

The same session also recorded the gradient trigger
\[
\chi_\lambda(r)=\lambda\,\bigl|\partial_r\ln V(r)\bigr|,
\]
and then used the steepest value on the active branch,
\[
\chi_{\rm peak}=\max \chi_\lambda(r),
\]
inside the collapse-time compiler.

So Stage 250 depends on the Stage-248/Session-II front end only through the pair
\[
\boxed{
\bigl(\lambda_{\rm eff},\chi_{\rm peak}\bigr).
}
\]
That is why the later width-sensitivity result is so important.

---

## 7. Session-III benchmark specialization

Use the Session-II/III reduced benchmark values
\[
V_{\rm peak}=3.42933112,
\qquad
V(r_0=5)=0.19999794,
\qquad
\lambda_{\rm eff}=0.42826825,
\]
\[
 g_{UV}=0.95,
 \qquad
 \chi_{\rm peak}=21.73204372,
 \qquad
 m_s=\mu_\eta=1836.15267343.
\]
Then:

### 7.1 Proton-proxy classical threshold speed

From the Stage-248 threshold-speed compiler,
\[
\boxed{
 v_{\rm crit,p}
 =
 \sqrt{\frac{2\bigl(V_{\rm peak}-V(r_0)\bigr)}{m_s}}
 \approx 0.05930851.
}
\]

### 7.2 Collapse time

The Stage-250 collapse compiler gives
\[
\boxed{
 t_{\rm collapse}
 =
 \sqrt{\frac{1836.15267343}{0.95\times 21.73204372}}
 \approx 9.43066476.
}
\]

### 7.3 Lower safe edge in energy

The heavy-throat edge is
\[
\boxed{
 E_{\rm safe,min}
 =
 3.42933112
 +
 \frac{0.95\times 21.73204372\times 0.42826825^2}{2}
 \approx 5.32265943.
}
\]

### 7.4 Lower safe edge in speed

Using the Stage-250 speed compiler,
\[
\boxed{
 v_{\rm safe,min}
 =
 \sqrt{\frac{2\bigl(E_{\rm safe,min}-0.19999794\bigr)}{1836.15267343}}
 \approx 0.07469791.
}
\]
Measured in proton-threshold units,
\[
\boxed{
\frac{v_{\rm safe,min}}{v_{\rm crit,p}}
\approx 1.25948037.
}
\]
So the Goldilocks edge sits about \(25.95\%\) above the proton-proxy classical threshold speed on this closure.

### 7.5 One-sided safe band on the sampled closure

The Session-III scan did not find an upper collapse edge before the top of the sampled band,
\[
E_{\rm max,scan}=80.93332737,
\]
so the reported safe window is the one-sided interval
\[
\boxed{
5.32265943\lesssim E_{\rm inc}\lesssim 80.93332737.
}
\]
The corresponding sampled speed band is
\[
\boxed{
0.07469791\lesssim v_0\lesssim 0.29654256.
}
\]

### 7.6 Dynamic cross-check with actual aligned trajectories

The same aligned branch carried actual barrier-region transit times between about
\[
0.204
\quad\text{and}\quad
4.054,
\]
all of which remained below
\[
 t_{\rm collapse}\approx 9.43066476.
\]
So the characteristic Stage-250 criterion is consistent with the sampled dynamic trajectories on the chosen aligned branch.

---

## 8. Trigger-width / steepness sensitivity theorem

Session III also repeated the same calculation using the raw model width
\[
\lambda=1
\]
instead of the trigger-width choice.
Then the branch recorded
\[
\chi_{\rm peak}^{\rm raw}\approx 50.74399964.
\]
The same Stage-250 formulas give
\[
\boxed{
 t_{\rm collapse}^{\rm raw}
 =
 \sqrt{\frac{1836.15267343}{0.95\times 50.74399964}}
 \approx 6.17163516,
}
\]
\[
\boxed{
 E_{\rm safe,min}^{\rm raw}
 =
 3.42933112
 +
 \frac{0.95\times 50.74399964\times 1^2}{2}
 \approx 27.53273095.
}
\]
So the safe edge shifts upward dramatically when the width/steepness packet is changed.

This is the exact reduced form of the Session-III sensitivity verdict:

- the dominant lever is not heavy mass by itself,
- it is the pair
  \[
  (\lambda_{\rm eff},\chi_{\rm peak}),
  \]
  i.e. the geometric steepness packet imported from the Stage-248 trigger-width analysis.

At the formula level this is obvious from
\[
E_{\rm edge}-V_{\rm peak}
\propto
\lambda_{\rm eff}^2\chi_{\rm peak}.
\]

---

## 9. What this stage accomplishes physically

Stage 250 turns the Session-III “Goldilocks” story into an explicit derivation-stack object.

### 9.1 It separates the exact event chain from the characteristic survival closure

The Stage-248 branch already determines the exact path-time integral
\(T_{\rm traj}\).
Stage 250 then adds the **declared characteristic-timescale reduction** used to stress-test survival on the relaxed branch:
\[
T_{\rm traj}\quad\leadsto\quad t_{\rm cross},
\qquad
V\text{-leg instability}\quad\leadsto\quad t_{\rm collapse}.
\]
This keeps the claim-status firewall explicit.

### 9.2 It proves the lower-edge structure algebraically

The Goldilocks edge is not a scan artifact.
Inside the declared closure it is forced by the exact algebra
\[
\mathcal S(E)<1
\iff
E>V_{\rm peak}+\frac{m_s}{\mu_\eta}\frac{g_{UV}\chi_{\rm peak}\lambda_{\rm eff}^2}{2}.
\]

### 9.3 It identifies the real control parameters

The lower survival edge is controlled by
\[
\lambda_{\rm eff},\qquad \chi_{\rm peak},\qquad g_{UV},\qquad \mu_\eta/m_s,
\]
not by heavy scaling alone.
That is exactly why the trigger-width choice mattered so much in the Session-III sweep.

### 9.4 It explains why the benchmark window was one-sided

As long as
\[
 t_{\rm collapse}=\text{constant in }E,
 \qquad
 t_{\rm cross}(E)\downarrow,
\]
there is only one lower edge.
Any upper edge requires new physics beyond the present closure.

---

## 10. What is still missing

Stage 250 is still **not** a completed moving-throat survival theorem.
Several items remain outside its scope.

1. The collapse-time law is a reduced characteristic closure, not yet a solved geometry-time output of the full moving-throat PDE.
2. The aligned-vs-anti-aligned selection enters only through the chosen benchmark branch; the scalar timescale compiler itself is orientation-blind.
3. The stage does not yet include the damping/export kernel that can stabilize slower cold events.  That belongs to Stage 251.
4. The physical-unit mapping of \(\lambda_{\rm eff},\chi_{\rm peak},\mu_\eta\) into a material-specific device remains downstream.

So the correct reading is:

- Stage 248 gave the dynamic crossing branch,
- Stage 249 tagged the preferred hidden-channel export branch,
- Stage 250 now gives the structural survival inequality on that relaxed branch,
- and Stage 251 must replace the purely characteristic collapse law by a microscopic damping/export kernel.

---

## 11. SymPy-backed status

The accompanying audit script verifies all of the following.

1. The exact event-chain transit integral derived from the Stage-248 energy law.
2. The characteristic crossing-time compiler from the active width \(\lambda_{\rm eff}\).
3. The unstable-leg growth rate and collapse-time compiler.
4. The exact stability ratio \(\mathcal S(E)\) and the general lower-edge formula.
5. The heavy-throat scaling theorem and cancellation of explicit mass dependence when \(\mu_\eta=\alpha m_s\).
6. The exact speed-space compiler
   \[
   v_{\rm safe,min}^2=v_{\rm crit,new}^2+\lambda_{\rm eff}^2 g_{UV}\chi_{\rm peak}/\mu_\eta.
   \]
7. The full Session-III benchmark specialization, including
   \(t_{\rm collapse}\), \(E_{\rm safe,min}\), \(v_{\rm safe,min}\), the threshold-speed ratio, the sampled speed band, and the raw-width sensitivity values.
8. The inequality check that the sampled aligned transit times remain below the characteristic collapse time on the chosen benchmark branch.

---

## 12. Immediate next step

The next clean derivation move is now sharply defined.

Do **not** reopen the whole barrier chain.
Instead:

1. keep the Stage-250 lower-edge compiler,
2. replace the phenomenological characteristic-collapse law by a microscopic wall/BdG/localized-Maxwell/mixed export kernel,
3. derive the actual damping/export timescale from that kernel,
4. and then recompile the cold-survival threshold without using an imposed envelope law.

That is exactly the Stage-251 theorem gate.
