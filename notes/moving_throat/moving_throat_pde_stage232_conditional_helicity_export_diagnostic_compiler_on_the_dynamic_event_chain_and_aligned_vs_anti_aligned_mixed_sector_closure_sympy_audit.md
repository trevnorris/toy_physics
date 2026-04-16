# Moving-Throat PDE — Stage 232: Conditional Helicity-Export Diagnostic Compiler on the Dynamic Event Chain and the Aligned-vs-Anti-Aligned Mixed-Sector Closure

## Status

**Exact within**

1. the exact projected helicity identities already carried by the localized `4+1` Maxwell / plasma extension,
2. the Stage-231 dynamic event chain
   \[
   r=r_E(t),
   \qquad
   E=\frac{m_s}{2}\dot r^{\,2}+V_{\rm eff}^{(230)}(r),
   \]
3. the declared Session-II reduced mixed/vortical orientation closure,
4. and the Session-II benchmark specialization recorded in the barrier-session write-up.

This stage is **diagnostic rather than constitutive**.  It does not change the scalar event chain from Stage 231, and it does not upgrade the aligned-vs-anti-aligned comparison into a theorem about microscopic spin support.  It attaches an unresolved-helicity export ledger to the already-compiled crossing branch.

---

## Purpose

Stage 231 ended with the scalar event chain

\[
(r_{\rm peak},V_{\rm peak}),
\qquad
v_{\rm crit,new},
\qquad
r_\pm(E),
\qquad
I_{\rm new}(E),
\qquad
\Xi_{\rm turn}(E),
\qquad
\lambda_{\rm th}(E).
\]

What it deliberately did **not** yet include was the additional Session-II diagnostic that compared the two reduced magnetic/vortical orientation branches.

The next missing object was therefore:

> on the already-compiled dynamic crossing branch, how much unresolved / transverse-structure helicity is exported into the hidden sector, and how does that export differ between the aligned and anti-aligned reduced orientation closures?

That is exactly what Stage 232 does.

---

## Provenance

This stage is the correct continuation of Stage 231 for three reasons.

1. The parent/plasma theory already carries an exact projected-vs-resolved helicity ledger with a precise subscale-helicity transfer equation.
2. Stage 231 already supplied the dynamic event chain and the same turning-point diagnostics used in the Session-II run.
3. The barrier-session write-up explicitly states that, inside the chosen reduced mixed-sector closure, both aligned and anti-aligned branches reached contact, while the difference between them was the efficiency of unresolved helicity export rather than the existence of scalar crossing itself.

So Stage 232 is not a new scattering compiler.  It is a conditional **hidden-channel export diagnostic** layered on top of the Stage-231 path.

---

## 0. Why this stage is needed

Before this step, the derivation stack could say all of the following:

- the lowered branch has a finite threshold-speed window,
- the subbarrier turning points move relative to Coulomb,
- the WKB exponent is reduced,
- and the turning-point branch carries the diagnostics
  \(\Xi_{\rm turn}\) and \(\lambda_{\rm th}\).

But the stack still had no explicit derivation object for the Session-II statement

> aligned spins export more subscale helicity than anti-aligned ones on the chosen mixed-sector closure.

Stage 232 turns that statement into an auditable reduced compiler.

---

## 1. Exact subscale-helicity ledger carried by the parent theory

The projected helicity machinery is already exact in the plasma extension.
At fixed `w`,

\[
\partial_t(\mathbf A\!\cdot\!\mathbf B)
+\nabla_3\!\cdot\!(\Phi\,\mathbf B+\mathbf E\times\mathbf A)
=-2\,\mathbf E\!\cdot\!\mathbf B.
\]

Projecting in `w` gives

\[
\partial_t\overline h+\nabla_3\cdot\overline{\mathbf F}_h
=-2\,\overline{\mathbf E\cdot\mathbf B},
\qquad
\overline h\equiv\overline{\mathbf A\cdot\mathbf B}.
\]

The resolved brane helicity built from projected fields is

\[
h_{\rm res}=\overline{\mathbf A}\cdot\overline{\mathbf B},
\qquad
\partial_t h_{\rm res}+\nabla_3\cdot\mathbf F_{h,{\rm res}}
=-2\,\overline{\mathbf E}\cdot\overline{\mathbf B}.
\]

So the exact **subscale / unresolved helicity** is

\[
\boxed{
h_{\rm sub}
:=
\overline h-h_{\rm res}
=
\overline{\mathbf A'\cdot\mathbf B'}.
}
\]

Subtracting the two identities gives the exact transfer law

\[
\boxed{
\partial_t h_{\rm sub}+\nabla_3\cdot\mathbf F_{h,{\rm sub}}
=
-2\Bigl(\overline{\mathbf E\cdot\mathbf B}
-\overline{\mathbf E}\cdot\overline{\mathbf B}\Bigr)
=
-2\,\overline{\mathbf E'\cdot\mathbf B'}.
}
\]

So the hidden-sector export diagnostic is not ad hoc: it is exactly the projected covariance helicity ledger of the parent theory.

---

## 2. Volume-integrated event-chain diagnostic attached to Stage 231

Let `\mathcal V_{\rm br}` be the reduced brane control region used to monitor the event chain.  Define the integrated unresolved helicity

\[
H_{\rm sub}(t;E)
:=
\int_{\mathcal V_{\rm br}} h_{\rm sub}(\mathbf x,t;E)\,d^3x.
\]

Integrating the exact local equation gives

\[
\boxed{
\frac{dH_{\rm sub}}{dt}
=
-\Phi_h(t;E)-2\,\mathcal C_h(t;E),
}
\]

where

\[
\Phi_h(t;E)
:=
\oint_{\partial\mathcal V_{\rm br}} \mathbf F_{h,{\rm sub}}\cdot d\mathbf S,
\qquad
\mathcal C_h(t;E)
:=
\int_{\mathcal V_{\rm br}}\overline{\mathbf E'\cdot\mathbf B'}\,d^3x.
\]

Stage 231 already supplies the reduced dynamic path
\(r=r_E(t)\),
so Stage 232 simply evaluates the helicity diagnostic **along that same path**.

The Stage-231 turning-point tags

\[
\Xi_{\rm turn}(E)=\Xi_1\bigl(r_+(E)\bigr),
\qquad
\lambda_{\rm th}(E)=\left|\frac{E}{V'\bigl(r_+(E)\bigr)}\right|
\]

are therefore carried here only as path labels for the benchmark branch.
They do not change the helicity algebra.

---

## 3. Minimal aligned-vs-anti-aligned mixed-sector closure

Introduce the reduced orientation label

\[
\sigma\in\{+1,-1\},
\qquad
\sigma=+1\;\text{(aligned)},
\qquad
\sigma=-1\;\text{(anti-aligned)}.
\]

The minimal closure assumption is that the hidden-sector helicity export is linear in this two-state label:

\[
\Phi_{h,\sigma}(t;E)=\Phi_{h,0}(t;E)+\sigma\,\Phi_{h,1}(t;E),
\]
\[
\mathcal C_{h,\sigma}(t;E)=\mathcal C_{h,0}(t;E)+\sigma\,\mathcal C_{h,1}(t;E).
\]

Then the integrated export rate is exactly of the form

\[
\boxed{
\dot H_{{\rm sub},\sigma}(t;E)
=
\Gamma_0(t;E)+\sigma\,\Gamma_1(t;E),
}
\]

with

\[
\Gamma_0:=-\Phi_{h,0}-2\mathcal C_{h,0},
\qquad
\Gamma_1:=-\Phi_{h,1}-2\mathcal C_{h,1}.
\]

Define the normalized orientation asymmetry

\[
\boxed{
\alpha_h(t;E):=\frac{\Gamma_1(t;E)}{\Gamma_0(t;E)}.
}
\]

Whenever \(\Gamma_0>0\), the branch takes the canonical form

\[
\boxed{
\dot H_{{\rm sub},\sigma}(t;E)
=
\Gamma_0(t;E)\bigl[1+\sigma\,\alpha_h(t;E)\bigr].
}
\]

### 3.1 Positivity and preference conditions

For both branches to export positively, it is enough that

\[
\boxed{
\Gamma_0(t;E)>\bigl|\Gamma_1(t;E)\bigr|
\iff
\bigl|\alpha_h(t;E)\bigr|<1.
}
\]

For the aligned branch to dominate instantaneously, it is enough that

\[
\boxed{
\Gamma_1(t;E)>0
\iff
\alpha_h(t;E)>0.
}
\]

So the entire aligned-vs-anti-aligned comparison collapses to the sign and size of one reduced scalar \(\alpha_h\).

### 3.2 Separation from the scalar event chain

The scalar Stage-231 crossing path is still controlled only by

\[
E,
\qquad
V_{\rm eff}^{(230)}(r),
\qquad
r_\pm(E),
\qquad
I_{\rm new}(E).
\]

The orientation label `\sigma` enters only through the hidden-channel export functional \(H_{{\rm sub},\sigma}\).  That is why this stage is a diagnostic layered on top of Stage 231 rather than a replacement for it.

---

## 4. Exact peak and integrated ratio compilers

### 4.1 Instantaneous / peak ratio

At any common reference point on the same dynamic path,

\[
\boxed{
R_{\rm inst}(t;E)
:=
\frac{\dot H_{{\rm sub},+}(t;E)}{\dot H_{{\rm sub},-}(t;E)}
=
\frac{1+\alpha_h(t;E)}{1-\alpha_h(t;E)}.
}
\]

Therefore the orientation asymmetry is recovered exactly from the observed ratio by

\[
\boxed{
\alpha_h(t;E)
=
\frac{R_{\rm inst}(t;E)-1}{R_{\rm inst}(t;E)+1}.
}
\]

In particular, at the export peak,

\[
R_{\rm pk}
=
\frac{1+\alpha_{\rm pk}}{1-\alpha_{\rm pk}},
\qquad
\alpha_{\rm pk}:=\alpha_h(t_{\rm pk};E).
\]

### 4.2 Integrated ratio on a common event interval

Define the integrated exported unresolved helicity over a common event interval `I_E` on the same crossing branch,

\[
\mathcal H_\sigma(E)
:=
\int_{I_E}\dot H_{{\rm sub},\sigma}(t;E)\,dt.
\]

Then

\[
\mathcal H_+(E)-\mathcal H_-(E)
=
2\int_{I_E}\Gamma_1(t;E)\,dt,
\]
\[
\mathcal H_+(E)+\mathcal H_-(E)
=
2\int_{I_E}\Gamma_0(t;E)\,dt.
\]

Define the `\Gamma_0`-weighted mean asymmetry

\[
\boxed{
\bar\alpha_h(E)
:=
\frac{\int_{I_E}\Gamma_1(t;E)\,dt}{\int_{I_E}\Gamma_0(t;E)\,dt}.
}
\]

Then the integrated ratio obeys the exact same Möbius form

\[
\boxed{
R_{\rm int}(E)
:=
\frac{\mathcal H_+(E)}{\mathcal H_-(E)}
=
\frac{1+\bar\alpha_h(E)}{1-\bar\alpha_h(E)}.
}
\]

Hence

\[
\boxed{
\bar\alpha_h(E)
=
\frac{R_{\rm int}(E)-1}{R_{\rm int}(E)+1}.
}
\]

Any common overall export scale, including the Session-II parameter \(\eta_h\), multiplies both \(\mathcal H_+\) and \(\mathcal H_-\) equally and therefore cancels out of \(R_{\rm int}\).  So the aligned-vs-anti-aligned preference is a pure asymmetry diagnostic rather than a statement about the absolute size of the hidden-channel source.

### 4.3 Final-helicity ratio on equal-start runs

If the two runs start with the same initial unresolved helicity, and in particular if

\[
H_{{\rm sub},+}(t_{\rm in};E)=H_{{\rm sub},-}(t_{\rm in};E)=0,
\]

then the final values satisfy

\[
\boxed{
\frac{H_{{\rm sub},+}(t_{\rm out};E)}{H_{{\rm sub},-}(t_{\rm out};E)}
=
R_{\rm int}(E).
}
\]

So the benchmark may be checked either by integrated export or by final unresolved-helicity load.

---

## 5. Session-II benchmark specialization

The Session-II benchmark used the same reduced path family already compiled in Stage 231, with turning-point tags

\[
\Xi_{\rm turn}\approx 0.34437471,
\qquad
\lambda_{\rm th}\approx 0.42826825,
\]

and the above-threshold contact demonstration at

\[
v_0\approx 2.59221845.
\]

The reported helicity-export outputs were

\[
\max\dot H_{{\rm sub},+}
\approx 281.79830789,
\qquad
\max\dot H_{{\rm sub},-}
\approx 56.96878122,
\]

\[
H_{{\rm sub},+}^{\rm final}
\approx 20.58070146,
\qquad
H_{{\rm sub},-}^{\rm final}
\approx 5.00843357.
\]

### 5.1 Peak-ratio reconstruction

From the reported peak rates,

\[
R_{\rm pk}
=
\frac{281.79830789}{56.96878122}
\approx 4.94653917.
\]

So the reconstructed peak asymmetry is

\[
\boxed{
\alpha_{\rm pk}
=
\frac{R_{\rm pk}-1}{R_{\rm pk}+1}
\approx 0.6636699192.
}
\]

Thus the aligned branch carries about two-thirds of the instantaneous export asymmetry scale at the benchmark peak, while the anti-aligned branch still remains positive because \(|\alpha_{\rm pk}|<1\).

### 5.2 Integrated / final-ratio reconstruction

From the reported integrated ratio,

\[
R_{\rm int}\approx 4.10920923,
\]

the weighted mean asymmetry is

\[
\boxed{
\bar\alpha_h
=
\frac{R_{\rm int}-1}{R_{\rm int}+1}
\approx 0.6085499908.
}
\]

Using the reported final unresolved-helicity loads,

\[
\frac{H_{{\rm sub},+}^{\rm final}}{H_{{\rm sub},-}^{\rm final}}
=
\frac{20.58070146}{5.00843357}
\approx 4.10920923,
\]

which matches the integrated ratio to displayed precision.
So the benchmark is consistent with equal-start bookkeeping for the two branches.

### 5.3 Exact interpretation of the benchmark inequalities

The benchmark therefore satisfies

\[
0<\bar\alpha_h\approx 0.60855<\alpha_{\rm pk}\approx 0.66367<1.
\]

This has three immediate consequences.

1. **Both** branches export positive unresolved helicity throughout the representative active region.
2. The aligned branch exports strictly more than the anti-aligned branch.
3. The asymmetry is strongest near the export peak and averages down over the full trajectory.

So the reduced comparison is not “aligned crosses, anti-aligned fails.”
It is

> both branches can share the same scalar crossing corridor, but the aligned label unloads unresolved repulsive structure more efficiently into the hidden transverse sector.

That is exactly the diagnostic content the session write-up attributed to the run.

---

## 6. Why this stage is diagnostic rather than a spin theorem

This derivation deliberately stops at the reduced export ledger.
The label

\[
\sigma\in\{+1,-1\}
\]

is only the Session-II **orientation closure index** carried by the mixed/vortical diagnostic.  It is not yet a derived microscopic spin quantum number.

So Stage 232 proves only the following conditional statement:

> on the Stage-231 event chain and within the declared reduced mixed-sector closure, the aligned orientation label produces a larger unresolved-helicity export than the anti-aligned label whenever the reduced asymmetry scalar satisfies \(0<\alpha_h<1\).

It does **not** yet prove intrinsic spin-`1/2`, a microscopic Pauli sector, or a completed spin-support theorem of the parent defect ontology.

---

## 7. What this stage adds to the audit trail

After Stage 231 the stack could already say where the lowered same-charge branch turned, tunneled, and crossed.

After Stage 232 the stack can additionally say:

1. the unresolved-helicity export ledger is anchored to an exact parent-theory projection identity,
2. the aligned-vs-anti-aligned diagnostic is a one-scalar mixed/vortical asymmetry problem,
3. the instantaneous and integrated preference ratios are exactly invertible,
4. the Session-II benchmark lies strictly inside the “both branches export, aligned dominates” regime,
5. and the comparison is a **hidden-channel unloading diagnostic**, not yet a theorem of microscopic spin support.

That is precisely the right status for the next branch in the derivation stack.

---

## 8. Benchmark summary packet

For the representative Session-II branch,

\[
\boxed{
R_{\rm pk}\approx 4.94653917,
\qquad
\alpha_{\rm pk}\approx 0.66366992,
}
\]

\[
\boxed{
R_{\rm int}\approx 4.10920923,
\qquad
\bar\alpha_h\approx 0.60854999,
}
\]

with matching final-load ratio

\[
\boxed{
\frac{H_{{\rm sub},+}^{\rm final}}{H_{{\rm sub},-}^{\rm final}}
\approx 4.10920923.
}
\]

So the Stage-232 benchmark packet is

\[
\boxed{
\bigl(\Xi_{\rm turn},\lambda_{\rm th},R_{\rm pk},R_{\rm int},\alpha_{\rm pk},\bar\alpha_h\bigr)
=
\bigl(0.34437471,\;0.42826825,\;4.94653917,\;4.10920923,\;0.66366992,\;0.60854999\bigr).
}
\]

---

## 9. Immediate next step

The next stage should **not** reinterpret this diagnostic as a spin theorem.
The right continuation is instead to ask how this hidden-channel export diagnostic interacts with the structural survival problem:

- does the branch that exports more unresolved helicity also survive long enough to use that unloading advantage,
- and how does that preference feed into the crossing-vs-collapse / Goldilocks compiler that compares transit against dressing-leg failure?

That is the correct bridge into the next stage.
