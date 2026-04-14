
# Same-Charge Barrier Audit — Stage 004: Resonance/Linewidth Tradeoff, the Dispersive No-Free-Lunch Theorem, and the Linear Survival Window

## 0. Purpose

Stage 003 killed the naive dynamic story but left one narrow escape hatch alive:

> perhaps a **resonant dispersive** mixed-sector window can amplify the already-known short-range attractive families enough to matter before the same pole simply turns into absorptive pumping / leakage.

This stage turns that possibility into an exact line-shape audit.

The job is not to solve the full driven two-throat PDE.
It is much sharper:

1. isolate the local simple-pole normal form of the dynamic one-port mixed bundle,
2. derive the exact relation between the conservative dispersive gain and the absorptive load,
3. specialize that relation to the wall-like pole of the Stage-003 bundle,
4. and state the first honest linear survival criterion.

The result is again restrictive:

> at linear order, the **largest possible conservative barrier reshaping** occurs exactly where the absorptive load is of the same order.
> If one insists on a cleaner low-loss window, the maximum conservative enhancement is suppressed by a simple exact factor.

So after Stage 004, the dynamic corridor is no longer “maybe resonance helps.”
It is:

> only if the actual PDE branch supplies a pole with a sufficiently large residue-to-linewidth ratio can the linear mixed-sector route survive.

---

## 1. Frozen input carried forward

### 1.1 Dynamic one-port bundle from Stage 003

The Stage-003 dynamic reduced correction has the exact collinear-source form

\[
\mathfrak V_{\rm mix}(x,\omega)
=
-\frac12\,\chi_s(\omega)\,\mathcal S(x,\omega)^2,
\qquad
\chi_s(\omega)=\frac{\mathcal N_s(\omega)}{\Delta_\Pi(\omega)\,D_\Pi(\omega)}.
\]

For the primitive same-charge source families
\[
\mathcal S_Q(x)=\frac1{x^3},
\qquad
\mathcal S_Y(x)=\frac{e^{-2\kappa x}}{x},
\]
the spatial kernel class is already frozen:

\[
\mathfrak V_{\rm mix}(x,\omega)
=
-\frac12\left[
\frac{\mathcal C_6(\omega)}{x^6}
+
2\mathcal C_4(\omega)\frac{e^{-2\kappa x}}{x^4}
+
\mathcal C_2(\omega)\frac{e^{-4\kappa x}}{x^2}
\right].
\]

So the dynamic audit is **not** about inventing a new spatial family.
It is about what the pole structure can do to the coefficients of the existing short-range families.

### 1.2 Why the pole analysis is the right next move

The 4D/plasma ontology keeps the mixed-sector variables

\[
A_w,\qquad F_{\mu w},\qquad J^w
\]

alive beyond the strict far-field brane reduction, and the moving-throat / 5PN stack has already narrowed the remaining weak-axisymmetric grouped payload to the outgoing prefactor side. So once the static and first dynamic no-gos are accepted, the only linear corridor left is a mixed-sector **resonant** one.

---

## 2. Exact simple-pole normal form

Let

\[
F_\Pi(\omega):=\Delta_\Pi(\omega)\,D_\Pi(\omega)
\]

be the full dynamic denominator entering the reduced susceptibility.
Suppose the conservative branch has a simple real pole at \(\omega=\omega_*\):

\[
F_0(\omega_* )=0,
\qquad
F_0'(\omega_*)\neq 0,
\qquad
F_0:=F_\Pi\big|_{\Pi=0}.
\]

Define the outgoing-port sensitivity
\[
Z_* := -\partial_\Pi F_\Pi(\omega)\big|_{(\omega_*,\,\Pi=0)}.
\]

Then to first order in detuning
\[
\delta:=\omega-\omega_*,
\]
and to first order in the passive/outgoing port,
\[
F_\Pi(\omega)
=
F_0'(\omega_*)\,\delta
-
\Pi(\omega_*)\,Z_*
+
O(\delta^2,\Pi\delta,\Pi^2).
\]

If the port is passive,
\[
\Pi(\omega_*) = i\,\Gamma_*,
\qquad
\Gamma_*\ge 0,
\]
then every reduced coefficient near that pole has the universal Breit–Wigner form

\[
\chi_s(\omega)
\approx
\frac{A_*}{\delta-i\gamma_*},
\qquad
A_*:=\frac{\mathcal N_s(\omega_*)}{F_0'(\omega_*)},
\qquad
\gamma_*:=\frac{\Gamma_* Z_*}{|F_0'(\omega_*)|}.
\]

So the whole linear resonance problem collapses to two real scalars:

- the signed residue scale \(A_*\),
- the positive linewidth \(\gamma_*\).

---

## 3. Specialization to the wall-like pole

The Stage-003 one-port bundle is controlled by the reduced wall denominator
\[
D_\Pi(\omega)
=
K_B(\omega) - \frac{Q_\Pi(\omega)}{\Delta_\Pi(\omega)}.
\]

The exact Stage-003 derivative identity is
\[
\partial_\Pi D_\Pi(\omega) = -N(\omega),
\]
where
\[
N(\omega)=\frac{P(\omega)^2}{\Delta_\Pi(\omega)^2}.
\]

So at a **wall-like** conservative pole
\[
D_0(\omega_*)=0,
\qquad
\Delta_0(\omega_*)\neq 0,
\]
one has the local normal form

\[
D_\Pi(\omega)
=
D_0'(\omega_*)\,\delta
-
\Pi(\omega_*)\,N_*
+
O(\delta^2,\Pi\delta,\Pi^2),
\qquad
N_*:=N(\omega_*).
\]

Therefore the wall susceptibility itself is

\[
\chi_{qq}(\omega)
=
\frac1{D_\Pi(\omega)}
\approx
\frac{1}{D_0'(\omega_*)}\,
\frac1{\delta-i\gamma_{\rm wall}},
\qquad
\gamma_{\rm wall}:=\frac{\Gamma_*\,N_*}{|D_0'(\omega_*)|}.
\]

This is the first exact dynamic self-limiting statement of the audit:

> the same transfer factor \(N_*\) that helps the outgoing quadrupole bridge also widens the wall pole once the passive/outgoing port is restored.

So a larger transfer strength is not a free win. It simultaneously increases the linewidth that weakens the dispersive amplification.

---

## 4. Exact dispersive/absorptive tradeoff theorem

Now write the universal local line shape as

\[
\chi_*(\omega)=\frac{A_*}{\delta-i\gamma_*},
\qquad
\gamma_*>0.
\]

Rationalizing the denominator gives

\[
\chi_*(\omega)
=
A_*\frac{\delta+i\gamma_*}{\delta^2+\gamma_*^2}.
\]

So the conservative and absorptive pieces are exactly

\[
\Re \chi_*(\omega)=A_*\frac{\delta}{\delta^2+\gamma_*^2},
\qquad
\Im \chi_*(\omega)=A_*\frac{\gamma_*}{\delta^2+\gamma_*^2}.
\]

Introduce the dimensionless detuning ratio

\[
r:=\frac{|\delta|}{\gamma_*}.
\]

Then the line-shape magnitudes collapse to

\[
|\Re\chi_*|
=
\frac{|A_*|}{\gamma_*}\,
\frac{r}{1+r^2},
\qquad
|\Im\chi_*|
=
\frac{|A_*|}{\gamma_*}\,
\frac{1}{1+r^2},
\]
and therefore

\[
\boxed{
\frac{|\Re\chi_*|}{|\Im\chi_*|}=r.
}
\]

So the ratio of useful conservative reshaping to absorptive loading is **nothing but the detuning in linewidth units**.

### 4.1 Exact maximum of the conservative line shape

The dispersive factor
\[
f(r):=\frac{r}{1+r^2}
\]
has derivative
\[
f'(r)=\frac{1-r^2}{(1+r^2)^2}.
\]

So the exact maximum occurs at
\[
r=1,
\]
that is,
\[
|\delta|=\gamma_*.
\]

At that point
\[
\boxed{
\max_r |\Re\chi_*|
=
\frac{|A_*|}{2\gamma_*},
}
\]
and simultaneously
\[
\boxed{
|\Re\chi_*|=|\Im\chi_*|.
}
\]

This is the cleanest no-free-lunch theorem in the linear audit:

> the largest possible conservative dispersive enhancement appears precisely where the absorptive load is of the same size.

### 4.2 Exact low-loss bound

Suppose one demands that the absorptive part be at most a fraction \(\eta\) of the conservative part:

\[
|\Im\chi_*| \le \eta\,|\Re\chi_*|,
\qquad
0<\eta\le 1.
\]

Using \(|\Re|/|\Im|=r\), this is equivalent to
\[
r\ge \frac1\eta.
\]

Since \(f(r)\) decreases for \(r\ge 1\), the largest conservative magnitude allowed by that low-loss condition occurs at the boundary \(r=1/\eta\). Therefore

\[
\boxed{
\sup_{\,|\Im\chi_*|\le \eta|\Re\chi_*|}\,|\Re\chi_*|
=
\frac{|A_*|}{\gamma_*}\,
\frac{\eta}{1+\eta^2}.
}
\]

For small \(\eta\),
\[
\frac{\eta}{1+\eta^2}
=
\eta-\eta^3+O(\eta^5).
\]

So in a genuinely low-loss window, the best linear conservative enhancement scales only **linearly** with the allowed loss fraction.

This is the central Stage-004 theorem.

---

## 5. Barrier language and absorbed-power language

Take one of the already-frozen spatial families \(S_j(x)\) from Stage 003.
Near a simple pole its coefficient has the local form
\[
\chi_j(\omega)\approx \frac{A_j}{\delta-i\gamma_*}.
\]

Then the in-phase barrier reshaping contribution is

\[
U_j^{\rm disp}(x,\omega)
=
-\frac12\,\Re\chi_j(\omega)\,S_j(x)^2,
\]

while the out-of-phase absorbed-power diagnostic is

\[
\overline P_j(x,\omega)
=
-\omega\,\Im\mathfrak V_j(x,\omega)
=
\frac{\omega}{2}\,|\Im\chi_j(\omega)|\,S_j(x)^2.
\]

So the exact tradeoff becomes

\[
\boxed{
\frac{\overline P_j(x,\omega)}{\omega\,|U_j^{\rm disp}(x,\omega)|}
=
\frac{|\Im\chi_j|}{|\Re\chi_j|}
=
\frac1r
=
\frac{\gamma_*}{|\delta|}.
}
\]

This makes the physical reading completely transparent:

- **on resonance** \((\delta=0)\), the conservative reshaping vanishes and the response is purely absorptive;
- **at the dispersive optimum** \((|\delta|=\gamma_*)\), the power load and the conservative barrier term are locked one-to-one;
- **in a low-loss window** \((|\delta|\gg\gamma_*)\), the absorptive fraction is smaller, but the conservative gain is proportionally suppressed.

---

## 6. Quality-factor form of the same theorem

Introduce the local quality factor
\[
Q_*:=\frac{\omega_*}{2\gamma_*}.
\]

Then
\[
\frac{|\delta|}{\omega_*}=\frac{r}{2Q_*}.
\]

So if one imposes the same low-loss condition
\[
|\Im\chi_*| \le \eta |\Re\chi_*|,
\]
the exact detuning requirement becomes

\[
\boxed{
\frac{|\omega-\omega_*|}{\omega_*}
\ge
\frac{1}{2Q_*\eta}.
}
\]

This is useful because it tells us how demanding the “clean dispersive” corridor really is.

- If \(Q_*\) is modest, the required low-loss detuning is not a tiny perturbation of the resonance.
- If \(Q_*\) is very large, the corridor is narrow in absolute frequency, but the residue-to-linewidth scale \(|A_*|/\gamma_*\) can also become large enough to matter.

So the surviving linear question is not “does resonance exist?”
It is:

> does the actual moving-throat branch support a pole with sufficiently large \(Q_*\) and sufficiently large residue scale \(|A_*|\) that the bound above can beat the same-charge barrier threshold?

---

## 7. What survives after Stage 004

Stage 003 already proved that linear monochromatic driving never creates a new spatial family.
Stage 004 now adds the exact resonance theorem.

### 7.1 What is dead

The following picture is now strongly disfavored:

> tune a monochromatic mixed-sector drive to resonance and obtain a large conservative barrier reduction without paying a comparable pumping / leakage price.

That does not happen at linear order.

- On resonance the effect is purely absorptive.
- At the point of maximal conservative gain the absorptive and conservative pieces are equal in magnitude.
- In lower-loss windows the best achievable conservative gain is bounded by
  \[
  \frac{|A_*|}{\gamma_*}\frac{\eta}{1+\eta^2}.
  \]

So there is no linear resonant “free lunch.”

### 7.2 What is still alive

A linear dynamic corridor still survives, but only in the narrow form:

1. one of the already-known short-range families must already be spatially relevant in the tunneling region,
2. the actual PDE branch must supply a pole with sufficiently large \(|A_*|/\gamma_*|\),
3. and that branch must remain admissible under the same 5PN / outgoing constraints already frozen earlier.

So after Stage 004, the dynamic same-charge route lives or dies on a **residue-to-linewidth** question, not on a generic resonance slogan.

---

## 8. Best current verdict after Stage 004

The idea is still alive, but it has narrowed again.

- Stage 001 killed the naive Coulomb/KK story.
- Stage 002 killed the hope for a brand-new static mixed-sector spatial law.
- Stage 003 killed the linear outgoing-phase shortcut.
- Stage 004 now kills the naive “resonance fixes it” shortcut.

What survives is only this:

> a high-quality mixed-sector pole could still amplify one of the already-existing short-range attractive families enough to matter, but only if its residue-to-linewidth ratio is large enough to overcome the exact low-loss tradeoff bound.

That is a real theorem gate, not a slogan.

---

## 9. Immediate next step

The continuation point is now very sharp.

1. Keep the exact Stage-003 dynamic one-port bundle.
2. Choose an explicit primitive branch family for
   \[
   K,\ M,\ C,\ \varpi,\ \Omega_U,\ \Omega_W,\ G_U,\ G_W,\ R.
   \]
3. Compute the first actual pole locations \(\omega_*\), residues \(A_*\), and linewidths \(\gamma_*\).
4. Compare
   \[
   \frac{|A_*|}{\gamma_*}\frac{\eta}{1+\eta^2}
   \]
   against the barrier-softening size demanded by the reduced same-charge problem.

That is the smallest honest continuation point after the resonance/linewidth audit.
