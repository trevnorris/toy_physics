# Moving-Throat PDE — Stage 221: Resonance/Linewidth Tradeoff, the Dispersive No-Free-Lunch Theorem, and the Linear Survival Window

## Status

**Exact within the carried isotropic one-port wall/BdG/Maxwell/mixed closure, once one expands near a simple conservative pole**
\[
F_0(\omega_*)=0,
\qquad
F_0'(\omega_*)\neq 0,
\]
and, for the wall-like specialization,
\[
D_0(\omega_*)=0,
\qquad
\Delta_0(\omega_*)\neq 0.
\]

The local line-shape statements are exact to first order in detuning and passive/outgoing dressing after this normal form is taken. The accompanying SymPy audit verifies the line-shape algebra on a sign-fixed local slice
\[
F_0'(\omega_*)>0,
\qquad
D_0'(\omega_*)>0,
\]
and the absolute-value forms quoted below follow immediately by restoring magnitudes at the end.

This stage is the first post-Stage-220 resonance/linewidth insertion of actual mixed-port pole logic into the now-closed local mixed-ray ledger.
It does **not** solve the full driven two-throat moving PDE.
It computes the exact local resonance theorem implied by the one-port mixed bundle and reduces the linear dynamic survival question to one residue-to-linewidth inequality.

---

## Purpose

Stage 220 killed the naive dynamic story but left one narrow escape hatch alive:

> perhaps a **resonant dispersive** mixed-sector window can amplify the already-known short-range attractive families enough to matter before the same pole simply turns into absorptive pumping / leakage.

The job of Stage 221 is not to solve the full driven PDE.
It is much sharper:

1. isolate the local simple-pole normal form of the dynamic one-port mixed bundle,
2. derive the exact relation between conservative dispersive gain and absorptive load,
3. specialize that relation to the wall-like pole of the Stage-220 bundle,
4. and state the first honest **linear survival window**.

The main outputs are:

1. the exact local Breit–Wigner normal form of the reduced susceptibility near a simple conservative pole,
2. the exact wall-like specialization in which the same outgoing transfer factor that helps the quadrupole bridge also widens the wall pole,
3. the exact resonance tradeoff theorem
   \[
   \frac{|\Re\chi_*|}{|\Im\chi_*|}=\frac{|\delta|}{\gamma_*},
   \]
4. the exact no-free-lunch optimum
   \[
   \max |\Re\chi_*|=\frac{|A_*|}{2\gamma_*}
   \]
   at equal conservative and absorptive strength,
5. the exact low-loss bound
   \[
   \sup_{|\Im\chi_*|\le \eta|\Re\chi_*|}|\Re\chi_*|
   =
   \frac{|A_*|}{\gamma_*}\frac{\eta}{1+\eta^2},
   \]
6. and the resulting linear survival window
   \[
   \frac{|A_j|}{\gamma_*}\frac{\eta}{1+\eta^2}S_j(x)^2
   \ge
   2\,\Delta V_{\rm req}(x)
   \]
   as the first honest low-loss gate for any linear resonant same-charge claim.

So Stage 221 keeps the dynamic same-charge idea alive, but only as a **residue-to-linewidth** question on one of the already-known short-range families.

Script-backed status:
- `scripts/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py`
  checks the simple-pole normal form, Stage-220 derivative identity, wall-like
  specialization, exact tradeoff identities, and linear survival-window algebra;
  its constructive numeric slice is probe-only and not part of the proof path.
- `mathematica/moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_mathematica_audit.wl`
  mirrors the same symbolic theorem path in the second CAS and keeps the same
  numeric slice explicitly probe-only.

---

## 1. Frozen input carried forward

### 1.1 Dynamic one-port bundle from Stage 220

The Stage-220 dynamic reduced correction has the exact collinear-source form
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

So Stage 221 is **not** about inventing a new spatial family.
It is about what the local pole structure can do to the coefficients of the existing short-range families.

### 1.2 Why the pole analysis is the right next move

The parent `4+1` / plasma ontology keeps the mixed-sector variables
\[
A_w,
\qquad
F_{\mu w},
\qquad
J^w
\]
alive beyond the strict far-field brane reduction, and the moving-throat / 5PN stack has already narrowed the remaining weak-axisymmetric grouped payload to the outgoing-prefactor side.
So once the static and first dynamic no-gos are accepted, the only linear corridor left is a mixed-sector **resonant** one.

---

## 2. Exact simple-pole normal form

Let
\[
F_\Pi(\omega):=\Delta_\Pi(\omega)\,D_\Pi(\omega)
\]
be the full dynamic denominator entering the reduced susceptibility.
Suppose the conservative branch has a simple real pole at \(\omega=\omega_*\):
\[
F_0(\omega_*)=0,
\qquad
F_0'(\omega_*)\neq 0,
\qquad
F_0:=F_\Pi\big|_{\Pi=0}.
\]

Define the outgoing-port sensitivity
\[
Z_*:=-\partial_\Pi F_\Pi(\omega)\big|_{(\omega_*,\,\Pi=0)}.
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
\Pi(\omega_*)=i\,\Gamma_*,
\qquad
\Gamma_*\ge 0,
\]
then every reduced coefficient near that pole has the universal Breit–Wigner form
\[
\boxed{
\chi_s(\omega)
\approx
\frac{A_*}{\delta-i\gamma_*},
\qquad
A_*:=\frac{\mathcal N_s(\omega_*)}{F_0'(\omega_*)},
\qquad
\gamma_*:=\frac{\Gamma_* Z_*}{|F_0'(\omega_*)|}.
}
\]

So the whole linear resonance problem collapses to two real scalars:

- the signed residue scale \(A_*\),
- the positive linewidth \(\gamma_*\).

---

## 3. Specialization to the wall-like pole

The Stage-220 one-port bundle is controlled by the reduced wall denominator
\[
D_\Pi(\omega)
=
K_B(\omega)-\frac{Q_\Pi(\omega)}{\Delta_\Pi(\omega)}.
\]

The exact Stage-220 derivative identity is
\[
\partial_\Pi D_\Pi(\omega)=-N(\omega),
\qquad
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
\gamma_{\rm wall}:=\frac{\Gamma_* N_*}{|D_0'(\omega_*)|}.
\]

This is the first exact dynamic self-limiting statement of the audit:

> the same transfer factor \(N_*\) that helps the outgoing quadrupole bridge also widens the wall pole once the passive/outgoing port is restored.

So a larger transfer strength is not a free win.
It simultaneously increases the linewidth that weakens the dispersive amplification.

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
\Re\chi_*(\omega)=A_*\frac{\delta}{\delta^2+\gamma_*^2},
\qquad
\Im\chi_*(\omega)=A_*\frac{\gamma_*}{\delta^2+\gamma_*^2}.
\]

Introduce the dimensionless detuning ratio
\[
r:=\frac{|\delta|}{\gamma_*}.
\]
Then the line-shape magnitudes collapse to
\[
|\Re\chi_*|
=
\frac{|A_*|}{\gamma_*}\,\frac{r}{1+r^2},
\qquad
|\Im\chi_*|
=
\frac{|A_*|}{\gamma_*}\,\frac{1}{1+r^2},
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
obeys the exact identity
\[
f(r)-\frac12
=
-\frac{(r-1)^2}{2(1+r^2)}.
\]
So the exact maximum occurs at
\[
r=1,
\qquad
\text{i.e.}
\qquad
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
|\Im\chi_*|\le \eta\,|\Re\chi_*|,
\qquad
0<\eta\le 1.
\]
Using \(|\Re|/|\Im|=r\), this is equivalent to
\[
r\ge \frac1\eta.
\]
The exact comparison identity is
\[
\frac{r}{1+r^2}-\frac{\eta}{1+\eta^2}
=
-\frac{(r-\eta)(\eta r-1)}{(1+r^2)(1+\eta^2)}.
\]
So on the low-loss domain \(r\ge 1/\eta\) with \(0<\eta\le 1\), the right-hand side is nonpositive and the largest conservative magnitude allowed by that low-loss condition occurs at the boundary \(r=1/\eta\). Therefore
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
\frac{\eta}{1+\eta^2}=
\eta-\eta^3+O(\eta^5).
\]
So in a genuinely low-loss window, the best linear conservative enhancement scales only **linearly** with the allowed loss fraction.

This is the central Stage-221 theorem.

---

## 5. Barrier language and absorbed-power language

Take one of the already-frozen spatial families \(S_j(x)\) from Stage 220.
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

## 6. Quality-factor form and the exact linear survival window

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
|\Im\chi_*|\le \eta|\Re\chi_*|,
\]
the exact detuning requirement becomes
\[
\boxed{
\frac{|\omega-\omega_*|}{\omega_*}
\ge
\frac{1}{2Q_*\eta}.
}
\]

This says the clean dispersive corridor is only parametrically easy when the actual branch provides a large enough quality factor.

Now let \(\Delta V_{\rm req}(x)\) be the local barrier reduction required at radius \(x\) for one spatial family \(S_j(x)\) to matter.
Inside the same low-loss window, the exact Stage-221 bound implies
\[
|U_j^{\rm disp}(x,\omega)|
\le
\frac12\,\frac{|A_j|}{\gamma_*}\,
\frac{\eta}{1+\eta^2}
S_j(x)^2.
\]
Therefore a necessary low-loss **linear survival condition** is
\[
\boxed{
\frac{|A_j|}{\gamma_*}
\frac{\eta}{1+\eta^2}
S_j(x)^2
\ge
2\,\Delta V_{\rm req}(x)
}
\]
for at least one admissible pole and at least one detuning with the sign chosen so that \(-\Re\chi_j\,S_j(x)^2\) actually lowers the same-charge barrier.
Equivalently,
\[
\boxed{
\frac{|A_j|}{\gamma_*}
\ge
\frac{2\,\Delta V_{\rm req}(x)}{S_j(x)^2}
\frac{1+\eta^2}{\eta}.
}
\]

So the linear resonant corridor is now completely explicit: it lives or dies on whether the actual branch can provide a sufficiently large **residue-to-linewidth ratio** on one of the already-known short-range families.

---

## 7. What survives after Stage 221

Stage 220 already proved that linear monochromatic driving never creates a new spatial family.
Stage 221 now adds the exact resonance theorem.

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
3. the detuning must remain inside an admissible low-loss window set by the same quality-factor bound above,
4. and that branch must remain admissible under the same outgoing / grouped / branch constraints already frozen earlier.

So after Stage 221, the dynamic same-charge route lives or dies on a **residue-to-linewidth** question, not on a generic resonance slogan.

---

## 8. Best current verdict after Stage 221

The idea is still alive, but it has narrowed again.

- Stage 252 closed the local mixed-ray search sieve.
- Stage 253 killed the hope for a brand-new static mixed-sector spatial law.
- Stage 220 killed the linear outgoing-phase shortcut.
- Stage 221 now kills the naive “resonance fixes it” shortcut.

What survives is only this:

> a high-quality mixed-sector pole could still amplify one of the already-existing short-range attractive families enough to matter, but only if its residue-to-linewidth ratio is large enough to overcome the exact low-loss tradeoff bound.

That is a real theorem gate, not a slogan.

---

## 9. Script-backed status

The accompanying SymPy audit verifies:

1. the generic simple-pole reduction
   \[
   \chi_s(\omega)\approx \frac{A_*}{\delta-i\gamma_*},
   \]
   from the local denominator normal form
   \[
   F_\Pi(\omega)=F_0'(\omega_*)\,\delta-\Pi(\omega_*)Z_*,
   \]
   on a sign-fixed slice,
2. the exact Stage-220 derivative identity
   \[
   \partial_\Pi D_\Pi(\omega)=-N(\omega),
   \]
3. the wall-like pole specialization
   \[
   \chi_{qq}(\omega)\approx \frac{1}{D_0'(\omega_*)}\frac{1}{\delta-i\gamma_{\rm wall}},
   \]
4. the exact real/imaginary line-shape formulas,
5. the exact resonance optimum
   \[
   \max |\Re\chi_*|=\frac{|A_*|}{2\gamma_*}
   \]
   at \(r=1\),
6. the exact low-loss factorization identity
   \[
   \frac{r}{1+r^2}-\frac{\eta}{1+\eta^2}
   =-
   \frac{(r-\eta)(\eta r-1)}{(1+r^2)(1+\eta^2)},
   \]
7. the barrier/absorbed-power ratio
   \[
   \frac{\overline P_j}{\omega |U_j^{\rm disp}|}=\frac{\gamma_*}{|\delta|},
   \]
8. the quality-factor detuning bound
   \[
   \frac{|\omega-\omega_*|}{\omega_*}\ge \frac{1}{2Q_*\eta},
   \]
9. and the linear survival-window inequality in its residue-to-linewidth form.

Supporting file:
- `moving_throat_pde_stage221_resonance_linewidth_tradeoff_dispersive_no_free_lunch_theorem_and_linear_survival_window_sympy_audit.py`

---

## 10. Immediate next step

The continuation point is now very sharp.

1. Keep the exact Stage-220 dynamic one-port bundle.
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
