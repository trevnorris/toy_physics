# Moving-Throat PDE — Stage 234: Microscopic Damping/Export Kernel Replacing the Phenomenological `V`-Leg Envelope Law

## Status

**Exact within**

1. the Stage-228 non-rigid mouth/dressing packet with active dressing variable `V`,
2. the Stage-233 local unstable-leg coefficient
   \[
   \kappa_V:=g_{UV}\chi_{\rm peak},
   \qquad
   s_0:=\sqrt{\kappa_V/\mu_\eta}=1/t_{\rm collapse,0},
   \]
3. the moving-throat scalar derivative-coupled outgoing lane,
4. the Stage-12/13 selected mixed outgoing quadrupole coefficient,
5. linear projection of the active `V` leg onto those outgoing lanes,
6. and the causal Laplace-domain growth test.

This stage is a **microscopic export-kernel compiler**, not a materials theorem and not a full nonlinear damping solve.
It replaces the Session-IV phenomenological envelope law
\[
\gamma_{\rm tot}\,\dot V
\]
by the first actual odd export kernel implied by the wall/BdG/localized-Maxwell/mixed reduction.

---

## Purpose

Stage 233 turned the relaxed branch into a Goldilocks-window criterion,
\[
\mathcal S(E)=\frac{t_{\rm cross}(E)}{t_{\rm collapse}},
\]
but it still used the Session-IV effective shedding law
\[
\gamma_{\rm tot}=\gamma_{\rm vac}+\gamma_{\rm lattice},
\qquad
 t_{\rm collapse}^{\rm damped}=\frac{1}{\gamma_{\rm crit}-\gamma_{\rm tot}}
 \ \ (\gamma_{\rm tot}<\gamma_{\rm crit}).
\]
The session write-up was explicit that this was only an **envelope closure**, not a microscopic theorem of the completed PDE.

So the next honest question is:

> what is the first actual passive/export law seen by the active dressing leg once the wall/BdG/localized-Maxwell/mixed outgoing channels are projected back onto `V`?

That is exactly what Stage 234 compiles.

The main results are:

1. the active `V` leg inherits a **cubic-plus-quintic** odd export kernel,
2. the selected mixed quadrupole lane contributes the exact quintic coefficient
   \[
   \Gamma_5=\Pi_{V-}^2\,\Gamma_{5,-},
   \qquad
   \Gamma_{5,-}=\frac{a^5}{27c_s^5}P_{0,-},
   \]
3. the derivative-coupled scalar outlet contributes the exact cubic coefficient
   \[
   \Gamma_3=\Pi_{V0}^2\,\Gamma_{3,0},
   \qquad
   \Gamma_{3,0}=\gamma_1\eta_0^2\frac{\Omega_{U,0}^4}{\Delta_0^2},
   \]
4. the resulting microscopic collapse law is not Ohmic and does **not** admit a finite analogue of the Session-IV unconditional threshold `\gamma_{\rm crit}`,
5. but it does give an exact event-level safe surface replacing `\gamma_{\rm safe}`.

So Stage 234 is the right bridge between the reduced barrier session and the actual moving-throat outgoing machinery.

---

## Provenance

This stage pulls together four earlier threads.

- Stage 228 activated the non-rigid `U/V` lane and made `V` a live reduced dressing variable.
- Stage 233 converted that live leg into the Goldilocks-window question
  \[
  \text{crossing outruns collapse}.
  \]
- The moving-throat outgoing program already proved that a derivative-coupled scalar outlet first appears at odd order `\omega^3`, while the selected mixed quadrupole outlet appears at odd order `\omega^5`.
- The session write-up then stated explicitly that the remaining gap was a **higher-fidelity damping/export closure** replacing the phenomenological envelope law.

So Stage 234 is not a speculative add-on. It is the direct derivation-stack replacement for the Session-IV `\gamma_{\rm tot}` closure.

---

## 0. Why this stage is needed

Before this step, the stack could already say:

- the relaxed branch lowers the near-contact barrier,
- a dynamic crossing corridor exists,
- aligned orientation exports hidden structure more efficiently on the chosen closure,
- and a phenomenological shedding rate can rescue a slow cold event.

But it still could **not** say:

- what odd export operator the active `V` leg actually sees,
- how that operator relates to the selected mixed outgoing branch,
- whether the finite Session-IV `\gamma_{\rm crit}` survives microscopic replacement,
- or how the later heat-partition stage should be written without inserting a by-hand damping law.

Stage 234 closes that gap.

---

## 1. Microscopic odd export channels seen by the active `V` leg

Write the active dressing variable as its projection onto the first outgoing scalar and selected outgoing quadrupole channels,
\[
V(t)=\Pi_{V0}\,q_0(t)+\Pi_{V-}\,q_-(t)+\cdots.
\]
Here:

- `q_0` is the derivative-coupled scalar outlet,
- `q_-` is the selected mixed quadrupole wall/outgoing mode,
- `\Pi_{V0}` and `\Pi_{V-}` are the active-leg projection amplitudes.

### 1.1 Derivative-coupled scalar outlet

In the reduced scalar mixed-port lane, take the first derivative coupling
\[
 g_{W,0}(\omega)=\eta_0\,\omega,
 \qquad g_{A,0}(\omega)=0,
\]
with outgoing scalar port law
\[
 \Pi_0^{\rm out}(\omega)=i\gamma_1\omega+O(\omega^3).
\]
From the Stage-4 transfer factor,
\[
 N_0(\omega)=
 \frac{\bigl(A_0(\omega)g_{W,0}(\omega)\bigr)^2}{\bigl(A_0(\omega)W_0(\omega)-R_0^2\bigr)^2},
\]
where
\[
A_0(\omega)=\Omega_{U,0}^2-\omega^2,
\qquad
W_0(\omega)=\Omega_{W,0}^2-\omega^2,
\qquad
\Delta_0:=\Omega_{U,0}^2\Omega_{W,0}^2-R_0^2,
\]
one gets the exact low-frequency expansion
\[
N_0(\omega)=\eta_0^2\frac{\Omega_{U,0}^4}{\Delta_0^2}\,\omega^2+O(\omega^4).
\]
So the wall-operator odd correction begins at cubic order,
\[
\delta D_0^{\rm odd}(\omega)
=
-i\Gamma_{3,0}\,\omega^3+O(\omega^5),
\]
with exact coefficient
\[
\boxed{
\Gamma_{3,0}
=
\gamma_1\eta_0^2\frac{\Omega_{U,0}^4}{\Delta_0^2}.
}
\]
Projecting this outlet onto the active dressing leg gives
\[
\boxed{
\Gamma_3=\Pi_{V0}^2\,\Gamma_{3,0}.
}
\]

### 1.2 Selected mixed quadrupole outlet

From the selected-mode outgoing branch,
\[
\delta D_-^{\rm odd}(\omega)
=
-i\Gamma_{5,-}\,\omega^5+O(\omega^7),
\]
with exact selected-mode coefficient
\[
\boxed{
\Gamma_{5,-}
=
\frac{a^5}{27c_s^5}P_{0,-},
\qquad
P_{0,-}=\frac{\beta_0 s_-}{\lambda_-}.
}
\]
Projecting onto the active `V` leg gives the exact quintic export coefficient
\[
\boxed{
\Gamma_5=\Pi_{V-}^2\,\Gamma_{5,-}
=
\Pi_{V-}^2\frac{a^5}{27c_s^5}P_{0,-}.
}
\]

### 1.3 First microscopic export kernel

So the first active-leg odd export kernel already implied by the moving-throat reduction is
\[
\boxed{
\Sigma_{\rm exp}(\omega)
=
-i\Gamma_3\omega^3-i\Gamma_5\omega^5+O(\omega^7),
\qquad
\Gamma_3,\Gamma_5\ge 0.
}
\]
This is the stage’s first main theorem.

It says that the phenomenological Session-IV `\gamma_{\rm tot}\dot V` law is not the first microscopic export law of the branch.
The first microscopic law is **super-Ohmic**: cubic if the derivative-coupled scalar lane is active, quintic if the selected quadrupole lane is the first live outlet.

---

## 2. Causal equation for the active `V` leg

Let
\[
\kappa_V:=g_{UV}\chi_{\rm peak},
\qquad
s_0:=\sqrt{\kappa_V/\mu_\eta}=\frac{1}{t_{\rm collapse,0}},
\]
be the undamped Stage-233 unstable-leg coefficient and growth rate.
Absorb the even low-frequency conservative renormalizations into the effective `\mu_\eta` and `\kappa_V` notation.
Then the first microscopic active-leg equation is
\[
\boxed{
\mu_\eta\,\ddot V-\kappa_V V
+\Gamma_3 V^{(3)}
-\Gamma_5 V^{(5)}
=\mathcal S_V(t).
}
\]
Here `\mathcal S_V(t)` denotes any retained forcing from leakage/work or branch transport, but the stability analysis below uses the homogeneous equation.

Equivalently, in Laplace space with growth variable `s`,
\[
\boxed{
\mathcal D_V(s)
=
\mu_\eta s^2-\kappa_V+\Gamma_3 s^3+\Gamma_5 s^5.
}
\]
So the exact microscopic export kernel is
\[
\boxed{
K_{\rm exp}(s)=\Gamma_3 s^3+\Gamma_5 s^5.
}
\]

This replaces the Session-IV envelope scalar `\gamma_{\rm tot}`.

---

## 3. Exact passive power identities

The active-leg odd force is
\[
F_{\rm odd}=\Gamma_3 V^{(3)}-\Gamma_5 V^{(5)}.
\]
Its power input is
\[
\dot V F_{\rm odd}
=
\Gamma_3\dot V V^{(3)}-\Gamma_5\dot V V^{(5)}.
\]
Using the exact Schott identities
\[
\dot q\,q^{(3)}=\frac{d}{dt}(\dot q\,\ddot q)-(\ddot q)^2,
\]
\[
\dot q\,(-q^{(5)})
=-\frac{d}{dt}(\dot q\,q^{(4)}-\ddot q\,q^{(3)})-(q^{(3)})^2,
\]
one gets
\[
\boxed{
\dot V F_{\rm odd}
=
\frac{d}{dt}\mathcal S_{\rm odd}
-\Gamma_3\,\ddot V^{\,2}
-\Gamma_5\,\dddot V^{\,2},
}
\]
with Schott storage
\[
\boxed{
\mathcal S_{\rm odd}
=
\Gamma_3\dot V\ddot V
-\Gamma_5\bigl(\dot V V^{(4)}-\ddot V\dddot V\bigr).
}
\]
So after the total derivative is removed, the microscopic export power is positive-definite:
\[
\boxed{
\mathcal P_{\rm exp}
=
\Gamma_3\,\ddot V^{\,2}+
\Gamma_5\,\dddot V^{\,2}\ge 0.
}
\]

This is the second main theorem of the stage.

It gives the first actual wall/support/mixed export ledger for the active dressing leg.

---

## 4. Collapse law from the microscopic kernel

For homogeneous growth, set
\[
V(t)=V_0 e^{st},
\qquad s>0.
\]
Substituting into the microscopic equation gives the exact characteristic polynomial
\[
\boxed{
F(s):=\Gamma_5 s^5+\Gamma_3 s^3+\mu_\eta s^2-\kappa_V=0.
}
\]

### 4.1 Exact uniqueness of the positive growth root

For
\[
\Gamma_3,\Gamma_5,\mu_\eta,\kappa_V>0,
\]
one has
\[
F(0)=-\kappa_V<0,
\qquad
\lim_{s\to\infty}F(s)=+\infty,
\]
and
\[
F'(s)=5\Gamma_5 s^4+3\Gamma_3 s^2+2\mu_\eta s>0
\qquad (s>0).
\]
Therefore:
\[
\boxed{
\text{for every finite }\Gamma_3,\Gamma_5\ge0,
\text{ there is exactly one positive real growth root }s_+.
}
\]

This gives the key structural verdict:

> the minimal microscopic export kernel does **not** admit a finite analogue of the Session-IV unconditional threshold `\gamma_{\rm crit}`.

The positive growth root can be reduced, but not removed, by any finite passive cubic/quintic coefficient on this minimal closure.

### 4.2 Small-kernel slowdown formula

For weak export,
\[
\Gamma_3,\Gamma_5\ll 1,
\]
expand around
\[
 s_0=\sqrt{\kappa_V/\mu_\eta}.
\]
Then the unique positive root is
\[
\boxed{
 s_+
 =
 s_0-
 \frac{\Gamma_3 s_0^2+\Gamma_5 s_0^4}{2\mu_\eta}
 +O(\Gamma^2).
 }
\]
So both microscopic channels slow the collapse, and the quintic channel is weighted by the extra factor `s_0^2` relative to the cubic one.

---

## 5. Exact event-level safe surface replacing `\gamma_{\rm safe}`

Let
\[
 s_c:=\frac{1}{t_{\rm cross}}
\]
be the Stage-233 crossing-rate target for a chosen event.
Since `F(s)` is strictly increasing for `s>0`, the event is safe iff
\[
 s_+\le s_c
 \iff
 F(s_c)\ge 0.
\]
Therefore the exact microscopic safe condition is
\[
\boxed{
\Gamma_3 s_c^3+\Gamma_5 s_c^5
\ge
\mu_\eta\bigl(s_0^2-s_c^2\bigr).
}
\]
This is the stage’s third main theorem.

It is the direct microscopic replacement for the Session-IV `\gamma_{\rm safe}` inequality.

### 5.1 Normalized form

Define the mass-normalized coefficients
\[
\widehat\Gamma_3:=\frac{\Gamma_3}{\mu_\eta},
\qquad
\widehat\Gamma_5:=\frac{\Gamma_5}{\mu_\eta}.
\]
Then the safe surface becomes the exact half-plane
\[
\boxed{
\widehat\Gamma_3+s_c^2\widehat\Gamma_5
\ge
\frac{s_0^2-s_c^2}{s_c^3}.
}
\]
So the microscopic event-safe region is linear in the normalized cubic/quintic export strengths.

### 5.2 One-channel limits

If only the cubic lane is active,
\[
\boxed{
\widehat\Gamma_{3,\rm safe}
=
\frac{s_0^2-s_c^2}{s_c^3}.
}
\]
If only the selected quadrupole quintic lane is active,
\[
\boxed{
\widehat\Gamma_{5,\rm safe}
=
\frac{s_0^2-s_c^2}{s_c^5}.
}
\]
Hence
\[
\frac{\widehat\Gamma_{5,\rm safe}}{\widehat\Gamma_{3,\rm safe}}=\frac{1}{s_c^2}.
\]
So on slow cold events the quintic branch is less efficient than the cubic branch by the exact factor `1/s_c^2`.

---

## 6. Session-IV benchmark specialization

The Session-IV cold aligned event used the characteristic data
\[
 t_{\rm cross}\approx 1.82169718,
 \qquad
 t_{\rm collapse,0}\approx 0.14402764,
 \qquad
 \gamma_{\rm crit}\approx 6.94311167,
\]
so, using the displayed Session-IV growth-rate benchmark,
\[
 s_c=\frac{1}{t_{\rm cross}}\approx 0.5489386551,
\qquad
 s_0\approx \gamma_{\rm crit}\approx 6.94311167.
\]
(The inverse of the rounded printed `t_{\rm collapse,0}` agrees only to the displayed report precision.)
Therefore
\[
 s_c^2\approx 0.3013336471,
\qquad
 \frac{s_0^2-s_c^2}{s_c^3}\approx 289.61004918.
\]
So the exact benchmark safe surface becomes
\[
\boxed{
\widehat\Gamma_3+0.3013336471\,\widehat\Gamma_5
\ge 289.61004918.
}
\]
In particular,
\[
\boxed{
\widehat\Gamma_{3,\rm safe}\approx 289.61004918,
\qquad
\widehat\Gamma_{5,\rm safe}\approx 961.09429528.
}
\]
So for the default cold event, a pure quintic selected-mode outlet would need about
\[
\frac{961.09429528}{289.61004918}\approx 3.318
\]
times the normalized strength of a pure cubic derivative-coupled scalar outlet to achieve the same event-level rescue.

This gives the first clean microscopic reading of the Session-IV damping story.
The phenomenological envelope closure did not distinguish these super-Ohmic outlet classes; the microscopic kernel does.

---

## 7. Channel-resolved vacuum/lattice export decomposition

The microscopic kernel can now be split by channel rather than by a by-hand `3:1` damping partition.
Write
\[
\Gamma_3=\Gamma_3^{\rm vac}+\Gamma_3^{\rm lat},
\qquad
\Gamma_5=\Gamma_5^{\rm vac}+\Gamma_5^{\rm lat},
\]
with all four coefficients nonnegative.
Then the exact exported-power split is
\[
\boxed{
\mathcal P_{\rm vac}
=
\Gamma_3^{\rm vac}\ddot V^{\,2}+
\Gamma_5^{\rm vac}\dddot V^{\,2},
}
\]
\[
\boxed{
\mathcal P_{\rm lat}
=
\Gamma_3^{\rm lat}\ddot V^{\,2}+
\Gamma_5^{\rm lat}\dddot V^{\,2}.
}
\]
Therefore the integrated exported energies are
\[
\boxed{
E_{\rm vac}
=
\int dt\,\bigl(\Gamma_3^{\rm vac}\ddot V^{\,2}+
\Gamma_5^{\rm vac}\dddot V^{\,2}\bigr),
}
\]
\[
\boxed{
E_{\rm lat}
=
\int dt\,\bigl(\Gamma_3^{\rm lat}\ddot V^{\,2}+
\Gamma_5^{\rm lat}\dddot V^{\,2}\bigr).
}
\]
So Stage 235 no longer needs a phenomenological heat partition. It can start from an exact kernel-level one.

---

## 8. What this stage achieves physically

Stage 234 changes the barrier thread in four concrete ways.

### 8.1 It replaces `\gamma_{\rm tot}` by an actual microscopic export operator

The first honest `V`-leg export law is not Ohmic. It is
\[
K_{\rm exp}(s)=\Gamma_3 s^3+\Gamma_5 s^5,
\]
with the cubic coefficient coming from derivative-coupled scalar export and the quintic coefficient coming from the selected mixed quadrupole outlet.

### 8.2 It explains why the Session-IV `\gamma_{\rm crit}` was closure-specific

Because the microscopic characteristic polynomial always has one positive real growth root for finite passive coefficients, the finite unconditional threshold
\[
\gamma_{\rm crit}\approx 6.94311167
\]
is not a theorem of the minimal microscopic kernel. It is a property of the Session-IV envelope closure.

### 8.3 It preserves an exact event-safe criterion

Although there is no finite microscopic unconditional threshold, there is still an exact event-level safe half-plane,
\[
\widehat\Gamma_3+s_c^2\widehat\Gamma_5
\ge
\frac{s_0^2-s_c^2}{s_c^3},
\]
which is the correct microscopic replacement for `\gamma_{\rm safe}`.

### 8.4 It opens the exact heat-partition path needed by Stage 235

The exported power now splits channel by channel into vacuum and lattice pieces without inserting an ad hoc partition ratio.

---

## 9. What is still missing

This stage is still not the full damping theorem of the completed PDE.
The remaining open objects are now sharply identified:

1. the actual projection amplitudes `\Pi_{V0}` and `\Pi_{V-}` on the realized branch,
2. the actual scalar derivative-coupled port coefficient `\gamma_1` on that branch,
3. the actual selected-mode prefactor `P_{0,-}` and therefore `\Gamma_{5,-}`,
4. and the actual vacuum-vs-lattice split of those microscopic coefficients.

So the active bottleneck has moved from “invent a damping law” to “compute the realized branch coefficients of the microscopic export kernel.”

---

## 10. Immediate next step

The next stage is now well defined.

1. Keep the exact Stage-234 kernel
   \[
   K_{\rm exp}(s)=\Gamma_3 s^3+\Gamma_5 s^5.
   \]
2. Use the exact power formulas above to compute
   \[
   E_{\rm vac},\qquad E_{\rm lat},
   \]
   along a chosen cold event.
3. Pull the microscopic coefficients into the condensed-matter map rather than using the Session-IV `3:1` partition.

That is precisely the Stage-235 heat-partition / cold-survival compiler.
