# Moving-Throat PDE — Stage 235: Vacuum-vs-Lattice Heat Partition and Cold-Survival Compiler from the Microscopic Export Kernel

## Status

**Exact within**

1. the Stage-234 microscopic export kernel
   \[
   K_{\rm exp}(s)=\Gamma_3 s^3+\Gamma_5 s^5,
   \]
2. the exact channel split
   \[
   \Gamma_3=\Gamma_3^{\rm vac}+\Gamma_3^{\rm lat},
   \qquad
   \Gamma_5=\Gamma_5^{\rm vac}+\Gamma_5^{\rm lat},
   \qquad
   \Gamma_n^{\rm vac},\Gamma_n^{\rm lat}\ge 0,
   \]
3. the exact exported-power identities from Stage 234,
4. an arbitrary finite event window \([0,T]\) with square-integrable \(\ddot V\) and \(\dddot V\),
5. the single-growth cold-event specialization
   \[
   V(t)=V_{\rm in}e^{s t},
   \]
6. and the Stage-234 event-safe equality on the cold edge
   \[
   \Gamma_3 s_c^3+\Gamma_5 s_c^5=\mu_\eta(s_0^2-s_c^2).
   \]

This stage is a **channel-resolved energy compiler**. It does not yet map the result to \(\lambda_{\rm ep}\omega_D\), \(k_{\rm eff}\), or \(T_{\max}\); that is the Stage-236 job.

---

## Purpose

Stage 234 replaced the Session-IV phenomenological envelope law
\[
\gamma_{\rm tot}\,\dot V
\]
by the first microscopic passive/export kernel carried by the moving-throat wall/BdG/localized-Maxwell/mixed stack,
\[
\Sigma_{\rm exp}(\omega)=-i\Gamma_3\omega^3-i\Gamma_5\omega^5+O(\omega^7).
\]
It also gave the exact channel-resolved power split
\[
\mathcal P_{\rm vac}
=
\Gamma_3^{\rm vac}\ddot V^{\,2}+\Gamma_5^{\rm vac}\dddot V^{\,2},
\qquad
\mathcal P_{\rm lat}
=
\Gamma_3^{\rm lat}\ddot V^{\,2}+\Gamma_5^{\rm lat}\dddot V^{\,2},
\]
but it deliberately stopped one step before the actual heat-partition problem.

So the next honest question is:

> for a chosen cold event, what fraction of the exported energy goes into the vacuum outlet, what fraction goes into the lattice outlet, and what total export is minimally required for the cold event to survive?

That is exactly what Stage 235 compiles.

The main results are:

1. the exact vacuum-vs-lattice split depends on one scalar event-shape quotient,
2. on a single-growth event that quotient collapses to \(s^2\),
3. the Stage-234 safe surface implies an exact **safe-edge exported-energy theorem**,
4. the Session-IV `3:1` heat split becomes a microscopic coefficient surface rather than a by-hand partition,
5. and the microscopic event-equivalent damping rates that Stage 236 must map into condensed-matter variables are now explicit.

---

## Provenance

This stage is the direct continuation Stage 234 itself asked for.

- Stage 233 converted the relaxed branch into a cold-event survival problem.
- Stage 234 replaced the phenomenological damping scalar by the microscopic cubic-plus-quintic export kernel and derived exact channel-resolved exported-power formulas.
- Session IV used a by-hand vacuum/lattice split to estimate exported energies and cold-event rescue.
- Session V then mapped that phenomenological lattice component into condensed-matter quantities.

Stage 235 is therefore the missing bridge:

- it removes the ad hoc heat partition,
- it computes the exact microscopic vacuum/lattice fractions on a chosen event,
- and it packages the result in the form that Stage 236 can translate into \(\lambda_{\rm ep}\omega_D\), trap stiffness, and thermal-survival variables.

---

## 0. Why this stage is needed

Before this step the stack could already say:

- which microscopic channels export energy from the active `V` leg,
- that the total safe-event condition is
  \[
  \Gamma_3 s_c^3+\Gamma_5 s_c^5
  \ge
  \mu_\eta(s_0^2-s_c^2),
  \]
- and that Session IV’s default split numerically placed most of the dissipated energy in the lattice channel.

But it still could **not** say:

- what the exact vacuum-vs-lattice fraction is on a general cold event,
- how that fraction depends on the relative cubic vs quintic loading of the two channels,
- how the phenomenological `3:1` split is written microscopically,
- or what exact exported energy is minimally required on the cold-survival edge.

Stage 235 closes those gaps.

---

## 1. Exact channel-resolved exported-energy ledger

Split the Stage-234 microscopic coefficients into vacuum and lattice pieces:
\[
\Gamma_3=\Gamma_3^{\rm vac}+\Gamma_3^{\rm lat},
\qquad
\Gamma_5=\Gamma_5^{\rm vac}+\Gamma_5^{\rm lat}.
\]
For an arbitrary finite event window \([0,T]\), define the exact shape integrals
\[
\mathcal I_2(T):=\int_0^T \ddot V^{\,2}\,dt,
\qquad
\mathcal I_3(T):=\int_0^T \dddot V^{\,2}\,dt.
\]
Then the integrated exported energies are exactly
\[
\boxed{
E_{\rm vac}(T)=\Gamma_3^{\rm vac}\,\mathcal I_2(T)+\Gamma_5^{\rm vac}\,\mathcal I_3(T),
}
\]
\[
\boxed{
E_{\rm lat}(T)=\Gamma_3^{\rm lat}\,\mathcal I_2(T)+\Gamma_5^{\rm lat}\,\mathcal I_3(T),
}
\]
so the total exported energy is
\[
\boxed{
E_{\rm exp}(T)=\Gamma_3\,\mathcal I_2(T)+\Gamma_5\,\mathcal I_3(T).
}
\]

This is the exact channel-resolved replacement for the Session-IV `3:1` heat split.

### 1.1 One-scalar shape quotient

Assume \(\mathcal I_2(T)>0\) and define the exact shape quotient
\[
\boxed{
\mathfrak r_V(T):=\frac{\mathcal I_3(T)}{\mathcal I_2(T)}.
}
\]
Then the partition fractions are
\[
\boxed{
 f_{\rm vac}(\mathfrak r_V)
 =
 \frac{\Gamma_3^{\rm vac}+\Gamma_5^{\rm vac}\mathfrak r_V}
 {\Gamma_3+\Gamma_5\mathfrak r_V},
}
\]
\[
\boxed{
 f_{\rm lat}(\mathfrak r_V)
 =
 \frac{\Gamma_3^{\rm lat}+\Gamma_5^{\rm lat}\mathfrak r_V}
 {\Gamma_3+\Gamma_5\mathfrak r_V},
 \qquad
 f_{\rm vac}+f_{\rm lat}=1.
}
\]
So the entire heat-partition problem collapses to one scalar event-shape quotient.

### 1.2 Exact speed-drift law of the partition

Treat \(\mathfrak r_V\) as the event-shape variable. Then
\[
\boxed{
\frac{d f_{\rm lat}}{d\mathfrak r_V}
=
\frac{\Gamma_5^{\rm lat}\Gamma_3^{\rm vac}-\Gamma_3^{\rm lat}\Gamma_5^{\rm vac}}
     {(\Gamma_3+\Gamma_5\mathfrak r_V)^2}.
}
\]
So the partition drifts upward toward the lattice side iff
\[
\Gamma_5^{\rm lat}\Gamma_3^{\rm vac}>
\Gamma_3^{\rm lat}\Gamma_5^{\rm vac},
\]
i.e. iff the lattice outlet is relatively more quintic-loaded than the vacuum outlet.

The endpoint limits are exact:
\[
\boxed{
 f_{\rm lat}(0)=\frac{\Gamma_3^{\rm lat}}{\Gamma_3},
 \qquad
 \lim_{\mathfrak r_V\to\infty} f_{\rm lat}(\mathfrak r_V)=\frac{\Gamma_5^{\rm lat}}{\Gamma_5}.
}
\]
So slow events are controlled by the cubic split, while rapidly varying events are controlled by the quintic split.

---

## 2. Single-growth cold-event specialization

For the reduced cold-event model, use the one-rate growth orbit
\[
\boxed{
V(t)=V_{\rm in}e^{st},
\qquad 0\le t\le T,
\qquad s>0.
}
\]
Then
\[
\dot V=sV,
\qquad
\ddot V=s^2V,
\qquad
\dddot V=s^3V.
\]
Therefore the exact shape integrals become
\[
\boxed{
\mathcal I_2(T)
=
\frac{V_{\rm in}^2 s^3}{2}\bigl(e^{2sT}-1\bigr),
}
\]
\[
\boxed{
\mathcal I_3(T)=s^2\mathcal I_2(T).
}
\]
So the Stage-235 shape quotient collapses exactly to
\[
\boxed{
\mathfrak r_V=s^2.
}
\]

### 2.1 Event-equivalent damping rates

Define the exact velocity-weighted event integral
\[
\mathcal I_1(T):=\int_0^T \dot V^{\,2}\,dt
=
\frac{V_{\rm in}^2 s}{2}\bigl(e^{2sT}-1\bigr).
\]
Then the exported energies can be written as
\[
\boxed{
E_{\rm vac}(T)=\gamma_{\rm vac}^{\rm eq}(s)\,\mathcal I_1(T),
\qquad
\gamma_{\rm vac}^{\rm eq}(s):=\Gamma_3^{\rm vac}s^2+\Gamma_5^{\rm vac}s^4,
}
\]
\[
\boxed{
E_{\rm lat}(T)=\gamma_{\rm lat}^{\rm eq}(s)\,\mathcal I_1(T),
\qquad
\gamma_{\rm lat}^{\rm eq}(s):=\Gamma_3^{\rm lat}s^2+\Gamma_5^{\rm lat}s^4,
}
\]
with total event-equivalent export rate
\[
\boxed{
\gamma_{\rm eff}^{\rm eq}(s)=\Gamma_3 s^2+\Gamma_5 s^4.
}
\]
So the microscopic replacement for the Session-IV phenomenological pair
\[
(\gamma_{\rm vac},\gamma_{\rm lattice})
\]
is the rate pair
\[
(\gamma_{\rm vac}^{\rm eq}(s),\gamma_{\rm lat}^{\rm eq}(s))
\]
evaluated on the actual event shape.

The partition fractions become
\[
\boxed{
 f_{\rm vac}(s)=\frac{\gamma_{\rm vac}^{\rm eq}(s)}{\gamma_{\rm eff}^{\rm eq}(s)},
 \qquad
 f_{\rm lat}(s)=\frac{\gamma_{\rm lat}^{\rm eq}(s)}{\gamma_{\rm eff}^{\rm eq}(s)}.
}
\]

---

## 3. Exact cold-survival compiler on the safe edge

Let
\[
 s_c:=\frac{1}{t_{\rm cross}},
 \qquad
 s_0:=\frac{1}{t_{\rm collapse,0}}=\sqrt{\kappa_V/\mu_\eta},
\]
and specialize to the **safe-edge** equality from Stage 234,
\[
\boxed{
\Gamma_3 s_c^3+\Gamma_5 s_c^5
=
\mu_\eta(s_0^2-s_c^2).
}
\]
Take the event window to be the characteristic crossing time,
\[
T=t_{\rm cross}=\frac{1}{s_c}.
\]
Then
\[
\mathcal I_1\Bigl(\frac{1}{s_c}\Bigr)
=
\frac{V_{\rm in}^2 s_c}{2}(e^2-1),
\]
so the safe-edge total exported energy is exactly
\[
\boxed{
E_{\rm exp,min}^{\rm safe}
=
\frac{V_{\rm in}^2}{2}(e^2-1)\,\mu_\eta(s_0^2-s_c^2).
}
\]
This is the stage’s first main theorem.

It says that once the event is pinned to the cold-survival edge, the **minimum total exported energy** is completely fixed by

- the initial active-leg amplitude \(V_{\rm in}\),
- the Stage-233/234 rate deficit \(s_0^2-s_c^2\),
- and the effective dressing inertia \(\mu_\eta\),

and is independent of how the microscopic vacuum/lattice split is distributed between the cubic and quintic channels.

### 3.1 Channel-resolved safe-edge energies

The safe-edge channel energies are therefore
\[
\boxed{
E_{\rm vac,min}^{\rm safe}=f_{\rm vac}(s_c)
\frac{V_{\rm in}^2}{2}(e^2-1)\,\mu_\eta(s_0^2-s_c^2),
}
\]
\[
\boxed{
E_{\rm lat,min}^{\rm safe}=f_{\rm lat}(s_c)
\frac{V_{\rm in}^2}{2}(e^2-1)\,\mu_\eta(s_0^2-s_c^2).
}
\]
So the Stage-234 cold-survival surface and the Stage-235 heat partition now combine into one exact compiler.

### 3.2 Safe-edge event-equivalent rate

At the same edge, the total event-equivalent export rate obeys
\[
\boxed{
\gamma_{\rm eff,safe}^{\rm eq}
:=
\Gamma_3 s_c^2+\Gamma_5 s_c^4
=
\mu_\eta\frac{s_0^2-s_c^2}{s_c}.
}
\]
Therefore
\[
\boxed{
\gamma_{\rm vac,safe}^{\rm eq}=f_{\rm vac}(s_c)\,\mu_\eta\frac{s_0^2-s_c^2}{s_c},
\qquad
\gamma_{\rm lat,safe}^{\rm eq}=f_{\rm lat}(s_c)\,\mu_\eta\frac{s_0^2-s_c^2}{s_c}.
}
\]
This is the exact microscopic object that later condensed-matter maps must use. It is **not** numerically the same object as Session-IV’s envelope parameter \(\gamma_{\rm lattice}\).

---

## 4. The Session-IV `3:1` split as a microscopic coefficient surface

The Session-IV default split was
\[
 f_{\rm vac}=\frac14,
 \qquad
 f_{\rm lat}=\frac34.
\]
On the Stage-235 microscopic compiler, evaluated on the cold-survival edge, this becomes
\[
\boxed{
 f_{\rm lat}(s_c)=\frac34
 \iff
 \Gamma_3^{\rm lat}+\Gamma_5^{\rm lat}s_c^2
 =
 3\Bigl(\Gamma_3^{\rm vac}+\Gamma_5^{\rm vac}s_c^2\Bigr).
}
\]
So the old by-hand partition is now a linear surface in microscopic coefficient space.

### 4.1 Speed-independent special family

If the vacuum/lattice split is the same in both channels, i.e.
\[
\Gamma_3^{\rm lat}=\phi\,\Gamma_3,
\qquad
\Gamma_5^{\rm lat}=\phi\,\Gamma_5,
\]
then
\[
\boxed{
 f_{\rm lat}(s)=\phi
 \qquad\text{for every event speed }s.
}
\]
So the Session-IV `3:1` rule is recovered exactly by the microscopic family
\[
\boxed{
\phi=\frac34.
}
\]
That is the simplest reduced continuation of the old phenomenological split.

### 4.2 Speed-dependent split

If the cubic and quintic channels do **not** share the same vacuum/lattice fraction, then the partition becomes event-speed dependent through the exact law in Section 1.2. In that case:

- slower events sample mostly the cubic split,
- faster events sample mostly the quintic split,
- and the Session-IV `3:1` fraction can hold at one chosen cold edge without holding elsewhere.

That is the second main theorem of the stage.

---

## 5. Session-IV benchmark specialization

Use the Session-IV cold-event data carried into Stage 234:
\[
 t_{\rm cross}\approx 1.82169718,
 \qquad
 s_c=1/t_{\rm cross}\approx 0.5489386551,
 \qquad
 s_c^2\approx 0.3013336471,
\]
\[
 s_0\approx \gamma_{\rm crit}\approx 6.94311167,
 \qquad
 \mu_\eta=1.
\]
Then the exact safe-edge event-equivalent export rate is
\[
\boxed{
\gamma_{\rm eff,safe}^{\rm eq}
=
\frac{s_0^2-s_c^2}{s_c}
\approx 87.26925235.
}
\]
On the speed-independent `3:1` microscopic family \((\phi=3/4)\), this gives
\[
\boxed{
\gamma_{\rm vac,safe}^{\rm eq}\approx 21.81731309,
\qquad
\gamma_{\rm lat,safe}^{\rm eq}\approx 65.45193926.
}
\]
These are the exact microscopic event-equivalent rates replacing the phenomenological Session-IV pair.

### 5.1 Safe-edge energy match to Session IV

Session IV reported the safe-edge dissipated energy
\[
E_{\rm diss,total}^{\rm sess}\approx 0.01033460.
\]
If one uses this only as a **benchmark calibration** of the initial active-leg amplitude, the safe-edge energy theorem gives
\[
V_{\rm in,match}
=
\sqrt{\frac{E_{\rm diss,total}^{\rm sess}}
{\tfrac12(e^2-1)\mu_\eta(s_0^2-s_c^2)}}
\approx 8.21771260\times 10^{-3}.
\]
Then the Stage-235 channel formulas reproduce the Session-IV split exactly:
\[
\boxed{
E_{\rm vac,min}^{\rm safe}\approx 0.00258365,
\qquad
E_{\rm lat,min}^{\rm safe}\approx 0.00775095,
}
\]
with
\[
E_{\rm exp,min}^{\rm safe}\approx 0.01033460.
\]
So the old benchmark numbers are absorbed by the new microscopic compiler without keeping the heat partition phenomenological.

This is a calibration consistency check, not a theorem that the realized branch must use the speed-independent `3:1` family.

---

## 6. What this stage achieves physically

Stage 235 changes the barrier thread in four concrete ways.

### 6.1 It removes the ad hoc vacuum/lattice heat split

The exported energies are now exact functions of the microscopic coefficients
\[
(\Gamma_3^{\rm vac},\Gamma_3^{\rm lat},\Gamma_5^{\rm vac},\Gamma_5^{\rm lat})
\]
and one scalar event-shape quotient.

### 6.2 It isolates the true event variable controlling the partition

For a general event the partition is controlled by
\[
\mathfrak r_V=\mathcal I_3/\mathcal I_2,
\]
while on the cold single-growth branch it reduces exactly to
\[
\mathfrak r_V=s^2.
\]
So the partition is sensitive not only to the outlet coefficients, but also to the dynamical roughness of the event.

### 6.3 It produces an exact safe-edge exported-energy theorem

At cold survival threshold, the minimum exported energy is fixed by
\[
\frac{V_{\rm in}^2}{2}(e^2-1)\mu_\eta(s_0^2-s_c^2),
\]
and the vacuum/lattice shares are obtained by a pure fraction multiplier.

### 6.4 It gives Stage 236 the right microscopic object to map

The quantity that later condensed-matter equations must map is not the old envelope scalar \(\gamma_{\rm lattice}\), but the event-equivalent microscopic rate
\[
\gamma_{\rm lat}^{\rm eq}(s_c)=\Gamma_3^{\rm lat}s_c^2+\Gamma_5^{\rm lat}s_c^4.
\]
That is the exact bridge needed for the materials stage.

---

## 7. What is still missing

Stage 235 is still not the full materials theorem of the completed PDE.
The remaining open objects are now sharply narrowed to:

1. the realized branch values of
   \[
   \Gamma_3^{\rm vac},\Gamma_3^{\rm lat},\Gamma_5^{\rm vac},\Gamma_5^{\rm lat},
   \]
2. the realized event shape quotient \(\mathfrak r_V\) on the cold branch,
3. the physical-unit conversion linking
   \[
   \gamma_{\rm lat}^{\rm eq}(s_c)
   \]
   to condensed-matter turnover variables,
4. and any branch-dependent calibration of the initial active-leg amplitude if one wants absolute per-event energies rather than fractions.

So the live bottleneck has moved from “invent a damping partition” to “compute the realized microscopic coefficient split and then map it physically.”

---

## 8. Immediate next step

Stage 236 is now sharply defined.

1. Keep the Stage-235 microscopic event-equivalent lattice rate
   \[
   \gamma_{\rm lat}^{\rm eq}(s_c)=\Gamma_3^{\rm lat}s_c^2+\Gamma_5^{\rm lat}s_c^4.
   \]
2. Map that quantity into the electron-phonon turnover condition replacing the old phenomenological \(\gamma_{\rm lattice}\).
3. Carry the same microscopic compiler into the harmonic-trap stiffness and Korringa-limited spin-survival maps.

That is the right Stage-236 condensed-matter companion.
