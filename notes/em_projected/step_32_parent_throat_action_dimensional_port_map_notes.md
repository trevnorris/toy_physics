# Parent Throat Action — Dimensional Port Map

## Purpose

The next missing piece is not another reduced branch search. It is the exact
bridge between the SI-side targets and the reduced PDE variables actually used
by the runtime patch.

This step compiles that bridge explicitly.

---

## 1. Frozen carry-forward constants

From the 4D summaries, the runtime patch inherits:

\[
n=5,\qquad
\alpha_{\rm opt}=\frac{n-1}{2}=2,
\qquad
\kappa_{\rm add}=\frac12,
\]
\[
\alpha^2=\frac34,
\qquad
a_H=0,
\qquad
K_{\rm vec}=\frac{2}{\pi^2}.
\]

These are no longer free knobs.

---

## 2. Exact source-map normalization

The compact PDE note’s exact normalization surface is

\[
\widehat m_0^{\,2}\,S_{\rm port}\,P_0
=
\frac{54Gc_s^5}{5a^5c^5}.
\]

So the reduced prefactor target is

\[
P_0^{\rm target}
=
\frac{54Gc_s^5}{5S_{\rm port}a^5c^5\widehat m_0^{\,2}}.
\]

Equivalently,

\[
S_{\rm port}
=
\frac{54Gc_s^5}{5a^5c^5\widehat m_0^{\,2}P_0},
\]

and

\[
c_s
=
\left(
\frac{5P_0S_{\rm port}a^5c^5\widehat m_0^{\,2}}{54G}
\right)^{1/5}.
\]

So the exact dimensional bridge is carried by \(S_{\rm port}\), not by
\(\widehat m_0\) or \(\lambda_{\rm out}\).

---

## 3. What the hardcoded Branch-B patch actually fixes

The reduced Branch-B scripts enforce

\[
\widehat m_0^{\,2}P_0=\frac{54}{5}.
\]

Substituting that into the exact normalization law gives

\[
S_{\rm port}
=
\frac{Gc_s^5}{a^5c^5}.
\]

This is the key result.

On the hardcoded target sheet:

- \(\widehat m_0\) is dimensionless,
- \(P_0\) is dimensionless,
- \(\lambda_{\rm out}\) is dimensionless,
- \(S_{\rm port}\) carries the full unresolved SI bridge.

So the direct-SI electron screen from `step_31` was really the special source-map
convention

\[
S_{\rm port}=1.
\]

Under that convention,

\[
c_s=\left(\frac{a^5c^5}{G}\right)^{1/5}.
\]

That is exactly the closure `step_31` falsified.

---

## 4. Frequency and amplitude maps

The exact pole scale is

\[
\Omega_Q=\frac{3c_s}{2a},
\qquad
\omega_{\rm red}=\frac{\omega_{\rm phys}}{\Omega_Q}
=
\frac{2a\,\omega_{\rm phys}}{3c_s}.
\]

So any physical-spectrum comparison must first pass through this frequency map.

For the hardcoded Branch-B amplitude patch,

\[
\lambda_{\rm out}=\frac{P_0}{P_{0,\rm base}}
=
e^{\delta\ln T_{\rm eff}^2}
=
1-\sigma.
\]

So \(\lambda_{\rm out}\) is a pure dimensionless transfer-shape amplitude ratio.

---

## 5. Transfer-shape / prefactor bridge

From the reduced compiler,

\[
P_0=\frac{K_{\rm bl}}{D_0}T_{\rm eff}^2.
\]

From the exact selected-branch transfer shape,

\[
T_{\rm eff}^2
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}\frac{1-\epsilon_\eta}{R_{\rm target}}.
\]

Combining these gives

\[
\frac{P_0^{\rm target}}{T_{\rm eff}^2}
=
\frac{8}{\pi^2}
\frac{R_{\rm target}}{S_{\rm port}\widehat m_0^{\,2}(1-\epsilon_\eta)}.
\]

Equivalently,

\[
R_{\rm target}
=
\frac{27\pi^2Gc_s^5}{20a^5c^5}\frac{1-\epsilon_\eta}{T_{\rm eff}^2}.
\]

So the SI-side port map can be written either through \(P_0\) or through
\(T_{\rm eff}^2\). They are exactly linked.

---

## 6. Interpretation

This step identifies the real unit bridge.

1. The reduced Branch-B patch does **not** by itself fix a physical \(c_s\).
2. It fixes only the dimensionless combination \(\widehat m_0^{\,2}P_0=54/5\).
3. The entire unresolved SI mapping lives in \(S_{\rm port}\).
4. Therefore the fastest falsification strategy is to test explicit
   source-map conventions, not to pretend the reduced patch already comes with a
   unique dimensionalization.

That is why `step_31` is a valid fast falsifier of the naive direct-SI map, but
not yet of the whole hardcoded Branch-B patch.
