# Parent Throat Action — Electron Fast Falsification

## Purpose

The shortest useful falsification screen is not to solve one free equation for
\(c_s\). That remains underdetermined until a dimensional map is fixed.

The sharper test is:

1. pick the most obvious direct physical map,
2. plug in real electron data,
3. demand one independent electron observable.

This step does exactly that for the moderate branch-B patch.

---

## 1. Conditional map being tested

The screen now writes the dimensional bridge with an explicit port scale:

\[
a = r_e,
\qquad
S_{\rm port}\,\widehat m_0^{\,2}\lambda_{\rm out}P_0^{\rm base}
=\frac{54Gc_s^5}{5a^5c^5}.
\]

The most naive direct physical map is the special case \(S_{\rm port}=1\), with
no additional dimensional prefactor inserted between the reduced branch and the
physical target.

This is the simplest map you would try if you wanted to turn the moderate
branch-B patch directly into an electron analog.

The moderate patch data are imported from `step_29`'s exported corridor records
rather than duplicated inside this script. The values are

\[
P_0^{\rm base}=0.00023523241237876055.
\]

For the moderate corridor:

\[
\lambda_{\rm out}=20
\quad\Rightarrow\quad
P_0=0.004704648247575211,
\quad
\widehat m_0\approx 47.91244113628207,
\]

\[
\lambda_{\rm out}=50
\quad\Rightarrow\quad
P_0=0.011761620618938028,
\quad
\widehat m_0\approx 30.302488449879455.
\]

Both satisfy

\[
\widehat m_0^{\,2}\lambda_{\rm out}P_0^{\rm base}
=\widehat m_0^{\,2}P_0
=\frac{54}{5}=10.8
\]

in reduced units.

The script treats these `step_29` values as load-bearing data. As a mutation
guard, omitting the \(\lambda_{\rm out}=50\) factor changes the required
\(S_{\rm port}\) by a factor of exactly `50`.

The script also mutates the port scale itself. Replacing the direct
\(S_{\rm port}=1\) by \(2\) leaves a residual of about `-10.8`, and replacing
the Compton-calibrated \(S_{\rm port}\) by `1.01*S_port` leaves a residual of
about `-1.1054545016851435e+50`. These mutations make the dimensional port
factor load-bearing rather than a decorative label.

---

## 2. Direct-SI closure prediction

Using

- \(G=6.67430\times 10^{-11}\),
- \(c=299792458\),
- \(\hbar=1.054571817\times 10^{-34}\),
- \(m_e=9.1093837015\times 10^{-31}\),
- \(r_e=2.8179403262\times 10^{-15}\),

the direct-SI closure \(S_{\rm port}=1\) gives

\[
c_s=
\left(\frac{5(54/5)\,r_e^5 c^5}{54G}\right)^{1/5}
\approx 9.159491211330665\times 10^{-5}.
\]

So the analog pole scale is

\[
\Omega_Q=\frac{3c_s}{2a}
\approx 4.8756308603324455\times 10^{10}.
\]

The electron Compton angular frequency is

\[
\omega_C=\frac{m_ec^2}{\hbar}
\approx 7.76344071105011\times 10^{20}.
\]

Therefore

\[
\frac{\Omega_Q}{\omega_C}
\approx 6.280244857660478\times 10^{-11}.
\]

Equivalently,

\[
\frac{\hbar\Omega_Q}{m_ec^2}
\approx 6.280244857660478\times 10^{-11}.
\]

So under this direct-SI map, the analog pole scale is about

\[
1.59\times 10^{10}
\]

times too small to sit at the electron rest-energy scale.

---

## 3. Reverse calibration

If instead one forces the analog pole to the electron Compton scale,

\[
\Omega_Q=\omega_C,
\]

then

\[
c_s=\frac{2a\omega_C}{3}
\approx 1458460.8433153937.
\]

That in turn requires

\[
K_{\rm required}
=\frac{54Gc_s^5}{5a^5c^5}
\approx 1.1054545016851398\times 10^{52}.
\]

Compared to the reduced branch value \(54/5\), the required port scale is

\[
S_{\rm port}
=\frac{K_{\rm required}}{54/5}
\approx 1.023568983041796\times 10^{51}.
\]

So if you force the pole scale to match the electron, the dimensional port map
must insert about fifty-one orders of magnitude.

---

## 4. Interpretation

This is the fastest clean falsification result so far.

Under the naive direct-SI electron map:

1. either you accept the reduced normalization target, in which case the pole
   scale lands about ten orders too low;
2. or you force the pole scale to the electron Compton value, in which case the
   required port factor is about \(1.02\times10^{51}\).

So the naive direct electron identification is falsified immediately.

That is a **conditional** falsification:

it rejects the simplest direct map from the reduced branch-B patch to electron
data.

It does **not** yet reject the whole analog program, because the dimensional
port map and any additional physical normalization factors have not been fixed.
