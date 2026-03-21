# 2PN inner-throat freeze + execution plan

## Frozen baseline for the next derivation stage

1. **Conservative, low-frequency, local closure first.**
   Work in the conservative near-zone limit and treat the throat response through the low-\(\omega\) local expansion of the mouth operator before touching any nonlocal or dissipative channel.

2. **Keep the existing 2PN bookkeeping split.**
   Preserve the harness separation:
   - self / optical sector,
   - static pair-counted sector,
   - response correction,
   - wake / comparable-mass cross sector.

3. **Freeze the geometry model at the reduced level.**
   Use
   - Family-1 confinement,
   - W1 walls-as-potential ontology,
   - reduced geometry DOFs \((a,L)\),
   - and for the first response prototype keep \(L=\Lambda a\).

4. **Freeze the throat microphysics to the already-calibrated \(n=5\) vacuum.**
   Use the same stiff barotrope as the 1PN program:
   \[
   P(\rho)=K\rho^5,
   \qquad
   c_s^2 \propto \rho^4,
   \qquad
   \rho/\rho_0 = (1-4u)^{1/4},\quad u\equiv U/c^2.
   \]

5. **Freeze the port convention for the one-body response prototype.**
   Use the monopole mouth channel only, with
   - effort variable = brane enthalpy perturbation \(u=\delta h\),
   - flux variable = normal brane mass flux,
   - sign convention chosen so passive response is positive-definite.

6. **Freeze the current static target.**
   Keep \(\lambda_\rho=1/2\) from the isotropic test-mass static match. The remaining open one-body problem is then purely dynamic.

## What the SymPy prototype established

The current self+static test-mass candidate differs from the exact isotropic Schwarzschild target only by
\[
\Delta L = 4 u^2 \nu,
\qquad
u\equiv U/c^2,\quad \nu\equiv v^2/c^2.
\]

A single quadratic correction to the Bernoulli denominator,
\[
D(u)=1-4u+\chi u^2,
\]
changes only the \(u^2\nu\) slot at 2PN. The coefficient becomes
\[
6-\chi/2.
\]
Matching the isotropic target fixes
\[
\chi=8.
\]

So the exact isotropic one-body target through 2PN can be written as
\[
L_{\rm iso}^{\rm 2PN}
=
-(1-u)
\sqrt{1-\frac{\nu}{1-4u+8u^2}}
-\frac12 u^2 + \frac14 u^3
+ O(\text{3PN}).
\]

This is the cleanest current target for the missing response sector.

## What the second-order 1DOF throat closure already gives

Using the normalized \(11:2:5\) operating-point partition,
\[
F(a,\rho)=11\frac{\rho^2}{a}+2\frac{1}{\rho a^2}+5\rho^5 a^3,
\]
the equilibrium series around \((\rho,a)=(1,1)\) is
\[
a(\rho)=1-\frac{57}{64}\,\epsilon + \frac{123717}{131072}\,\epsilon^2 + O(\epsilon^3),
\qquad
\epsilon\equiv \rho/\rho_0-1.
\]
This reproduces the known breathing slope.

After dividing out the baseline density scaling \(\kappa_\rho=1\), the pure internal response factor is
\[
R_{\rm PV}(\epsilon)=1+\frac32\epsilon+\frac{151}{256}\epsilon^2+O(\epsilon^3).
\]
Composing with exact \(n=5\) Bernoulli density gives
\[
a(u)=1+\frac{57}{64}u+\frac{298821}{131072}u^2+O(u^3),
\]
\[
R_{\rm PV}(u)=1-\frac32 u-\frac{425}{256}u^2+O(u^3).
\]

## Immediate next derivation target

The missing one-body 2PN problem is now narrow:

> derive the quadratic mouth / port / DtN correction that upgrades the raw Bernoulli denominator
> \(1-4u\) to the response-corrected denominator \(1-4u+8u^2\).

A minimal way to phrase that in the reduced throat model is
\[
D_{\rm eff}(u)=(1-4u)\big[1+\mu\,\delta a(u)^2\big]+O(u^3),
\qquad
\delta a(u)=a(u)-1.
\]
Then the required coefficient is
\[
\mu = \frac{8}{(57/64)^2} = \frac{32768}{3249}\approx 10.0856.
\]
This is the simplest concrete reduced-model target to send into the next Mathematica script.

## Script roadmap from here

1. **Refactor the prototype into a Mathematica-ready referee script.**
   Goal: verify the \(\chi=8\) denominator result and the second-order throat closure algebra symbolically in WL.

2. **Add a symbolic mouth-operator reduction block.**
   Parameterize the low-\(\omega\) monopole operator as
   \[
   Z_{00}(\omega;u)=Z_2(u)\,\omega^2 + Z_4(u)\,\omega^4+\cdots,
   \]
   and derive the condition on \(Z_2(u)\) that reproduces \(\chi=8\).

3. **Try two concrete closure routes for \(Z_2(u)\).**
   - pure breathing / geometry normalization,
   - breathing + port normalization / projection factor.

4. **Only after the one-body slot closes, return to the comparable-mass 2PN wake sector.**

