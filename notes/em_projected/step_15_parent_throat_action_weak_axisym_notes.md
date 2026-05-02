# Parent Throat Action — Weak-Axisymmetric `P2` Law and the Wall-Only Gate Result

## Purpose

This note asks what the promoted throat action does to the existing
weak-axisymmetric grouped-`P2` program.

The key question is not whether the promoted action changes the known weak
axisymmetric pattern. The key question is whether the promotion gives a
**parent-level origin** for that pattern, and whether wall anisotropy alone can
close the live `5PN`-style gates.

The answer is sharp:

1. the promoted action reproduces the same grouped signature
   `(20,21,22) ~ (1,1/2,-1)`,
2. therefore every first-order axisymmetric wall anisotropy still lies on the
   weak-axisymmetric line `b = 3 a`,
3. but a **pure wall-only** anisotropy cannot close the even gates nontrivially.

So promotion of `S_\Sigma` is necessary for parent completeness, but not by
itself sufficient to solve the remaining weak-axisymmetric bundle.

---

## 1. Weak axisymmetric perturbation of the promoted action

Take the isotropic quadratic data from the parent-complete action and perturb
it by a small axisymmetric `Y_{20}` component:

\[
\mu_\eta \to \mu_\eta + \varepsilon\,\delta\mu(w)Y_{20}(\Omega),
\]
\[
T_w \to T_w + \varepsilon\,\delta T_w(w)Y_{20}(\Omega),
\]
\[
T_\Omega \to T_\Omega + \varepsilon\,\delta T_\Omega(w)Y_{20}(\Omega),
\]
\[
K_\eta \to K_\eta + \varepsilon\,\delta K_\eta(w)Y_{20}(\Omega).
\]

The grouped `l=2` triple overlap then supplies the exact lane signature

\[
\lambda_{20}=1,
\qquad
\lambda_{21}=\frac12,
\qquad
\lambda_{22}=-1.
\]

The accompanying SymPy audit derives this ratio from the Wigner/Gaunt overlap
\(\int Y_{20}Y_{2m}Y_{2,-m}\,d\Omega\), together with the real-harmonic factor
\((-1)^m\) coming from \(Y_{2m}^*=(-1)^mY_{2,-m}\), rather than taking the
three lane weights as independent inputs.

So any wall-sector first-order shift has the form

\[
X_{20}=X^{(0)}+\varepsilon X^{(1)},
\qquad
X_{21}=X^{(0)}+\frac{\varepsilon}{2}X^{(1)},
\qquad
X_{22}=X^{(0)}-\varepsilon X^{(1)}.
\]

---

## 2. Exact grouped anomaly law

Using the grouped trace/anomaly variables

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\qquad
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}{2},
\]

a weak wall anisotropy obeys

\[
\bar x = x^{(0)},
\qquad
a_x = \frac{\varepsilon}{4}x^{(1)},
\qquad
b_x = \frac{3\varepsilon}{4}x^{(1)}.
\]

Hence

\[
\boxed{b_x = 3 a_x.}
\]

So the promoted throat action preserves the same weak-axisymmetric line already
isolated elsewhere in the grouped-`P2` program.

---

## 3. Wall contributions to the live weak-axisymmetric gates

Let the wall-only first-order shifts in the grouped `P2` lane be

\[
\delta M_A = \varepsilon\lambda_A M_1,
\qquad
\delta K_A = \varepsilon\lambda_A K_1^{\rm wall}.
\]

Then the wall contribution to the denominator slopes is simply

\[
D_{01}^{(A)} = \delta K_A,
\qquad
D_{21}^{(A)} = -\delta M_A,
\qquad
D_{41}^{(A)} = 0.
\]

So the two live even gates become

\[
K_{1,\rm wall} = D_{21}+\frac{D_{01}}{9} = -\delta M + \frac{\delta K}{9},
\]
\[
H_{{\rm even},\rm wall} = D_{41}-\frac{2}{3}D_{21}-\frac{D_{01}}{27} = \frac{2}{3}\delta M - \frac{\delta K}{27}.
\]

Solving

\[
K_{1,\rm wall}=0,
\qquad
H_{{\rm even},\rm wall}=0,
\]

gives only

\[
\boxed{\delta K=0,\qquad \delta M=0.}
\]

So a pure wall-only axisymmetric deformation cannot close the even gates
nontrivially.

This is best read as a linear obstruction inside the reduced wall-only gate
system. The updated script now obtains the wall-only formulas by specializing
the full weak-axisymmetric gate system to zero support/mixed slopes and then
checks that the resulting \(2\times 2\) coefficient matrix has determinant
\(1/27\), so only the trivial wall-only slope packet solves both at once. It
is not being claimed here as a deeper universal no-go theorem detached from the
reduced model.

That recovers, now from the parent-promoted wall action itself, the same
qualitative verdict that the later reduced weak-axisymmetric notes found from
the full grouped packet.

---

## 4. Consequence for the prefactor slope

If the outgoing-transfer sector stays isotropic at first order, then the
prefactor slope receives only the conservative denominator shift:

\[
\Xi_{\rm load}^{(A)} = -\frac{D_{01}^{(A)}}{D_0} = -\varepsilon\lambda_A\frac{K_1^{\rm wall}}{D_0}.
\]

So the promoted wall action again gives a grouped weak-axisymmetric defect on
exactly the same line

\[
\boxed{b = 3 a.}
\]

Likewise the static prefactor perturbation

\[
\delta P_0^{(A)} = -\frac{N_0 D_{01}^{(A)}}{D_0^2}
\]

has the same grouped signature.

So the promotion of `S_\Sigma` gives a parent-level origin for the grouped
transport pattern, but it does **not** remove the need for Maxwell/mixed and
support participation in the full prefactor and even-gate closure.

---

## 5. Physical reading

The promoted action clarifies the weak-axisymmetric story in three ways.

First, the grouped real `P2` transport line is not a mysterious algebraic
artifact. It is the direct consequence of how a `Y_{20}` deformation talks to
`l=2` wall modes.

Second, the wall sector now has a genuine parent-level meaning:
its slope data are not abstract fit parameters but first variations of the
parent throat action around the isotropic branch.

Third, the result is also a narrow obstruction:

> parent-level wall promotion alone does not solve the weak-axisymmetric `5PN`
> gate problem. The support and Maxwell/mixed sectors still have to participate.

That is exactly the kind of structurally useful result one wants from the
promotion.

---

## 6. Immediate next theorem gate

The next sharp move is now smaller than “solve the full nonlinear PDE.”

It is:

1. take the promoted throat action as the parent wall block,
2. combine its wall-side slope data with the already-derived support and
   Maxwell/mixed slope data,
3. test the full weak-axisymmetric packet
   `\{\Xi_1, K_1, H_{\rm even}\}`
   on an actual branch export,
4. and see whether the realized branch stays near the monomial-preserving /
   prefactor-compatible surface.

So the promotion moves the program forward in a very concrete way:
`K` and `M` become branch outputs of a parent action, but the branch still has
to earn the full packet. 
