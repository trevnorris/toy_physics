# Step 3 — Fixed-current Maxwell/current-loop closure

## Purpose

This step adds one explicit closure: each finite throat mouth is approximated at large separation by a small circular current loop. This is the magnetism-like branch.

The result of this step is not a universal theorem of the fluxoid sector. It is the leading force law under a fixed-current localized Maxwell/current closure.

## Geometry

Take two coaxial circular loops with radii \(R_1,R_2\), separation \(d\), and currents \(I_1,I_2\), with \(d\gg R_A\).

The Neumann mutual inductance is
\[
M(d)=\frac{\mu_0}{4\pi}\oint\!\oint
\frac{d\boldsymbol\ell_1\cdot d\boldsymbol\ell_2}{|\mathbf r_1-\mathbf r_2|}.
\]

Parameterizing the two loops gives
\[
\mathbf r_1(t_1)=(R_1\cos t_1,R_1\sin t_1,0),
\qquad
\mathbf r_2(t_2)=(R_2\cos t_2,R_2\sin t_2,d).
\]
Then
\[
d\boldsymbol\ell_1\cdot d\boldsymbol\ell_2
=R_1R_2\cos(t_1-t_2)\,dt_1dt_2,
\]
and
\[
|\mathbf r_1-\mathbf r_2|^2
=d^2+R_1^2+R_2^2-2R_1R_2\cos(t_1-t_2).
\]
Translation symmetry reduces the double integral to the one-variable integral
\[
M(d)
=\frac{\mu_0 R_1R_2}{2}
\int_0^{2\pi}\frac{\cos u\,du}
{\sqrt{d^2+R_1^2+R_2^2-2R_1R_2\cos u}}.
\]

## Far-field finite-mouth expansion

Expanding for \(d\gg R_A\),
\[
\boxed{
M(d)
=\frac{\mu_0\pi R_1^2R_2^2}{2d^3}
\left[1-\frac{3}{2}\frac{R_1^2+R_2^2}{d^2}
+\frac{15}{8}\frac{R_1^4+3R_1^2R_2^2+R_2^4}{d^4}
+O\!\left(\frac{R^6}{d^6}\right)
\right].
}
\]

The leading term is exactly the magnetic dipole-dipole result for moments
\[
m_A=I_A\pi R_A^2.
\]

## Fixed-current mechanical potential

At fixed currents, the mechanical interaction potential is
\[
U_I(d)=-I_1I_2M(d).
\]
Therefore
\[
F_d=-\frac{\partial U_I}{\partial d}=I_1I_2M'(d),
\]
where \(F_d<0\) means the separation decreases, i.e. attraction.

Using the expansion above,
\[
\boxed{
F_d
= -\frac{3\mu_0\pi I_1I_2R_1^2R_2^2}{2d^4}
+\frac{15\mu_0\pi I_1I_2R_1^2R_2^2(R_1^2+R_2^2)}{4d^6}
+O\!\left(\frac{R^8}{d^8}\right).
}
\]
Equivalently,
\[
F_d
= -\frac{3\mu_0\pi I_1I_2R_1^2R_2^2}{2d^4}
\left[1-\frac{5}{2}\frac{R_1^2+R_2^2}{d^2}+\cdots\right].
\]
The script also prints the next displayed asymptotic term,
\[
\frac{35}{8}\frac{R_1^4+3R_1^2R_2^2+R_2^4}{d^4},
\]
inside the force ratio.

## Sign rule in this closure

For \(d\gg R_A\), the leading sign is fixed:

\[
I_1I_2>0 \quad\Longrightarrow\quad F_d<0\quad\text{attraction},
\]
\[
I_1I_2<0 \quad\Longrightarrow\quad F_d>0\quad\text{repulsion}.
\]

This is the ordinary parallel-current attraction / antiparallel-current repulsion law.

## Step verdict

Under a fixed-current Maxwell/current-loop closure,
\[
\boxed{\text{same global current sense attracts, opposite global current sense repels.}}
\]

The next step translates this global-current rule into the local swirl labels of two facing throat mouths.

## SymPy audit

The script derives the reduced Neumann integral from the double-line parameterization, then verifies the asymptotic mutual-inductance expansion, the fixed-current potential, and the finite-size force correction.

Run:

```bash
python step_03_current_loop_closure_mutual_inductance_sympy.py
```
