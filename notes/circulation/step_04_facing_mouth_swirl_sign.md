# Step 4 — Global current sign versus local swirl labels for two facing mouths

## Purpose

Step 3 derived the fixed-current magnetic rule in a **global** coordinate convention:

\[
\text{same global current sense attracts.}
\]

This step translates that rule into the **local swirl labels** of two finite throat mouths.

This is the step that resolves the apparent tension between “parallel currents attract” and the physical intuition that two facing swirls may need to look opposite to attract.

## Orientation definitions

Let \(\hat d\) point from mouth 1 to mouth 2. Write the effective global current sign as
\[
I_A=I_0\,\sigma_A\,n_A,
\qquad I_0>0,
\qquad \sigma_A=\pm1.
\]
Here \(n_A\) is the local fluxoid/swirl label and \(\sigma_A\) tells us how positive local swirl maps into a common global loop direction.

In the audit, positive local swirl is represented by the right-handed loop basis around the local normal. For a loop of radius \(R\),
\[
\mathbf r(\phi)=R(\cos\phi\,\hat e+\sin\phi\,\hat f),
\qquad
\hat e\times\hat f=\hat n.
\]
The oriented current-axis vector is
\[
\frac12\int_0^{2\pi}\mathbf r\times \frac{d\mathbf r}{d\phi}\,d\phi
=\pi R^2\hat n.
\]
The local swirl field used in the audit is tangential,
\[
\mathbf v_\phi(\phi)=\frac{\Gamma_0}{2\pi R}\,\hat t(\phi),
\qquad
\hat t=\frac{1}{R}\frac{d\mathbf r}{d\phi},
\]
so the circulation is
\[
\oint\mathbf v_\phi\cdot d\boldsymbol\ell=\Gamma_0.
\]
With the global axis \(\hat d\), the coaxial orientation sign is
\[
\sigma_A=\hat n_A\cdot\hat d.
\]

From Step 3, the leading force is
\[
F_d=-\frac{K}{d^4}\,\sigma_1\sigma_2 n_1n_2,
\qquad K>0,
\]
with \(F_d<0\) meaning attraction.

Therefore the attraction condition is
\[
\boxed{\sigma_1\sigma_2 n_1n_2>0.}
\]

## Case 1: parallel local normals

If the two local positive-swirl conventions map to the same global axis,
\[
\sigma_1\sigma_2=+1,
\]
then attraction requires
\[
n_1n_2>0.
\]
So same local swirl attracts.

## Case 2: facing mouths

For two mouths facing each other, the local outward/face normals are opposite. The oriented loop integrals then give
\[
\sigma_1\sigma_2=-1.
\]
The attraction condition becomes
\[
-n_1n_2>0,
\]
or
\[
\boxed{n_1n_2<0.}
\]

So in the **facing-mouth local convention**:

\[
\boxed{\text{opposite local swirl labels attract, same local swirl labels repel.}}
\]

This is exactly the magnetism-compatible interpretation: the attractive state still has the same **global current/dipole alignment**, but because the mouths face one another, that alignment appears as opposite local circulation labels.

## Important caveat

This result is conditional on the fixed-current Maxwell/current-like closure. It does not follow from fluxoid quantization alone. Step 6 states the corresponding mixed-sector plumbing condition.

## SymPy audit

The script computes the tangential swirl circulation and oriented area/current-axis vectors for local normals \(+\hat z\) and \(-\hat z\), then verifies the sign table for parallel-normal and facing-mouth conventions.

Run:

```bash
python step_04_facing_mouth_swirl_sign_sympy.py
```
