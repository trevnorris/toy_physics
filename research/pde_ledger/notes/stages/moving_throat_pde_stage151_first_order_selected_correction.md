# Moving-Throat PDE — Stage 151: First-Order Source Correction Selected by the Full Mouth Profile

## Goal

Project the exact full-profile mouth residual \(R_*(x)\) onto the previously derived
Family-1 rigidity formulas.

This produces the **actual** first-order non-exponential correction selected by the
full coupled GNLS + localized-Maxwell mouth layer.

---

## 1. First-order self-consistent source law

The exact self-consistent positive mouth source is
\[
\Sigma_{\rm full}(x)\propto e^{-\Phi_*(x)}
=
e^{-\Pi_* x - R_*(x)}.
\]

Expanding about the canonical exponential source
\[
\Sigma_*(x)=\frac{\Pi_*e^{-\Pi_*x}}{1-e^{-\Pi_*}}
\]
gives the normalized first-order correction
\[
\boxed{
\Sigma_{\rm act}(x)
=
\Sigma_*(x)\Big[1-\widetilde R_*(x)\Big]
+O(R_*^2),
\qquad
\widetilde R_*(x):=R_*(x)-\langle R_*\rangle_*,
}
\]
where
\[
\langle f\rangle_*:=\int_0^1 \Sigma_*(x)f(x)\,dx.
\]

So the actual selected deformation is not free:
it is exactly the centered residual of the full mouth potential.

---

## 2. Only two moment shifts matter

As in Stages 197–198, define
\[
c(x)=\cos\!\left(\frac{\pi x}{2}\right),
\qquad
K_q(x)=\frac{\cosh\!\left(\frac{\pi}{2}(1-x)\right)}{\cosh(\pi/2)}.
\]

Then the actual first-order shifts are
\[
\boxed{
\delta \mathfrak g_{\rm act}
=
\int_0^1 c(x)\,\delta\Sigma_*(x)\,dx
=
-\operatorname{Cov}_*(c,R_*),
}
\]
\[
\boxed{
\delta \mathcal S_{\rm act}
=
\int_0^1 K_q(x)\,\delta\Sigma_*(x)\,dx
=
-\operatorname{Cov}_*(K_q,R_*),
}
\]
with
\[
\operatorname{Cov}_*(f,h)
=
\langle fh\rangle_*-\langle f\rangle_*\langle h\rangle_*.
\]

So the full coupled mouth layer selects a unique pair of moment shifts; no new
branch ambiguity is introduced.

---

## 3. Projection onto the rigidity kernel

The canonical lower branch remains defined by
\[
\mathfrak g=\mathfrak g_*,
\]
so the electrochemical bias retunes by
\[
\boxed{
\delta\Pi_{\rm act}
=
-\frac{\delta\mathfrak g_{\rm act}}{\mathfrak g_*'}
=
\frac{\operatorname{Cov}_*(c,R_*)}{\mathfrak g_*'}.
}
\]

The normalized mouth traction shift is
\[
\boxed{
\delta \widehat T_{m,{\rm act}}
=
A_T\,\delta\mathfrak g_{\rm act}
+
B_T\,\delta\mathcal S_{\rm act}
=
- A_T\,\operatorname{Cov}_*(c,R_*)
- B_T\,\operatorname{Cov}_*(K_q,R_*).
}
\]

So the full-mouth selected correction is now completely explicit:
it is the projection of the exact residual \(R_*(x)\) against the same rigidity
kernel already derived in Stage 249.

---

## 4. Meaning

The mouth-side ambiguity is now reduced to one actual physical object:

\[
\boxed{
R_*(x)
}
\]

Once \(R_*(x)\) is known from the channel solve, the induced source correction,
bias retuning, and traction retuning follow automatically.

So the next question is no longer “which positive family?”
It is simply: what are the actual numerical covariances of \(R_*(x)\) on the
explicit Family-1 branch?