# Moving-Throat PDE — Stage 149: Exact Bundle Inversion of the Last Four Irreducible Branch Drifts

## Purpose

Stage 148 reduced the explicit lower compensated Family-1 branch to four irreducible microscopic drifts,
\[
\delta\ln \mathcal Z_q,
\qquad
\delta\ln \rho_w,
\qquad
\delta\ln c_s,
\qquad
\delta\ln a.
\]

The right next move is **not** to guess them one by one.  The grouped wall/BdG/Maxwell/mixed bundle already exposes four natural observables that determine them algebraically:
\[
\Theta_w,
\qquad
K_s,
\qquad
K_q,
\qquad
P_0.
\]

Here:

- \(\Theta_w\) is the explicit Family-1 wall-depth datum,
- \(K_s\) is the shell stiffness,
- \(K_q\) is the mixed side-channel stiffness,
- \(P_0\) is the isotropic static outgoing prefactor,
  \[
  P_0=\frac{N_0}{D_0}.
  \]

So the open problem is now narrower than “solve for four unrelated drifts.”
It is:

> derive the four bundle observables \((\Theta_w,K_s,K_q,P_0)\) on the actual branch.

Once those are known, the last four microscopic drifts follow exactly.

---

## 1. The four exact branch laws to invert

On the explicit Family-1 wall branch, the wall-depth datum is
\[
\Theta_w = 25\,\lambda_\mu^2\rho_w^2.
\]
So at fixed \(\lambda_\mu\),
\[
\boxed{
\delta\ln \Theta_w = 2\,\delta\ln \rho_w.
}
\]

On the healing-locked shell branch,
\[
K_s \propto a^2\rho_w,
\]
so
\[
\boxed{
\delta\ln K_s = 2\,\delta\ln a + \delta\ln \rho_w.
}
\]

On the exact lower compensated D/N branch,
\[
K_q \propto \mathcal Z_q\,\frac{c_s^2}{L_W^2},
\qquad
\delta\ln L_W = \delta\ln a,
\]
so
\[
\boxed{
\delta\ln K_q = \delta\ln \mathcal Z_q + 2\,\delta\ln c_s - 2\,\delta\ln a.
}
\]

And on the isotropic outgoing quadrupole branch,
\[
P_0 = \frac{54G}{5c^5}\,\frac{c_s^5}{a^5},
\]
so
\[
\boxed{
\delta\ln P_0 = 5\bigl(\delta\ln c_s - \delta\ln a\bigr).
}
\]

These four equations are the complete inversion system.

---

## 2. Exact solution for the four remaining microscopic drifts

Solving the system above gives
\[
\boxed{
\delta\ln \rho_w
=
\frac12\,\delta\ln \Theta_w,
}
\]
\[
\boxed{
\delta\ln a
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln \Theta_w,
}
\]
\[
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln \Theta_w
+
\frac15\,\delta\ln P_0,
}
\]
\[
\boxed{
\delta\ln \mathcal Z_q
=
\delta\ln K_q
-
\frac25\,\delta\ln P_0.
}
\]

So the last four branch drifts are no longer open in any diffuse sense.
They are exact algebraic images of \((\Theta_w,K_s,K_q,P_0)\).

---

## 3. Equivalent full-bundle form using \(P_0=N_0/D_0\)

Because the grouped-bundle Stage-5/6 normalization uses
\[
P_0=\frac{N_0}{D_0},
\qquad
\delta\ln P_0 = \delta\ln N_0 - \delta\ln D_0,
\]
we may rewrite the last two drifts directly in full-bundle language:
\[
\boxed{
\delta\ln c_s
=
\frac12\,\delta\ln K_s
-
\frac14\,\delta\ln \Theta_w
+
\frac15\bigl(\delta\ln N_0-\delta\ln D_0\bigr),
}
\]
\[
\boxed{
\delta\ln \mathcal Z_q
=
\delta\ln K_q
-
\frac25\bigl(\delta\ln N_0-\delta\ln D_0\bigr).
}
\]

So the true grouped wall/BdG/Maxwell normalization data appear only through the isotropic static ratio \(N_0/D_0\), exactly as the 2.5PN normalization package suggested.

---

## 4. Frozen-wall corollary on the explicit Family-1 branch

Stage 60 already extracted the natural explicit wall datum
\[
\Theta_w^{(\chi)} \approx 4.06863235008162\,\lambda_\mu^2.
\]
So if the explicit Family-1 wall profile and \(\lambda_\mu\) are held fixed at first order, then
\[
\delta\ln \Theta_w = 0.
\]
In that restricted branch-preserving case the inversion simplifies to
\[
\boxed{
\delta\ln \rho_w = 0,
\qquad
\delta\ln a = \frac12\,\delta\ln K_s,
}
\]
\[
\boxed{
\delta\ln c_s = \frac12\,\delta\ln K_s + \frac15\,\delta\ln P_0,
\qquad
\delta\ln \mathcal Z_q = \delta\ln K_q - \frac25\,\delta\ln P_0.
}
\]

Numerically, the corresponding wall density on the explicit branch is
\[
\rho_w^{(\chi)}
=
\sqrt{\Theta_w^{(\chi)}/25}\,\lambda_\mu^{-1}
\approx
0.403417022451042\,\lambda_\mu^{-1}.
\]

So even this explicit wall branch already removes one of the four would-be free drifts unless the wall profile itself is allowed to deform.

---

## 5. What Stage 149 changes

Before this step, the remaining unresolved statement was:

> compute \((\mathcal Z_q,\rho_w,c_s,a)\) from the grouped real wall/BdG/Maxwell bundle.

After this step, the exact theorem gate is smaller:

> compute \((\Theta_w,K_s,K_q,P_0)\), or equivalently \((\Theta_w,K_s,K_q,N_0/D_0)\), on the actual isotropic branch.

That is enough to recover all four microscopic drifts exactly.

So the remaining PDE-facing bottleneck is no longer the direct computation of four separate microscopic quantities.  It is the derivation of:

1. the explicit wall-depth drift \(\delta\ln\Theta_w\),
2. the shell stiffness drift \(\delta\ln K_s\),
3. the mixed-channel stiffness drift \(\delta\ln K_q\),
4. and the isotropic outgoing normalization drift \(\delta\ln P_0\).

Once those are known, the Stage-148 lower-branch transport is fully closed at first order.
