
# Moving-Throat PDE — Stage 118: Outlet-Consistent Mouth Closure

## Goal

Use the concrete two-channel throat-core compensation law from Stages 97–99 as a
first explicit closure for the Family-1 mouth-layer gains.

The point is to see what the coupled mouth-layer solve predicts when the
mouth-bias channels inherit the same shell-versus-mixed weighting that preserved
the canonical outgoing quadrupole branch.

---

## 1. Nontrivial compensated outlet ratio

The exact compensated core outlet found in Stage 98 was
\[
\delta\Lambda_{\rm core}(z)
=
4\sigma_*-\frac{\sigma_*}{1-z^2/3-i z^5/9}.
\]

So the nontrivial compensated outlet weights the static shell lane and the mixed
pole lane in the ratio
\[
\boxed{4 : -1.}
\]

As a first **outlet-consistent mouth closure**, impose the same ratio on the
dimensionless mouth-layer gains:
\[
\boxed{
M_s = 4\Sigma_m,
\qquad
M_q = -\Sigma_m,
}
\]
with \(\Sigma_m>0\) a single residual mouth-gain amplitude.

This is not yet a theorem of the full PDE. It is the cleanest way to make the
mouth-layer solve and the outlet-core compensation law talk to each other.

---

## 2. Reduced one-parameter fixed-point law

With
\[
M_s = 4\Sigma_m,
\qquad
M_q = -\Sigma_m,
\]
the Family-1 mouth-layer equation becomes
\[
\boxed{
\Pi
=
\Sigma_m\left[4-\mathcal S_q(\Pi)\right].
}
\]

So the actual source-shape bias is now selected by a **single** dimensionless
mouth gain \(\Sigma_m\).

At the canonical point computed below,
\[
\mathcal S_q(\Pi_*)\approx 0.6581<1,
\]
so the right-hand side there lies strictly between \(3\Sigma_m\) and
\(4\Sigma_m\). In particular, the canonical outlet-consistent branch lives at a
moderate \(\Pi\) rather than an extreme mouth-localization limit.

---

## 3. Canonical Family-1 value

Requiring the canonical compensation point
\[
\Pi=\Pi_*
\approx 1.50882951349316
\]
gives the exact outlet-consistent gain
\[
\boxed{
\Sigma_m^*
=
\frac{\Pi_*}{4-\mathcal S_q(\Pi_*)}
\approx 0.451485277739090.
}
\]

Therefore the corresponding shell and mixed gains are
\[
\boxed{
M_s^* = 4\Sigma_m^*
\approx 1.80594111095636,
\qquad
M_q^* = -\Sigma_m^*
\approx -0.451485277739090.
}
\]

So the canonical mouth-bias point is realized by a **moderate** one-parameter
gain rather than a delicate large-cancellation limit.

---

## 4. Direct comparison with the Stage-114 pure-mechanical threshold

At Stage 114 the pure-mechanical threshold was
\[
\Pi_* \approx 1.5088
\]
with no resolved mixed-lane correction.

The outlet-consistent coupled solve refines that picture:

- the shell lane must supply
  \[
  M_s^*\approx 1.80594111095636;
  \]
- but the mixed lane contributes back with
  \[
  M_q^*\mathcal S_q(\Pi_*)
  \approx -0.297111597463199151;
  \]
- leaving the exact net bias
  \[
  \Pi_*.
  \]

So the mixed D/N lane is not a tiny perturbation. It supplies an \(O(0.3)\)
correction to the shell demand while keeping the total bias in the moderate
canonical range.

---

## Result

Under the first outlet-consistent mouth closure, the full coupled mouth-layer
selection law has collapsed to
\[
\boxed{
\Pi
=
\Sigma_m\left[4-\mathcal S_q(\Pi)\right],
}
\]
and the canonical Family-1 branch is selected at
\[
\boxed{
\Sigma_m^*\approx 0.451485277739090.
}
\]

So the remaining mouth-layer ambiguity is now one dimensionless gain amplitude,
not an arbitrary source profile and not an arbitrary parent slope combination.
