# Moving-Throat PDE — Stage 139: Renormalized Canonical Branch Under Full Core–Mouth Co-Evolution

## Goal

Determine the actual Family-1 traction required to restore the exact lower
compensated branch once the mouth profile and the core loading ratio are allowed
to co-evolve together.

Equivalently, solve the exact self-consistent condition
\[
\mathfrak g[\Sigma_{\Sigma_0}]=\mathfrak g_*,
\]
rather than imposing the old canonical traction by hand.

---

## 1. Unique renormalized canonical gain

Scanning the exact self-consistent Family-1 map shows that the fixed-point moment
\(\mathfrak g_{\rm fp}(\Sigma_0)\) increases monotonically over the physical
interval used in the solve. Therefore the compensation-restoration equation
\[
\mathfrak g_{\rm fp}(\Sigma_0)=\mathfrak g_*
\]
has a unique positive root on that analyzed interval.

That numerically located root is
\[
\boxed{
\Sigma_0^{\rm can}
\approx 4.651033550168876.
}
\]

Using the self-matched susceptibility closure,
\[
\Sigma_0=\frac{20}{9}\widehat T_m^2,
\]
this corresponds to the renormalized canonical traction
\[
\boxed{
\widehat T_{m,\rm can}
\approx 1.446708366456762.
}
\]

So the full co-evolving reduced closure on the analyzed positive branch window does
preserve the lower compensated branch, but only after a substantial upward
renormalization of the mouth traction.

---

## 2. The restored canonical fixed point

At that renormalized traction, the exact self-consistent fixed point satisfies
\[
\boxed{
\mathfrak g_{\rm can}=\mathfrak g_*
\approx 0.758035078944663,
\qquad
\mathcal R_{\rm can}=\frac14.
}
\]

Its selected mixed response and mouth bias are
\[
\boxed{
\mathcal S_{\rm can}
\approx 0.6703621156734617,
\qquad
\Pi_{\rm can}
\approx 3.8715643774790087.
}
\]

So the co-evolving compensated Family-1 branch is still perfectly regular and still
finite — it just sits at a higher traction/bias point than the earlier
fixed-core correction suggested.

---

## 3. Comparison with earlier mouth-only corrections

The earlier fixed-core mouth analysis gave approximately
\[
(\Pi_{\rm corr},\widehat T_{m,\rm corr})
\approx (2.4159,1.1731),
\]
with a one-step nonlinear check closer to
\[
(\Pi_1,\widehat T_{m,1})
\approx (2.5391,1.2104).
\]

The full co-evolving canonical solve instead gives
\[
\boxed{
(\Pi_{\rm can},\widehat T_{m,\rm can})
\approx
(3.8716,1.4467).
}
\]

Relative to the original canonical point
\[
(\Pi_*,\widehat T_{m,*})
\approx (1.5088,0.9015),
\]
the exact co-evolving compensation costs

\[
\boxed{
\frac{\Sigma_0^{\rm can}}{\Sigma_0^*}-1
\approx 1.5754070949223031,
}
\]

\[
\boxed{
\frac{\Pi_{\rm can}}{\Pi_*}-1
\approx 1.5659389234213572,
}
\]

\[
\boxed{
\frac{\widehat T_{m,\rm can}}{\widehat T_{m,*}}-1
\approx 0.6048074946616844.
}
\]

So exact preservation of the canonical outgoing quadrupole branch under the analyzed
full core–mouth co-evolving closure requires roughly a

- \(60.48\%\) traction increase,
- and a \(156.59\%\) mouth-bias increase,

relative to the original lower compensated point.

---

## 4. Meaning

This is the sharpest Family-1 mouth/core result so far.

Inside the explicit Family-1 co-evolving closure on the analyzed positive branch window:

1. the lower compensated branch **survives** full co-evolution,
2. it remains the unique physically admissible compensated branch,
3. but the old canonical point is not the final self-consistent point,
4. and exact compensation is recovered only on a **renormalized** finite-bias,
   finite-traction branch.

So the mouth-side problem is no longer branch selection and no longer a broad
profile ambiguity. It is now a quantitative renormalization problem for the
unique regular canonical branch.
