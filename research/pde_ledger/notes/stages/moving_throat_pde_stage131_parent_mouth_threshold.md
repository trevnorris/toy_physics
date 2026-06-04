
# Moving-Throat PDE — Stage 131: Parent Micro-Threshold for Canonical Mouth Compensation

## Goal

Translate the explicit Family-1 compensation point
\[
\Pi_* \approx 1.50882951349316
\]
into a direct parent-level threshold on the localized Maxwell + confinement
mouth data.

---

## 1. The parent bias parameter

From Stage 129,
\[
\Pi_m=\frac{V_1L}{\Theta_\sigma},
\qquad
V_1=
\left.\partial_z\delta V_{\rm conf}\right|_{\rm m}
-
q_*\left.\partial_zA_0\right|_{\rm m}.
\]

So the explicit canonical branch condition is simply
\[
\boxed{
\Pi_m=\Pi_*.
}
\]

Equivalently,
\[
\boxed{
V_1
=
\frac{\Pi_*\,\Theta_\sigma}{L}
\approx
1.50882951349316\,\frac{\Theta_\sigma}{L}.
}
\]

This is the first direct parent-level threshold for the mouth-source law.

---

## 2. Localized-Maxwell / mechanical split

Writing
\[
V_1=T_m-q_*A_0',
\qquad
T_m:=\left.\partial_z\delta V_{\rm conf}\right|_{\rm m},
\qquad
A_0':=\left.\partial_zA_0\right|_{\rm m},
\]
the canonical compensated branch requires
\[
\boxed{
T_m-q_*A_0'
=
1.50882951349316\,\frac{\Theta_\sigma}{L}.
}
\]

So the explicit outlet-core problem has collapsed again:

- \(T_m\) and \(A_0'\) may each vary,
- but only the **single linear combination**
  \[
  T_m-q_*A_0'
  \]
  matters for the mouth-source bias on this boundary-layer branch.

---

## 3. Linearized deviation law

Near the compensated point,
\[
\mathfrak g_{\Pi}
=
\mathfrak g_-^{F1}
+
\mathfrak g'(\Pi_*)(\Pi-\Pi_*)
+O((\Pi-\Pi_*)^2),
\]
with
\[
\mathfrak g'(\Pi_*)\approx 0.0714453558083195.
\]

Therefore
\[
\boxed{
\mathfrak g_{\Pi}-\mathfrak g_-^{F1}
\approx
0.0714453558083195\,(\Pi_m-\Pi_*).
}
\]

So the remaining source-shape uncertainty is now first-order equivalent to the
single parent bias mismatch
\[
\Pi_m-\Pi_*=
\frac{L}{\Theta_\sigma}\left(T_m-q_*A_0'\right)-\Pi_*.
\]

---

## 4. Meaning

The mouth-source problem is no longer:

- which sign branch,
- which family,
- or whether compensation is even reachable.

It is now just:

\[
\boxed{
\text{does the real mouth layer satisfy }
T_m-q_*A_0'
\approx
1.509\,\Theta_\sigma/L\ ?
}
\]

That is a concrete microscopic threshold, not a broad qualitative criterion.
