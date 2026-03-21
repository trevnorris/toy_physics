# 2PN Family-1 soft-wall strict boundary-layer no-go: preliminary result

## 1. What was tested

This step asks whether the passed static support-channel data can already be reproduced by the **strictest** steep-wall interpretation of the Family-1 mouth geometry:

- one isotropic normal penetration moment \(A\),
- one isotropic tangential-gradient moment \(B\),
- and the full geometric pullback of a flared mouth sphere
  \[
  \sigma(\mu)=1-q\mu^2+r\mu^4.
  \]

The resulting quadratic-flare truncation of the pullback factors is
\[
J(\mu) = -7\mu^8 r^2 + 6\mu^6 q r + 8\mu^6 r^2 - \mu^4 q^2 - 8\mu^4 q r + 2\mu^4 r + 2\mu^2 q^2 - 2\mu^2 q + 1,
\]
\[
F_{\theta}(\mu) = 8\mu^8 r^2 - 8\mu^6 q r - 8\mu^6 r^2 + 2\mu^4 q^2 + 8\mu^4 q r - 2\mu^2 q^2 + 1,
\]
\[
F_{\phi}(\mu) = -8\mu^8 r^2 + 8\mu^6 q r + 8\mu^6 r^2 - 2\mu^4 q^2 - 8\mu^4 q r + 2\mu^2 q^2 + 1.
\]

The strict surface action is then
\[
E^{\rm strict}_{\rm BL}
=
\frac12\int d\Omega
\left[
A J(\mu)\,\Psi^2
+
B\left(
F_{\theta}(\mu)(\partial_{\theta}\Psi)^2
+
F_{\phi}(\mu)\frac{(\partial_{\phi}\Psi)^2}{\sin^2\theta}
\right)
\right].
\]

---

## 2. Exact channel formulas

Projecting this strict model onto the \(\ell=0,1,2\) mouth harmonics gives
\[
K_{00}=\frac{7}{15}Aq^2-\frac{26}{35}Aqr-\frac23 Aq+\frac{23}{63}Ar^2+\frac25 Ar+A,
\]
\[
K_{1\perp}=\frac{11}{35}Aq^2-\frac25 Aqr-\frac25 Aq+\frac{13}{77}Ar^2+\frac{6}{35}Ar+A
+\frac{8}{35}Bq^2-\frac{32}{105}Bqr+\frac{32}{231}Br^2+2B,
\]
\[
K_{10}=\frac{27}{35}Aq^2-\frac{10}{7}Aqr-\frac65 Aq+\frac{25}{33}Ar^2+\frac67 Ar+A
-\frac{16}{35}Bq^2+\frac{64}{105}Bqr-\frac{64}{231}Br^2+2B,
\]
\[
K_{20}=\frac{13}{21}Aq^2-\frac{94}{77}Aqr-\frac{22}{21}Aq+\frac{6185}{9009}Ar^2+\frac67 Ar+A
-\frac{16}{7}Bq^2+\frac{320}{77}Bqr-\frac{320}{143}Br^2+6B,
\]
\[
K_{21}=\frac{13}{21}Aq^2-\frac{230}{231}Aqr-\frac67 Aq+\frac{205}{429}Ar^2+\frac{10}{21}Ar+A
+\frac{8}{21}Bq^2-\frac{96}{77}Bqr+\frac{800}{1001}Br^2+6B,
\]
\[
K_{22}=\frac{5}{21}Aq^2-\frac{58}{231}Aqr-\frac27 Aq+\frac{25}{273}Ar^2+\frac{2}{21}Ar+A
+\frac{16}{21}Bq^2-\frac{64}{77}Bqr+\frac{320}{1001}Br^2+6B.
\]

These are exact for the quadratic-flare truncation of the strict isotropic layer.

---

## 3. Best-fit test against the passed support data

The target support values are
\[
K_{00}=\frac4{45},\qquad
K_{1\perp}=\frac27,\qquad
K_{10}=\frac14,
\]
\[
K_{20}=\frac49,\qquad
K_{21}=\frac23,\qquad
K_{22}=\frac83.
\]

A multi-seed least-squares search over \((A,B,q,r)\) gives the robust best fit
\[
A\approx -0.0177094,
\qquad
B\approx 0.2505433,
\qquad
q\approx -3.5661065,
\qquad
r\approx -3.0166765,
\]
with
\[
\sum_i \Delta_i^2 \approx 0.5536733194,
\qquad
\max_i |\Delta_i| \approx 0.4390744081.
\]

The channel residuals at that best fit are
\[
\Delta_{00}\approx -0.1497425,
\qquad
\Delta_{1\perp}\approx +0.3824726,
\qquad
\Delta_{10}\approx -0.2656747,
\]
\[
\Delta_{20}\approx -0.1804115,
\qquad
\Delta_{21}\approx +0.4390744,
\qquad
\Delta_{22}\approx -0.2984083.
\]

So the residual obstruction is not small and does not disappear under reseeding.

---

## 4. What this means

This is a useful negative result.

The passed static support sector is **not** just the strict isotropic steep-wall layer pulled back through Family-1 flare geometry.

So:

- pure isotropic penetration + pure geometry pullback is too small a model,
- at least one extra local wall-stress / traction / profile degree of freedom is required,
- and the reduced variational wall-Hessian derived in the companion step is therefore a **genuine new structure**, not a disguised rewrite of the strict isotropic layer.

That sharply narrows the next PDE target:

> derive the extra local wall degree of freedom from the actual soft-wall / traction-balance physics, rather than trying to squeeze the full support sector out of isotropic penetration alone.
