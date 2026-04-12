# 5PN Stage 350–351 — Exact Family-1 branch location on the refreshed `\Lambda_{EM}` geometry

This pass stops compressing the support/source side and actually solves the exact operator-selected support/source fixed-point equation on the first explicit Family-1 branch, using the **refreshed**
\[
\Lambda_{EM}=\frac{\sqrt2\,\pi}{x_{01}}
\]
geometry rather than the old `37/20` shorthand.

## What was numerically located

Using
\[
\Lambda_\ell = \frac{L}{\ell}=20\Lambda_{EM}\approx 36.94973154240256,
\qquad
\eta=\Lambda_\ell,
\qquad
\kappa=\frac{1440\pi^2}{x_{01}^2}\approx 2457.5087899001137,
\]
the exact Robin/support constants are

\[
y\tan y=\eta
\quad\Longrightarrow\quad
y\approx 1.5294278190457656,
\]
\[
A_K=\frac{\kappa+\pi^2/4}{\kappa+y^2}\approx 1.0000521380385143,
\]
\[
\Delta_0\approx 1.7377393923469950\times 10^{-4},
\qquad
\Delta_\infty\approx 2.0172162594593645\times 10^{-2},
\]
\[
\zeta_{\max}=A_K\frac{\pi^2}{4}\approx 2.4675297457259358.
\]

Then, using the two explicit wall-depth extractions already carried in the notes,
\[
\Theta_w^{(\chi)}\approx 4.06863235008162\,\lambda_\mu^2,
\qquad
\Theta_w^{(J)}\approx 0.927552032539308\,\lambda_\mu^2,
\]
and setting the benchmark \(\lambda_\mu=1\), the exact wall/source figure of merit is
\[
\Xi = W_{\rm wall}=100\,\Theta_w\,\Lambda_\ell^2.
\]

This gives:

### χ-weighted extraction
\[
\Xi_\chi \approx 5.5548332017764099\times 10^5,
\]
and the exact fixed-point equation
\[
Pe=\Xi_\chi\,\Delta(Pe;\kappa,\eta)
\]
has the numerically located branch point
\[
Pe_*^{(\chi)} \approx 11155.7265863205869.
\]

At that root,
\[
\zeta_{\rm phys}^{(\chi)} \approx 2.4675296478814376,
\qquad
\rho_{\alpha,\max}^{(\chi)} = 1+\zeta_{\rm phys}^{(\chi)} \approx 3.4675296478814376.
\]

### J-weighted extraction
\[
\Xi_J \approx 1.2663707072528143\times 10^5,
\]
and the exact fixed-point branch point is
\[
Pe_*^{(J)} \approx 2504.9703142859238.
\]

At that root,
\[
\zeta_{\rm phys}^{(J)} \approx 2.4675278051675084,
\qquad
\rho_{\alpha,\max}^{(J)} = 1+\zeta_{\rm phys}^{(J)} \approx 3.4675278051675084.
\]

So both explicit operator-selected Family-1 branches sit extremely close to the exact Family-1 ceiling
\[
\zeta_{\max}\approx 2.4675297457259358.
\]

## What this means for the actual isotropic branch

The natural isotropic passive/outgoing grouped-`P2` branch still requires only
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13
\]
in the unblocked case, and it stays lowest-twin-safe throughout the admissible blocked interval.

Comparing against the numerically located Family-1 roots gives very large margins:

### χ branch
\[
\zeta_{\rm phys}^{(\chi)}-\zeta_{\rm req} \approx 2.1341963145481043,
\]
\[
\rho_{\alpha,\max}^{(\chi)}-\frac43 \approx 2.1341963145481043.
\]

### J branch
\[
\zeta_{\rm phys}^{(J)}-\zeta_{\rm req} \approx 2.1341944718341751,
\]
\[
\rho_{\alpha,\max}^{(J)}-\frac43 \approx 2.1341944718341751.
\]

So the support/source side is no longer just “probably okay.” On the refreshed exact `\Lambda_{EM}` Family-1 branch it is numerically located and **deep inside** the lowest-symmetric-twin-safe region.

## Finish-line status after the numerical location

This pass closes the strongest remaining **numerically accessible** part of the problem from the current stack:

- the explicit support/source branch has now been located numerically;
- it is strongly non-bottlenecked;
- on the canonical compact passive/outgoing branch one still has
  \[
  \chi_Q=1,\qquad N_Q=1
  \]
  exactly.

What is **still not numerically present in the files** is a completed PDE-selected point for the orbit-lock packet:
\[
d\ln R_{\rm tr}=0,\qquad d\ln R_{\rm target}=0,\qquad d\ln\epsilon_\eta=0.
\]

So the honest finish-line verdict is now:

1. **support/source** — numerically located and safe;
2. **outgoing normalization** — exactly fixed on the canonical passive/outgoing branch;
3. **orbit lock** — still only given as exact drift surfaces, not as a numerically selected PDE point.

That is the sharpest landing available from the current notes and scripts without a stronger direct solve of the completed moving-throat operator.
