# Parent Throat Action — Transfer-Amplitude Interpretation

## Purpose

`step_25` and `step_26` showed that the `step_24` outgoing-amplitude coordinate
cannot be read naively as the PDE scalar \(N_Q\). On that interpretation the
Packet-A finish line is not improved, and the natural-source outgoing burden
stays huge.

That still leaves a different possibility:

> perhaps \(\lambda_{\rm out}\) is not \(N_Q\), but a genuinely independent
> transfer/outlet amplitude coordinate.

This step isolates that possibility exactly.

---

## 1. Uniform outgoing-packet scaling is algebraically coherent

From the exact grouped-prefactor formulas,

\[
P_0=\frac{N_0}{D_0},
\qquad
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]
\[
P_4=
\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^{\,2}N_0}{D_0^3}.
\]

So a uniform outgoing-packet scaling

\[
(N_0,N_2,N_4)\mapsto
\lambda_{\rm out}(N_0,N_2,N_4)
\]

forces

\[
(P_0,P_2,P_4)\mapsto
\lambda_{\rm out}(P_0,P_2,P_4).
\]

The grouped outgoing compiler then gives

\[
K_0=P_0,
\qquad
K_2=P_2+A P_0,
\qquad
K_4=P_4+A P_2+B P_0,
\qquad
\Gamma_5=G_5P_0,
\]

so the full compiled branch packet also scales uniformly:

\[
(K_0,K_2,K_4,\Gamma_5)\mapsto
\lambda_{\rm out}(K_0,K_2,K_4,\Gamma_5).
\]

At the same time the normalized shape ratios

\[
K_2/K_0,\qquad K_4/K_0
\]

are unchanged.

So `step_24` did not invent an algebraically inconsistent direction. It picked
an exact amplitude-type scaling direction for the outgoing packet.

---

## 2. Exact PDE-side amplitude channel

The PDE ledger’s exact transfer-shape theorem gives

\[
T_{\rm eff}^2=\frac{N_0}{K}.
\]

Therefore a uniform outgoing amplitude scaling at fixed \(K\) produces

\[
T_{\rm eff}^2\mapsto \lambda_{\rm out}T_{\rm eff}^2.
\]

This is important because it shows there is already an exact PDE-side scalar
amplitude channel that is separate from the outgoing-normalization scalar
\(N_Q\).

So after `step_26`, the right statement is not

> “`step_24` is impossible.”

The right statement is

> “`step_24` is impossible **as \(N_Q\)**, but it may still be possible as an
> independent transfer/outlet amplitude coordinate.”

---

## 3. Exact hybrid outlet family

The Stage-95 hybrid Robin/mixed branch gives the exact branch-B law

\[
\chi_B(\sigma,\gamma)=\frac{1-9\sigma\gamma}{1-\sigma}.
\]

On the canonical quintic branch \(\gamma=\frac19\), this collapses to

\[
\chi_B=1.
\]

At the same time, the exact Stage-95 scaled identity says the outgoing branch
scales by

\[
1-\sigma.
\]

So the minimal exact outlet family already contains an independent amplitude
coordinate:

\[
\lambda_{\rm out}=1-\sigma,
\qquad
\chi_Q=1.
\]

This is the clean escape hatch that survives `step_25` and `step_26`.

---

## 4. Quantitative burden

The script maps representative `step_24` points into the minimal Stage-95
branch-B amplitude parameter:

- best sampled point under \(\widehat m_0^{\rm req}\le 50\):
  \[
  \lambda_{\rm out}=20
  \quad\Rightarrow\quad
  \sigma=-19,
  \]
- best sampled point with \(Q_{\rm iso}\le 1\):
  \[
  \lambda_{\rm out}=2000
  \quad\Rightarrow\quad
  \sigma=-1999,
  \]
- best sampled point with \(Q_{\rm iso}\le 0.5\):
  \[
  \lambda_{\rm out}=2000
  \quad\Rightarrow\quad
  \sigma=-1999.
  \]

So the remaining burden is no longer algebraic admissibility. It is the size of
the required amplitude deformation in the smallest exact PDE family we know how
to write down.

---

## 5. Interpretation

This step changes the diagnosis again, but in a narrower way than `step_24`
did.

1. The large-`\(\lambda_{\rm out}\)` improvement is not admissible as the PDE
   scalar \(N_Q\).
2. But exact PDE algebra does admit a different kind of independent outgoing
   amplitude coordinate.
3. The minimal explicit family is the Stage-95 branch-B hybrid outlet family,
   where \(\chi_Q\) stays canonical and the amplitude scales like
   \(1-\sigma\).
4. The current useful reduced points still sit at large negative \(\sigma\),
   especially \(\sigma=-1999\) for the strongest frontier points.

So the next realization question is now very specific:

> can a realized branch generate a moderately sized transfer/outlet amplitude
> coordinate of this type, or are the currently useful reduced points only
> accessible in a strongly nonperturbative deformation regime?
