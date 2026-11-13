# S17 source-fidelity review of the first record

The reviewed record is preserved at `528d4e64`, SHA-256 `1e4802b348ab34fb9baf1fafdffecc1babbf9900c47d131f1909c67b821e8c1a`. The local reviewer checked the seven manuscripts and ledger sources at `438d270c`. Grok received the same record text with file and web tools disabled. Neither reviewer received the other's report. [Local report and retrieved passages](local.md) · [Grok report](grok.md)

The local verdict was “revisions required for source fidelity.” Grok opened with “Scoped clear (record text only).” Its report also listed five substantive wording issues and requested local verification of the sources. The original record required correction.

| Finding | Verified decision |
|---|---|
| Local 1, 2.5PN overlap | Carry the theorem's natural-branch condition `m̂₀ = 1 + O(a²/r²)`, with `m̂₀ ≠ 0`, beside the small-body assumptions. Absolute canonical normalization is a separate test. |
| Local 2, 3PN static coefficient | The imported reference fixes the **total** static coefficient. Subtract the contribution already fixed by the grouped quadrupole fit to identify the residual assigned to geometry. |
| Local 3, 4PN conservative response | Add the three grouped components on the physical branch and the required isotropy check `a₂ = b₂ = 0` before the amplitude comparison. The explicit toy tail has a separate transport condition. |
| Local 4, S15 item 4 | Restrict the deferred all-separation classification to the closed-form geometric departures on the power-law profile family. The physical response to changing the separation remains open. |
| Grok 1, use of GR targets | In the results for 1PN, 2PN, 3PN and local 4PN identify which reference data went into the reconstruction where the agreement is first stated. |
| Grok 2–4, physical benchmark and normalization | A free drain can be tested for acceleration and precession. Testing the first paper's inertia mechanism is a further comparison. Attribute equality of the static potentials to the theorem, the Keplerian result to the assigned central acceleration, and the adopted source strength to the reduced monopole calibration. |
| Grok 5, inheritance | Limit the statement about later 4D papers to the lower-order particle dynamics they explicitly carry forward. Cite each import beside that statement. |

The local report gives pinned lookups for all four source findings. The earlier pinned source excerpts reproduce the seven title/DOI assignments and cited matches. Additional lines for the static-limit and lower-order-inheritance wording can be retrieved with the command below. [Source collection](../../../directives/_measurements/S17_record_list_review_p1/sources.md)

```bash
git show 438d270c:research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex | sed -n '359,367p;405,410p'
git show 438d270c:research/4d_2pn/paper/4d_2pn.tex | sed -n '101,105p'
git show 438d270c:research/4d_2_5pn/paper/4d_2_5pn.tex | sed -n '35,38p'
git show 438d270c:research/4d_3pn/paper/4d_3pn.tex | sed -n '34p'
git show 438d270c:research/4d_4pn/paper/4d_4pn.tex | sed -n '870,872p'
```

```text
\medskip\noindent
\textbf{Static Limit Theorem.}
\emph{For a strictly time-independent mass density $\rho(\mathbf{x})$, and for boundary conditions that drive the lag potential $\PhiL$ to zero at spatial infinity, the late-time solution of the toy model field equations reduces to the Newtonian Poisson solution. In particular, the total potential satisfies}
\begin{equation}
  \nabla^2 \Phi(\mathbf{x}) = 4\pi G\,\rho(\mathbf{x}),
  \qquad
  \Phi(\mathbf{x}) = \PhiP(\mathbf{x}),
\end{equation}
\emph{and the scalar lag sector does not contribute to the gravitational field at 0PN order.}

In the static regime described above, test bodies move under the central potential~\eqref{eq:newton_potential_mu}. The equation of motion for a test defect of mass $m$ in this potential is
\begin{equation}
  m \frac{\dd^2 \mathbf{x}}{\dd t^2} = -m\,\bm{\nabla} \Phi(\mathbf{x})
  = -m\,\frac{\mu}{r^3}\,\mathbf{x}.
\end{equation}
We present a controlled derivation of the full conservative two-body second
post-Newtonian (\(\TwoPN\)) sector from the unified \(4\)D toy model.  Carrying
forward the exact projection / Poisson-hook / worldtube reduction logic and the
already-fixed Newtonian+\(\OnePN\) ledger, we close the missing one-body
\(\TwoPN\) response slot by a minimal \(\DtN\)-invariant denominator that
unified $4$D toy model, namely the nonspinning two-body $2.5$PN sector.
Carrying forward the exact projection / leakage / Poisson-hook / worldtube
reduction framework together with the already-fixed Newtonian+$1$PN+$2$PN
conservative hierarchy, we first show that the frozen conservative kernels alone
We present the full conservative two-body $3$PN assembly for the $4$D toy-model program within the same declared closure hierarchy used at $1$PN and $2$PN. Starting from the exact $4$D reduction backbone, the frozen Newtonian/$1$PN/$2$PN ledger, and the $2.5$PN narrowing to the grouped real $P_{20}\oplus P_{21}\oplus P_{22}$ quadrupole bundle, we close the complete conservative near-zone ledger through order $c^{-6}$ in the fixed ADM chart. The strict Schwarzschild one-body gate fixes a minimal three-parameter repair ledger
The present paper imports the conservative Newtonian+$1$PN+$2$PN+$3$PN ledger as
frozen input, not as open fit parameters.  The surviving lower-order viability
conditions are
```
