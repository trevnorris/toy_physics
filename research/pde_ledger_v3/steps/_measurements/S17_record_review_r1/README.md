# S17 review of the corrected orbital record

The reviewed version is preserved at `a9db3f1d`, SHA-256 `8c027a1a5cd61e71bafdab58b0460ef6e2d8bf8de3e63fe664fc87adf62c0bec`. The local reviewer's verdict was “Three substantive qualifications are missing in the first two entries.” Grok wrote “The record is clear as a cite-only import for the Newtonian mouth acceleration, the matched conservative sectors, and both conditional higher-order results, except two passages that widen what those entries may claim.” Grok reviewed only the record text. The local reviewer checked the seven primary sources, the bibliography and the supporting ledger at `438d270c`. [Local report](local.md) · [Grok report](grok.md)

| Verified finding | Disposition |
|---|---|
| Local 1. Local flow in the orbital paper | Its 1PN calculation subtracts smooth background flow in a comoving patch, assumes the residual local mouth fluxes sum to zero and invokes uniform-drift invariance. State those assumptions and their untested applicability to two inward drains. |
| Local 2–3 and Grok 1. Inputs to the coefficient bridge | Attribute each result to its own calculation. The half-displaced-mass result concerns a spherical cross-section, uniformity along the extra dimension, and a low-Mach, incompressible, irrotational exterior with a no-penetration boundary for translation. The quartic energy coefficient uses a trapped mode with isotropic linear dispersion and group velocity equal to the defect velocity. The wake coefficients follow within the separate refractive, EIH and energy-partition matches. |
| Grok 2. Breathing protocol | Imposing the orbital target selects the slope and the flow-to-wave energy ratio. The model's equilibrium virial relation then fixes the pressure–volume share. Cite the two steps separately. |

While these reviews were running, the user clarified the current scientific priority on 2026-11-13. The seven papers supply historical comparisons. Work on gravity beyond 0PN requires the brane's response to gravity and the effect of that response on the attraction. Carry this direction explicitly into the record and retain the 0PN free-drain acceleration as the current force target.

The local command below retrieves the assumptions and virial relation from the pinned sources.

```bash
git show 438d270c:research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex | sed -n '217p;238,242p;963p'
git show 438d270c:research/4d_1pn_bridge/paper/4d_1pn_bridge.tex | sed -n '707,713p;776p;1224,1252p;3557,3566p'
git show 438d270c:research/4d_1pn_full/paper/4d_1pn_full.tex | sed -n '1712,1726p;5972,5992p'
```

```text
The effective $1$PN dynamics used in this paper, however, are constructed to be \emph{uniform--drift invariant} (UDI). In the flux--neutral, irrotational far--field regime, one can show that the net hydrodynamic force on a bounded set of mouths is invariant under a uniform boost of the entire configuration,
By contrast, the local $1$PN derivation in this paper is performed in a \emph{flux--neutral frame}. We work in a comoving patch where the smooth background flow has been factored out and the localized mouths satisfy
\begin{equation}
  \sum_i Q_i^\mathrm{local} = 0 .
\end{equation}
This removes the monopole of the local configuration so that the far field of the \emph{perturbation} decays as $\delta \mathbf{v} = O(r^{-3})$. In this regime boost--dependent surface terms in the traction integral vanish at infinity, and the effective $1$PN forces respect uniform--drift invariance.
Second, the uniform--drift invariance (UDI) postulated in Section~\ref{sec:udi} warrants direct numerical verification. While the analytic $1$PN derivation relies on boost invariance to discard aether--wind terms, fully nonlinear numerical simulations in a boosted frame could quantify the threshold at which UDI breaks down. This would allow us to map the precise limits of the effective Lorentz invariance arising from the superfluid hydrodynamics.
The no-penetration boundary condition on the surface \(r=a\) is
\begin{equation}
  \partial_n \phi \big|_{r=a} = 0
  \quad\Longrightarrow\quad
  \partial_n \phi' \big|_{r=a} = -U\,(\hat{\mathbf e}_1\cdot \hat{\mathbf n})
  = -U\cos\theta,
\end{equation}
In the unified 4D picture, the defect of interest is a corridor extended along the extra coordinate \(w\), i.e.\ a hyper-cylinder with cross-section \(S^2\) in \((x,y,z)\) and approximate \(w\)-extent \(L\).  In the idealized limit where the corridor is uniform in \(w\), the exterior-flow potential can be chosen \(w\)-independent:
\subsection{Rest energy as trapped-mode energy}

We idealize the defect as admitting at least one stable transverse standing mode with characteristic transverse wavenumber $k_\perp$ (set by the throat geometry and boundary/confinement physics). Denote the characteristic propagation speed of this mode by $c$.\footnote{In the broader model, $c$ is the limiting wave speed relevant to operational rulers/clocks (the same constant that appears in the 1PN bookkeeping). The derivations below require only that the trapped mode has an isotropic linear dispersion at the scales of interest.}

For a massless dispersion,
\begin{equation}
\omega = c\,|\mathbf k| ,
\end{equation}
a purely transverse standing mode in the defect rest frame has
\begin{equation}
|\mathbf k| = k_\perp, \qquad \omega_0 = c\,k_\perp,
\end{equation}
so the defect's rest energy is identified with the mode energy,
\begin{equation}
E_0 \equiv \hbar \omega_0 = \hbar c k_\perp,
\qquad
m_0 \equiv \frac{E_0}{c^2}.
\label{eq:E0_m0_def}
\end{equation}
The microphysical relation between $k_\perp$ and geometric parameters (e.g.\ $k_\perp\simeq \pi/a$ for a lowest transverse mode in a simple cavity) is not needed for the relativistic coefficients derived below; only the existence of a stable trapped mode matters.

\subsection{Boosted trapped mode and the \texorpdfstring{$\gamma$}{gamma} energy factor}

Consider uniform motion of the defect at speed $v$ along (say) the $x$-direction. For the trapped wave pattern to be carried along with the moving defect (no secular slippage), its \emph{group velocity} along $x$ must match the defect velocity:
\begin{equation}
v \;=\; v_g \;=\; \frac{\partial \omega}{\partial k_x}.
\label{eq:vg_condition}
\end{equation}
With $\omega = c\sqrt{k_\perp^2 + k_x^2}$, this gives
\begin{enumerate}
\item the far field is low-Mach and effectively incompressible,
\item the exterior flow is irrotational (no vorticity shed to infinity),
\item the defect is sufficiently smooth that a no-penetration condition is valid,
\item the throat is approximately uniform in $w$ over the region that dominates the
      exterior kinetic energy (so that $\partial_w\phi\approx 0$ in that domain).
\end{enumerate}

Departures from these conditions generate corrections (compressibility, radiation
of vorticity, $w$-nonuniform structure, boundary-condition leakage, etc.), but the
Since \(\krho=1\), the target value \(\kPV=3/2\) is equivalent to
\begin{equation}
\frac{d\ln F_{\rm eq}}{d\ln \rho}
=
\krho+\kPV
=
\frac52.
\label{eq:dlnF-target}
\end{equation}
Applying this to \eqref{eq:dlnF-n5} gives
\begin{equation}
\frac{11+7x}{4+5x}=\frac52
\qquad\Longrightarrow\qquad
x=\frac{2}{11}.
\label{eq:x-two-elevenths}
Once \(x=2/11\) is known, the virial relation fixes the internal energy partition
uniquely. Since
\begin{equation}
x=\frac{E_f}{E_w}=\frac{2}{11},
\end{equation}
we have
\begin{equation}
\frac{E_{\rm PV}}{E_w}
=
\frac{1+2x}{3}
=
\frac{1+4/11}{3}
=
\frac{5}{11}.
\end{equation}
Therefore
\begin{equation}
\boxed{
E_w:E_f:E_{\rm PV}=11:2:5.
}
\label{eq:app-11-2-5}
```
