# S17 orbital source excerpts for the review disposition and step record

The command reads the README, the seven main manuscripts, the step plan, the earlier citation note and the accepted S15 and S16 records from commit `438d270c`. Run from the repository root. The numbered lines are excerpts from that revision. The DOI count is a literal search of each paper’s main TeX source for its own assigned DOI.

## Command

```bash
python3 - <<'PY'
import subprocess

ref = "438d270c0fa5576167ab237c61301f1435ef7800"

def source(path):
    return subprocess.check_output(["git", "show", f"{ref}:{path}"], text=True).splitlines()

def excerpt(path, ranges):
    lines = source(path)
    print(path)
    for lo, hi in ranges:
        for number in range(lo, hi + 1):
            print(f"{number} {lines[number - 1]}".rstrip())

papers = [
    ("1pn_orbital_dynamics", "19449058", [(38, 40), (349, 355), (376, 399), (709, 714), (733, 741), (746, 752), (870, 883)]),
    ("4d_1pn_bridge", "19449653", [(91, 91), (1065, 1067), (1113, 1121), (1163, 1169), (1184, 1193), (1278, 1284)]),
    ("4d_1pn_full", "19450102", [(97, 97), (203, 204), (550, 564), (574, 586), (673, 689), (697, 719), (812, 820), (825, 836), (862, 891), (939, 946), (1712, 1748), (2074, 2085), (2120, 2123), (2146, 2154), (3094, 3107), (5972, 5976)]),
    ("4d_2pn", "19450284", [(101, 118), (3901, 3907), (6011, 6015)]),
    ("4d_2_5pn", "19492270", [(35, 60), (791, 797), (815, 825), (2456, 2475), (2537, 2542), (2552, 2561)]),
    ("4d_3pn", "19501724", [(34, 48), (8331, 8338), (8533, 8537), (8969, 8971)]),
    ("4d_4pn", "19561056", [(43, 43), (2169, 2175), (2368, 2397), (2406, 2417), (2662, 2674), (3420, 3425), (3484, 3502), (3520, 3528)]),
]
readme = source("README.md")
for name, suffix, ranges in papers:
    doi = "10.5281/zenodo." + suffix
    [(number, row)] = [(i, line) for i, line in enumerate(readme, 1) if doi in line]
    path = f"research/{name}/paper/{name}.tex"
    print(f"README.md {number} {row}")
    print(f"Assigned DOI occurrences in main TeX: {sum(line.count(doi) for line in source(path))}")
    excerpt(path, ranges)
    print()
excerpt("research/pde_ledger_v3/V3_STEP_PLAN.md", [(747, 763), (1303, 1308)])
excerpt("research/pde_ledger_v2/notes/stages/ledger_stage029_pn_corpus_doi_cite.md", [(19, 21), (50, 53)])
excerpt("research/pde_ledger_v3/steps/S15_two_drain_region.md", [(3, 5), (47, 47), (65, 67)])
excerpt("research/pde_ledger_v3/steps/S16_moving_region.md", [(394, 402)])
PY
```

## Literal stdout

```text
README.md 156 - *Newtonian and 1PN Orbital Dynamics from a Superfluid Defect Toy Model.* [10.5281/zenodo.19449058](https://doi.org/10.5281/zenodo.19449058)
Assigned DOI occurrences in main TeX: 0
research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex
38 We construct and analyze a superfluid defect toy model that reproduces both Newtonian (0PN) gravity and the leading post-Newtonian (1PN) perihelion precession of a test body in a central field. The model consists of a homogeneous superfluid background and sink-like defects whose effective gravitational potential splits into two scalar pieces: an instantaneous Poisson sector $\PhiP$ and a finite-speed ``lag'' sector $\PhiL$ governed by a wave equation with propagation speed $\cS$. For a static central defect in the test-mass limit we show that the retarded scalar solution collapses exactly to the Poisson solution, so that $\PhiL$ vanishes and the near-zone potential is strictly Newtonian, $\Phi_{\rm tot}(r) = -\mu/r$ with no $1/\cS^2$ corrections. As a result the scalar sector generates no 1PN perihelion precession: the entire 1PN correction is encoded kinematically in a position-dependent effective mass $m_{\rm eff}(r) = m[1 + \sigma(r)]$ with $\sigma(r) = \beta \mu/(\cS^2 r)$. Matching the Schwarzschild 1PN precession for $\cS = c$ and $\mu = GM$ requires $\beta = 3$.
39
40 We interpret this coefficient as a sum of three hydrodynamic contributions, $\beta = \kappa_\rho + \kappa_{\rm add} + \kappa_{\rm PV}$, coming respectively from density depletion in the cavitation region, classical added mass of entrained fluid, and internal pressure--volume inertia of a compressible throat. The first two pieces are supported by simplified hydrodynamic derivations and companion Mathematica calculations within the toy-model parameterization, yielding $\kappa_\rho = 1$ and $\kappa_{\rm add} = 1/2$. The remaining piece is not derived here; instead it is fixed by the 1PN matching condition to be $\kappa_{\rm PV} = 3/2$, which we view as an effective field theory constraint on the bulk equation of state and throat compressibility. This provides a concrete hydrodynamic target for future microphysical derivations.
349 The difference $\PhiL$ solves the homogeneous Laplace equation with boundary condition $\PhiL \to 0$ at spatial infinity.\footnote{More generally, one can impose any common boundary condition on $\PhiP$ and $\PhiL$ at a large but finite radius and then take the limit as the boundary is pushed to infinity.} The only regular solution of this homogeneous equation is a constant. Any constant offset in the total potential is physically irrelevant, since only gradients of $\Phi$ enter the acceleration. We are therefore free to choose this constant to be zero, which yields
350 \begin{equation}
351   \PhiL(\mathbf{x}) \to 0,
352   \qquad
353   \Phi(\mathbf{x}) \equiv \PhiP(\mathbf{x}) + \PhiL(\mathbf{x}) \to \PhiP(\mathbf{x}),
354 \end{equation}
355 in the strict static limit.
376 To make the connection with standard Newtonian gravity explicit, consider a single isolated sink defect of mass $M$ localized at the origin. In the point-particle limit the mass density is
377 \begin{equation}
378   \rho(\mathbf{x}) = M\,\delta^{(3)}(\mathbf{x}).
379 \end{equation}
380 The Static Limit Theorem implies that the total potential satisfies the Poisson equation
381 \begin{equation}
382   \nabla^2 \Phi(\mathbf{x}) = 4\pi G M\,\delta^{(3)}(\mathbf{x}),
383 \end{equation}
384 with boundary condition $\Phi(\mathbf{x}) \to 0$ as $|\mathbf{x}| \to \infty$.
385
386 The spherically symmetric solution is the familiar Newtonian potential,
387 \begin{equation}
388   \Phi(r) = -\frac{\mu}{r},
389   \qquad
390   \mu \equiv G M,
391   \label{eq:newton_potential_mu}
392 \end{equation}
393 where $r = |\mathbf{x}|$ and $\mu$ is the gravitational parameter. The corresponding radial acceleration experienced by a test defect at radius $r$ is
394 \begin{equation}
395   a_r(r) = -\frac{\dd \Phi}{\dd r}
396           = -\frac{\mu}{r^2},
397   \label{eq:newton_accel}
398 \end{equation}
399 which is precisely the inverse-square law. In this way, the effective parameter $\mu$ in the toy model is directly identified with $GM$ in the usual Newtonian and GR notation.
709 In the present context, the defect throat plays a similar role to a compact body moving through the superfluid: its motion induces a flow pattern that carries kinetic energy. If the throat maintains an approximately fixed effective volume $V_{\mathrm{cav}}$ and the flow around it is approximately potential on orbital scales, it is natural to assign it an added mass
710 \begin{equation}
711   m_{\mathrm{add}}
712   \simeq \frac{1}{2}\,\rho(r) V_{\mathrm{cav}},
713 \end{equation}
714 where $\rho(r)$ is the local background density. In terms of the reference mass $m = \rho_0 V_{\mathrm{cav}}$, this corresponds to a fractional contribution
733 once all factors are expressed in terms of the PN parameter $\mu/[\cS^2 a(1-e^2)]$. This effective identification is what enters the derivation of the total coefficient $2\beta$ in Eq.~\eqref{eq:Delta_phi_total_beta}: the added-mass contribution combines with the inertia modulation to shift the coefficient of $u$ in the orbit equation.
734
735 A detailed calculation in the companion notebook \texttt{kappa\_add\_derivation.wl} shows that the cross-terms vanish by symmetry and the dipole self-energy of the entrained flow yields an effective inertia equal to one half of the displaced fluid mass, supporting the identification $\kappa_{\mathrm{add}} = 1/2$ in our parametrization.
736
737 We therefore associate
738 \begin{equation}
739   \kappa_{\mathrm{add}} = \frac{1}{2}
740 \end{equation}
741 with the classical added mass of the defect throat moving through the superfluid.
746 The remaining piece, $\kappa_{\mathrm{PV}}$, is associated with unsteady pressure--volume work and compressible response in the vicinity of the defect throat. In the language of Eq.~\eqref{eq:beta_decomp}, this contribution would encode any additional inertia that arises when the throat accelerates and the surrounding fluid must be compressed and rarefied in time.
747
748 In the corrected 1PN accounting the scalar lag sector contributes nothing, so matching the full GR precession with $\beta = 3$ requires
749 \begin{equation}
750   \kappa_{\mathrm{PV}}
751   \simeq 3 - 1 - \frac{1}{2}
752   = \frac{3}{2},
870 is measured in the same way as in the scalar-only case. Over the range of orbits where the PN parameter $\mu/[\cS^2 a(1-e^2)]$ is small, the measured precession agrees with the analytic expression
871 \begin{equation}
872   \Delta \varphi_{\mathrm{tot}}
873   \simeq 2\pi\beta\,\frac{\mu}{\cS^2 a (1-e^2)}
874   = \frac{6\pi \mu}{\cS^2 a (1-e^2)},
875 \end{equation}
876 to within the same numerical accuracy. When $\mu = GM$ and $\cS = c$, this matches the standard GR 1PN formula,
877 \begin{equation}
878   \Delta \varphi_{\mathrm{GR}}
879   = \frac{6\pi GM}{c^2 a (1-e^2)}.
880 \end{equation}
881 In other words, once the inertia modulation with $\beta = 3$ is included, the toy model reproduces the full GR 1PN precession in the test-mass, central-field limit at the level probed by the reduced-orbit numerics.
882
883 We have also run fully three-dimensional ``dynamic source'' experiments in which both the central mass and orbiting test bodies are represented as defects in the grid-based PDE solver, with the density field $\rho(\mathbf{x},t)$ recomputed from the ensemble at each time step. In these runs the central mass is not pinned at the origin but moves slightly in response to the orbital back-reaction, ensuring that the source density is genuinely time-dependent and that the lag equation remains active. In a ``slow-light'' stress-test regime with artificially small $\cS$ (so that the PN parameter is enhanced), the measured precessions are large compared to the grid-induced noise and exhibit the same scaling with semi-major axis and eccentricity as the analytic 1PN prediction. These fully PDE-based experiments are consistent with the reduced-orbit results and provide an independent confirmation that the finite-speed scalar lag plus position-dependent kinetic prefactor can together reproduce GR-like 1PN dynamics in the toy model.

README.md 145 - *Deriving Key Post-Newtonian Coefficients from the Unified 4D Toy Model.* [10.5281/zenodo.19449653](https://doi.org/10.5281/zenodo.19449653)
Assigned DOI occurrences in main TeX: 0
research/4d_1pn_bridge/paper/4d_1pn_bridge.tex
91 This paper is a derivation supplement to the unified $4$D throat model of Ref.~\cite{norris2026_4d_model}, written to make explicit which coefficients used earlier in the series are fixed by the $4$D structure (and which must remain symbolic pending additional closure). We show that a defect uniform in the extra coordinate behaves as a throat with hypercylindrical topology, for which the exterior potential--flow problem yields the classical added-mass coefficient $\kappa_{\mathrm{add}}=\tfrac12$ (whereas a $4$D bubble would instead give $\tfrac13$), providing a sharp topology discriminator. We then show that matching the weak-field refractive coefficient uniquely selects the stiff polytropic exponent $n=5$. In the $1$PN interaction sector, EIH matching produces a one-parameter family of couplings, which collapses uniquely to $\alpha^2=\tfrac34$ and $K_{\rm vec}=2/\pi^2$ once the $n=5$ thermodynamic relation is imposed. Treating rest mass as wave-supported energy reproduces the special-relativistic expansion $E=\gamma E_0$, including the universal $v^4$ coefficient $3/8$. We also summarize a controlled reduction to brane-effective equations, making leakage/source terms explicit, and we outline the additional response modeling required to determine the remaining undetermined coefficients.
1065   \label{eq:EIH_target_cross}
1066 \end{equation}
1067 In this section we use \eqref{eq:EIH_target_cross} as the target ``data'' that the toy model must reproduce.
1113 Equations \eqref{eq:invariant1}--\eqref{eq:invariant2} exhibit the EIH-matching structure clearly:
1114 \emph{EIH matching alone} fixes two invariants and leaves a one-parameter family.
1115 Solving for $(\alpha^2,K_{\rm vec})$ in terms of $a_H^2$ gives
1116 \begin{equation}
1117   \alpha^2(a_H^2)=\frac{3}{4}(1-a_H^2),
1118   \qquad
1119   K_{\rm vec}(a_H^2)=\frac{2}{\pi^2(1-a_H^2)},
1120   \qquad
1121   0\le a_H^2<1,
1163 Using \eqref{eq:U_over_P}, the simplest thermodynamic closure consistent with that definition is
1164 \begin{equation}
1165   \alpha^2_{\rm thermo}
1166   = 1-\frac{U}{P}
1167   =1-\frac{1}{n-1}.
1168   \label{eq:alpha_thermo_general}
1169 \end{equation}
1184 Substituting back into \eqref{eq:family_alphaK} yields the unique point
1185 \begin{equation}
1186   \boxed{
1187   \alpha^2=\frac{3}{4},\qquad
1188   a_H=0,\qquad
1189   K_{\rm vec}=\frac{2}{\pi^2}.
1190   }
1191   \label{eq:alphaK_final}
1192 \end{equation}
1193 As a direct check, inserting \eqref{eq:alphaK_final} into \eqref{eq:CL_Cpar_forms} reproduces the EIH targets:
1278 &= 1 + \frac{1}{2}\frac{v^2}{c^2} + \frac{3}{8}\frac{v^4}{c^4} + \mathcal O\!\left(\frac{v^6}{c^6}\right),
1279 \\[4pt]
1280 E(v)
1281 &= E_0
1282 + \frac{1}{2}\left(\frac{E_0}{c^2}\right)v^2
1283 + \frac{3}{8}\left(\frac{E_0}{c^2}\right)\frac{v^4}{c^2}
1284 + \mathcal O\!\left(\frac{v^6}{c^4}\right)

README.md 148 - *A Controlled Derivation of the Full First Post-Newtonian Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19450102](https://doi.org/10.5281/zenodo.19450102)
Assigned DOI occurrences in main TeX: 0
research/4d_1pn_full/paper/4d_1pn_full.tex
97 We present a controlled derivation of the full conservative two-body first post-Newtonian (1PN) sector from the four-spatial-dimensional brane--bulk toy model. Starting from the exact bulk action, projection map, open-system balance laws, and quasi-static near-zone Poisson hook on the brane, we use a coherent-defect / small-body worldtube reduction to obtain the Newtonian point-particle limit and then assemble the conservative 1PN ledger sector by sector. Newtonian matching fixes the scalar mass-dressing coefficient \(\krho=1\) (historically written \(q=1\) in earlier gravity-side bookkeeping), weak-field optical consistency fixes the polytropic exponent \(n=5\), and wave-supported rest energy yields the universal free-particle \(v^4\) coefficient. Local mass scaling plus pair counting generate the static \(G^2/r^2\) term, exterior potential-flow topology fixes the added-mass coefficient \(\kappa_{\rm add}=1/2\) for a \(w\)-uniform throat, and a constructive isotropic Fourier-space wake basis reproduces the Einstein--Infeld--Hoffmann coefficients \(C_\parallel=-7/2\) and \(C_L=-1/2\) with \(\alpha^2=3/4\), \(a_H=0\), and \(K_{\rm vec}=2/\pi^2\). Under an explicit adiabatic one-degree-of-freedom breathing closure we further obtain \(\kappa_{\rm PV}=3/2\), so that \(\beta_{1\rm PN}=\kappa_\rho+\kappa_{\rm add}+\kappa_{\rm PV}=3\). The assembled reduced two-body Lagrangian is then algebraically identical to the standard conservative EIH form, and the test-mass limit reproduces the standard perihelion shift. The result is a full conservative 1PN derivation within a declared closure hierarchy that makes explicit which ingredients are fixed, protocol-fixed, and still symbolic, with remaining freedom confined to generalized internal response, finite-size structure, and beyond-1PN open-system effects.
203     Newtonian scalar coupling & \(\krho=1\) & Newtonian matching of the scalar-dressed particle action (historically \(q=1\)) & \Controlled \\
204     EOS exponent & \(n=5\) & weak-field optical/refraction coefficient & \Controlled \\
550   \item \(\Proj{\rho}\approx \rho_0\) is slowly varying over the region of interest,
551   \item \(\partial_t\Proj{\rho}\) is quasi-static at the order retained,
552   \item the advective correction \((\gradThree \Proj{\rho})\!\cdot\!(\gradThree\phibr+\mathbf{v}_T)\) is subleading, and
553   \item the transverse sector is negligible or perturbative,
554 \end{enumerate}
555 then \eqref{eq:dictionary-longitudinal-identity} reduces to the Poisson hook
556 \begin{equation}
557 \lapThree \phibr
558 \;\simeq\;
559 \frac{S_{\rm eff}}{\rho_0},
560 \qquad
561 S_{\rm eff}\equiv \Sleak\ \text{plus any retained slow-source terms}.
562 \label{eq:dictionary-poisson-hook}
563 \end{equation}
564 For a localized effective source,
574 r=|\x|.
575 \label{eq:dictionary-inverse-square}
576 \end{equation}
577 This is the route by which inverse-square longitudinal behavior enters the present
578 program. In the next section, a controlled worldtube reduction turns that kinematic
579 Poisson-sector structure into the effective Newtonian potential \(\PhiN\) used in the
580 particle Lagrangian.
581
582 It is important not to overstate the claim at this stage. Equation
583 \eqref{eq:dictionary-longitudinal-identity} is exact. Equation
584 \eqref{eq:dictionary-poisson-hook} is a regime statement. The identification of the
585 resulting near-zone scalar with the Newtonian potential used by a reduced point particle
586 is a further controlled constitutive step carried out later in the paper.
673 Operationally, we therefore define the measured Newtonian potential \(\PhiN\) as the
674 properly normalized near-zone scalar channel on the brane,
675 \begin{equation}
676 \PhiN(\x,t)\equiv \lambda_\Phi\,\phibr(\x,t),
677 \label{eq:newtonian-potential-definition}
678 \end{equation}
679 with \(\lambda_\Phi\) fixed by monopole calibration. For an isolated compact source \(B\),
680 the near-zone solution outside the source support is then
681 \begin{equation}
682 \Phi_B(\x,t)
683 =
684 -\frac{Gm_B}{|\x-\mathbf X_B(t)|}
685 + O\!\left(\frac{Gm_B a_B^2}{|\x-\mathbf X_B|^3}\right),
686 \label{eq:isolated-source-newtonian-potential}
687 \end{equation}
688 where \(G\) is the measured brane normalization of the scalar channel and \(a_B\) is the
689 source size. The remaining task in this section is to show that a second compact defect
697   \item \textbf{Compact support / small body.}
698   The defect occupies a region of characteristic size \(a\) that is small compared with
699   the inter-body separation and with the scale over which the external field varies:
700   \(a\ll r\).
701
702   \item \textbf{Smooth external field across the support.}
703   Over the defect worldtube, the external potential may be Taylor-expanded about the
704   center-of-mass position.
705
706   \item \textbf{Negligible worldtube boundary fluxes at the retained order.}
707   Any residual leakage, entrainment, or stress flux through the defect boundary is beyond
708   the Newtonian order kept here.
709
710   \item \textbf{Defect-centered coordinates.}
711   The spatial origin is chosen so that the defect dipole moment vanishes.
712
713   \item \textbf{Quasi-static near zone.}
714   The dropped terms in the exact longitudinal identity remain parametrically small, so
715   the Poisson hook of \cref{eq:dictionary-poisson-hook} is valid.
716
717   \item \textbf{Conservative sector only.}
718   Radiation, dissipative exchange, and higher-response effects are deferred to later
719   work and are not part of the Newtonian limit derived in this section.
812 At Newtonian order, the second term is an internal energy contribution. It can renormalize
813 the defect rest energy, but it does not alter the leading external center-of-mass
814 kinematics so long as the coherent-translation closure is maintained.
815
816 This is the first place where the reduced description already displays an important
817 feature that later becomes the Newtonian matching condition of the scalar sector:
818 the same total monopole \(M\) appears in the kinetic term and in the leading coupling to
819 the external potential below. In the reduced Newtonian description, inertial and
820 gravitational mass therefore coincide automatically.
825 Let the defect move in a smooth external potential \(\PhiN(\x,t)\). The worldtube force is
826 \begin{equation}
827 \mathbf F_{\rm ext}(t)
828 =
829 -\int_{W_t}\rho_{\rm def}(\boldsymbol{\xi},t)\,
830 \gradThree \PhiN\!\bigl(\mathbf X_{\rm cm}(t)+\boldsymbol{\xi},t\bigr)\,\dd^3\xi
831 +
832 \mathbf F_{\partial W},
833 \label{eq:worldtube-force}
834 \end{equation}
835 where \(\mathbf F_{\partial W}\) denotes boundary-flux terms. By assumption,
836 \(\mathbf F_{\partial W}\) is beyond the retained Newtonian order.
862 \eqref{eq:defect-dipole-zero}. Therefore the leading force is monopolar:
863 \begin{equation}
864 \mathbf F_{\rm ext}
865 =
866 - M\,\gradThree \PhiN(\mathbf X_{\rm cm},t)
867 +
868 O\!\left(Q\,\nabla^3\PhiN\right).
869 \label{eq:monopole-force-leading}
870 \end{equation}
871
872 For an external field sourced by a distant compact body,
873 \(\PhiN\sim Gm/r\), the monopole force scales as \(MGm/r^2\), whereas the leading
874 finite-size correction scales as \(MGm\,a^2/r^4\). Their ratio is therefore
875 \begin{equation}
876 \frac{|\mathbf F_{\rm quad}|}{|\mathbf F_{\rm mono}|}
877 =
878 O\!\left(\frac{a^2}{r^2}\right).
879 \label{eq:quadrupole-suppression}
880 \end{equation}
881 This is the explicit small-body control parameter of the Newtonian reduction.
882
883 Combining \eqref{eq:defect-total-momentum} with \eqref{eq:monopole-force-leading}, we
884 obtain the reduced center-of-mass equation of motion
885 \begin{equation}
886 M\ddot{\mathbf X}_{\rm cm}
887 =
888 - M\,\gradThree\PhiN(\mathbf X_{\rm cm},t)
889 +
890 O\!\left(M\,\frac{a^2}{r^2}\,\frac{\PhiN}{r}\right).
891 \label{eq:newtonian-center-of-mass-eom}
939 Equation \eqref{eq:newtonian-two-body-lagrangian} is the baseline pair Lagrangian on top
940 of which the subsequent \(1\)PN self, cross, and static nonlinear sectors will be added.
941
942 \claimclass{\Controlled for the identification of the Poisson-hook scalar with the
943 effective Newtonian potential, the coherent-defect / small-body worldtube closure, and
944 the neglect of boundary-flux terms at the retained order. Exact for the Euler--Lagrange
945 equivalence between \cref{eq:newtonian-particle-lagrangian} and the leading monopole
946 equation of motion once the reduced Lagrangian has been stated.}
1712 Since \(\krho=1\), the target value \(\kPV=3/2\) is equivalent to
1713 \begin{equation}
1714 \frac{d\ln F_{\rm eq}}{d\ln \rho}
1715 =
1716 \krho+\kPV
1717 =
1718 \frac52.
1719 \label{eq:dlnF-target}
1720 \end{equation}
1721 Applying this to \eqref{eq:dlnF-n5} gives
1722 \begin{equation}
1723 \frac{11+7x}{4+5x}=\frac52
1724 \qquad\Longrightarrow\qquad
1725 x=\frac{2}{11}.
1726 \label{eq:x-two-elevenths}
1727 \end{equation}
1728 Therefore the adiabatic one-degree-of-freedom closure predicts
1729 \begin{equation}
1730 \boxed{\kPV=\frac32.}
1731 \label{eq:kPV-three-halves}
1732 \end{equation}
1733
1734 The same solution fixes the internal energy partition:
1735 \begin{equation}
1736 E_w:E_f:E_{\rm PV}=11:2:5,
1737 \label{eq:partition-11-2-5}
1738 \end{equation}
1739 or, equivalently,
1740 \begin{equation}
1741 (f_w,f_f,f_{\rm PV})
1742 =
1743 \left(\frac{11}{18},\frac{2}{18},\frac{5}{18}\right).
1744 \label{eq:partition-fractions}
1745 \end{equation}
1746 This is a useful point to emphasize. Even though \(\kPV\) is a protocol-level quantity,
1747 once the protocol is declared its value is not merely fitted; it comes with a definite
1748 predicted partition of internal energy among wave, flow, and pressure-volume sectors.
2074 Imposing the EIH targets \eqref{eq:wake-eih-target} on
2075 \eqref{eq:wake-C-parallel-long} gives
2076 \begin{align}
2077 -\frac72
2078 &=
2079 K_{\rm vec}\pi^2\bigl(-1+a_H^2-\alpha^2\bigr),
2080 \label{eq:wake-match-para}
2081 \\[3pt]
2082 -\frac12
2083 &=
2084 K_{\rm vec}\pi^2\bigl(-1+a_H^2+\alpha^2\bigr).
2085 \label{eq:wake-match-long}
2120 Thus EIH matching alone does \emph{not} yet pick a unique point. It fixes two invariants
2121 and leaves one apparent freedom, which may be parameterized by \(a_H^2\). Put
2122 differently, the cross-tensor structure is already highly constrained by the wake basis,
2123 but one more physical input is needed to collapse the family.
2146 The minimal thermodynamic closure used here identifies the longitudinal wake fraction with
2147 the non-internal part of the barotropic budget:
2148 \begin{equation}
2149 \alpha^2_{\rm thermo}
2150 \equiv
2151 1-\frac{U}{P}
2152 =
2153 1-\frac{1}{n-1}.
2154 \label{eq:wake-thermo-alpha}
3094     \(\krho\) & fixed & Newtonian matching gives \(\krho=1\) (historically \(q=1\)) \\
3095     \(n\) & fixed & optical consistency gives \(n=5\) \\
3096     \(1/8\) quartic Lagrangian coefficient & fixed & universal relativistic/wave-supported kinematics \\
3097     Static \(G^2/r^2\) coefficient & fixed & local mass scaling + pair counting gives \(-1/2\) \\
3098     \(\kadd\) & fixed & throat topology gives \(\kadd=1/2\) \\
3099     \(\alpha^2\), \(a_H\), \(K_{\rm vec}\) & fixed & wake + thermodynamic closure give \(\alpha^2=3/4\), \(a_H=0\), \(K_{\rm vec}=2/\pi^2\) \\
3100     \(C_\parallel\), \(C_L\) & fixed & \((-7/2,-1/2)\) from the wake construction \\
3101     \(L_{\rm 1PN}^{\rm derived}\) & fixed & exact equality to the conservative EIH target \\
3102     \(\kPV\) & protocol-fixed & equals \(3/2\) in the adiabatic 1-DOF closure \\
3103     \(\betaPN\) & protocol-fixed & equals \(3\) once \(\kPV=3/2\) is adopted \\
3104     \(\KEOS\) & symbolic & dimensional EOS normalization remains free \\
3105     \(\Wproj(w)\), \(\Zloc(w)\) & symbolic & only structural properties are used here \\
3106     finite-size moments / higher multipoles & symbolic & suppressed in the monopole small-body limit \\
3107     general response functions & symbolic & \(\kPV(\omega,\text{protocol})\), damping, compliance, etc. \\
5972 Once \(x=2/11\) is known, the virial relation fixes the internal energy partition
5973 uniquely. Since
5974 \begin{equation}
5975 x=\frac{E_f}{E_w}=\frac{2}{11},
5976 \end{equation}

README.md 149 - *A Controlled Derivation of the Full Conservative Second Post-Newtonian Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19450284](https://doi.org/10.5281/zenodo.19450284)
Assigned DOI occurrences in main TeX: 0
research/4d_2pn/paper/4d_2pn.tex
101 We present a controlled derivation of the full conservative two-body second
102 post-Newtonian (\(\TwoPN\)) sector from the unified \(4\)D toy model.  Carrying
103 forward the exact projection / Poisson-hook / worldtube reduction logic and the
104 already-fixed Newtonian+\(\OnePN\) ledger, we close the missing one-body
105 \(\TwoPN\) response slot by a minimal \(\DtN\)-invariant denominator that
106 preserves the frozen \(\OnePN\) coefficient and reproduces the isotropic
107 Schwarzschild test-mass target through \(\TwoPN\).  Minimal Bernoulli
108 exactification and quadratic local mass scaling then supply the self/static
109 \(\TwoPN\) seed, with the pure-static isotropic gate fixing
110 \(\lambda_\rho=\tfrac12\).  The remaining comparable-mass mismatch against the
111 generic-frame \(\ADM\) target collapses to a finite algebraic residual with a
112 unique solution for the added cross-sector coefficients \(\{q_i,t_i,s_1\}\).
113 After assembly, the full reduced two-body Lagrangian Legendre-transforms
114 exactly to the standard generic-frame conservative \(\ADM\) Hamiltonian through
115 \(\TwoPN\).  The paper therefore gives a full conservative \(\TwoPN\)
116 derivation within a declared closure hierarchy, while the appendices preserve
117 the broader throat-response constructions and the remaining dynamic
118 \(\DtN\)/PDE observables.
3901 What remains open is the deeper derivation of this closure from the fully dynamical
3902 inner-throat PDE and its frequency-dependent \(\DtN\) operator.
3903 In particular, the present appendix does not derive the analytic structure of
3904 \(F(G)\) beyond the minimal quadratic term needed at \(\TwoPN\), nor does it determine
3905 the dynamic pole scales of the appendix-side mouth-response operator.  Those are genuine
3906 inner-throat observables and belong to the later PDE / \(\DtN\) program rather than to
3907 the conservative algebraic proof completed in the main text.
6011 PDE-side no-go benchmark and for the one-pole-per-channel axisymmetric
6012 \(\DtN\) completion as the minimal dynamic scaffolding consistent with the solved
6013 conservative \(\TwoPN\) cross target.  This appendix does not claim a full moving-throat
6014 PDE solution; it identifies the exact zero-frequency target and the first open dynamic
6015 observables that the inner-throat program must still derive.}

README.md 150 - *A Conditional Derivation of the Point-Particle 2.5PN Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19492270](https://doi.org/10.5281/zenodo.19492270)
Assigned DOI occurrences in main TeX: 0
research/4d_2_5pn/paper/4d_2_5pn.tex
35 unified $4$D toy model, namely the nonspinning two-body $2.5$PN sector.
36 Carrying forward the exact projection / leakage / Poisson-hook / worldtube
37 reduction framework together with the already-fixed Newtonian+$1$PN+$2$PN
38 conservative hierarchy, we first show that the frozen conservative kernels alone
39 cannot generate a local $2.5$PN reaction term.  We then use a minimal retarded
40 quadrupole prototype to recover the Burke--Thorne member of the standard
41 Iyer--Will family, which fixes the structural target for any successful
42 completion.  The main body of the paper is a channel audit.  On the natural
43 compact/passive/outgoing small-body branch, the dangerous lower odd scalar
44 channels are derivative-mediated or super-Ohmic and are pushed above universal
45 point-particle $2.5$PN; the carried dipole wake is likewise demoted by the
46 outgoing $\ell=1$ threshold law and strict small-body scaling.  The only
47 surviving universal point-particle route is the orbital/worldtube STF
48 quadrupole.  In particular, the solved $2$PN $P_2$ bundle already exhausts the
49 real STF $\ell=2$ content, the low-frequency $P_2\to\mathrm{STF}$ matching is
50 scalar under isotropy, the passive outgoing $\ell=2$ branch carries the correct
51 positive $+i\omega^5$ parity/sign, and the static overlap gate passes on the
52 natural branch.  The strongest honest result is therefore a \emph{conditional}
53 $2.5$PN theorem within a declared closure hierarchy: if the completed
54 moving-throat PDE realizes this passive/outgoing quadrupole branch, then the
55 reduced nonspinning two-body dynamics contains the standard local quadrupolar
56 Burke--Thorne / Iyer--Will reaction term, with all other currently derived odd
57 channels entering only as higher-order finite-size corrections.  The remaining
58 serious open problem is no longer generic dissipative physics, but the final
59 passive/outgoing quadrupole normalization,
60 \(\hat m_0^2\Gamma_5 = 2G/(5c^5)\), on the actual moving-throat branch.
791 Comparison with the standard Burke--Thorne force,
792 \[
793 a_i^{\rm BT}=-\frac{2G}{5c^5}\,x^j\bigl(I_{ij}^{\rm orb}\bigr)^{(5)},
794 \]
795 fixes the canonical GR normalization target
796 \[
797 \gamma_{\rm GR}=\frac{2G}{5c^5}.
815 then the canonically normalized odd coefficient is controlled only by the product
816 \[
817 \gamma_{\rm toy}=\hat m_0^2\Gamma_5.
818 \]
819 Matching to Burke--Thorne is therefore equivalent to the scalar equality
820 \[
821 \boxed{
822 \hat m_0^2\Gamma_5=\frac{2G}{5c^5}.
823 }
824 \]
825 This is the sharpest normalization target used in the present paper.
2456 The cleanest theorem statement is now conditional rather than unconditional.
2457 Assume that the completed moving-throat PDE satisfies all of the following on
2458 its physical low-frequency branch:
2459 \begin{enumerate}
2460   \item there is no new direct scalar source beyond the continuity-mediated and
2461   reciprocal-boundary ontology isolated in
2462   Sections~\ref{sec:scalar-audit-demotion} and
2463   \ref{sec:benchmark-framework},
2464   \item the scalar and vector outlets are compact and outgoing at low
2465   frequency, so that the surviving scalar and dipole channels begin at the
2466   derivative/high-order thresholds established earlier,
2467   \item the low-frequency defect/worldtube sector is rotationally invariant,
2468   \item the surviving STF quadrupole branch is passive and outgoing,
2469   \item the strict small-body hierarchy remains valid,
2470   \item and the static orbital/source overlap remains on the natural branch,
2471   \[
2472   \hat m_0 = 1 + O\!\left(\frac{a^2}{r^2}\right),
2473   \qquad
2474   \hat m_0\neq 0.
2475   \]
2537 and that the outgoing $\ell=2$ branch has the correct positive
2538 $+i\omega^5$ parity/sign.  In other words, the surviving universal branch is
2539 already fixed as a source-content, parity, and sign statement.  What is no
2540 longer open is whether the model can support \emph{some} quadrupolar odd term.
2541 What remains open is whether the actual moving-throat PDE realizes the branch
2542 with the correct normalization.
2552 and the retarded STF matching coefficient on the same branch is $\Gamma_5$, then
2553 the canonically normalized odd coefficient is controlled only by the product
2554 \[
2555 \gamma_{\rm quad}^{\rm eff}=\hat m_0^2\Gamma_5.
2556 \]
2557 Matching to Burke--Thorne is therefore equivalent to the scalar condition
2558 \[
2559 \boxed{
2560 \hat m_0^2\Gamma_5=\frac{2G}{5c^5}.
2561 }

README.md 151 - *A Full Conservative Derivation of the Two-Body 3PN Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19501724](https://doi.org/10.5281/zenodo.19501724)
Assigned DOI occurrences in main TeX: 0
research/4d_3pn/paper/4d_3pn.tex
34 We present the full conservative two-body $3$PN assembly for the $4$D toy-model program within the same declared closure hierarchy used at $1$PN and $2$PN. Starting from the exact $4$D reduction backbone, the frozen Newtonian/$1$PN/$2$PN ledger, and the $2.5$PN narrowing to the grouped real $P_{20}\oplus P_{21}\oplus P_{22}$ quadrupole bundle, we close the complete conservative near-zone ledger through order $c^{-6}$ in the fixed ADM chart. The strict Schwarzschild one-body gate fixes a minimal three-parameter repair ledger
35 \[
36 (\mu_{\rho,3},d_3,s_{24})=\left(\frac14,-\frac{45}{4},-\frac1{16}\right).
37 \]
38 The comparable-mass sector then splits exactly into three lanes: the ordinary-Lagrangian counterimage of the universal free two-body relativistic Hamiltonian, an exact grouped real $P_2$ middle-block closure, and a unique geometry-side static completion. Equivalently,
39 \[
40 \Delta \hat L_3^{\rm GR}
41 =
42 \Delta l_1 v^8
43 +
44 L_{P_2}^{\rm mid}
45 +
46 \Delta l_{15}^{(g)}U^4,
47 \]
48 with $\Delta H_3=-\Delta L_3$ in the fixed ADM chart. The apparent pure-kinetic remainder collapses to the compiler image of the free $3$PN Hamiltonian, while the remaining static slot is uniquely geometry-side after sigma-collapse. Thus the conservative $3$PN bottleneck is no longer algebraic: within the declared hierarchy, the two-body $3$PN ledger is fully assigned, and the live unresolved problem is the deeper moving-throat derivation of the grouped-$P_2$ and geometry constitutive lanes together with the already-narrow outgoing quadrupole-normalization gap inherited from the $2.5$PN program.
8331 \subsection{Exact fit to the solved GR middle block}
8332 \label{app:richer_grouped_p2_compiler_inverse_map:gr_fit}
8333
8334 Now substitute the solved GR middle-block residual
8335 \[
8336 (L_6,\dots,L_{14})=(\Delta l_6,\dots,\Delta l_{14}).
8337 \]
8338 The inverse map then gives one unique grouped coefficient set:
8533 The same conclusion can be read directly on the target coefficients. The imported generic-frame static target coefficient is
8534 \[
8535 c_U^{\rm target}=-\frac{227}{24}+\frac{21\pi^2}{32},
8536 \]
8537 while the richer grouped-$P_2$ middle-block fit already forces the static companion
8969 In that precise sense, the conservative $3$PN two-body ledger is now fully assigned in the fixed ADM chart. What remains open is no longer algebraic coefficient matching.
8970
8971 The open problems are now dynamical and constitutive. First, the richer grouped real $P_2$ constitutive map should be derived directly from the moving-throat problem rather than by inverse matching. Second, the unique geometry completion should be derived directly from the throat-side geometry/compressibility dynamics. Third, those two conservative lanes must be carried into the already narrowed $2.5$PN bridge, where the live universal route has already been reduced to the grouped real quadrupole bundle and its passive/outgoing normalization problem.

README.md 152 - *A Conditional Full Conservative Derivation of the Two-Body 4PN Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19561056](https://doi.org/10.5281/zenodo.19561056)
Assigned DOI occurrences in main TeX: 0
research/4d_4pn/paper/4d_4pn.tex
43 We present the full conservative fourth post-Newtonian (4PN) assembly for the 4D toy-model program within the same declared closure hierarchy used at lower orders. The strict Schwarzschild one-body gate, exact quartic Legendre compiler, Hamiltonian-first generic-frame lift, canonical comparable-mass local residual, aligned-seed correction, and referee reconstruction together show that the entire local instantaneous 4PN sector already has an exact generic-frame ordinary representative. We then isolate the hereditary sector and show that its coefficient is not a new independent datum, but is fixed by the same canonically normalized orbital/worldtube STF quadrupole branch that controls the 2.5PN Burke--Thorne channel, through $C_{\rm tail}=(GM/2c^3)\,\gamma_{\rm quad}^{\rm eff}$. Hence 4PN opens no new normalization constant: if $\gamma_{\rm quad}^{\rm eff}=2G/(5c^5)$, then $C_{\rm tail}=G^2M/(5c^8)$ and the full conservative 4PN local+tail structure matches the standard GR target. The remaining gap is therefore the same passive/outgoing quadrupole-normalization theorem on the moving-throat branch already isolated by the 2.5PN program.
2169 side is therefore sharply narrower than the local side ever was.  The purpose of
2170 this section is to make that narrowing explicit.  First, we freeze the standard
2171 conservative GR tail functional and its universal local logarithmic shadow.
2172 Second, we show that the hereditary coefficient is tied exactly to the same
2173 canonically normalized STF quadrupole coefficient that already controls the
2174 $2.5$PN Burke--Thorne channel.  Third, we isolate the only remaining model-side
2175 ambiguity as a single scalar transport factor rather than a new tensor channel,
2368 be the canonically normalized STF quadrupole coefficient inherited from the
2369 $2.5$PN audit.  Then the hereditary side may be parameterized by one scalar
2370 transport factor $\Theta_{\rm tail}$ through
2371 \begin{equation}
2372 \boxed{
2373 \alpha_{\rm tail}^{\rm toy}
2374 =
2375 \Theta_{\rm tail}
2376 \frac{GM}{2c_s^3}
2377 \gamma_{\rm toy}.
2378 }
2379 \label{eq:4pn-tail-toy-transport}
2380 \end{equation}
2381 On the branch with $c_s=c$ and
2382 \begin{equation}
2383 \gamma_{\rm toy}=
2384 \gamma_{\rm quad}^{\rm eff},
2385 \label{eq:4pn-tail-toy-canonical-branch}
2386 \end{equation}
2387 Eq.~\eqref{eq:4pn-tail-toy-transport} reduces to the exact bridge
2388 \eqref{eq:4pn-tail-exact-bridge} when
2389 \begin{equation}
2390 \Theta_{\rm tail}=1.
2391 \label{eq:4pn-tail-theta-one}
2392 \end{equation}
2393 This is the sharpest way to say what is and is not still open.  The
2394 representation/source problem has already been solved by the $2.5$PN quadrupole
2395 audit.  The local $4$PN existence problem has already been solved by the local
2396 referee chain.  The only possible remaining hereditary mismatch is one scalar
2397 monopole-scattering transport factor.
2406   \item the unique compatible conservative hereditary coefficient is the bridge
2407   coefficient \eqref{eq:4pn-tail-exact-bridge},
2408   \item therefore the $4$PN tail sector introduces no new independent
2409   quadrupole-normalization datum,
2410   \item and full GR matching of the hereditary coefficient is equivalent to the
2411   already-known $2.5$PN normalization target
2412   \begin{equation}
2413   \gamma_{\rm quad}^{\rm eff}
2414   =
2415   \frac{2G}{5c^5}.
2416   \label{eq:4pn-tail-gr-recovery-condition}
2417   \end{equation}
2662 \begin{equation}
2663 \Gamma_5
2664 =
2665 P_0\,\frac{a_{\rm th}^5}{27c_s^5},
2666 \qquad
2667 \gamma_{\rm quad}^{\rm eff}
2668 =
2669 \hat m_0^{\,2}\Gamma_5
2670 =
2671 \hat m_0^{\,2}P_0\,\frac{a_{\rm th}^5}{27c_s^5},
2672 \label{eq:interface-gamma-from-P0}
2673 \end{equation}
2674 with $\hat m_0$ the source-map normalization factor on the natural branch.
3420   \item In that Hamiltonian chart there is no local span obstruction.  The raw
3421   exchange-symmetric generic-frame scaffold is already complete against the full
3422   fixed-chart local interaction target, the canonical comparable-mass residual
3423   is exactly reachable, and the only ordinary-chart obstruction is an
3424   aligned-seed issue rather than a failure of the local generic-frame lift
3425   itself.
3484 L_{4\mathrm{PN}}^{\mathrm{tail}},
3485 \qquad
3486 C_{\rm tail}
3487 =
3488 \frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
3489 }
3490 \label{eq:conclusions-final-split-4pn}
3491 \end{equation}
3492 Full conservative GR recovery is therefore equivalent to the single scalar
3493 condition
3494 \begin{equation}
3495 \gamma_{\rm quad}^{\rm eff}=\frac{2G}{5c^5},
3496 \label{eq:conclusions-gamma-target-4pn}
3497 \end{equation}
3498 which immediately yields
3499 \begin{equation}
3500 C_{\rm tail}=\frac{G^2M}{5c^8}=C_{\rm tail}^{\rm GR}.
3501 \label{eq:conclusions-tail-gr-4pn}
3502 \end{equation}
3520 First, the paper does not solve the fully dynamical moving-throat PDE and then
3521 recover the $4$PN theorem as an automatic bulk consequence.  The present result
3522 remains a full conservative derivation only within the declared reduction and
3523 closure hierarchy already used by the lower-order papers.  That status is not a
3524 weakness of the theorem; it is the correct statement of its scope.
3525
3526 Second, the paper does not yet derive from first principles the passive/
3527 outgoing STF quadrupole normalization on the actual moving-throat branch.  The
3528 remaining universal target may be stated in either of the equivalent forms

research/pde_ledger_v3/V3_STEP_PLAN.md
747 ### S17 · The PN ladder — cite-only
748 Seven separately-published papers, Zenodo-DOI'd (`ledger_stage029_pn_corpus_doi_cite.md:32-40`).
749 ⛔ **Do not re-derive.** Record what each establishes and its inputs.
750
751 ⚠ **Two precisions the word "audited" hides — checked 2026-07-31:**
752 - **`4d_2_5pn` and `4d_4pn` are *conditional* derivations** (their own titles say so). ⛔ Do not carry
753   them forward as unconditional; the condition is part of the imported result.
754 - **The DOIs are README-authoritative, not source-cross-validated.** *"none of the seven papers declares
755   its own Zenodo DOI in source"* (`:51`). ⛔ Repeat that caveat rather than dressing them as
756   source-verified.
757 - stage029 itself is **cite-only**: no scripts, no dual-engine, no re-audit. *"029 does NOT re-audit the
758   PN physics."* ⇒ The audit lives in the papers, not in the ledger.
759
760 ### S18 · 2.5PN Burke–Thorne matchback
761 ⚠ **Record honestly**: `INV2` is `a`-free — `K̄₀ = (54/5)Gc_s⁵/(a⁵c⁵)` times `a⁵/(27c_s⁵)` is
762 `2G/(5c⁵)` identically. It is the **rational identity** `54/(5·27) = 2/5`, i.e. calibrated closure
763 consistency, ⛔ **not** a first-principles derivation of `G`.
1303 | `a_WT` (worldtube profile width) | controls the `O(a_WT^2/r^2)` correction |
1304 | the `a_WT` / `a_mouth` BRIDGE | DOES NOT EXIST (A11). Never close it by substitution |
1305 | `M` and the supplied multipoles | inputs of the S16 reduction, not outputs |
1306 | `G` / `N_inf,3` normalization, `I_F,12` | the force magnitude is calibrated through them |
1307 | `Theta_Q` and the `J -> Q` map | S15 must not substitute silently |
1308 | `kappa_add`, `kappa_PV`, `kappa_rho`, optical `n = 5` | the PN response packet S17 imports |
research/pde_ledger_v2/notes/stages/ledger_stage029_pn_corpus_doi_cite.md
19 **⚠⚠ Scope invariant (load-bearing honesty): Stage 029 does NOT re-derive, re-audit, or re-verify any PN result.** The
20 PN physics audit lives in the seven corpus papers. 029 is provenance only. The ledger's OWN gravity earning is Stages
21 001–028; 029 records the external audited corpus those stages calibrate-predict against.
50 The top-level `README.md` "Papers List" (16 DOI-bearing bullets) is the **sole authoritative source** for these DOIs.
51 A read of all seven `.tex` sources confirms that **none of the seven papers declares its own Zenodo DOI in source** —
52 the two 4D 1PN papers cite only the foundational Action DOI `19449589`, and the other five carry no Zenodo/DOI
53 citation. Consequently:
research/pde_ledger_v3/steps/S15_two_drain_region.md
3 **Result.** S15 supplies a test a two-drain flow must pass **if its signed net shell crossing is to scale as the inverse square of separation**. All three components still contain unspecified fields, even on the power-law profile family. The two engines agree exactly on the general crossing and its geometric departure. Their expressions for the additional field-response departure agree at 30 tested points per component. Which response an evolving pair of drains produces remains open. [Measurements M1–M3](_measurements/S15_record_measurements.md) · [Accepted comparator](../directives/_measurements/S15_reading_comparator_run_record.md)
4
5 The Newtonian target is the companion-induced acceleration of a freely moving mouth. A finite shell can store momentum while material crosses its boundary. S15 neither identifies how much reaches the core nor supplies the core's relation between momentum and motion. Whether inverse-square shell crossing is necessary or sufficient for Newtonian acceleration remains unresolved. S21 receives the conditional shell test. Throat acceleration, its direction and its normalization remain open. [Moving-region account](S16_moving_region.md) · [Requirements register](../SUBSTRATE_REQUIREMENTS.md)
47 **Item 4 remains deferred.** Wolfram reports 96 Reduce returns, 42 exponent-sensitive returns and six sufficient-only cases, with 192 general-field cases inapplicable. The six concern the two-drain-labelled region momentum and its rate, in each Cartesian component. Their printed sufficient conditions and corresponding restricted quantities survive as one-engine evidence. SymPy attempted no item-4 case in the accepted build. The stopped independent trial assigns no final classification, and the 144-object independent survey never started. The all-separation profile classification and the actual field responses remain open. [Per-object Wolfram returns and restrictions, M3–M4](_measurements/S15_record_measurements.md) · [Accepted reading build](../directives/_measurements/S15_reading_build_review_r5_disposition.md) · [Deferral](../directives/S15_reading_comparator_decision_list.md)
65 The first population route files the user's flow-comparison decision against the existing `R-S21-01`, alongside its distinct light-sector conditions. The comparison stays OPEN. If S21 tests inverse-square signed shell crossing on the stated family, it must carry the geometric departure, live field response, pressure and viscosity inputs, both surfaces and their motion. If it tests Newtonian gravity, it needs companion-induced **free acceleration** relative to an isolated drain under the same ambient conditions. Separate the surroundings from self-field (`R-S16-01`), account for the momentum stored in the shell and exchanged at the sleeve/core and at each slab face (`R-S16-03`), and identify the core contents (S22) and inertial response (`R-S16-04`, S8). The compactness assumption, where used, stays with `R-S16-02`. The physical suppliers of the passive surroundings field, face-resolved momentum and compactness have no owner assigned in S16. [Register entries and S15 pass](../SUBSTRATE_REQUIREMENTS.md) · [S16 record](S16_moving_region.md)
66
67 The slow limit and far behaviour remain conditional on `R-S14a-01`, jointly owed by S12, S8/S22 and the unnamed face/interior laws. S22 carries the core, the source map `{J, Θ_Q, N_∞,3, Q}` and the coefficient `I_F,12`. C3 remains the gap between drain rate and mass. The full sign retains quantum, confinement, Maxwell and profile residuals. The forms of `Π_o` and `f` have no assigned owner here. Extending steady P6 and resolving possible stress-carried momentum under P5 remain S8 material questions. A Newtonian orbit, its corrections, brane feedback, local thinning, and mouth support lie beyond this step. [Retained questions](../directives/S15_decision_list.md#retained-open) · [S22 debt rows](../V3_STEP_PLAN.md) · [Defect C3](../DEFECT_REGISTER.md)
research/pde_ledger_v3/steps/S16_moving_region.md
394 The plan's response-side approximation, supplied mass/multipoles, compactness, smooth
395 external field, calibrated potential, discarded boundary fluxes and reduced-force equation
396 are **not printed as an S16 reduction** here. The list implements the later user decisions
397 with live count, count-moment and boundary-flux objects. The list author's filing leaves
398 `M`, `a`, multipoles of `M`, compactness, source normalization and any flux suppression
399 as debts for S22's table.
400 `a_WT` remains distinct from `a_mouth`; neither control radius nor a count second moment
401 delivers their A11 bridge. S16 supplies no smooth surroundings field to expand, no
402 Gaussian profile width and no source normalization. (Retrieve M01/M04/M07.)
```
