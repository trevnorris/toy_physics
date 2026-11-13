# S17 decision-list source retrieval

The command reads the README, the seven main manuscripts, the step plan, the earlier citation note and the accepted S15 record from commit `438d270c`. Run from the repository root. The numbered lines are excerpts from that revision. The DOI count is a literal search of each paper’s main TeX source for its own assigned DOI.

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
    ("1pn_orbital_dynamics", "19449058", [(38, 40), (376, 396)]),
    ("4d_1pn_bridge", "19449653", [(91, 91)]),
    ("4d_1pn_full", "19450102", [(97, 97), (697, 719), (818, 820), (825, 837), (939, 946), (3094, 3107)]),
    ("4d_2pn", "19450284", [(101, 118)]),
    ("4d_2_5pn", "19492270", [(35, 60)]),
    ("4d_3pn", "19501724", [(34, 48)]),
    ("4d_4pn", "19561056", [(43, 43)]),
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
excerpt("research/pde_ledger_v3/V3_STEP_PLAN.md", [(747, 763), (1308, 1308)])
excerpt("research/pde_ledger_v2/notes/stages/ledger_stage029_pn_corpus_doi_cite.md", [(19, 21), (50, 53)])
excerpt("research/pde_ledger_v3/steps/S15_two_drain_region.md", [(5, 5), (47, 47), (65, 65)])
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

README.md 145 - *Deriving Key Post-Newtonian Coefficients from the Unified 4D Toy Model.* [10.5281/zenodo.19449653](https://doi.org/10.5281/zenodo.19449653)
Assigned DOI occurrences in main TeX: 0
research/4d_1pn_bridge/paper/4d_1pn_bridge.tex
91 This paper is a derivation supplement to the unified $4$D throat model of Ref.~\cite{norris2026_4d_model}, written to make explicit which coefficients used earlier in the series are fixed by the $4$D structure (and which must remain symbolic pending additional closure). We show that a defect uniform in the extra coordinate behaves as a throat with hypercylindrical topology, for which the exterior potential--flow problem yields the classical added-mass coefficient $\kappa_{\mathrm{add}}=\tfrac12$ (whereas a $4$D bubble would instead give $\tfrac13$), providing a sharp topology discriminator. We then show that matching the weak-field refractive coefficient uniquely selects the stiff polytropic exponent $n=5$. In the $1$PN interaction sector, EIH matching produces a one-parameter family of couplings, which collapses uniquely to $\alpha^2=\tfrac34$ and $K_{\rm vec}=2/\pi^2$ once the $n=5$ thermodynamic relation is imposed. Treating rest mass as wave-supported energy reproduces the special-relativistic expansion $E=\gamma E_0$, including the universal $v^4$ coefficient $3/8$. We also summarize a controlled reduction to brane-effective equations, making leakage/source terms explicit, and we outline the additional response modeling required to determine the remaining undetermined coefficients.

README.md 148 - *A Controlled Derivation of the Full First Post-Newtonian Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19450102](https://doi.org/10.5281/zenodo.19450102)
Assigned DOI occurrences in main TeX: 0
research/4d_1pn_full/paper/4d_1pn_full.tex
97 We present a controlled derivation of the full conservative two-body first post-Newtonian (1PN) sector from the four-spatial-dimensional brane--bulk toy model. Starting from the exact bulk action, projection map, open-system balance laws, and quasi-static near-zone Poisson hook on the brane, we use a coherent-defect / small-body worldtube reduction to obtain the Newtonian point-particle limit and then assemble the conservative 1PN ledger sector by sector. Newtonian matching fixes the scalar mass-dressing coefficient \(\krho=1\) (historically written \(q=1\) in earlier gravity-side bookkeeping), weak-field optical consistency fixes the polytropic exponent \(n=5\), and wave-supported rest energy yields the universal free-particle \(v^4\) coefficient. Local mass scaling plus pair counting generate the static \(G^2/r^2\) term, exterior potential-flow topology fixes the added-mass coefficient \(\kappa_{\rm add}=1/2\) for a \(w\)-uniform throat, and a constructive isotropic Fourier-space wake basis reproduces the Einstein--Infeld--Hoffmann coefficients \(C_\parallel=-7/2\) and \(C_L=-1/2\) with \(\alpha^2=3/4\), \(a_H=0\), and \(K_{\rm vec}=2/\pi^2\). Under an explicit adiabatic one-degree-of-freedom breathing closure we further obtain \(\kappa_{\rm PV}=3/2\), so that \(\beta_{1\rm PN}=\kappa_\rho+\kappa_{\rm add}+\kappa_{\rm PV}=3\). The assembled reduced two-body Lagrangian is then algebraically identical to the standard conservative EIH form, and the test-mass limit reproduces the standard perihelion shift. The result is a full conservative 1PN derivation within a declared closure hierarchy that makes explicit which ingredients are fixed, protocol-fixed, and still symbolic, with remaining freedom confined to generalized internal response, finite-size structure, and beyond-1PN open-system effects.
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
837
939 Equation \eqref{eq:newtonian-two-body-lagrangian} is the baseline pair Lagrangian on top
940 of which the subsequent \(1\)PN self, cross, and static nonlinear sectors will be added.
941
942 \claimclass{\Controlled for the identification of the Poisson-hook scalar with the
943 effective Newtonian potential, the coherent-defect / small-body worldtube closure, and
944 the neglect of boundary-flux terms at the retained order. Exact for the Euler--Lagrange
945 equivalence between \cref{eq:newtonian-particle-lagrangian} and the leading monopole
946 equation of motion once the reduced Lagrangian has been stated.}
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

README.md 152 - *A Conditional Full Conservative Derivation of the Two-Body 4PN Sector from the Unified 4D Toy Model.* [10.5281/zenodo.19561056](https://doi.org/10.5281/zenodo.19561056)
Assigned DOI occurrences in main TeX: 0
research/4d_4pn/paper/4d_4pn.tex
43 We present the full conservative fourth post-Newtonian (4PN) assembly for the 4D toy-model program within the same declared closure hierarchy used at lower orders. The strict Schwarzschild one-body gate, exact quartic Legendre compiler, Hamiltonian-first generic-frame lift, canonical comparable-mass local residual, aligned-seed correction, and referee reconstruction together show that the entire local instantaneous 4PN sector already has an exact generic-frame ordinary representative. We then isolate the hereditary sector and show that its coefficient is not a new independent datum, but is fixed by the same canonically normalized orbital/worldtube STF quadrupole branch that controls the 2.5PN Burke--Thorne channel, through $C_{\rm tail}=(GM/2c^3)\,\gamma_{\rm quad}^{\rm eff}$. Hence 4PN opens no new normalization constant: if $\gamma_{\rm quad}^{\rm eff}=2G/(5c^5)$, then $C_{\rm tail}=G^2M/(5c^8)$ and the full conservative 4PN local+tail structure matches the standard GR target. The remaining gap is therefore the same passive/outgoing quadrupole-normalization theorem on the moving-throat branch already isolated by the 2.5PN program.

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
5 The Newtonian target is the companion-induced acceleration of a freely moving mouth. A finite shell can store momentum while material crosses its boundary. S15 neither identifies how much reaches the core nor supplies the core's relation between momentum and motion. Whether inverse-square shell crossing is necessary or sufficient for Newtonian acceleration remains unresolved. S21 receives the conditional shell test. Throat acceleration, its direction and its normalization remain open. [Moving-region account](S16_moving_region.md) · [Requirements register](../SUBSTRATE_REQUIREMENTS.md)
47 **Item 4 remains deferred.** Wolfram reports 96 Reduce returns, 42 exponent-sensitive returns and six sufficient-only cases, with 192 general-field cases inapplicable. The six concern the two-drain-labelled region momentum and its rate, in each Cartesian component. Their printed sufficient conditions and corresponding restricted quantities survive as one-engine evidence. SymPy attempted no item-4 case in the accepted build. The stopped independent trial assigns no final classification, and the 144-object independent survey never started. The all-separation profile classification and the actual field responses remain open. [Per-object Wolfram returns and restrictions, M3–M4](_measurements/S15_record_measurements.md) · [Accepted reading build](../directives/_measurements/S15_reading_build_review_r5_disposition.md) · [Deferral](../directives/S15_reading_comparator_decision_list.md)
65 The first population route files the user's flow-comparison decision against the existing `R-S21-01`, alongside its distinct light-sector conditions. The comparison stays OPEN. If S21 tests inverse-square signed shell crossing on the stated family, it must carry the geometric departure, live field response, pressure and viscosity inputs, both surfaces and their motion. If it tests Newtonian gravity, it needs companion-induced **free acceleration** relative to an isolated drain under the same ambient conditions. Separate the surroundings from self-field (`R-S16-01`), account for the momentum stored in the shell and exchanged at the sleeve/core and at each slab face (`R-S16-03`), and identify the core contents (S22) and inertial response (`R-S16-04`, S8). The compactness assumption, where used, stays with `R-S16-02`. The physical suppliers of the passive surroundings field, face-resolved momentum and compactness have no owner assigned in S16. [Register entries and S15 pass](../SUBSTRATE_REQUIREMENTS.md) · [S16 record](S16_moving_region.md)
```
