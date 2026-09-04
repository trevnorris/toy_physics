---
title: Paper Moving Throat Pde
type: source
status: current
sources:
- research/pde/paper/pde.tex
last_updated: '2026-09-02'
---

## Purpose and scope

### source-paper-moving-throat-pde--framework-scope — Moving-throat response framework

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper presents a framework for treating a finite brane–bulk defect throat as a distributed response system. It connects a \((4+1)\)-dimensional matter-and-gauge parent theory to a lifted throat geometry, a linearized wall–support–gauge PDE, grouped real \(\ell=2\) reductions, and branch-facing quadrupole diagnostics. It does not claim a solved nonlinear moving-throat branch or prove that the desired passive/outgoing branch is physically realized.

Sources:

- `research/pde/paper/pde.tex` — `\label{sec:introduction}`
- `research/pde/paper/pde.tex` — `\label{sec:discussion-gap}`

## Source-unit map

- Entry point: `research/pde/paper/pde.tex`
- Role: `primary`
- Read mode: `semantic`
- Shape: `file`
- Source kind: `paper`

## Key statements

### source-paper-moving-throat-pde--exact-parent-system — Parent equations, identities, and mixed observables

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

Variation of the declared action gives the gauged GNLS and localized-Maxwell equations. The paper separately states exact number-current continuity, requires total-current conservation for gauge consistency, and gives the vorticity–gauge identity away from phase singularities. It identifies \(E_w=F_{w0}\) and \(C_a=F_{aw}\) as exact gauge-invariant parent-theory observables. These clauses include displayed identities and consistency requirements as well as variational equations, so the collection is not attributed wholly to action variation.

Sources:

- `research/pde/paper/pde.tex` — `\label{eq:exact-gnls}`
- `research/pde/paper/pde.tex` — `\label{eq:exact-continuity}`
- `research/pde/paper/pde.tex` — `\label{eq:exact-maxwell}`
- `research/pde/paper/pde.tex` — `\label{eq:current-conservation}`
- `research/pde/paper/pde.tex` — `\label{eq:vorticity-gauge}`
- `research/pde/paper/pde.tex` — `\label{eq:mixed-invariants}`

### source-paper-moving-throat-pde--projection-reduction-firewall — Exact projected continuity

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Once the normalized brane projection kernel is declared, projecting exact bulk continuity gives the open-system brane law
\[
\partial_t\rho_{\mathrm{brane}}+\nabla_3\!\cdot\mathbf j_{\mathrm{brane}}
=S_{\mathrm{leak}},
\]
with \(S_{\mathrm{leak}}\) equal to the stated boundary-and-weight leakage term.

Sources:

- `research/pde/paper/pde.tex` — `\label{eq:projected-continuity}`
- `research/pde/paper/pde.tex` — `\label{eq:leakage}`

### source-paper-moving-throat-pde--charge-notation-firewall — Charge, lane, profile, and reduction ontology

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper assigns electric-charge sign to \(\eta_Q\), not to circulation or geometry; reinterprets the historical bare \(q=1\) as the mass-dressing coefficient \(\kappa_\rho=1\); and defines \(20/21/22\) as grouped real \(P_2\) lanes rather than spacetime indices. It distinguishes the Maxwell localization profile \(Z(w)\) from the observation kernel \(W(w)\). Suppression of \(A_w\), \(J^w\), and \(F_{\mu w}\) is confined to a controlled far-field zero-mode reduction and does not erase those channels from the microscopic ontology.

Sources:

- `research/pde/paper/pde.tex` — `\label{sec:status-parent}`
- `research/pde/paper/pde.tex` — `\label{eq:zero-mode-assumptions}`
- `research/pde/paper/pde.tex` — `\label{eq:effective-maxwell}`

### source-paper-moving-throat-pde--geometry-lift — Distributed geometry and linear response problem

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The framework promotes the throat to the level set \(\Sigma(\mathbf X,t)=r-R(\Omega,w,t)\), retains the earlier radius and depth variables as collective moments, and uses real \(\ell=2\) harmonics as wall-displacement channels. After adopting a passive quadratic wall action as an effective closure, it constructs coupled linearized wall, BdG-matter, localized-Maxwell, and geometry equations and defines—but does not solve—the mouth/worldtube response operator.

Sources:

- `research/pde/paper/pde.tex` — `\label{eq:geometry-lift-levelset}`
- `research/pde/paper/pde.tex` — `\label{eq:geometry-lift-a-moment}`
- `research/pde/paper/pde.tex` — `\label{eq:geometry-lift-L-moment}`
- `research/pde/paper/pde.tex` — `\label{eq:geometry-lift-harmonic-decomp}`
- `research/pde/paper/pde.tex` — `\label{sec:linearized-pde}`
- `research/pde/paper/pde.tex` — `\label{eq:linpde-wall-action}`
- `research/pde/paper/pde.tex` — `\label{eq:linpde-bdg-matrix}`
- `research/pde/paper/pde.tex` — `\label{eq:linpde-maxwell}`
- `research/pde/paper/pde.tex` — `\label{eq:linpde-geometry}`
- `research/pde/paper/pde.tex` — `\label{eq:linpde-response-operator}`

### source-paper-moving-throat-pde--reduced-kernels — Conservative response moments

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Within a declared stable-mode truncation of the linearized hierarchy, Schur-complement elimination of matter-support and internal gauge coordinates produces
\[
D_A^{\mathrm{cons}}(\omega)
=D_{A0}+D_{A2}\omega^2+D_{A4}\omega^4+\mathcal O(\omega^6).
\]
The coefficients separate wall data from BdG moments \(B_{A,n}\) and localized-Maxwell/mixed moments \(Z_{A,n}\). The elimination is exact for the selected reduced mode list, not for the unsolved nonlinear PDE.

Sources:

- `research/pde/paper/pde.tex` — `\label{sec:reduced-system}`
- `research/pde/paper/pde.tex` — `\label{eq:reduced-p2-bdg-kernel}`
- `research/pde/paper/pde.tex` — `\label{eq:reduced-gauge-self-energy}`
- `research/pde/paper/pde.tex` — `\label{eq:reduced-full-conservative-kernel}`
- `research/pde/paper/pde.tex` — `\label{eq:reduced-conservative-moments}`

### source-paper-moving-throat-pde--isotropic-collapse — Isotropic collapse and weak-axisymmetric fingerprint

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For an \(O(3)\)-invariant reference throat and reduced kernels, the real \(20/21/22\) lanes share common conservative and outgoing-transfer coefficients, so their grouped anisotropy coordinates vanish. For a first-order pure axisymmetric quadrupole perturbation of that carrier, the lane signature is constrained to
\[
(\lambda_{20},\lambda_{21},\lambda_{22})=(1,\tfrac12,-1),
\qquad b=3a.
\]
The associated prefactor deformation is governed by the single dimensionless slope
\[
\Xi_1=\frac{P_1}{\bar P_0}
=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
\]
These results classify the reduced isotropic and weak-axisymmetric regimes; they neither establish that a realized branch is isotropic nor identify a unique microscopic source of anisotropy.

Sources:

- `research/pde/paper/pde.tex` — `\label{eq:grouped-STF-orthonormality}`
- `research/pde/paper/pde.tex` — `\label{eq:isotropic-collapse}`
- `research/pde/paper/pde.tex` — `\label{eq:isotropic-grouped-anomalies}`
- `research/pde/paper/pde.tex` — `\label{eq:weakax-lambda-signature}`
- `research/pde/paper/pde.tex` — `\label{eq:weakax-b-equals-3a}`
- `research/pde/paper/pde.tex` — `\label{eq:weakax-Xi1-def}`
- `research/pde/paper/pde.tex` — `\label{eq:weakax-Xi1-ratio}`

### source-paper-moving-throat-pde--outgoing-normalization — Isotropic normalization identity

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Inside the isotropic reduced bridge, the leading odd \(\omega^5\) coefficient depends on throat data through \(P_0=N_0/D_0\), the dimensionful scale carrier \(\hat m_0\), and the optional outgoing factor \(\chi_Q\). Matching the stated Burke–Thorne target gives
\[
\hat m_0^{\,2}\chi_Q P_0
=\frac{54Gc_s^5}{5a^5c^5}.
\]
The separate expression \(\hat m=1+\mathcal O(a^2/r^2)\) denotes a dimensionless profile and must not be conflated with \(\hat m_0\). On the canonical compact outgoing branch \(\chi_Q=1\). This is an algebraic gate within the reduced bridge, not evidence that an actual throat realizes the required \(P_0\).

Sources:

- `research/pde/paper/pde.tex` — `\label{eq:outgoing-Gamma5}`
- `research/pde/paper/pde.tex` — `\label{eq:outgoing-BT-target}`
- `research/pde/paper/pde.tex` — `\label{eq:outgoing-factorized-normalization}`
- `research/pde/paper/pde.tex` — `\label{eq:outgoing-natural-source-map}`
- `research/pde/paper/pde.tex` — `\label{eq:canonical-chiQ}`
- `research/pde/paper/pde.tex` — `\label{eq:outgoing-canonical-normalization}`

## Computed evidence represented by the source

### source-paper-moving-throat-pde--computed-evidence-boundary — Verification artifacts are not prepared evidence

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper describes symbolic audits for projector algebra, overlap models, Schur complements, normalization identities, and weak-axisymmetric transport. The executable files, exact invocations, literal outputs, and separate interpretation records are not supplied in the prepared unit. Therefore, the command/output/interpretation chain is **not supplied in the prepared unit**, and no statement in this capsule is classified as measured.

Sources:

- `research/pde/paper/pde.tex` — `\label{app:verification-map}`

## Assumptions, exclusions, and open questions

### source-paper-moving-throat-pde--closure-boundaries — Closure and reduction assumptions

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The distributed wall action, finite stable-mode selection, finite-throat support branch, low-frequency expansion, isotropic carrier, and weak-axisymmetric perturbation are closure or reduction choices. Conclusions depending on them apply only within their declared regimes and do not establish unique nonlinear defect dynamics. The paper characterizes the work as exploratory and calls for subject-matter review.

Sources:

- `research/pde/paper/pde.tex` — `\label{eq:linpde-wall-action}`
- `research/pde/paper/pde.tex` — `\label{sec:reduced-system}`
- `research/pde/paper/pde.tex` — `\label{eq:reduced-full-conservative-kernel}`
- `research/pde/paper/pde.tex` — `\label{eq:outgoing-actual-one-pole-branch}`
- `research/pde/paper/pde.tex` — `\label{eq:weakax-axisymmetric-perturbation}`
- `research/pde/paper/pde.tex` — heading `\section*{Acknowledgments}`

### source-paper-moving-throat-pde--branch-realization-gap — Physical realization remains open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The paper leaves unresolved whether an actual stationary moving-throat solution supplies admissible overlap data, grouped isotropy or a controlled weak-axisymmetric departure, positive static stability, the required outgoing normalization, and an acceptable value of \(\Xi_1\). Deriving the normalization identity and fingerprint does not show that a realized passive/outgoing branch provides the required \(P_0\) and anisotropic response. The proposed next step is branch construction and testing.

Sources:

- `research/pde/paper/pde.tex` — `\label{tab:claim-status}`
- `research/pde/paper/pde.tex` — `\label{eq:howto-branch-data}`
- `research/pde/paper/pde.tex` — `\label{eq:howto-stability}`
- `research/pde/paper/pde.tex` — `\label{eq:howto-main-normalization}`
- `research/pde/paper/pde.tex` — `\label{eq:howto-Xi1}`
- `research/pde/paper/pde.tex` — `\label{sec:discussion-gap}`

## Revision and supersession relationships

### source-paper-moving-throat-pde--no-declared-supersession — No whole-framework replacement declared

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper assigns a companion work the role of derivation and proof archive and describes specialized continuation packets as downstream of the present framework. It does not declare that either replaces this framework paper as a whole.

Sources:

- `research/pde/paper/pde.tex` — `\label{sec:discussion-gap}`

## Related topics and scripts

- `memory/sources/pde-audit.md`
- `memory/sources/software/stage1-solver.md`
- `memory/sources/software/stage1-verdict.md`