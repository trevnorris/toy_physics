---
title: Paper Unified 4D
type: source
status: current
sources:
- research/4d/paper/4d.tex
last_updated: '2026-08-25'
---

## Purpose and scope

### source-paper-unified-4d--purpose-action — Four-spatial-dimensional action framework

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

The paper specifies a \(4+1\)-dimensional action containing a confined gauged nonlinear Schrödinger matter field, a localized Maxwell sector, and two throat-geometry coordinates \(a(t)\) and \(L(t)\). From this declared action, variation yields the GNLS and localized Maxwell equations; global \(U(1)\) phase invariance yields the exact bulk continuity identity.

Sources:

- `research/4d/paper/4d.tex` — `\label{sec:action}`
- `research/4d/paper/4d.tex` — `\label{sec:eom-matter}`
- `research/4d/paper/4d.tex` — `\label{sec:eom-maxwell}`
- `research/4d/paper/4d.tex` — `\label{sec:eom-continuity}`

## Source-unit map

- Entry point and primary member: `research/4d/paper/4d.tex`
- Shape: single-file paper
- Read mode: `semantic`

## Key statements

### source-paper-unified-4d--projection-reduction-firewall — Projection is not dimensional reduction

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper provisionally declares a firewall between projection, which defines brane observables using a normalized measurement kernel \(W(w)\), and dimensional reduction, which integrates over \(w\) under an explicit ansatz to obtain effective couplings. These operations have different purposes and cannot be substituted for one another.

Sources:

- `research/4d/paper/4d.tex` — `\label{sec:braneobs}`
- `research/4d/paper/4d.tex` — `\label{sec:emreduction}`

### source-paper-unified-4d--charge-ontology — Fixed charge and localization-dependent coupling

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

For a unit-normalized charged-defect branch, the paper declares \(Q_{\rm tot}=q_\star=\eta_Q e_\star\). This total electric charge is independent of \(a\), \(L\), damping parameters \(\Gamma\), and continuous throat breathing. The canonically normalized brane coupling \(q_{\rm eff}=q_\star/\sqrt{Z_{\rm int}}\), by contrast, depends on electromagnetic localization, while circulation belongs to a separate vortical/fluxoid sector.

Sources:

- `research/4d/paper/4d.tex` — `\label{sec:fields}`
- `research/4d/paper/4d.tex` — `\label{eq:qeff-canonical}`
- `research/4d/paper/4d.tex` — `\label{sec:eom-winding}`

### source-paper-unified-4d--vorticity-fluxoid — Derived vortical and fluxoid identities

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Away from phase singularities, the Madelung relation gives \(\Omega_{ij}=-(q_\star/m)F_{ij}\). Fluxoid quantization is further scoped to a closed brane loop encircling the defect mouth in a region where \(\psi\) is single-valued; there the gauge-invariant phase winding is integral and yields the stated circulation–flux relation. Neither identity defines electric charge.

Sources:

- `research/4d/paper/4d.tex` — `\label{eq:vorticity-gauge}`
- `research/4d/paper/4d.tex` — `\label{eq:fluxoid-quantization}`
- `research/4d/paper/4d.tex` — `\label{eq:circulation-fluxoid}`

### source-paper-unified-4d--projected-continuity — Exact projected continuity and leakage

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Projecting the exact bulk continuity identity with \(W(w)\) yields an exact brane continuity equation with leakage \(S_{\rm leak}\), including both the asymptotic boundary term and the kernel-gradient term involving \(W'(w)j^w\). The simplified expression without the boundary term additionally assumes \(Wj^w\to0\) as \(|w|\to\infty\).

Sources:

- `research/4d/paper/4d.tex` — `\label{eq:proj-continuity}`
- `research/4d/paper/4d.tex` — `\label{eq:Sleak-general}`
- `research/4d/paper/4d.tex` — `\label{eq:boundary-vanish}`
- `research/4d/paper/4d.tex` — `\label{eq:Sleak-simplified}`

### source-paper-unified-4d--poisson-regime — Exact longitudinal identity and conditional Poisson limit

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

A Helmholtz decomposition of the projected brane velocity and the projected continuity equation produce the exact longitudinal identity for \(\varphi\). Quasi-static density, a suppressed density-gradient/advection correction, and approximately constant \(\rho_{\rm brane}\) are sufficient regime conditions for the stated Poisson approximation; longitudinal dominance is optional but common, not necessary. The \(1/r^2\) longitudinal-velocity scaling additionally requires a localized source and suitable three-dimensional Green-function boundary behavior, and its identification with gravitational acceleration would require a separate constitutive step.

Sources:

- `research/4d/paper/4d.tex` — `\label{eq:longitudinal-identity}`
- `research/4d/paper/4d.tex` — `\label{sec:poisson-regime}`
- `research/4d/paper/4d.tex` — `\label{eq:poisson-approx}`
- `research/4d/paper/4d.tex` — `\label{sec:poisson-invsq}`

### source-paper-unified-4d--maxwell-zero-mode — Controlled Maxwell zero-mode reduction

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Under \(A_w=0\), \(\partial_wA_\mu=0\), and the presentation’s Lorenz gauge, integrating the localized Maxwell equation gives \(\partial_\mu F^{\mu\nu}=(\mu_0/Z_{\rm int})J_{\rm eff}^{\nu}\). The \(N=w\) equation also requires the total current \(J_\psi^w+J_{\rm ext}^w\) to be negligible or suppressed under the pure zero-mode ansatz. A compatible \(\nu\)-source profile proportional to \(Z(w)\) is an additional pointwise consistency condition; mismatched profiles excite higher modes.

Sources:

- `research/4d/paper/4d.tex` — `\label{eq:zero-mode}`
- `research/4d/paper/4d.tex` — `\label{sec:emreduction-ansatz}`
- `research/4d/paper/4d.tex` — `\label{eq:brane-maxwell}`
- `research/4d/paper/4d.tex` — `\label{eq:source_profile_matching}`

### source-paper-unified-4d--geometry-force-ledger — Derived generalized-force ledger

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For the declared conservative energy, the geometry forces are defined by \(F_a=-\partial_aH_{\rm tot}\) and \(F_L=-\partial_LH_{\rm tot}\). The paper derives the confinement-loaded matter contributions and the optional ratio-penalty forces within that ledger.

Sources:

- `research/4d/paper/4d.tex` — `\label{eq:ledger-forces}`
- `research/4d/paper/4d.tex` — `\label{eq:Fa-matter}`
- `research/4d/paper/4d.tex` — `\label{eq:FL-matter}`
- `research/4d/paper/4d.tex` — `\label{eq:ratio-forces}`

### source-paper-unified-4d--geometry-closure — Adopted phenomenological geometry dynamics

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The curvature-completed geometry energy, effective inertias, damping coefficients, optional ratio penalty, and damped two-coordinate evolution laws are adopted phenomenological choices rather than a microscopic wall derivation. Overdamped relaxation and instantaneous minimization are controlled limits of this baseline closure.

Sources:

- `research/4d/paper/4d.tex` — `\label{sec:geometry-Egeom}`
- `research/4d/paper/4d.tex` — `\label{eq:Egeom}`
- `research/4d/paper/4d.tex` — `\label{sec:geometry-dynamics}`
- `research/4d/paper/4d.tex` — `\label{eq:geom-eom-again}`

## Computed evidence represented by the source

### source-paper-unified-4d--symbolic-evidence-boundary — Symbolic checks lack a sealed measurement chain

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper names two Mathematica commands and reports that their harnesses check equation forms and internal symbolic consistency, but the executable sources, literal outputs with operands and residuals, and a separate interpretation record are not members of this prepared unit: not supplied in the prepared unit. The reported checks therefore do not establish external physical validity or validate the Poisson regime for a concrete scenario.

Sources:

- `research/4d/paper/4d.tex` — `\label{app:wl-repro}`
- `research/4d/paper/4d.tex` — `\label{app:wl-interpret}`

## Assumptions, exclusions, and open questions

### source-paper-unified-4d--open-validation — Exploratory status and outstanding validation

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The proposal is explicitly exploratory. Symbolic checks address internal consistency only, and external subject-matter validation remains outstanding. Open work includes choosing and testing confinement and projection profiles, quantifying the Poisson-regime corrections, incorporating higher and mixed electromagnetic modes, deriving effective mechanical parameters microscopically, and developing richer or multi-throat geometry closures.

Sources:

- `research/4d/paper/4d.tex` — heading `\section*{Acknowledgments}`
- `research/4d/paper/4d.tex` — `\label{sec:discussion-limitations}`
- `research/4d/paper/4d.tex` — `\label{sec:discussion-outlook}`

## Revision and supersession relationships

### source-paper-unified-4d--revision-scope — Corrections and extensions without formal supersession

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper describes a corrected charge/current ontology, extends the geometry energy with a curvature/bending term, and connects dynamic geometry to quasi-static limits. These are scoped corrections and extensions within the paper; the prepared unit does not declare formal supersession of another source or capsule.

Sources:

- `research/4d/paper/4d.tex` — `\label{sec:action-em}`
- `research/4d/paper/4d.tex` — `\label{sec:discussion-explains}`
- `research/4d/paper/4d.tex` — `\label{sec:geometry-Egeom}`
- `research/4d/paper/4d.tex` — `\label{sec:geometry-dynamics}`

## Related topics and scripts

No related memory pages or script domains were supplied for this task.
