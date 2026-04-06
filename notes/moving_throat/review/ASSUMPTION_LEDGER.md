# Moving-Throat PDE — Assumption Ledger

## Scope

This ledger records the working hypotheses that materially matter for the
present proof chain after the full review and the complete SymPy/Mathematica
audit buildout.

This is the first hardening tranche only. It covers the choke-point stages
identified in
`notes/moving_throat/review/PROOF_HARDENING_PLAN.md`:

- `003-004`
- `020-031`
- `041`
- `089`
- `108`
- `125-140`
- `153`
- `168-170`

The goal here is not to restate every algebraic assumption appearing in every
CAS file. The goal is to make the mathematically relevant hypotheses explicit:

- positivity or stability conditions,
- nonvanishing denominators and selected minors,
- branch/root choices,
- perturbative or asymptotic limits,
- and places where the statement is knowingly more local than global.

## Conventions

- `Explicit` means stated in the note or theorem body.
- `CAS` means used directly by the SymPy/Mathematica audits.
- `Branch` means a chosen root, positivity branch, or stable-side restriction.
- `Open edge` means a regime the present stage does not settle.

## Stage `003`

- `Explicit`: stable BdG normal-mode reduction on an isotropic reference throat.
- `CAS`: positive/internal stable BdG frequencies and regular inverse
  `(\Omega_0^2 - \omega^2 I)^{-1}` away from resonance.
- `Branch`: weak-coupling perturbative pole shift is taken only on the separated
  branch `varpi^2 > Omega_eta^2`.
- `Open edge`: exact resonance and full non-isotropic grouped closure are not
  handled here.

## Stage `004`

- `Explicit`: reduced linear/quadratic wall + Maxwell + mixed model with
  outgoing `l=2` dressing.
- `CAS`: positive reduced parameters such as `M`, `K`, `Omega_A`, `Omega_W`,
  `a`, `c_s`, `gamma_1`; conservative stability effectively requires
  `Omega_A^2 Omega_W^2 - R^2 > 0`.
- `Branch`: outgoing response is used only to first order in `Pi_out`; the
  active quadrupole port is passive/outgoing.
- `Branch`: the scalar rescue statement assumes derivative-only coupling into
  the port-active block, with no direct non-derivative wall coupling there.
- `Open edge`: full grouped-`P2` Maxwell channelization and multimode outgoing
  coupling are deferred.

## Stage `020`

- `Explicit`: first finite-throat linear continuum operator on `s in [0,L]`
  with N/N wall-`U` modes and D/N half-wave `W/phi` modes.
- `CAS`: `s, L > 0`; positive masses and stiffnesses; real couplings.
- `CAS`: conservative stability conditions
  `K_U K_eta^eff > c_etaU^2` and `K_U K_W^eff > c_UW^2 sigma`.
- `Branch`: static/conservative branch only; no outgoing dressing yet.
- `Open edge`: this stage gives explicit continuum data, but not a theorem that
  the completed PDE necessarily lands in the admissible/stable region.

## Stage `021`

- `Explicit`: uses the dimensionless map built from Stage `020`.
- `CAS`: `eps_eta, eps_W, rho, Z_W, delta0, Lambda > 0`.
- `Branch`: physical/stable placement uses `0 < eps_eta, eps_W < 1`.
- `Branch`: the natural transfer branch is the nonvanishing one `1 + rho > 0`.
- `Open edge`: the map is exact only after the actual PDE kernel ratios are
  supplied; this stage does not derive them.

## Stage `022`

- `Explicit`: first axial `U` splitting of the Stage-`020` operator.
- `CAS`: positive effective stiffnesses and masses; `delta0, deltaU > 0`;
  real couplings; nonzero split parameters where denominators require it.
- `Branch`: flat limit `delta_U = 0` is treated as a special degeneration.
- `Branch`: directional conclusions are taken on the constructive branch
  `rho0 > 0`.
- `Open edge`: source/loading collinearity is lost for generic
  `delta_U != 0`, so later normalization stages must admit noncollinear
  geometry.

## Stage `023`

- `Explicit`: generalized rank-1 selected branch parameterized by the loading
  ratio `q = z1 / z0`.
- `CAS`: `A0, delta, alpha, z0, xi, R_U, eps > 0`, `q` real and nonzero,
  `eta` real.
- `Branch`: stable selected branch is restricted to `0 <= xi < 1`.
- `Branch`: Stage-`018/019` flat-`U` formulas are recovered only on the
  `R_U = 1` limit.
- `Open edge`: the actual support direction is still unknown here.

## Stage `024`

- `Explicit`: two rank-1 loading directions `z` and `y`.
- `CAS`: `A0, delta, xi > 0`, `R_U > 0`, and real `m, n, q, r`.
- `Branch`: stable branch still uses `0 <= xi < 1`.
- `Branch`: the support theorem requires the denominator
  `delta + (1 + r^2) xi - m (q - r)^2` to stay nonzero.
- `Open edge`: the physically relevant alignment of `y` remains unresolved
  (`y || z` versus `y || v`).

## Stage `025`

- `Explicit`: general three-direction geometry `(z, y, s)` inserted into the
  selected-mode normalization law.
- `CAS`: `A0, delta, xi > 0`, `R_U, eps > 0`, and real `m, q, r, t`.
- `Branch`: sign conclusions are interpreted on the constructive split-`U`
  branch `R_U < 1`.
- `Branch`: tracking and source-tied formulas are specializations of the
  general normalization law, not generic equalities.
- `Open edge`: this stage still waits on the actual continuum support direction.

## Stage `026`

- `Explicit`: first symmetry-allowed `U/phi` coupling while keeping a rank-1
  support lane after elimination of the split `U` doublet.
- `CAS`: positive `K_U`, `K_eta^eff`, `K_phi^eff`, `delta_U`, `mu_eta`,
  `mu_phi`, `sigma`; real/nonzero couplings where required.
- `Branch`: source-tied closure occurs only in the minimal-kernel limit
  `sigma0 = 0` or unsplit limit `delta_U = 0`.
- `Branch`: exact tracking occurs only on the codimension-one interference
  match `sigma0 = rho0`.
- `Open edge`: this stage extracts the actual support lane but does not yet
  close the physical selected branch.

## Stage `027`

- `Explicit`: actual continuum inputs `(M_mix, R_U, M_supp, R_phi)` are inserted
  into the rank-2 selected-branch formulas.
- `CAS`: `xi, delta > 0`, nonnegative `M_mix, M_supp`, and `R_U, R_phi > 0`.
- `Branch`: the physical quadratic root is the `+` branch selected by
  continuity to `xi = 0` at zero load.
- `Branch`: source-tied and tracking are recognized as exact special surfaces
  `R_phi = 1` and `R_phi = R_U`.
- `Open edge`: without a concrete kernel, this stage does not decide whether
  the physical PDE lands on source-tied, tracking, or intermediate closure.

## Stage `028`

- `Explicit`: coherent local-kernel hypothesis that `W` and `phi` couple to the
  same local wall/`U` density.
- `CAS`: positive `lambda_W`, `lambda_phi`, `gamma`, masses, `K_U`, `chi_0`,
  `delta_U`, `Z_W`, `Z_phi`; later formulas also assume `xi, delta, lambda0 > 0`.
- `Branch`: constructive coherent branch uses `chi_0 > 0`, `delta_U > 0`.
- `Branch`: the key structural specialization is exact tracking
  `rho0 = sigma0`, hence `R_tr = R_U = R_phi`.
- `Open edge`: the note’s exact normalization-law collapse is now supported by
  Mathematica, but the physical kernel hypothesis itself is still an input.

## Stage `029`

- `Explicit`: physical tracking branch is compared against the flat branch.
- `CAS`: `xi, delta, R > 0`.
- `Branch`: the physically relevant range is `0 < R < 1`, `0 < xi < 1`,
  `delta > 0`; the note also records bounds for `0 <= R <= 1`.
- `Branch`: monotonicity conclusions depend on positivity of the factored
  polynomials in this domain.
- `Open edge`: coefficient positivity is clear algebraically but only partly
  asserted programmatically in the original SymPy check.

## Stage `030`

- `Explicit`: coherent local-kernel tracking branch with explicit support
  amplitude ratio `zeta`.
- `CAS`: positive `Keta_eff`, `KU`, `KW_eff`, `Kphi_eff`, `lamW`, `lamphi`,
  `gamma`, `c_etaU`, `Tw`, `L`, `TU`, `G`, `cs`, `a`, `c`, `muW`, `zeta`.
- `Branch`: coherent support formulas are used on `0 < eps < 1`,
  `0 <= zeta < 1 / eps`.
- `Branch`: `R_target` independence of `zeta` belongs to the coherent tracking
  branch, not the unrestricted rank-2 problem.
- `Open edge`: the stage defines the enhancement law, but not the actual
  physical `zeta`.

## Stage `031`

- `Explicit`: support compensation theorem on the coherent tracking branch.
- `CAS`: `xi, delta, R > 0`, `eps, zeta > 0`, `M_mix > 0`, `S_req, S_crit > 0`,
  `nu > 0`.
- `Branch`: theorem domain is `0 < eps < 1`, `0 <= zeta < 1 / eps`,
  `0 < xi < 1`, and mixed-only baseline below the critical load.
- `Branch`: the theorem is stable-side only; `zeta_req` is defined before
  softening.
- `Open edge`: no reduced-level obstruction remains, but the true PDE may still
  fail physically if its realized `zeta` never reaches `zeta_req`.

## Stage `041`

- `Explicit`: coupled support/source drift-diffusion fixed point for `Pe`.
- `CAS`: `Pe, alpha, eta, Xi > 0`.
- `Branch`: constructive branch uses `Pe >= 0` and stationary zero-flux
  reduction.
- `Branch`: weak-coupling asymptotic `Pe_* = Xi Delta_0 + O(Xi^2)` is only a
  small-`Xi` statement.
- `Open edge`: the stage proves existence and bracketing, not a unique
  closed-form root for all `Xi`.

## Stage `089`

- `Explicit`: canonical outgoing reduced closure for the quadrupole channel.
- `CAS`: positive `G`, `c`, `c_s`, `a`, `m0hat`, `chi_Q`, `N_Q`.
- `Branch`: canonical outgoing DtN branch is `chi_Q = 1`.
- `Branch`: natural source-map branch uses `mhat_0 = 1 + O(a^2 / r^2)`.
- `Branch`: GR-target matching is taken in the strict point-particle limit.
- `Open edge`: the reduction is conditional on the PDE actually realizing the
  canonical outgoing branch.

## Stage `108`

- `Explicit`: positive normalized localized mouth source on the first D/N
  interval.
- `CAS`: `z, L > 0`, `x` real, together with positivity/normalization of the
  source density.
- `Branch`: shell and mixed channels are driven by the same positive density, so
  the first cosine moment satisfies `g[sigma] in [0,1]`.
- `Branch`: the theorem excludes the upper Family-1 root `g_+^F1 > 1` and keeps
  only the lower positive-source branch.
- `Open edge`: this is still a theorem inside the positive-source setup; it does
  not rederive the full branch-value provenance locally.

## Stage `125`

- `Explicit`: explicit exponential positive mouth layer `g_c = g_Pi` and
  Family-1 self-matched closure `Sigma0 = (20/9) Tmhat^2`.
- `Branch`: Family-1 core gain law is fixed as
  `R_q(Pi) = (g_Pi - r_F1)^2 / (1 + r_F1^2)`.
- `Branch`: the self-consistent branch is the finite positive solution that
  yields the canonical point `(Pi_*, Tmhat_*)`.
- `Open edge`: the singular locus `1 - R_q S_q = 0` is not globally analyzed in
  this stage.

## Stage `126`

- `Explicit`: same explicit exponential mouth family as Stage `125`.
- `CAS`: `Pi > 0`.
- `Branch`: finite positive `Pi` always satisfies `0 < g_Pi < 1`.
- `Branch`: the equal-normalized branch `g_c = 1` occurs only in the singular
  limit `Pi -> infinity`, with `Tmhat ~ 0.72567 sqrt(Pi)`.
- `Open edge`: exclusion is only inside the explicit exponential mouth family.

## Stage `127`

- `Explicit`: positive-source theorem plus the self-consistent exponential
  mouth/core map.
- `Branch`: monotone exponential positive mouth family and positive-source bound
  are assumed.
- `Branch`: upper branch is excluded; equal-normalized branch is singular; the
  unique regular finite branch is the lower compensated branch.
- `Open edge`: uniqueness is proved only inside the explicit Family-1
  positive-mouth closure, not the full PDE.

## Stage `128`

- `Explicit`: status/consolidation stage for the mouth branch.
- `Branch`: carries forward the explicit exponential family, upper-branch
  exclusion, and equal-normalized singularity from Stages `113`, `125`, `127`.
- `Open edge`: still assumes the real mouth layer is close enough to the
  explicit exponential family to make this reduction meaningful.

## Stage `129`

- `Explicit`: small positive convex deformation
  `Sigma_eps = (1 - eps) Sigma_* + eps varpi`.
- `CAS`: `Pi > 0`, `x` real, `gbar, Sbar, eps` real.
- `Branch`: first-order regime only, `0 <= eps << 1`.
- `Branch`: the deformation is retuned to stay on the canonical lower
  compensated branch `g = g_*`.
- `Open edge`: all formulas are perturbative to first order; no finite-`eps`
  theorem is claimed.

## Stage `130`

- `Explicit`: first-order traction-rigidity kernel on the compensated branch.
- `CAS`: `Pi > 0`, `x, eps, lam` real.
- `Branch`: lower compensated branch structure is frozen at `R_q = 1/4`.
- `Branch`: the traction law is a first-order consequence of the Stage-`129`
  retuning formulas, not a global exact law.
- `Open edge`: higher-order mouth corrections are not controlled here.

## Stage `131`

- `Explicit`: two explicit positive non-exponential mouth families plus affine
  interpolation.
- `CAS`: `Pi > 0`, `x, eps, lam` real.
- `Branch`: evaluations use the first-order rigidity kernel from Stage `130`.
- `Branch`: the “zero-shift” interpolation claim is comparative and approximate,
  not a new exact theorem.
- `Open edge`: stage-level issues were cosmetic/numeric; the mathematical result
  remains first-order and family-restricted.

## Stage `132`

- `Explicit`: status stage for mouth rigidity.
- `Branch`: simply inherits the `129-131` deformation/rigidity regime.
- `Open edge`: no stage-specific audit; the conclusion remains a summary of the
  explicit Family-1 mouth closure rather than a general PDE theorem.

## Stage `133`

- `Explicit`: explicit Family-1 branch with `kappa_s = 0`, `kappa_q = pi/2`,
  canonical exponential source, and compensated outlet
  `M_s^* = 4 Sigma_m^*`, `M_q^* = -Sigma_m^*`.
- `CAS`: `x, Pi, Sigma > 0`.
- `Branch`: everything is evaluated on the canonical compensated branch.
- `Branch`: the residual theorem is exact for the explicit branch, but it is
  still local to that branch data.
- `Open edge`: the stage proves tangent matching and broadening of the full
  profile, not the final branch retuning.

## Stage `134`

- `Explicit`: linearization of the actual full source around the canonical
  branch.
- `Branch`: assumes
  `Sigma_act = Sigma_* (1 - Rtilde_*) + O(R_*^2)`.
- `Branch`: uses frozen canonical derivatives and the Stage-`130` rigidity
  coefficients.
- `Open edge`: the actual numerical covariances still have to be computed; this
  stage is structural.

## Stage `135`

- `Explicit`: canonical residual from Stage `133` is evaluated numerically and
  fed through the Stage-`134` first-order correction law.
- `Branch`: still a mouth-only correction with the core outlet frozen.
- `Branch`: one-step nonlinear Picard update is used as a consistency check,
  not as a full convergence theorem.
- `Open edge`: this stage does not yet include core-mouth co-evolution.

## Stage `136`

- `Explicit`: status stage consolidating the mouth correction.
- `Branch`: inherits the mouth-only correction regime from Stages `133-135`.
- `Open edge`: the next unresolved variable is full core-mouth co-evolution.

## Stage `137`

- `Explicit`: any positive normalized source `Sigma`, exact Family-1 core law,
  self-consistent Boltzmann fixed point, and self-matched closure
  `Sigma0 = (20/9) Tmhat^2`.
- `CAS`: `g, r, dg` real; `Sigma0, dSigma0, Sstar, dS, Rstar, dR` real.
- `Branch`: canonical compensation is the branch condition `g = g_* iff R = 1/4`.
- `Branch`: defect transport `delta R = -delta g / sqrt(1 + r_F1^2) + O(delta g^2)`
  is a first-order statement around the canonical branch.
- `Open edge`: the exact traction needed to land back on `g_*` is not solved
  here.

## Stage `138`

- `Explicit`: exact co-evolving fixed-point map with traction frozen at the old
  canonical value `Sigma0 = Sigma0^*`.
- `Branch`: positive fixed point is found from the canonical seed on the fixed
  traction slice.
- `Branch`: the branch survives but is no longer exactly compensated.
- `Open edge`: uniqueness is supported by the solved branch and audits, but the
  note does not prove global uniqueness from arbitrary seeds.

## Stage `139`

- `Explicit`: exact co-evolving condition `g_fp(Sigma0) = g_*`.
- `Branch`: assumes monotonicity of `g_fp(Sigma0)` on the numerically solved
  interval used for the root-finding.
- `Branch`: the renormalized compensated branch is obtained numerically on that
  bracket.
- `Open edge`: the theorem is exact as a fixed-point condition, but the realized
  branch value is numerical rather than closed-form.

## Stage `140`

- `Explicit`: status stage for the full co-evolving fixed-point picture.
- `Branch`: inherits branch elimination from `108`, `126`, `127`, and the
  co-evolving fixed-point results from `138-139`.
- `Open edge`: summary numerics must continue to match the cited stage values
  exactly; the derivation itself is carried by prior stages, not proved anew
  here.

## Stage `153`

- `Explicit`: isotropic compensated canonical branch with pure linear grouped
  real `P2` anisotropy.
- `CAS`: nonzero real reduced coefficients such as `D0`, `N0`, `P0`, `sigma`,
  and real grouped perturbations.
- `Branch`: relies on Stage `152` constraint `delta C^(1,P2) = 0`.
- `Branch`: weak-axisymmetric specialization later uses the pattern
  `(1, 1/2, -1)` after the grouped derivation is done.
- `Open edge`: this gives the linear grouped outlet map, not the actual
  anisotropies realized on the physical branch.

## Stage `168`

- `Explicit`: coherent local D/N tracking branch at first grouped
  weak-axisymmetric/reference-branch order, freezing
  `chi0,*`, `deltaU,*`, `E_*`, `F_*`.
- `CAS`: `chi0s, deltaUs, epsWs, epss > 0`; microscopic reference quantities
  positive; remaining drift variables real.
- `Branch`: zero defect is identified with preservation of three direct
  microscopic monomials on the reference branch.
- `Branch`: this is a linearized rigidity statement, not yet a theorem for the
  full nonlinear microscopic branch.
- `Open edge`: the true branch preserving those monomials remains an open PDE
  theorem.

## Stage `169`

- `Explicit`: same coherent local D/N tracking and grouped weak-axisymmetric
  setting as Stage `168`.
- `CAS`: `chi, delta > 0`, `E, F` real, and the selected rank-3 minor is used.
- `Branch`: the scalar axial scale `L` is inert at first order in the present
  reduction.
- `Branch`: the constructive coherent branch requires `chi0,* > 0` so the
  selected minor is nonzero and the monomial-drift map has rank three.
- `Open edge`: this stage identifies the exact finite similarity orbit whose
  tangent space solves the linearized defect equations, but remains an
  infinitesimal/tangent-space theorem.

## Stage `170`

- `Explicit`: positive microscopic state space `M_+` and the coherent local D/N
  tracking branch.
- `CAS`: relevant microscopic variables are real, with
  `chi`, `deltaU`, `E`, `F > 0`.
- `Branch`: the finite quotient theorem is proved only on the positive-coupling
  sector and constructive coherent branch `chi0,* > 0`.
- `Branch`: equality of invariants is reduced to the same finite log-ratio
  system `M_* Delta x = 0`.
- `Open edge`: the theorem does not yet prove that the actual dynamical PDE
  branch stays on a single similarity orbit or preserves the three quotient
  invariants.

## Hardening Notes

- The highest-value unresolved assumption question is no longer algebraic; it is
  whether the real completed PDE actually lands on the constructive/coherent
  branches assumed in `028-031`, `089`, and `168-170`.
- The mouth-side chain `125-140` is mathematically solid inside the explicit
  Family-1 positive-mouth closure, but that closure itself is still a model
  choice, not a theorem about every positive mouth profile.
- The late orbit/quotient chain `168-170` is exact as an invariant-structure
  theorem inside the coherent positive sector, but it is still missing the final
  dynamical branch-selection theorem.
