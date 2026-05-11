# Moving-Throat PDE — Dependency Ledger

## Scope

This ledger makes the proof chain explicit for the first hardening tranche from
`notes/moving_throat/review/PROOF_HARDENING_PLAN.md`.

Covered stages:

- `003-004`
- `037-048`
- `058`
- `106`
- `125`
- `142-157`
- `170`
- `185-187`

For each stage, the ledger records:

- direct inputs actually used,
- exported formulas or theorem statements,
- output type,
- and the main downstream consumers inside the present choke-point set.

## Output-Type Legend

- `Exact identity`: algebraic identity or reduction formula.
- `Exact theorem`: exact branch, existence, or equivalence result under the
  recorded assumptions.
- `Perturbative`: first-order, weak-coupling, or asymptotic statement.
- `Numerical fixed point`: exact defining equation solved numerically.
- `Status`: consolidation stage only.

## Stage `003`

- `Inputs`: Stage `002` geometry-only wall reduction and the confined BdG sector.
- `Exports`: exact eliminated wall kernel `D0_eff`, low-frequency
  `K/M/N` renormalization, two-pole shift formulas, grouped-`P2`
  isotropy/anomaly diagnostics.
- `Output type`: exact identity plus perturbative pole shift.
- `Consumed by`: Stage `004`; conceptually by later reduced wall-sector checks.

## Stage `004`

- `Inputs`: Stage `003` conservative wall/matter lane; prior localized Maxwell
  and mixed-sector reductions; earlier outgoing `l=2` fingerprint theorem.
- `Exports`: exact conservative EM+mixed self-energy, outgoing transfer factor,
  wall quadrupole odd term, scalar odd-term rescue from `i omega` to `i omega^3`.
- `Output type`: exact identity plus first-order outgoing dressing.
- `Consumed by`: no direct later choke-point stage in this tranche, but it is a
  foundational reduced-sector checkpoint for the outgoing machinery.

## Stage `037`

- `Inputs`: Stage `036` admissibility geometry and the `034-036` branch data.
- `Exports`: explicit continuum operator data
  `(A, DeltaK_ax, beta0, alpha_mix, M_mix, delta)` and Schur-factorized wall
  kernel `Sigma_wall = Xi I + alpha v v^T`.
- `Output type`: exact identity.
- `Consumed by`: Stages `038`, `039`.

## Stage `038`

- `Inputs`: Stage `037` continuum formulas.
- `Exports`: exact five-ratio placement map and product law
  `R_target M_mix = 8 Lambda (1 - eps_W) / pi^2`.
- `Output type`: exact identity.
- `Consumed by`: Stage `039`; conceptually by the later coherent-kernel chain.

## Stage `039`

- `Inputs`: Stages `037-038`.
- `Exports`: split-`U` quantities
  `delta_split`, `eps_W^split`, direction factor `R_U`, splitting invariant
  `D_dir`; proof that scalar placement survives while source/loading
  collinearity generally fails.
- `Output type`: exact identity/theorem.
- `Consumed by`: Stages `040`, `043`, `044`, `045`, `047`.

## Stage `040`

- `Inputs`: Stage `039` and flat-limit comparison to Stages `035-036`.
- `Exports`: generalized selected-branch functions `G_q` and `F_(q,eta)`,
  recovering the flat-`U` branch at `R_U = 1`.
- `Output type`: exact identity.
- `Consumed by`: Stages `041`, `042`.

## Stage `041`

- `Inputs`: Stage `040`.
- `Exports`: exact support-loading theorem giving `n_req(m)` and strict
  monotonicity `dn/dm < 0`; tracking and source-tied specializations.
- `Output type`: exact theorem.
- `Consumed by`: Stages `042`, `043`, `044`.

## Stage `042`

- `Inputs`: Stages `040-041`.
- `Exports`: exact three-direction normalization law `F_(q,r,t)(xi,delta;m)`,
  tracking collapse, source-tied specialization `F_src`.
- `Output type`: exact identity.
- `Consumed by`: Stages `043`, `044`.

## Stage `043`

- `Inputs`: Stage `039` split-`U` continuum plus the rank-2 bottleneck from
  Stages `041-042`.
- `Exports`: support direction factor `R_phi`, split support blocking
  `eps_phi^split`, support baseline `M_supp`, exact tracking condition
  `sigma0 = rho0`.
- `Output type`: exact theorem/identity.
- `Consumed by`: Stages `044`, `045`.

## Stage `044`

- `Inputs`: Stages `039`, `041`, `042`, `043`.
- `Exports`: exact quadratic physical branch equation
  `xi^2 + B_cont xi + C_cont = 0`, continuum normalization gate
  `R_target = F_cont(xi_phys)`, mismatch penalty
  `lambda0 M_mix M_supp (R_U - R_phi)^2`.
- `Output type`: exact theorem.
- `Consumed by`: Stage `045`.

## Stage `045`

- `Inputs`: Stages `039`, `043`, `044`.
- `Exports`: coherent local-kernel tracking law
  `rho0 = sigma0`, hence `R_tr = R_U = R_phi`, total baseline
  `M_tr = M_mix + M_supp`, and one-parameter collapse of the Stage-`044` branch.
- `Output type`: exact theorem.
- `Consumed by`: Stages `046`, `047`, `048`.

## Stage `046`

- `Inputs`: Stage `045`.
- `Exports`: exact ordering theorems
  `dG_tr/dR < 0`, `dF_tr/dR > 0`, loading excess over the flat branch,
  normalization deficit, residual ordering `E_tr > E_flat`.
- `Output type`: exact theorem.
- `Consumed by`: Stages `047`, `048`.

## Stage `047`

- `Inputs`: Stages `039`, `045`, `046`.
- `Exports`: exact coherent support-enhancement factor
  `S(zeta;eps) = 1 + zeta (1 - eps)/(1 - zeta eps)`, total baseline
  `M_tr = M_mix S`, and proof that `R_target` is independent of `zeta`.
- `Output type`: exact identity/theorem.
- `Consumed by`: Stage `048`.

## Stage `048`

- `Inputs`: Stage `047` plus the tracking laws from `045-046`.
- `Exports`: exact critical load `M_crit`, inverse required-support map
  `zeta_req`, stable-side support compensation theorem, monotone branch response
  `d xi_phys / d zeta > 0`.
- `Output type`: exact theorem.
- `Consumed by`: no direct later choke-point stage in this tranche, but it
  closes the early coherent support-feasibility chain.

## Stage `058`

- `Inputs`: Stages `056-057`.
- `Exports`: exact coupled support/source fixed-point equation
  `Pe = Xi Delta(Pe;kappa,eta)`, interval bounds
  `Xi Delta_0 <= Pe_* <= Xi Delta_inf`, weak-coupling law
  `Pe_* = Xi Delta_0 + O(Xi^2)`.
- `Output type`: exact existence/bracketing theorem plus perturbative expansion.
- `Consumed by`: later source/support chain beyond the present tranche.

## Stage `106`

- `Inputs`: Stage `101`, Stage `105`, and the reduced outgoing stack from
  `099/100/103/104`.
- `Exports`: canonical outgoing reduced closure with `N_Q = 1` and the GR-target
  quadrupole coefficients including the canonical `GammaBar_5`.
- `Output type`: exact theorem on the canonical outgoing branch.
- `Consumed by`: later outgoing/canonical closure stages beyond this tranche.

## Stage `125`

- `Inputs`: local positive-source kernel setup; Family-1 branch values carried
  from the earlier compensation chain.
- `Exports`: positive-source theorem excluding `g_+^F1` and retaining the unique
  physically admissible lower branch `g_-^F1`.
- `Output type`: exact theorem within the positive-source setup.
- `Consumed by`: Stages `144`, `157`.

## Stage `142`

- `Inputs`: Stages `130`, `138`, `140`.
- `Exports`: self-consistent Family-1 exponential mouth branch
  `Sigma0(Pi)`, `Tmhat(Pi)`, canonical point `(Pi_*, Tmhat_*)`.
- `Output type`: exact defining relations plus numerically identified canonical
  point.
- `Consumed by`: Stages `143`, `144`, `146`, `147`, `150`, `155`, `156`.

## Stage `143`

- `Inputs`: Stages `130`, `142`.
- `Exports`: theorem that `0 < g_Pi < 1` for every finite `Pi > 0`, with
  `g_c = 1` only in the singular limit `Pi -> infinity`.
- `Output type`: exact theorem plus asymptotic limit.
- `Consumed by`: Stages `144`, `157`.

## Stage `144`

- `Inputs`: Stages `125`, `130`, `142`, `143`.
- `Exports`: unique regular finite lower compensated branch inside the explicit
  Family-1 positive-mouth closure; exclusion of upper and equal-normalized
  branches.
- `Output type`: exact theorem within the explicit closure.
- `Consumed by`: Stages `145`, `146`, `154`, `155`, `156`, `157`.

## Stage `145`

- `Inputs`: Stages `130`, `142`, `144`.
- `Exports`: mouth-side status reduction to finite corrections around the unique
  regular branch `(Pi_*, Tmhat_*)`.
- `Output type`: status.
- `Consumed by`: interpretive only; sets context for `146-153`.

## Stage `146`

- `Inputs`: Stages `142`, `144`.
- `Exports`: first-order retuning law
  `delta Pi = -eps (gbar_varpi - g_*) / g'_*` and induced `delta S_q`, depending
  only on the two source moments `(gbar_varpi, Sbar_varpi)`.
- `Output type`: perturbative.
- `Consumed by`: Stages `147`, `149`.

## Stage `147`

- `Inputs`: Stages `142`, `146`.
- `Exports`: first-order rigidity kernel
  `delta Tmhat = eps [A_T (gbar_varpi - g_*) + B_T (Sbar_varpi - S_*)]` and the
  reduced effective kernel `W_*`.
- `Output type`: perturbative.
- `Consumed by`: Stages `148`, `149`, `151`.

## Stage `148`

- `Inputs`: Stage `147`, plus the earlier exact compensation fraction from
  Stage `126`.
- `Exports`: explicit evaluations of the rigidity kernel on two positive
  non-exponential families and their interpolation.
- `Output type`: perturbative/example evaluation.
- `Consumed by`: Stage `149`, Stage `152`.

## Stage `149`

- `Inputs`: Stages `146-148`.
- `Exports`: status conclusion that mouth-side ambiguity has reduced to one
  canonical branch plus one rigidity kernel and modest finite shifts.
- `Output type`: status.
- `Consumed by`: interpretive only.

## Stage `150`

- `Inputs`: Stages `135`, `142`.
- `Exports`: exact residual profile `R_*(x) = Phi_*(x) - Pi_* x` with tangent
  matching at `x = 0` and `R_*''(0) < 0`.
- `Output type`: exact identity/theorem on the canonical branch.
- `Consumed by`: Stages `151`, `152`, `153`.

## Stage `151`

- `Inputs`: Stages `147`, `150`.
- `Exports`: covariance formulas
  `delta g_act = -Cov_*(c, R_*)`, `delta S_act = -Cov_*(K_q, R_*)`, and induced
  first-order `delta Pi_act`, `delta Tmhat_act`.
- `Output type`: perturbative.
- `Consumed by`: Stage `152`.

## Stage `152`

- `Inputs`: Stages `148`, `150`, `151`.
- `Exports`: numerical canonical correction estimates
  `(Pi_corr, Tmhat_corr)`, one-step Picard update `(Pi_1, Tmhat_1)`, and
  effective broadening `lambda_eff`.
- `Output type`: perturbative numerical correction.
- `Consumed by`: Stage `153`, comparison in Stage `156`.

## Stage `153`

- `Inputs`: Stages `150-152`.
- `Exports`: status conclusion that the lower compensated branch survives but is
  shifted upward in both bias and traction.
- `Output type`: status.
- `Consumed by`: interpretive setup for `154-157`.

## Stage `154`

- `Inputs`: Stage `140` and the core-law data carried through `142/144`.
- `Exports`: exact co-evolving core-mouth map, canonical equivalence
  `g = g_* iff R = 1/4`, defect transport law, slope identity
  `Pi = Sigma0 (1 - R S)`.
- `Output type`: exact identity/theorem plus first-order defect transport.
- `Consumed by`: Stages `155`, `156`.

## Stage `155`

- `Inputs`: Stage `154`, canonical point from `142/144`.
- `Exports`: unique positive fixed point on frozen traction
  `(g_fp, R_fp, Pi_fp)`; theorem that the branch survives near the canonical
  point but is no longer exactly compensated.
- `Output type`: numerical fixed point.
- `Consumed by`: Stage `156`, Stage `157`.

## Stage `156`

- `Inputs`: Stages `154-155`, original canonical target from `142/144`,
  corrected comparison point from `152`.
- `Exports`: renormalized canonical co-evolving branch
  `(Sigma0_can, Tmhat_can, Pi_can)` with `R_can = 1/4`.
- `Output type`: numerical fixed point/root on an exact defining condition.
- `Consumed by`: Stage `157`.

## Stage `157`

- `Inputs`: Stages `125`, `143`, `144`, `155`, `156`.
- `Exports`: status summary that the explicit Family-1 closure is an honest
  nonlinear fixed point; frozen old traction gives a nearby non-compensated
  branch, while exact compensation requires the renormalized branch from `156`.
- `Output type`: status.
- `Consumed by`: interpretive only.

## Stage `170`

- `Inputs`: Stages `022`, `023`, `159`, `169`.
- `Exports`: exact linear grouped outlet map collapsing the bottleneck to lane
  variables `(K_A, G_A)` with the hidden-even relation
  `delta D_{A,4} = (2/3) delta D_{A,2} + delta D_{A,0}/27`.
- `Output type`: exact identity/theorem.
- `Consumed by`: later microscopic grouped/anisotropy chain leading to
  `185-187`.

## Stage `185`

- `Inputs`: Stages `182`, `183`, `184`.
- `Exports`: zero-defect equivalence with preservation of the three microscopic
  monomials `(C_tr, C_nt, epsilon_eta)` and the linear compatibility ledger for
  the grouped drift variables.
- `Output type`: exact linearized theorem.
- `Consumed by`: Stages `186`, `187`.

## Stage `186`

- `Inputs`: Stages `183`, `185`, and the invariant package from `184-185`.
- `Exports`: rank-3 monomial-drift matrix `M_*`, equivalence
  `Theta_1 = Xi_1 = R_1 = 0 iff M_* delta x = 0 iff delta x in T_id G_*`, and
  the exact five-parameter similarity orbit `G_*`.
- `Output type`: exact theorem with tangent-space closure.
- `Consumed by`: Stage `187`.

## Stage `187`

- `Inputs`: Stages `183`, `184-186`.
- `Exports`: finite orbit-fibre theorem for the invariant map
  `I = (C_tr, C_nt, epsilon_eta)`, solved finite log-ratio laws, and quotient
  identification `M_+ / G_* ~= (R_{>0})^3`.
- `Output type`: exact theorem in the positive coherent sector.
- `Consumed by`: current endpoint of the derivation chain.

## Hardening Notes

- The strongest remaining dependency risk is not local algebra; it is silent
  movement from “exact inside an explicit closure” to “exact for the full PDE.”
  That is most acute in Stages `106`, `125`, `142-157`, and `185-187`.
- The early coherent support chain `037-048` is structurally clean: every major
  object used later is defined before it is consumed, and the key branch
  specialization happens explicitly at Stage `045`.
- The late orbit chain is also clean at the dependency level: Stage `185`
  defines the invariant ledger, Stage `186` identifies its tangent-space orbit,
  and Stage `187` integrates that to the finite quotient statement.
