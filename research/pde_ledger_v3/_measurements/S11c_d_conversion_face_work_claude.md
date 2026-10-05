Continuing the review: I've read the proposal, the guide, and the staged sources for S11b, S12 and stage006. I also ran a paper check of the task_b2d identity and the integration by parts. Verdict and findings below.

**Verdict: COHERENT CONDITIONAL BRIDGE.** The corrections needed are bookkeeping and scope, not physical.

## 1. Moving-volume identities and the front sign

- **Transport identities hold.** Both Reynolds-type identities follow from the two balances: d/dt∫ρ = −∮ρ(v−b)·N, and the order identity adds −∮K·N and +∫ρΓ. Only the normal part of b matters.
- **Mass versus number.** Setting ρ=mn and K=mJ_χ with constant m is a pure unit change. Γ cancels in the sum of the order and disorder balances (A2), so there is no total-mass source. This is the S12 content, and the "frozen-wall sink plus remote return" shortcut is not reintroduced.
- **Front sign is right.** Steady-in-a-frame continuity gives constant j=ρ(v_n−V). Then d_z[χj+K_n]=ρΓ integrates to j=K_out−K_in−G with χ_in=1, χ_out=0. With z pointing from ordered to disordered, positive j means mass leaving the ordered side. This matches the S11b rule J_±=ρ_m(v_w−∂_tζ_±)(±1), where J>0 removes material from the slab. So G>0 gives J<0, which is accretion. The sign is consistent.
- **Scope is handled correctly.** The text says this is not an existence result and not an imposed background. It says curved or unsteady fronts return to the full balance. It also separates the DC integral G from the harmonic J_f by order of perturbation, and it does not claim the geometry join.
- **Minor fixes:**
  - State that "Γ is total" is a reading. P8 writes D_tχ=−M_χμ_χ+Γ_B, while A2 and the balance carry Γ_B and J_χ only. The proposal's reading is acceptable if the relaxation is absorbed in either Γ or K, never both.
  - The steady front also needs ∂_tn=0 in the moving frame and zero endpoint gradients, so the K_n=0 case is a hypothesis, not a free choice.

## 2. Integration by parts

- **Algebra is correct.** D_tχ=Γ−ρ⁻¹∇·K follows from the two balances, since n·D_tχ equals the left side of the order balance after continuity. Then −∫(μ/ρ)∇·K = +∫K·∇(μ/ρ) − ∮(μ/ρ)K·N.
- **Density placement is right.** The erratum's P_order=∫μ_χ D_tχ has no extra n, and μ_χ/ρ carries specific-energy units. The proposal does not use the erroneous handoff form μ n Γ.
- **Scope is right.** It claims only the order-work contribution. It correctly says this is neither a moving-volume energy law nor a mode-energy balance.
- **Gradient-energy boundary term.** μ_χ=δF/δχ includes −κ∇²χ, and the variation also leaves a boundary term κ∇χ·N D_tχ. The proposal names this as a missing term but does not write it. Writing it is cheap and removes any suspicion that it was dropped.
- **No double counting.** The kinetic adjunct is not added. The one risk is the relaxation piece −M_χμ_χ, which makes μΓ contain −Mμ², a dissipation. It must stay inside the single "total Γ" and not appear again.

## 3. Native face power

- **Identity checks on paper.** The task_b2d identity (p+Λ_XA)V̄+μ_sJ̄ = p·v̄_bulk+AJ̄+Λ_XAV̄ is an algebraic consequence of A=μ_s−p/ρ_m and v_bulk=V+J/ρ_m. The pJ̄/ρ_m and AJ̄ terms combine to μ_sJ̄. Re/2 is the harmonic mean convention. The proposal uses it faithfully and does not equate it to P_χ.
- **Missing joins before comparing μ_χ/ρ with μ_s or A_f.** Unit equality is not enough, because:
  - μ_s=μ_θ/ρ_br⁰ with μ_θ=δU/δθ at fixed u and e_W. That is a variation in the slab's density and thickness variable θ. μ_χ is a variation in χ, at fixed n or fixed some other variable. They are conjugate to different state variables, and the map between the (n,χ) description and (θ, ρ_br⁰) is not given.
  - A_f subtracts a bulk pressure perturbation δp_f/ρ_m. Phase conversion between ordered and disordered phases is driven by a difference of per-mass chemical potentials, or equivalently a grand-potential difference. Which pressure and which phase the ρ_m refers to is not fixed.
  - The S11b slab has two faces and wave perturbations on a rest bulk. A single DC front is a different object, and the rest-bulk wave operator does not carry background drain (a scope limit the guide notes).
  - Λ_X and the memory times τ_I are prescribed, not derived from F.

## 4. What the shear gate establishes

- **Storage.** Shear energy is stored only where χ>0, and ∂(χf_sh)/∂χ=f_sh≥0 at fixed u_d and fixed μ_R. At fixed displacement, ordering costs shear energy and disordering releases it. The proposal says exactly this and correctly does not infer that conversion pumps the wave.
- **Product rule.** d/dt∫_{Ω_m}χf = ∫[fD_tχ+χD_tf+χf∇·v] is correct for a material volume. The three terms must stay separate. Using u_d≠v, D_tf_sh contains ∂_t u_d dynamics, which no supplied equation fixes. If μ_R depends on n, D_tf also gets a compressional piece.
- **Left open:**
  - **Creation state of u_d.** Newly ordered material arrives with some displacement. The balances do not say whether u_d is carried with the material, continuous across the front, or reset. This decides how much shear energy appears or disappears at the front.
  - **Mode inertia.** The kinetic term is ½mn|u|² with material velocity. The relation of that to a trapped shear mode's inertia and to ∂_t u_d is not supplied. The proposal does not treat it as a solved inertia, which is correct.
  - **Mouth stress.** The shear stress ∂(χf)/∂(∂u_d) is derivable from F. But a front also feels a configurational force from f_sh∇χ and the jump of the energy-momentum flux, which is not written. The first is derived from F. Neither gives a sign or a magnitude without the throat core.

## 5. Smallest next physics deliverable

One short paper, not a solver, with two separated parts.

**Part A: derived, from the existing F and balances.**
- The normal energy-momentum (Eshelby) jump across a χ front, plus the shear-stress traction on the mouth, written for a given u_d and χ profile.
- The material-volume energy bookkeeping for the ordered region, listing every term: order work, kinetic/compressional work, gradient boundary work, shear storage, and boundary transport. This replaces the claim that the pieces cannot be joined with an explicit ledger that has named gaps.

**Part B: the genuinely missing functions, each named, none invented.**
1. The creation rule for u_d at the front, which I recommend asking the user to choose. Options are (i) u_d is a material-label displacement carried continuously across the front (simplest; shear energy appears automatically because of the χ gate), (ii) u_d resets to zero on re-ordering, or (iii) the stored shear energy is conserved across the front. I recommend (i) as the minimal choice, but it is a physical decision for the user.
2. The mode inertia linking u_d to the material kinetic term.
3. f_throat and f_mix (placeholders), Γ_return, Γ_drain, J_χ, and the orientation/core coupling.
4. The (n,χ)→(θ, μ_s, A_f) map needed for any comparison with the S11b face closure.

None of this predicts an electric sign or leakage, and the proposal correctly makes no such claim.

## Source coverage and limitations

- **Read in full:** `input/proposal.md` and `input/guide.txt`.
- **Read in part:**
  - `stage006-material-record.md`, lines 66–98 and 284–297 (balances, pins P1–P3, P6, P8, A2 grep hits, errata).
  - `s11b-shared-physics.md`, lines 190–240, plus grep hits for the outward-normal and affinity conventions.
  - `v3-step-plan.md`, lines 580–649 (dark-energy postulate and S12).
  - `s11b-source.py` task_b2d, lines 1725–1736, plus dimension and homogeneity hits.
- **Not inspected:** `native-interpretation.md`, `brane-bulk-handoff.md`, `s11b-result-record.md`, `interface-geometry.py`, and `stage006-material-source.py`, apart from grep hits. The claims about the geometry operands, the handoff erratum text, and the f_throat/f_mix placeholder status rest on the proposal and guide statements, not my reading.
- **Not done:** no CAS, code execution or replay. The algebra was checked by hand. I did not independently re-derive the pressure/area conventions of the native curved-geometry join, so that remains unverified.
- **Scope:** this assesses the paper step only. It is not a numerical-result acceptance, and no scientific worker, READY gate or leakage work is implied.