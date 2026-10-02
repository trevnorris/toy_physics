# Closed direct mixed kernel at grazing: proposed bounded method

2026-10-01. This is an author derivation for independent assessment, not an executed limit, evaluated convolution, cleared worker, or loss result. The next objective is to close the exact-match applicability gap of the saved both-face direct mixed increment without repeating its construction. The proposed output is a source-bound L1 limiting kernel and an explicit local bound. It is not a defect sweep.

## 1. Scope and actual saved starting point

Keep strict-rest-bulk LAB_HELD/RHO4_CONSTANT, real frequency 3, conserved edge momenta (1/5,1/10), W=1, L=10, rho_m=1/10, Lambda_A_0=1/100, tau_A=1/10, and independent eta/sigma_W. Physical permeability and memory remain live. Only the direct kernel grade (1,1) with source and consumer grade (0,0) is considered; pure eta^2 and sigma_W^2 and the other operator blocks are not completed here.

The saved raw diagnostic derives the lower boundary, joins both faces to native closure/reference/source/consumer records, and establishes the ordered contact cancellation for nonzero external depths. Its unfinished controls were completed using saved returns; no integral was evaluated. The packet supplies actual operands, exact residuals and controls rather than relying on an author's success label. These are source-joined author-run observations, not independent physics revalidation.

Use profile momenta k (input), l=k+Q (output), t (height transfer) and Q-t (slope transfer). The independent positive outgoing depths are qi=q(k), qh=q(k+t), qs=q(l-t), qo=q(l). The reflected route qs is not qh. At real frequency,

    q(p) = sqrt(9/cs^2 - 1/20 - p^2),               positive radicand,
           i sqrt(p^2 - (9/cs^2 - 1/20)),           negative radicand.

Zero radicands mean a limiting outgoing sheet, not pointwise division by zero. Both face extensions retain their native outward/lab signs; the saved scalar outward coefficient is equal on the two faces, while normal jets have opposite signs.

The saved factorization (physical-factorization-return.json) is C=t B with

    B = -i omega rho_m [
          k(2l-t)/(qi qo(qs+qo))
        + k(t+2k)/(qh qo(qh+qi))
        + qi/(qh qo) ].

It is an ordered coefficient, not the sum of two full-weight assignments. The saved full coefficient and physical dispersion identities establish this factorization on the regular domain. The complete contact coefficient at t=0 is zero there. With

    A(t) = L t / [4 sinh(pi L t/2)],  A(0)=1/(2 pi),
    h_hat(t) = (W/2)[delta(t)/2 + PV(A(t)/(i t))],
    j_hat(u) = L A(u)/2,

the saved ordinary raw density is (W L/(4i)) A(t) A(Q-t) B. No extra swapped term or 2 pi factor is introduced. The original contact/PV expression must remain saved alongside this representation.

The native face definition is I+a(omega)Z, where

    a(omega)=Lambda_A_0/[rho_m^2(1-i omega tau_A)],
    beta(omega)=a(omega) rho_m omega.

At omega=3, the saved values on BOTH faces are a=(100+30i)/109 and beta=(30+9i)/109. The two external factors multiply to qi qo/[(qi+beta)(qo+beta)]. The saved zero-grade reference trace is 1; the separate normal-jet map is i f qo for face f=+1 or -1. The first-shape iteration stays once and unchanged. An already integrated whole D never goes into another middle integral.

## 2. Proposed canceled expression and domain

Multiplying the saved B by the native two-leg factor gives the proposed identity

    Bc = -i omega rho_m / [(qi+beta)(qo+beta)] * [
           k(2l-t)/(qs+qo)
         + k(t+2k) qi/[qh(qh+qi)]
         + qi^2/qh ].                                      (1)

The closed reference density is G=(W L/(4i)) A(t) A(Q-t) Bc. This identity is a new hand rearrangement of saved operands, not a newly executed symbolic check. A later instrument must join it exactly to each saved closed raw/reference density before applying any limit. The normal jet is i f qo G; grade-zero native row/source multipliers are already saved and must be joined, not reconstructed as guesses.

Use a declared compact envelope cs in [1,2], k,l in [-3,3]. This is a local mathematical domain containing both saved modal matches, not a physical calibration or a new sweep. Those matches are cs=sqrt(3/2), k=+/-sqrt(595)/10 at LEFT and cs=sqrt(150/101), k=+/-sqrt(601)/10 at RIGHT. The domain also covers their forward/backward pairings. Its real acoustic radius kappa=sqrt(9/cs^2-1/20) is at least sqrt(11/5)>0. The case kappa=0 is explicitly excluded; there the two roots of one route merge, and the simple-root bound below does not apply. No assertion is made about that bulk-edge degeneracy or an impermeable beta=0 limit.

## 3. Outgoing regularization and uniform local bound

To assess the physical limiting prescription, let Omega=3+i delta, 0<delta<=1/10, and use q_delta(p)=sqrt(Omega^2/cs^2-1/20-p^2) in the first quadrant. Continue the SAME native a(Omega), beta(Omega) and the explicit omega factor in (1). This is the source's limiting-absorption prescription for a proof; it is not a replacement for physical permeability or a finite-regulator simulation. Its legitimacy and exact native joins require assessment. No conjugated power/current form is being analytically continued here.

Write D_delta=(1+tau_A delta)^2+9 tau_A^2. Directly from the saved native law,

    Re beta(Omega)=(Lambda_A_0/rho_m)*3/D_delta > 0,
    Im beta(Omega)=(Lambda_A_0/rho_m)*
                     [delta(1+tau_A delta)+9 tau_A]/D_delta > 0.

Thus beta and all four depths lie in the closed first quadrant, with beta strictly interior. For u,v in that quadrant, |u+v|>=max(|u|,|v|). In particular, |qi+beta| and |qo+beta| are at least beta_min=3000/11101; |qi/(qh+qi)|<=1 and 1/|qs+qo|<=1/|qs| away from zeros. These bounds extend almost everywhere to the boundary sheet. They do not require qi or qo to remain bounded away from zero.

On the declared compact domain, |Omega|<=4 and |qi|,|qo|<=5 suffice. Equation (1) therefore has the explicit envelope

    |Bc| <= (4 rho_m/beta_min^2) * [
                   (18+3|t|)/|qs| + (43+3|t|)/|qh| ].        (2)

Define kappa_delta=sqrt((9-delta^2)/cs^2-1/20), including delta=0. Throughout this domain kappa_delta>=sqrt(879)/20=:kappa_min>0. Since

    |q_delta(p)| >= sqrt(|kappa_delta^2-p^2|),
    |kappa_delta^2-p^2| >= kappa_min * dist(p,{+kappa_delta,-kappa_delta}),

each inverse depth is bounded by kappa_min^(-1/2) times the SUM of the inverse square-root distances to its two real comparison endpoints. In transfer coordinates the endpoints are -k+/-kappa_delta for qh and l+/-kappa_delta for qs. The comparison endpoints move with delta and cs; they are not claimed to be complex branch points at finite delta.

For any measurable real set E of measure m, the integral of |t-a|^(-1/2) over E is bounded by 2 sqrt(2m), uniformly in a. Equations (2) and this elementary rearrangement bound give uniform absolute continuity of the integral on every compact t interval. Different routes' endpoints may coincide: the bound remains a SUM of four inverse square roots, not a product of singular factors. This is the key collision check to assess.

The tanh factors have removable values at t=0 and t=Q and exponential tails. With |Q|<=6, A(t)A(Q-t) times the polynomial in (2) gives a uniform integrable tail outside a sufficiently large compact interval (all comparison endpoints remain bounded). Record an explicit tail estimate in the later instrument; no numerical integral is needed to establish it. Compact uniform integrability plus that tail and pointwise convergence off the finite limiting endpoint set should establish convergence of G in L1(R) along any parameter path in this domain, including delta->0 and qi/qo->0 together.

This is a uniform-integrability argument with moving singularities, not an unsupported assertion of one pointwise dominating function. The source-bound proof and its constants, limiting paths and closure domain are the substance of this review.

## 4. Limit, contact and strength of the proposed conclusion

Away from internal endpoints, if qi->0, the second and third terms in (1) tend to zero and the first tends to

    -i omega rho_m k(2l-t)/[beta(qo+beta)(qs+qo)].             (3)

If only qo->0, use (1) with qo=0 on its almost-everywhere domain. If both external depths vanish, (3) becomes

    -i omega rho_m k(2l-t)/(beta^2 qs).                      (4)

The two simultaneous real grazing configurations are l=k and l=-k, with |k|=kappa>0. Singular values at a finite number of t points are not assigned a finite pointwise value. Their almost-everywhere density and L1 limit define the object. The proposed result is finite, path-independent convolution against any bounded test function, not differentiability in cs, a pointwise regular kernel, a global operator norm bound, or smooth physical leakage.

At every delta>0, the common four-leg dispersion still gives C(t=0)=0, and t PV(1/t)=1 is valid with the smooth finite-delta coefficient. The contact is therefore zero before taking any limit. The uniform L1 argument would exclude a hidden concentrated contact term in the limiting ordinary density. Review this connection explicitly; a finite pointwise calculation or separately canceled summands cannot replace it. If analytic continuation/contact order fails, retain the narrower proven statement and mark the physical joint limit unresolved.

The factor i f qo is uniformly bounded and continuous in the domain, so the same L1 statement should pass to the separate native normal jet after exact joining. Saved grade-zero source/row factors depend on external momenta, not t; their finite polynomial symbols may multiply the bound on compact k,l. This needs actual coefficient/slot joins on both faces. A symbolic jet amplitude is not an incident light state. No conclusion transfers automatically to first-shape iteration, the full retained operator, the untruncated finite inverse, its current balance, a drain, or primitive calibration. Finite behavior of this one missing direct term does not prove that loss is smooth through speed matching.

## 5. Proposed bounded implementation and acceptance

After substantive method assessment, author one task-local instrument under the existing no-deadline pooled guard. It restores completed factorization, closure/reference and selected row operands without replaying their functions. It derives only the cancellation/limit/bound certificate absent from the saved work. Keep source-symbol assumptions, exact actual arguments, dimensions and Fourier normalization attached to each result.

Required outputs are: exact per-face joins for (1); the unchanged physical beta law and its sign/lower-bound certificate; four-leg root/domain certificates; simple-root and moving-endpoint bound with explicit tail; qi-only, qo-only and simultaneous same/opposite-momentum almost-everywhere limits; the contact/PV limiting argument; and separate pressure/normal-jet/grade-zero consumer applicability. Save operands and residuals before each guard. A sampled approach to grazing cannot substitute for the limit/bound. Unknown or nonzero joins stop with evidence. Do not numerically evaluate the convolution or build a finite matrix in this instrument.

Controls must expose a missing external factor, a wrong root sheet and incorrect lower normal-jet sign through the corresponding expression or certificate. They must be shown applicable and responsive, with the same source bindings; algebraic self-reconstruction alone is not a control. Keep completed old sensitivity evidence as provenance, not a claim of independent correctness. No completed result is rerun to manufacture controls.

If this bounded proof is supported, the next decision is what additional whole-operator/integration applicability is needed before a corrected near-unity finite pilot. If a nonintegrable term or unjoined physical prescription remains, report it as a method finding. Do not disguise it as a tolerance or formatting repair, expand into a general exterior theory, or launch the defect sweep automatically.
