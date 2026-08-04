# Relation candidates — review required
**REQUIRES REVIEW BEFORE INSERTION.** These prose/algebra records were not auto-compared, and no relation YAML rows were generated.
- Mathematica transcript: `/home/trevnorris/s11bB_build/wl_run5.txt`
- SymPy transcript: `/home/trevnorris/s11bB_build/py_run4.txt`

## Candidate 1: `WL_ENERGY_BASIS` ↔ `S11BB_ENERGY_BASIS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ENERGY_BASIS`

    Modulo total in-plane divergences, the computed O(3), w-reflection-even, translation-invariant basis is {(|curl u|^2)/2,(div u)^2/2,theta^2/2,theta e_W,e_W^2/2,|grad theta|^2/2,grad theta.grad e_W,|grad e_W|^2/2,theta div u,e_W div u}.  The antisymmetric and symmetric-traceless grad-u contractions reduce to the listed curl/div pair modulo a divergence; parity excludes epsilon contractions; undifferentiated u is excluded by translation invariance.

### SymPy — `S11BB_ENERGY_BASIS`

    fields/first gradients=(u,grad u,theta,grad theta,e_W,grad e_W). Modulo divergences the O(3), w-reflection, translation-invariant quadratic basis is {curl(u)^2,(div u)^2,theta div u,e_W div u,theta^2,theta e_W,e_W^2,|grad theta|^2,grad theta.grad e_W,|grad e_W|^2}

## Candidate 2: `WL_ENERGY_BASIS_INDEPENDENT_TERMS` ↔ `S11BB_ENERGY_BASIS_INDEPENDENT_TERMS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ENERGY_BASIS_INDEPENDENT_TERMS`

    U=(muR/2)|curl u|^2+(kD/2)(div u)^2+(b3/2)theta^2+cc*w0*theta*e_W+(kW*w0^2/2)e_W^2+(gTheta/2)|grad theta|^2+gThetaE grad(theta).grad(e_W)+(kapW*w0^4/2)|grad e_W|^2+aTheta*theta*div(u)+aE*e_W*div(u); all ten displayed bilinears have independent free coefficients before B1.

### SymPy — `S11BB_ENERGY_BASIS_INDEPENDENT_TERMS`

    all ten listed basis bilinears are independent before B1. The three possible grad-u contractions reduce modulo a divergence to curl^2 and div^2; O(3) parity forbids epsilon-pseudoscalar terms

## Candidate 3: `WL_ENERGY_BASIS_OMISSIONS` ↔ `S11BB_ENERGY_BASIS_OMISSIONS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ENERGY_BASIS_OMISSIONS`

    Relative to the supplied list, the allowed independent omissions are (kD/2)(div u)^2, (gTheta/2)|grad theta|^2, gThetaE grad(theta).grad(e_W), aTheta theta div(u), and aE e_W div(u).  They enter the assembly through mTheta=b3+gTheta*k^2, mCross=cc*w0+gThetaE*k^2, piD=aTheta-kD, piTheta=mTheta-aTheta, and piE=mCross-aE.

### SymPy — `S11BB_ENERGY_BASIS_OMISSIONS`

    the supplied list omits (div u)^2, theta div u, e_W div u, |grad theta|^2, and grad theta.grad e_W; they are carried as K_L,D_theta,D_e,A_theta,A_theta_e. A pinning |u|^2 term is forbidden by in-plane translations

## Candidate 4: `WL_BASIS_REDUNDANCY_UNDER_CONSTRAINT` ↔ `S11BB_BASIS_REDUNDANCY_UNDER_CONSTRAINT`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_BASIS_REDUNDANCY_UNDER_CONSTRAINT`

    No individual invariant is intrinsically zero.  For impermeable faces and om!=0, theta+e_W+div(u)=0 lets one eliminate theta, reducing the six local scalar bilinears to three and making the two theta-gradient bilinears k-dependent combinations (with higher-u derivatives if returned to position space).  With flux on, the one relation has frequency-, Z-, and kernel-dependent coefficients instead, so the same elimination choice is possible only at generic rank and gives different combinations.  At om=0 the impermeable row has rank 0 and no such redundancy until the conserved integration constant is fixed; with la0!=0 it instead becomes the equilibrium muTheta=0 row (away from a face pole).

### SymPy — `S11BB_BASIS_REDUNDANCY_UNDER_CONSTRAINT`

    impermeable, omega!=0, fixed zero integration constant: d+theta+e_W=0, so the six local quadratic forms in (d,theta,e_W) reduce to three. Flux-on: substitution contains the retarded face history and gives no redundancy among stored-energy field bilinears. At omega=0 the conserved integration constant must be fixed separately, so neither reduction is automatic

## Candidate 5: `WL_FACE_RESPONSE` ↔ `S11BB_FACE_RESPONSE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_FACE_RESPONSE`

    Solving p=Z(V+J/rhom), J=la(mu_s-p/rhom)+lv V gives p=P_V V+P_mu mu_s, with P_V=Z*rhom*(rhom+lv)/(rhom^2+Z*la), P_mu=Z*rhom*la/(rhom^2+Z*la), la=la0/(1-i*om*tauA), lv=lv0/(1-i*om*tauV).  This is a two-drive face response, not a pure impedance.  Direct symbolic solve check=True.

### SymPy — `S11BB_FACE_RESPONSE`

    delta_p=P_V V+P_mu mu_theta; P_V=omega*rho_m**2*(Lambda_V0/(-I*omega*tau_V + 1) + rho_m)/(q_out*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)); P_mu=Lambda_A0*omega*rho_m**2/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)); r_A and r_V remain distinct

## Candidate 6: `WL_FACE_RESPONSE_MU_COEFF` ↔ `S11BB_FACE_RESPONSE_MU_COEFF`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_FACE_RESPONSE_MU_COEFF`

    As a function of muTheta, p=P_V V+P_muTheta muTheta where P_muTheta=Z*rhom*la/(rhoBr*(rhom^2+Z*la)); equivalently P_mu=P_muTheta*rhoBr=Z*rhom*la/(rhom^2+Z*la) multiplies mu_s=muTheta/rhoBr.

### SymPy — `S11BB_FACE_RESPONSE_MU_COEFF`

    Lambda_A0*omega*rho_m**2/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) (coefficient of mu_theta; the response is therefore not a pure velocity impedance)

## Candidate 7: `WL_ZPERM_REDUCTION_CHECK` ↔ `S11BB_ZPERM_REDUCTION_CHECK`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ZPERM_REDUCTION_CHECK`

    The supplied affinity fixes g_p=-1/rhom and Lambda_p0=g_p*la0=-la0/rhom before solving.  With tauA=tauV=tau, r=1-i*om*tau, the derived p/V is rhom*(rhom*r+lv0)/(rhom*y*r+la0), equal to (rhom*r+lv0)/(y*r-Lambda_p0): True.  The general rA!=rV response is reported above and is not checked by this specialized acceptance value.

### SymPy — `S11BB_ZPERM_REDUCTION_CHECK`

    PASS; g_p=-1/rho_m and Lambda_p0=-Lambda_A0/rho_m; after tau_A=tau_V=tau, P_V-Z_perm=0; the unequal-time response is not checked by the supplied standard

## Candidate 8: `WL_TWO_PORT_POWER_IDENTITY` ↔ `S11BB_TWO_PORT_POWER_IDENTITY`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_TWO_PORT_POWER_IDENTITY`

    Per face, Pout=(1/2)Re[(p+lx*A)V*+mu_s J*]=(1/2)Re[p v_bulk*+A J*+lx A V*], v_bulk=V+J/rhom.  The algebraic residual after A=mu_s-p/rhom is 0; transferred-mass pressure work p J*/rhom and reciprocal-traction work lx A V* are both present.

### SymPy — `S11BB_TWO_PORT_POWER_IDENTITY`

    PASS; P_out=1/2 Re[(delta_p+Lambda_X*A)V*+mu_s J*]=1/2 Re[delta_p(V+J/rho_m)*+A J*+Lambda_X A V*]; transferred-mass pressure work delta_p J*/rho_m and reciprocal-traction work are both present

## Candidate 9: `WL_PORT_DISSIPATIVITY` ↔ `S11BB_PORT_DISSIPATIVITY`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_PORT_DISSIPATIVITY`

    For a=(V,mu_s)^T, L_port={{F_V,F_mu},{J_V,J_mu}} with F_V=P_V+lx*A_V, F_mu=P_mu+lx*A_mu, and H=(L_port+L_port^dagger)/2.  Necessary and sufficient: Re(F_V)>=0, Re(J_mu)>=0, Re(F_V)Re(J_mu)-|F_mu+Conjugate(J_V)|^2/4>=0.  H={{((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)) + Conjugate[(om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX))])/2, ((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)) + Conjugate[rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))]/(Conjugate[rhom]^2 + Conjugate[la0*om*rhom]/(Conjugate[qout]*(1 + I*Conjugate[om*tauA]))))/2}, {((rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV)))/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + Conjugate[(la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX))])/2, ((la0*rhom^2)/((1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (Conjugate[la0]*Conjugate[rhom]^2)/((Conjugate[rhom]^2 + Conjugate[la0*om*rhom]/(Conjugate[qout]*(1 + I*Conjugate[om*tauA])))*(1 + I*Conjugate[om*tauA])))/2}}.

### SymPy — `S11BB_PORT_DISSIPATIVITY`

    for a=(V,mu_s): H11=Re(T_V), H22=Re(rho_br0*J_mu), H12=(rho_br0*T_mu+conjugate(J_V))/2, H21=conjugate(H12); necessary-and-sufficient: H11>=0, H22>=0, det(H)>=0 (equivalently both eigenvalues nonnegative). T and J coefficients are the explicit face coefficients above

## Candidate 10: `WL_PORT_CONDITION_KIND` ↔ `S11BB_PORT_CONDITION_KIND`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_PORT_CONDITION_KIND`

    The port condition is not parameter-only: it also depends on real (om,k) through Z(om,k), including the outgoing/evanescent branch.

### SymPy — `S11BB_PORT_CONDITION_KIND`

    the principal-minor condition depends on coefficients and also on real (omega,k) through Z and q_out; it is not a parameter-only condition

## Candidate 11: `WL_COEFFICIENT_ADMISSIBILITY` ↔ `S11BB_COEFFICIENT_ADMISSIBILITY`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_COEFFICIENT_ADMISSIBILITY`

    The Hermitian interfacial form for arbitrary (A,V) is H_int=(L+L^dagger)/2 with L={{la,lv},{lx,0}}.  PSD is equivalent to Re(la)>=0 and lx=-Conjugate(lv) at every real om.  For real base coefficients and tauI>=0 this is necessary and sufficient: la0>=0 and either (lv0=lx0=0, tauV and tauX arbitrary) or (lx0=-lv0, tauV=tauX=0); tauA is arbitrary.

### SymPy — `S11BB_COEFFICIENT_ADMISSIBILITY`

    unconditional interfacial entropy condition for every real omega: Lambda_A0>=0 and [Lambda_V0=Lambda_X0=0 OR (Lambda_X0=-Lambda_V0 and tau_V=tau_X=0)]; with reciprocity the same admissible set results. Thus admissibility is a strict subset of the reciprocity region

## Candidate 12: `WL_ONSAGER_CONDITION` ↔ `S11BB_ONSAGER_CONDITION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ONSAGER_CONDITION`

    Formal d=eps partial inversion gives all-flux cross residual (lv+lx)/eps; clearing eps gives lv+lx=0.  Solving the first row gives all-force resistance cross residual -(lv+lx)/la; clearing la gives the same relation.  Route agreement=True.

### SymPy — `S11BB_ONSAGER_CONDITION`

    partial-inversion flux route and resistance route agree after clearing denominators: Lambda_V(omega)+Lambda_X(omega)=0; flux-route minus force-route=0

## Candidate 13: `WL_ONSAGER_RECIPROCITY` ↔ `S11BB_ONSAGER_RECIPROCITY`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ONSAGER_RECIPROCITY`

    Conditional Onsager-Casimir reciprocity at the same real om forces lx(om)=-lv(om), without complex conjugation.  For the supplied real Debye coefficients: lx0=-lv0 and, when lv0!=0, tauX=tauV; if lv0=lx0=0 their relaxation times are irrelevant.

### SymPy — `S11BB_ONSAGER_RECIPROCITY`

    conditional (time-reversal-symmetric) region: Lambda_X0=-Lambda_V0 and, if that common magnitude is nonzero, tau_X=tau_V; if both coefficients vanish their relaxation times are irrelevant. No conjugation is inserted

## Candidate 14: `WL_ONSAGER_DETERMINABLE` ↔ `S11BB_ONSAGER_DETERMINABLE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ONSAGER_DETERMINABLE`

    DETERMINABLE: both prescribed partial-inversion routes give lx=-lv; direct transposition of the mixed matrix was not used.

### SymPy — `S11BB_ONSAGER_DETERMINABLE`

    YES; the supplied partial-inversion conversion makes it determinate

## Candidate 15: `WL_RELAXATION_TIME_RELATIONS` ↔ `S11BB_RELAXATION_TIME_RELATIONS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_RELAXATION_TIME_RELATIONS`

    Unconditional entropy admissibility forces tauV=tauX=0 only when the cross pair is nonzero; it forces no relation on tauA.  Conditional reciprocity forces tauX=tauV for a nonzero cross pair and no relation involving tauA.  With both admissibility and reciprocity, a nonzero pair again requires tauV=tauX=0.

### SymPy — `S11BB_RELAXATION_TIME_RELATIONS`

    admissibility: tau_A unrestricted; if Lambda_V0 != 0, tau_V=tau_X=0, while for Lambda_V0=Lambda_X0=0 those times are irrelevant. Reciprocity alone: tau_V=tau_X when the cross pair is nonzero; tau_A unrestricted

## Candidate 16: `WL_GROWTH_INSIDE_ADMISSIBLE` ↔ `S11BB_GROWTH_INSIDE_ADMISSIBLE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_GROWTH_INSIDE_ADMISSIBLE`

    The explicit K0<0 growing root has la0=lv0=lx0=0 and hence is inside unconditional interfacial admissibility and conditional reciprocity (also inside their intersection).  A general root is inside exactly when la0>=0 and [lv0=lx0=0 or (lx0=-lv0 and tauV=tauX=0)], and inside reciprocity exactly when lx0=-lv0 with tauX=tauV for a nonzero pair; neither condition is used as a gate.

### SymPy — `S11BB_GROWTH_INSIDE_ADMISSIBLE`

    for any growing Omega_j, evaluate the above parameter-only interfacial condition and the root-frequency port minors H11,H22,det(H); membership is symbolic/parameter-dependent and is not used to remove the root

## Candidate 17: `WL_DECAY_INSIDE_ADMISSIBLE` ↔ `S11BB_DECAY_INSIDE_ADMISSIBLE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_DECAY_INSIDE_ADMISSIBLE`

    The explicit radiative decay roots also have the zero interface tuple and lie inside both regions.  General decay roots obey the same conditional membership tests as growing roots; radiation resistance is present even at the zero interface tuple.  Admissibility is strictly nested inside reciprocity, while admissibility intersect reciprocity equals admissibility for the supplied real Debye form.

### SymPy — `S11BB_DECAY_INSIDE_ADMISSIBLE`

    for any decaying Omega_j, evaluate the same unconditional and reciprocity-conditioned regions separately; membership is symbolic/parameter-dependent and is not a gate

## Candidate 18: `WL_CONSTRAINT` ↔ `S11BB_CONSTRAINT`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONSTRAINT`

    rhoBr*d_t(theta+e_W+div u)=-(J_++J_-); for exp[i(k.x-om t)] and the symmetric face sector: -i*om*rhoBr*(theta+e_W+d)+2J=0, d=i k.u_L.

### SymPy — `S11BB_CONSTRAINT`

    rho_br0*partial_t(theta+e_W+div u)=-(J_++J_-); Fourier: -i*omega*rho_br0*(theta+e_W+d)+J_++J_-=0

## Candidate 19: `WL_CONSTRAINT_TERM_ORIGINS` ↔ `S11BB_CONSTRAINT_TERM_ORIGINS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONSTRAINT_TERM_ORIGINS`

    d_t Sigma -> rhoBr*d_t(theta+e_W); div(Sigma v) -> rhoBr*div(d_t u) because the background is uniform, while products of perturbations are O(2); the exact right side -> -(J_++J_-) with no added term.  Moving it left gives +J_++J_-.

### SymPy — `S11BB_CONSTRAINT_TERM_ORIGINS`

    partial_t Sigma -> rho_br0 partial_t(theta+e_W); div(Sigma v) -> rho_br0 partial_t(div u) (background uniform; perturbation*velocity is O(2)); RHS -> -(J_++J_-); no other first-order terms

## Candidate 20: `WL_INTERNAL_DOF_COUNT` ↔ `S11BB_INTERNAL_DOF_COUNT`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_INTERNAL_DOF_COUNT`

    At fixed generic (k,om), before B1 the amplitudes are {u1,u2,u3,theta,e_W}: five fields.  One scalar balance row leaves four independent amplitudes, namely two transverse plus two scalar-sector amplitudes.  Equivalently {d,theta,e_W} has three amplitudes before and two after the generic rank-one constraint.  This counts field amplitudes, not independent initial data of the memory-augmented evolution.

### SymPy — `S11BB_INTERNAL_DOF_COUNT`

    at fixed nondegenerate (k,omega), one scalar balance relation reduces five field amplitudes to four: two transverse amplitudes plus two independent amplitudes in (u_L,theta,e_W). Memory kernels may require auxiliary initial data, but are not extra field amplitudes in this count

## Candidate 21: `WL_DOF_COUNTING_CONVENTION` ↔ `S11BB_DOF_COUNTING_CONVENTION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_DOF_COUNTING_CONVENTION`

    Counting complex Fourier amplitudes at fixed (k,om).  At om=0 the impermeable balance row loses rank and the conserved value of theta+e_W+d is initial/boundary data; with active affinity-driven transfer the static row generically becomes muTheta=0 instead.

### SymPy — `S11BB_DOF_COUNTING_CONVENTION`

    counted complex Fourier field amplitudes: before=(u_x,u_y,u_z,theta,e_W), five; after one nonzero-rank scalar constraint=four. At omega=0 this rank statement must include the stationary flux equation and a separately fixed conserved mass constant

## Candidate 22: `WL_INPLANE_EOM` ↔ `S11BB_INPLANE_EOM`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_INPLANE_EOM`

    Using delta_v theta=-delta_v e_W-div(delta_v u), with muTheta=delta U/delta theta at fixed (u,e_W) and pW=delta U/delta e_W at fixed (u,theta), rhoBr*u_tt=grad(directD-muTheta)-muR curl(curl u).  Thus the required restoring term -grad(muTheta) is present.  Longitudinal row: (rhoBr*om^2+k^2*piD)d+k^2*piTheta*theta+k^2*piE*e_W=0.

### SymPy — `S11BB_INPLANE_EOM`

    rho_br0*partial_t^2 u + mu_R curl(curl u)-grad[K_L div u+D_theta theta+D_e e_W]+grad(mu_theta)=0, with mu_theta=D_theta*d + e_W*(A_theta_e*k**2 + C*W0) + theta*(A_theta*k**2 + B_rho3); equivalently restoring force contains -grad(delta U/delta theta)

## Candidate 23: `WL_THICKNESS_EOM` ↔ `S11BB_THICKNESS_EOM`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_THICKNESS_EOM`

    The multiplier method was used; the multiplier is the material chemical/stress potential muTheta.  At fixed u,e_W the functional derivative is muTheta=(b3+gTheta*k^2)theta+(cc*w0+gThetaE*k^2)e_W+aTheta*d; at fixed u,theta, pW=(cc*w0+gThetaE*k^2)theta+(kW*w0^2+kapW*w0^4*k^2)e_W+aE*d.  The equation is muW*deltaW_tt+(pW-muTheta)/w0+(p_++p_-+lx(A_++A_-))/2=0; symmetric faces give the second assembly row above.

### SymPy — `S11BB_THICKNESS_EOM`

    -mu_W*omega^2*W0*e_W+(p_W-mu_theta)/W0+(delta_p+Lambda_X*A)=0 for the symmetric two-face amplitude, V_+=V_-=-i*omega*W0*e_W/2; row coefficients=Matrix([[D_theta*(Lambda_A0*omega*rho_m**2/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + Lambda_X0*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1)) + (D_e - D_theta)/W0, (A_theta*k**2 + B_rho3)*(Lambda_A0*omega*rho_m**2/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + Lambda_X0*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1)) + (-A_theta*k**2 + A_theta_e*k**2 - B_rho3 + C*W0)/W0, -W0*mu_W*omega**2 - I*W0*omega*(-Lambda_X0*omega*rho_m*(Lambda_V0/(-I*omega*tau_V + 1) + rho_m)/(q_out*(-I*omega*tau_X + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + omega*rho_m**2*(Lambda_V0/(-I*omega*tau_V + 1) + rho_m)/(q_out*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)))/2 + (A_theta_e*k**2 + C*W0)*(Lambda_A0*omega*rho_m**2/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + Lambda_X0*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1)) + (-A_theta_e*k**2 - C*W0 + W0**4*k**2*kappa_W + W0**2*k_W)/W0]])

## Candidate 24: `WL_BULK_FORCE_ON_THICKNESS` ↔ `S11BB_BULK_FORCE_ON_THICKNESS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_BULK_FORCE_ON_THICKNESS`

    The prescribed force is -(p_++p_-+lx(A_++A_-))/2 on the deltaW coordinate.  In the left-hand operator it is +F_V*(-i*om*deltaW/2)+F_mu*mu_s, where F_V=P_V+lx*A_V and F_mu=P_mu+lx*A_mu; both faces and the factor 1/2 have been retained.

### SymPy — `S11BB_BULK_FORCE_ON_THICKNESS`

    generalized external virtual work is -1/2 sum_s(delta_p_s+Lambda_X A_s) delta_v(deltaW), hence the symmetric thickness operator contains -i*omega*P_V/2 on deltaW and P_mu*mu_theta plus the analogous Lambda_X terms

## Candidate 25: `WL_RECIPROCAL_TRACTION_THICKNESS_EFFECT` ↔ `S11BB_RECIPROCAL_TRACTION_THICKNESS_EFFECT`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_RECIPROCAL_TRACTION_THICKNESS_EFFECT`

    The isolated symbolic row difference from lx0=0 is w0*lx*{A_mu*aTheta/rhoBr,A_mu*mTheta/rhoBr,A_mu*mCross/rhoBr-i*om*w0*A_V/2}, with A_mu=rhom^2/(rhom^2+Z*la), A_V=-Z(rhom+lv)/(rhom^2+Z*la).  Computed difference={-((aTheta*la0*om*rhom*w0)/(rhoBr*(la0*om + qout*rhom - I*om*qout*rhom*tauA))) + (aTheta*rhom*(I*la0*om + I*lx0*qout + lx0*om*qout*tauA + la0*om^2*tauX)*w0)/(rhoBr*(la0*om + qout*rhom - I*om*qout*rhom*tauA)*(I + om*tauX)), -(((b3 + gTheta*k^2)*la0*om*rhom*w0)/(rhoBr*(la0*om + qout*rhom - I*om*qout*rhom*tauA))) + ((b3 + gTheta*k^2)*rhom*(I*la0*om + I*lx0*qout + lx0*om*qout*tauA + la0*om^2*tauX)*w0)/(rhoBr*(la0*om + qout*rhom - I*om*qout*rhom*tauA)*(I + om*tauX)), ((I/2)*om^2*rhom*(I + om*tauA)*(I*lv0 + I*rhom + om*rhom*tauV)*w0^2)/((I*la0*om + I*qout*rhom + om*qout*rhom*tauA)*(I + om*tauV)) + (om^2*(I + om*tauA)*(lv0 + rhom - I*om*rhom*tauV)*((-I)*lx0 + I*rhom + om*rhom*tauX)*w0^2)/(2*(I*la0*om + I*qout*rhom + om*qout*rhom*tauA)*(I + om*tauV)*(I + om*tauX)) - (la0*om*rhom*w0*(gThetaE*k^2 + cc*w0))/(rhoBr*(la0*om + qout*rhom - I*om*qout*rhom*tauA)) + (rhom*(I*la0*om + I*lx0*qout + lx0*om*qout*tauA + la0*om^2*tauX)*w0*(gThetaE*k^2 + cc*w0))/(rhoBr*(la0*om + qout*rhom - I*om*qout*rhom*tauA)*(I + om*tauX))}.

### SymPy — `S11BB_RECIPROCAL_TRACTION_THICKNESS_EFFECT`

    Delta O_X acting on (d,theta,e_W)=Lambda_X*(A_M*mu_theta-i*omega*W0*A_V*e_W/2); exactly D-D|Lambda_X0=0=det(R1,R2,DeltaR3), DeltaR3=Matrix([[D_theta*Lambda_X0*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1), Lambda_X0*(A_theta*k**2 + B_rho3)*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1), I*Lambda_X0*W0*omega**2*rho_m*(Lambda_V0/(-I*omega*tau_V + 1) + rho_m)/(2*q_out*(-I*omega*tau_X + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + Lambda_X0*(A_theta_e*k**2 + C*W0)*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1)]])

## Candidate 26: `WL_THICKNESS_RESPONSE` ↔ `S11BB_THICKNESS_RESPONSE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_THICKNESS_RESPONSE`

    Introduce an added symmetric generalized face pressure fW on the right of the unmultiplied thickness equation and hold d=theta=0.  Then chiW=deltaW/fW=w0^2/T_e, where T_e=-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr and chiW=w0^2/(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr).

### SymPy — `S11BB_THICKNESS_RESPONSE`

    with F_W=-(R3_d*d+R3_theta*theta), output/input chi_W=deltaW/F_W=W0/R3_e; R3_e is the explicit third entry of S11BB_THICKNESS_EOM

## Candidate 27: `WL_RESPONSE_NORMALIZATION` ↔ `S11BB_RESPONSE_NORMALIZATION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_RESPONSE_NORMALIZATION`

    B3 output/input is deltaW (length) divided by an added face pressure fW; [chiW]=L^3 T^2 M^-1.  It is a susceptibility, not an effective inertia.

### SymPy — `S11BB_RESPONSE_NORMALIZATION`

    output/input=deltaW [L] divided by generalized thickness force density F_W [M L^-2 T^-2], so [chi_W]=[L^3 T^2 M^-1]

## Candidate 28: `WL_BULK_OPERATOR_BY_REGIME` ↔ `S11BB_BULK_OPERATOR_BY_REGIME`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_BULK_OPERATOR_BY_REGIME`

    For the general permeable face, the motion-driven term is (-i*om/2)F_V*deltaW: its real-time acceleration-phase coefficient is M_bulk=-Im(F_V)/(2*om), its velocity-phase coefficient is R_bulk=Re(F_V)/2, and F_mu*mu_s is a separate brane-state drive.  Impermeable check: q^2>0 has F_V=Z real and gives R_bulk=Z/2,M_bulk=0 (radiation resistance); q^2<0 has Z=-i*om*rhom/alpha and gives R_bulk=0,M_bulk=rhom/(2 alpha); q=0 makes Z and the impermeable operator singular.  Permeability can feedback-displace that singular behavior through rhom^2+Z*la.

### SymPy — `S11BB_BULK_OPERATOR_BY_REGIME`

    for real omega the pressure part is O_p=-i*omega*P_V/2=-omega^2 M_p-i*omega Gamma_p with M_p=-Im(P_V)/(2 omega), Gamma_p=Re(P_V)/2. Impermeable propagating: M_p=0,Gamma_p=Z/2>0. Impermeable evanescent: M_p=rho_m/(2 alpha),Gamma_p=0. Permeation generally makes both nonzero in either regime; reciprocal traction adds the same decomposition with P_V replaced by Lambda_X*A_V. q_out=0 is singular

## Candidate 29: `WL_MASS_INTERPRETATION_VALID_WHERE` ↔ `S11BB_MASS_INTERPRETATION_VALID_WHERE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_MASS_INTERPRETATION_VALID_WHERE`

    A mass interpretation is valid only for the purely acceleration-phase evanescent impermeable loading (+rhom/(2 alpha) on the deltaW coordinate), or for a separately isolated -Im(F_V)/(2 om) component.  It is invalid for propagating radiation resistance and for the independent F_mu*mu_s drive.

### SymPy — `S11BB_MASS_INTERPRETATION_VALID_WHERE`

    only a purely acceleration-phase load admits a mass interpretation: the impermeable evanescent load does (rho_m/(2 alpha) for deltaW); propagating radiation resistance and a generic permeable/reciprocal load do not

## Candidate 30: `WL_COMPRESSIONAL_RESPONSE` ↔ `S11BB_COMPRESSIONAL_RESPONSE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_COMPRESSIONAL_RESPONSE`

    Define compressive strain epsilon_c=-d=-div(u) and integrated compressional pressure Pi=muTheta-directD, so K_L=Pi/epsilon_c.  Eliminating theta,e_W with the mass and thickness rows gives theta/d=((aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) + (I*om*rhoBr - (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr))/(-((-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))) + ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr)), e_W/d=((I*om*rhoBr - (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr) + ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr))/(-((-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))) + ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr)), and K_L=-aTheta + kD - ((-aE + gThetaE*k^2 + cc*w0)*((I*om*rhoBr - (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr) + ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)))/(-((-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))) + ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr)) - ((-aTheta + b3 + gTheta*k^2)*((aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) + (I*om*rhoBr - (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr)))/(-((-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))) + ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr)).  C enters through mCross=cc*w0+gThetaE*k^2; kapW enters through mE=kW*w0^2+kapW*w0^4*k^2.

### SymPy — `S11BB_COMPRESSIONAL_RESPONSE`

    K_L^eff= sigma_L/d, sigma_L=K_L*d+D_theta*theta+D_e*e_W-mu_theta. Let Delta=R2_theta*R3_e-R2_e*R3_theta; theta/d=(-R2_d*R3_e+R2_e*R3_d)/Delta, e_W/d=(-R2_theta*R3_d+R2_d*R3_theta)/Delta, and K_L^eff=(K_L-D_theta)+(D_theta-b_theta)*theta/d+(D_e-b_theta_e)*e_W/d

## Candidate 31: `WL_LIMITS_AND_PATH` ↔ `S11BB_LIMITS_AND_PATH`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_LIMITS_AND_PATH`

    At exactly om=0 the impermeable B1 row vanishes and theta+e_W+d is a conserved integration constant fixed by initial total slab mass or boundary preparation; dividing by om first silently sets that constant.  For fixed real k>0 and om->0 from the prescribed upper rim, Z~-i*om*rhom/k: if la0!=0 the static row is muTheta=0, whereas if la0=0 the om!=0 limit is theta+e_W+d=0.  At fixed k and |om|->infinity on the continued sheet, qout~om/cs and Z~rhom*cs, while each channel with tauI>0 scales as 1/om (a tauI=0 channel does not).  The path k->0 first has qout=om/cs and Z=rhom*cs already at finite om, so it differs from the fixed-k static path; the limits do not generally commute.

### SymPy — `S11BB_LIMITS_AND_PATH`

    generic active slab affinity (Lambda_A0!=0), fixed real k!=0 and omega->0 from the prescribed upper rim: mu_theta=p_W=0 and K_eff->-D_e*(D_e*(A_theta*k**2 + B_rho3)/(A_theta*W0**4*k**4*kappa_W + A_theta*W0**2*k**2*k_W - A_theta_e**2*k**4 - 2*A_theta_e*C*W0*k**2 + B_rho3*W0**4*k**2*kappa_W + B_rho3*W0**2*k_W - C**2*W0**2) + D_theta*(-A_theta_e*k**2 - C*W0)/(A_theta*W0**4*k**4*kappa_W + A_theta*W0**2*k**2*k_W - A_theta_e**2*k**4 - 2*A_theta_e*C*W0*k**2 + B_rho3*W0**4*k**2*kappa_W + B_rho3*W0**2*k_W - C**2*W0**2)) - D_theta*(D_e*(-A_theta_e*k**2 - C*W0)/(A_theta*W0**4*k**4*kappa_W + A_theta*W0**2*k**2*k_W - A_theta_e**2*k**4 - 2*A_theta_e*C*W0*k**2 + B_rho3*W0**4*k**2*kappa_W + B_rho3*W0**2*k_W - C**2*W0**2) + D_theta*(W0**4*k**2*kappa_W + W0**2*k_W)/(A_theta*W0**4*k**4*kappa_W + A_theta*W0**2*k**2*k_W - A_theta_e**2*k**4 - 2*A_theta_e*C*W0*k**2 + B_rho3*W0**4*k**2*kappa_W + B_rho3*W0**2*k_W - C**2*W0**2)) + K_L. Fixed k with |omega|->infinity (mu_W!=0, finite Debye parameters): e_W->0, theta->-d and K_eff->A_theta*k**2 + B_rho3 - 2*D_theta + K_L. Taking impermeability before omega->0 instead preserves theta+e_W+d=constant and generally differs; limits/cuts do not commute

## Candidate 32: `WL_FROZEN_THICKNESS_IDENTIFICATION` ↔ `S11BB_FROZEN_THICKNESS_IDENTIFICATION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_FROZEN_THICKNESS_IDENTIFICATION`

    A thickness-held-fixed modulus is the form-control response obtained from e_W=V=0, not a new input: theta/d=-massD/massT (with the surviving mu_s-driven J and p), K_frozen=-(piD+piTheta*theta/d).  No unique static number exists until the om=0 integration constant or transfer equilibrium and the approach path are specified.

### SymPy — `S11BB_FROZEN_THICKNESS_IDENTIFICATION`

    holding e_W=0 gives K_frozen=-D_theta + K_L - (2*D_theta*Lambda_A0*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_A + 1) - I*omega*rho_br0)*(-A_theta*k**2 - B_rho3 + D_theta)/(2*Lambda_A0*(A_theta*k**2 + B_rho3)*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_A + 1) - I*omega*rho_br0) after the still-active mass/transfer row is solved. It is not a primitive single compression modulus; in the impermeable high-frequency reduction it becomes K_L-2D_theta+B_rho3+A_theta*k^2

## Candidate 33: `WL_LONGITUDINAL_DISPERSION` ↔ `S11BB_LONGITUDINAL_DISPERSION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_LONGITUDINAL_DISPERSION`

    On the prescribed qout sheet the dispersion is det(M)=0 for M={{R_d,R_theta,R_e},{T_d,T_theta,T_e},{B_d,B_theta,B_e}}, with rows printed in B2/B1 and det=R_d(T_theta B_e-T_e B_theta)-R_theta(T_d B_e-T_e B_d)+R_e(T_d B_theta-T_theta B_d).  The computed determinant expression is k^2*(-aE + gThetaE*k^2 + cc*w0)*(((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr) - ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)) - k^2*(-aTheta + b3 + gTheta*k^2)*((aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) - ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr)) + (k^2*(aTheta - kD) + om^2*rhoBr)*((-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) - ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr)).

### SymPy — `S11BB_LONGITUDINAL_DISPERSION`

    D(omega,k;q_out)=det(R_inplane,R_mass,R_thickness)=0 with unscaled rows R1=Matrix([[k**2*(-D_theta + K_L) - omega**2*rho_br0, k**2*(-A_theta*k**2 - B_rho3 + D_theta), k**2*(-A_theta_e*k**2 - C*W0 + D_e)]]), R2=Matrix([[2*D_theta*Lambda_A0*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_A + 1) - I*omega*rho_br0, 2*Lambda_A0*(A_theta*k**2 + B_rho3)*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_A + 1) - I*omega*rho_br0, 2*Lambda_A0*(A_theta_e*k**2 + C*W0)*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_A + 1) - I*W0*omega*(-Lambda_A0*omega*rho_m*(Lambda_V0/(-I*omega*tau_V + 1) + rho_m)/(q_out*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + Lambda_V0/(-I*omega*tau_V + 1)) - I*omega*rho_br0]]), R3=Matrix([[D_theta*(Lambda_A0*omega*rho_m**2/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + Lambda_X0*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1)) + (D_e - D_theta)/W0, (A_theta*k**2 + B_rho3)*(Lambda_A0*omega*rho_m**2/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + Lambda_X0*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1)) + (-A_theta*k**2 + A_theta_e*k**2 - B_rho3 + C*W0)/W0, -W0*mu_W*omega**2 - I*W0*omega*(-Lambda_X0*omega*rho_m*(Lambda_V0/(-I*omega*tau_V + 1) + rho_m)/(q_out*(-I*omega*tau_X + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + omega*rho_m**2*(Lambda_V0/(-I*omega*tau_V + 1) + rho_m)/(q_out*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)))/2 + (A_theta_e*k**2 + C*W0)*(Lambda_A0*omega*rho_m**2/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + Lambda_X0*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1)) + (-A_theta_e*k**2 - C*W0 + W0**4*k**2*kappa_W + W0**2*k_W)/W0]]); explicitly D=R11*(R22*R33-R23*R32)-R12*(R21*R33-R23*R31)+R13*(R21*R32-R22*R31)

## Candidate 34: `WL_ROOTS` ↔ `S11BB_ROOTS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ROOTS`

    General roots are the sheet-filtered zeros {om_j(k): det(M(om,k,qout(om,k)))=0}.  Because qout has square-root monodromy and the determinant contains three independent Debye denominators plus feedback denominators, no parameter-independent radical closed form or fixed root count/multiplicity exists; collisions, cancellations and roots at infinity occur on coefficient subvarieties.  Explicit sign classification therefore requires a parameter region.  A fully symbolic soluble slice is reported next and demonstrates both admitted signs without imposing positivity.

### SymPy — `S11BB_ROOTS`

    From the assembled determinant, substitution of k=0, Lambda_A0=Lambda_V0=Lambda_X0=0, and q_out=omega/c_s0 gives D_slice=(2*I*B_rho3*omega**3*rho_br0**2 - 4*I*C*W0*omega**3*rho_br0**2 + W0**2*c_s0*omega**4*rho_br0**2*rho_m + 2*I*W0**2*k_W*omega**3*rho_br0**2 - 2*I*W0**2*mu_W*omega**5*rho_br0**2)/(2*W0). The cancelled overall factor is -I*W0*omega**3*rho_br0**2 (including the omega^3 rank-loss factor); the surviving polynomial is (-2*B_rho3 + 4*C*W0 + I*W0**2*c_s0*omega*rho_m - 2*W0**2*k_W + 2*W0**2*mu_W*omega**2)/(2*W0**2), the previously-solved quadratic is I*c_s0*omega*rho_m/2 + mu_W*omega**2 - (B_rho3 - 2*C*W0 + W0**2*k_W)/W0**2, and their symbolic difference is 0. Solving the surviving derived polynomial gives roots=[(-I*W0*c_s0*rho_m - sqrt(16*B_rho3*mu_W - 32*C*W0*mu_W - W0**2*c_s0**2*rho_m**2 + 16*W0**2*k_W*mu_W))/(4*W0*mu_W), (-I*W0*c_s0*rho_m + sqrt(16*B_rho3*mu_W - 32*C*W0*mu_W - W0**2*c_s0**2*rho_m**2 + 16*W0**2*k_W*mu_W))/(4*W0*mu_W)]. The general dispersion remains the sheet-filtered determinant above and has no universal elementary closed form; exceptional multiplicities still satisfy D=partial_omega D=0

## Candidate 35: `WL_ROOT_STABILITY_CLASS` ↔ `S11BB_ROOT_STABILITY_CLASS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ROOT_STABILITY_CLASS`

    From the assembled determinant, substitution of k=0, la0=lv0=lx0=0, and qout=om/cs gives D_slice=(om^2*rhoBr*((2*I)*b3*om*rhoBr - (4*I)*cc*om*rhoBr*w0 + (2*I)*kW*om*rhoBr*w0^2 - (2*I)*muW*om^3*rhoBr*w0^2 + cs*om^2*rhoBr*rhom*w0^2))/2.  The cancelled overall factor is (-I)*om^3*rhoBr^2*w0^2 (including the om^3 rank-loss factor); the surviving polynomial is -kW + muW*om^2 + (I/2)*cs*om*rhom - (b3 - 2*cc*w0)/w0^2, the previously-reported quadratic is muW*om^2 + (I/2)*cs*om*rhom - (b3 - 2*cc*w0 + kW*w0^2)/w0^2, and their symbolic difference is 0.  Solving the surviving derived polynomial gives roots={((-1/4*I)*(cs*rhom - (I*Sqrt[16*b3*muW - 32*cc*muW*w0 + 16*kW*muW*w0^2 - cs^2*rhom^2*w0^2])/w0))/muW, ((-1/4*I)*(cs*rhom + (I*Sqrt[16*b3*muW - 32*cc*muW*w0 + 16*kW*muW*w0^2 - cs^2*rhom^2*w0^2])/w0))/muW}.  K0<0 gives one growing and one decaying root; K0=0 gives one static and one decaying root; K0>0 gives two decaying roots (overdamped or a damped oscillatory pair).

### SymPy — `S11BB_ROOT_STABILITY_CLASS`

    On the explicit slice omega_- is DECAYING for every K0. omega_+ is GROWING for K0<0, static at K0=0, and DECAYING for K0>0; for K0>0 the two decays are overdamped or a damped oscillatory pair. The separator is K0=B_rho3-2*C*W0+k_W*W0^2=0, equivalently C=(B_rho3+k_W*W0^2)/(2*W0): above this C omega_+ grows, below it both roots decay. On q_out=omega/c_s0 the growing root is a normal mode and each decaying root is an outgoing resonance; none is discarded

## Candidate 36: `WL_STABILITY_CONDITION` ↔ `S11BB_STABILITY_CONDITION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_STABILITY_CONDITION`

    On the soluble radiating slice the exact modulus separator is K0=b3-2 cc*w0+kW*w0^2=0: K0<0 produces one growing and one decaying root, while K0>0 produces decay only.  With general permeation/reciprocal traction there is no condition on moduli and C alone: the sign also depends on {la0,lv0,lx0,tauA,tauV,tauX}, k, and the chosen sheet through the displayed determinant.

### SymPy — `S11BB_STABILITY_CONDITION`

    The actual boundary on the explicit constrained slice is K0=B_rho3-2*C*W0+k_W*W0^2=0, because the k=0 in-plane row sets d=0 and mass balance sets theta=-e_W; for mu_W>0, K0<0 gives growth and K0>0 gives decay only. By contrast lambda_min(H_scalar(k))=0 for H_scalar=Matrix([
    [    K_L,               D_theta,                            D_e],
    [D_theta, A_theta*k**2 + B_rho3,          A_theta_e*k**2 + C*W0],
    [    D_e, A_theta_e*k**2 + C*W0, W0**4*k**2*kappa_W + W0**2*k_W]]) is the stronger unconstrained three-field Hessian boundary: its positive-definite side requires K_L>0, K_L*(B_rho3+A_theta*k^2)-D_theta^2>0, and det(H_scalar)>0. That is sufficient unconstrained positivity, not the slice stability boundary. The full free-interface boundary still depends on the interface coefficients and is D=0 with Im(omega)=0

## Candidate 37: `WL_IMAGINARY_PART` ↔ `S11BB_IMAGINARY_PART`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_IMAGINARY_PART`

    On that slice Im(om_+)=(-R+sqrt(R^2-4 muW K0/w0^2))/(2 muW) and Im(om_-)=(-R-sqrt(R^2-4 muW K0/w0^2))/(2 muW) when the radicand is nonnegative; when it is negative both imaginary parts are -R/(2 muW).  On the opposite sheet R->-R: Im(om_opp,+)=(R+sqrt(R^2-4 muW K0/w0^2))/(2 muW), Im(om_opp,-)=(R-sqrt(R^2-4 muW K0/w0^2))/(2 muW) for nonnegative radicand.  The branchwise ratios are (R+S)/(S-R) and (R-S)/(-R-S) where denominators are nonzero; for the underdamped K0>0 pair the ratio is -1.

### SymPy — `S11BB_IMAGINARY_PART`

    On the explicit slice let S=sqrt(R^2-4*mu_W*K0/W0^2). For nonnegative radicand, Im(omega_+)=(S-R)/(2*mu_W) and Im(omega_-)=-(R+S)/(2*mu_W): K0<0 gives signs (+,-), K0=0 gives (0,-), and 0<K0<=R^2*W0^2/(4*mu_W) gives (-,-). Above that threshold the roots are a damped oscillatory pair and both imaginary parts are -R/(2*mu_W)<0. On the opposite sheet the concrete roots are omega_opp,+=I*(c_s0*rho_m/2 + sqrt(c_s0**2*rho_m**2/4 - 4*mu_W*(B_rho3 - 2*C*W0 + W0**2*k_W)/W0**2))/(2*mu_W) and omega_opp,-=I*(c_s0*rho_m/2 - sqrt(c_s0**2*rho_m**2/4 - 4*mu_W*(B_rho3 - 2*C*W0 + W0**2*k_W)/W0**2))/(2*mu_W); for nonnegative radicand their imaginary parts are (R+S)/(2*mu_W) and (R-S)/(2*mu_W), with signs (+,-) for K0<0, (+,0) for K0=0, and (+,+) for K0>0, and their ratios to omega_+,omega_- are (R+S)/(S-R) and (R-S)/(-R-S) where defined. For the underdamped pair both opposite-sheet imaginary parts are R/(2*mu_W)>0 and both ratios are -1

## Candidate 38: `WL_DISSIPATION_ORIGIN` ↔ `S11BB_DISSIPATION_ORIGIN`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_DISSIPATION_ORIGIN`

    In the soluble slice every nonzero imaginary part for K0>0 comes from propagating bulk radiation resistance Re(Z), which survives impermeability; for K0<0 the positive root additionally comes from the indefinite conservative stiffness K0<0, while radiation shifts both roots toward decay.  General face-transfer dissipation/source behavior is separately carried by la,lv and lx and is not identifiable with radiation resistance.

### SymPy — `S11BB_DISSIPATION_ORIGIN`

    On the explicit slice all interface channels vanish. For K0>0 the negative imaginary parts come from propagating bulk radiation resistance Re(Z)=rho_m*c_s0; for K0<0 the growing/decaying pair is created by the indefinite constrained stiffness K0<0 and radiation shifts both toward decay. No slice root is made non-real by mass transfer or reciprocal traction. In the general determinant, interfacial transfer/reciprocal kernels can be a separate source or sink, distinguished by the impermeable and no-reciprocal form cuts

## Candidate 39: `WL_RECIPROCAL_TRACTION_ROOT_EFFECT` ↔ `S11BB_RECIPROCAL_TRACTION_ROOT_EFFECT`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_RECIPROCAL_TRACTION_ROOT_EFFECT`

    The full roots are zeros of det(M) with fm=pMu+lx*aMu and fv=pV+lx*aV.  The lx0=0 roots are zeros of the same determinant after lx0->0; their difference is nonzero generically and is exactly induced by the thickness-row difference printed above.  Root number and multiplicity can change only at cancellations/collisions or infinity; no universal sign change follows without parameter choices.

### SymPy — `S11BB_RECIPROCAL_TRACTION_ROOT_EFFECT`

    D-D_F=det(R1,R2,DeltaR3), DeltaR3=Matrix([[D_theta*Lambda_X0*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1), Lambda_X0*(A_theta*k**2 + B_rho3)*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1), I*Lambda_X0*W0*omega**2*rho_m*(Lambda_V0/(-I*omega*tau_V + 1) + rho_m)/(2*q_out*(-I*omega*tau_X + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + Lambda_X0*(A_theta_e*k**2 + C*W0)*(-Lambda_A0*omega*rho_m/(q_out*rho_br0*(-I*omega*tau_A + 1)*(Lambda_A0*omega*rho_m/(q_out*(-I*omega*tau_A + 1)) + rho_m**2)) + 1/rho_br0)/(-I*omega*tau_X + 1)]]); Lambda_X generically shifts roots and can split/merge them, while unchanged roots satisfy this difference=0. Multiplicity changes only on D=partial_omega D=0; the F cut sets tau_X irrelevant

## Candidate 40: `WL_BRANCH_REALAXIS_CHECK` ↔ `S11BB_BRANCH_REALAXIS_CHECK`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_BRANCH_REALAXIS_CHECK`

    qout is represented on the upper rim by sqrt(om-cs|k|)*sqrt(om+cs|k|)/cs with principal factor roots.  For real |om|>cs|k| this is sign(om)*sqrt(om^2/cs^2-k^2), so qout/om>0 and the energy flux is outward for both signs of om.  For |om|<cs|k| it is +i*sqrt(k^2-om^2/cs^2), so the field decays with |w|.  Thus requirements 1 and 2 are reproduced.

### SymPy — `S11BB_BRANCH_REALAXIS_CHECK`

    PASS; q_out=sgn(omega)*sqrt(omega^2/c_s0^2-k^2) for q^2>0 gives Re(Z)>0 for both signs of omega; q_out=i*alpha (alpha>0) for q^2<0 gives decay and Z=-i*omega*rho_m/alpha

## Candidate 41: `WL_BRANCH_DEGENERATE_POINT` ↔ `S11BB_BRANCH_DEGENERATE_POINT`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_BRANCH_DEGENERATE_POINT`

    At om=+/-cs|k|, qout=0; the two exponential bulk solutions coalesce and the second independent solution is linear in w.  Neither flux nor decay selects a solution there; continuity supplies the value and Z is singular.

### SymPy — `S11BB_BRANCH_DEGENERATE_POINT`

    omega=+/-c_s0*|k|: q_out=0, Z is singular, the two exponential bulk solutions coalesce and the second solution is linear in w; neither flux nor decay selects it, so continuity supplies the value

## Candidate 42: `WL_BRANCH_SENSITIVITY` ↔ `S11BB_BRANCH_SENSITIVITY`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_BRANCH_SENSITIVITY`

    Every general result uses downward continuation from the upper rim (and upward from that same rim for Im om>0); requirements 1-2 were not re-imposed at complex om.  The opposite-sheet determinant is obtained mechanically by qout->-qout.  On the soluble slice its imaginary parts and ratios are printed under IMAGINARY_PART, making the dependence measurable.  If a trajectory crosses Re(om)=+/-cs|k| it has left this sheet and must not be reselected.

### SymPy — `S11BB_BRANCH_SENSITIVITY`

    physical-sheet roots Omega_j^(+) solve D(Omega,k,q_out)=0; opposite-sheet roots Omega_j^(-) solve D(Omega,k,-q_out)=0; Im-parts are Im(Omega_j^(+)) and Im(Omega_j^(-)), with measurable ratio Im(Omega_j^(-))/Im(Omega_j^(+)) wherever the denominator is nonzero; free parameters prevent a universal numeric ratio

## Candidate 43: `WL_SHEET_OF_EACH_ROOT` ↔ `S11BB_SHEET_OF_EACH_ROOT`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_SHEET_OF_EACH_ROOT`

    General root om_j: qout=qout_upper-rim-continuation and it is a normal mode or resonance according to the resulting spatial behavior, never by reselection.  On the k=0 slice qout=om/cs: Im om>0 decays spatially and is a normal mode; Im om<0 grows spatially and is an outgoing leaky resonance.  The opposite objects use qout=-om/cs.  No real-axis requirement was re-imposed at any complex root.

### SymPy — `S11BB_SHEET_OF_EACH_ROOT`

    each root is labelled (Omega_j^(+),q_out continued from the upper rim along the prescribed vertical path); no real-axis decay/flux condition is re-imposed at complex Omega. It is a normal mode iff the resulting continued bulk factors decay, otherwise a resonance. Crossing Re(Omega)=+/-c_s0|k| leaves this sheet rather than triggering reselection

## Candidate 44: `WL_KERNEL_ORIENTATION_IDENTITIES` ↔ `S11BB_KERNEL_ORIENTATION_IDENTITIES`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_KERNEL_ORIENTATION_IDENTITIES`

    Primitive coefficient residuals {K_A-la,K_V-lv,K_X-lx}={0, 0, 0}; rational zero tests={True, True, True}.  Coverage is decisive only when the corresponding coefficient is nonzero and tauI>0; tauI=0 and coefficient=0 cases are correctly reported as orientation-indistinguishable.

### SymPy — `S11BB_KERNEL_ORIENTATION_IDENTITIES`

    A,V,X normalized numerators=(0,0,0) for K_I-Lambda_I before specialization; an active tau_I>0 channel has its bare retarded pole at omega=-i/tau_I

## Candidate 45: `WL_KERNEL_PROPAGATION_RESIDUALS` ↔ `S11BB_KERNEL_PROPAGATION_RESIDUALS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_KERNEL_PROPAGATION_RESIDUALS`

    Carrying inert {ellA,ellV,ellX} through face elimination gives residuals for {P_V,P_mu,A_V,A_mu,J_V,J_mu,F_V,F_mu}={0, 0, 0, 0, 0, 0, 0, 0}; assembly entry zero tests={{True, True, True}, {True, True, True}, {True, True, True}}; determinant residual under identical row normalization is zero=True.  No conjugation or om->-om was used.

### SymPy — `S11BB_KERNEL_PROPAGATION_RESIDUALS`

    face(pV,pMu,jV,jMu,tV,tMu), mass row, thickness row, determinant residuals=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

## Candidate 46: `WL_KERNEL_POLE_LOCATIONS` ↔ `S11BB_KERNEL_POLE_LOCATIONS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_KERNEL_POLE_LOCATIONS`

    Bare retarded poles are om=-i/tauA,-i/tauV,-i/tauX for active positive tauI, verified by the required pole test 1/(1/rI)==0 at rI=0.  In the face response the A bare pole cancels and is feedback-displaced to zeros of rhom^2*(1-i om tauA)+Z*la0; the V pole is retained in V-driven coefficients unless its numerator/channel vanishes; the X pole is retained in F=p+lx*A unless its channel or A vanishes.  The assembled equations inherit those retained/displaced poles; removable coefficient-zero factors are cancelled.  No upper-half-plane bare pole occurs.  Feedback-displaced half-plane is sign-indeterminate and is not an orientation verdict.

### SymPy — `S11BB_KERNEL_POLE_LOCATIONS`

    after cancellation: A bare pole -i/tau_A is feedback-displaced to roots of rho_m^2*(1-i*omega*tau_A)+Z*Lambda_A0=0; V bare pole -i/tau_V is generically retained in velocity-driven face/row objects; X bare pole -i/tau_X is generically retained in reciprocal traction/thickness objects. They cancel when the corresponding residue/channel vanishes. All bare retained poles are lower-half-plane for tau_I>0; feedback-displaced half-plane is sign-indeterminate and is not an orientation verdict. Dispersion zeros and Hermitian-form conjugates are excluded

## Candidate 47: `WL_CAUSALITY_CHECK` ↔ `S11BB_CAUSALITY_CHECK`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CAUSALITY_CHECK`

    PASS for every active memory channel: primitive rational identities and placeholder propagation residuals vanish, and every retained bare pole is at -i/tauI.  Explicit 1/expr==0 pole tests at {-i/tauA,-i/tauV,-i/tauX}={True, True, True}.  Absent or tauI=0 channels are indistinguishable rather than claimed as covered.

### SymPy — `S11BB_CAUSALITY_CHECK`

    PASS for all active nondegenerate channels: both the primitive rational identities and inert-placeholder propagation checks vanish; zero-time and absent-channel cases are indistinguishable rather than orientation-tested

## Candidate 48: `WL_GROWTH_ARTIFACT_DIAGNOSTICS` ↔ `S11BB_GROWTH_ARTIFACT_DIAGNOSTICS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_GROWTH_ARTIFACT_DIAGNOSTICS`

    For the explicit K0<0 growing slice root all interface memory channels are absent (orientation tests therefore indistinguishable, with no finite response-kernel poles); the pressure traction sign and full two-port checks below both vanish.  It lies on qout=om/cs, requirements 1-2 were not re-imposed, is spatially normalizable, and the zero interface tuple lies inside both admissibility and reciprocity regions.  Growth is therefore the indefinite K0<0 conservative direction, not a transposed kernel or sheet reselection artifact.

### SymPy — `S11BB_GROWTH_ARTIFACT_DIAGNOSTICS`

    For the concrete omega_+ root when K0<0: Lambda_A0=Lambda_V0=Lambda_X0=0, so all memory channels are absent, orientation is indistinguishable, all kernel-propagation residuals are 0, and there are no finite response-kernel poles. It lies on q_out=omega/c_s0 reached by prescribed upward continuation; real-axis flux/decay requirements were not re-imposed, and its bulk field is spatially normalizable. The zero interface tuple lies inside both unconditional admissibility and conditional reciprocity (the former is a strict subset of the latter). Thus its growth is the indefinite K0<0 constrained direction, not a kernel-orientation or sheet-reselection artifact

## Candidate 49: `WL_DECAY_ARTIFACT_DIAGNOSTICS` ↔ `S11BB_DECAY_ARTIFACT_DIAGNOSTICS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_DECAY_ARTIFACT_DIAGNOSTICS`

    For every explicit slice decay root the same absent-channel orientation/pole inventory applies and both balance checks vanish.  It lies on qout=om/cs without complex-frequency reselection and is an outgoing resonance; the zero interface tuple lies inside both regions.  Its decay shift is propagating radiation resistance, not mass transfer.  General non-real roots inherit the symbolic diagnostics printed above and are classified conditionally by their coefficient tuple.

### SymPy — `S11BB_DECAY_ARTIFACT_DIAGNOSTICS`

    For the concrete omega_- root for every K0, and omega_+ when K0>0: all interface channels are absent, orientation is indistinguishable, every kernel-propagation residual is 0, and there are no finite response-kernel poles. Each lies on q_out=omega/c_s0 reached by prescribed downward continuation without re-imposing spatial decay and is an outgoing resonance. The zero interface tuple lies inside both admissibility regions. Radiation resistance supplies the decay shift; for K0<0 omega_- is also the decaying member of the indefinite conservative pair, not an artifact

## Candidate 50: `WL_CONVENTION_CHECK_INPLANE` ↔ `S11BB_CONVENTION_CHECK_INPLANE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONVENTION_CHECK_INPLANE`

    PASS: eliminating delta_v theta with delta_v theta=-delta_v e_W-div(delta_v u) produces rhoBr*u_tt=grad(directD-muTheta)-muR curl curl u, so the restoring force -grad(delta U/delta theta) is present.  The common multiplier is muTheta; coefficient check=True.

### SymPy — `S11BB_CONVENTION_CHECK_INPLANE`

    PASS; eliminating delta_v theta=-delta_v e_W-div(delta_v u) gives constrained derivative deltaU/deltau|constraint=deltaU/deltau|theta,e+grad(mu_theta), so momentum carries restoring force -grad(mu_theta). Used elimination; the equivalent Lagrange multiplier is -mu_theta (the material chemical-potential/generalized pressure multiplier)

## Candidate 51: `WL_CONVENTION_CHECK_CONSERVATIVE` ↔ `S11BB_CONVENTION_CHECK_CONSERVATIVE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONVENTION_CHECK_CONSERVATIVE`

    With no bulk/interface, kapW=0,k=0 and theta+e_W=constant, K_check=d^2 U/de_W^2=b3-2 cc*w0+kW*w0^2.  The thickness equation gives om^2=K_check/(muW*w0^2), exactly equal to the independently reduced energy stiffness: True.  b3 appears explicitly.

### SymPy — `S11BB_CONVENTION_CHECK_CONSERVATIVE`

    PASS; no bulk/interfaces, kappa_W=0,k=0,J=0 gives theta+e_W=constant and omega^2=(B_rho3 - 2*C*W0 + W0**2*k_W)/(W0**2*mu_W); equation stiffness equals K_check=B_rho3 - 2*C*W0 + W0**2*k_W; B_rho3 does appear

## Candidate 52: `WL_CONSERVATIVE_POSITIVITY_INEQUALITY` ↔ `S11BB_CONSERVATIVE_POSITIVITY_INEQUALITY`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONSERVATIVE_POSITIVITY_INEQUALITY`

    For muW>0 the thickness om^2 is positive iff b3-2 cc*w0+kW*w0^2>0.  Positive definiteness of the (theta,e_W) Hessian, b3>0, kW*w0^2>0, b3*kW*w0^2-(cc*w0)^2>0, implies this inequality because it is the quadratic form on (-1,1); hence the check passes for every positive-definite supplied local U.

### SymPy — `S11BB_CONSERVATIVE_POSITIVITY_INEQUALITY`

    omega^2>0 iff (B_rho3-2*C*W0+k_W*W0^2)/mu_W>0. If mu_W>0 and the (theta,e_W) stored-energy Hessian is positive definite [B_rho3>0 and B_rho3*k_W>C^2], this inequality holds for every such U

## Candidate 53: `WL_ENERGY_SINKS` ↔ `S11BB_ENERGY_SINKS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ENERGY_SINKS`

    d(T+U)/dt=-Sum_s[(p_s+lx A_s)V_s+mu_s J_s] (real instantaneous pairing; take one-half real part harmonically).  Equivalently the named outgoing channels are bulk acoustic transport p_s v_bulk,s, interfacial conversion A_s J_s, and reciprocal traction lx A_s V_s.  Positive values are sinks; propagating outgoing Re(Z)>0 is a sink even when la0=lv0=0.

### SymPy — `S11BB_ENERGY_SINKS`

    d(T+U)/dt=-sum_s[(delta_p_s+Lambda_X*A_s)V_s+mu_s J_s] (real-time convolution pairing). Positive real time-averages are sinks: outgoing bulk pressure transport, interfacial conversion, transferred slab mass chemical work, and reciprocal traction work; propagating impermeable Re(Z)>0 is a definite radiation sink

## Candidate 54: `WL_ENERGY_SOURCES` ↔ `S11BB_ENERGY_SOURCES`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_ENERGY_SOURCES`

    The same three signed exchange terms are sources when their real power is negative; free coefficients can make interfacial conversion or reciprocal traction source-valued.  No sign was imposed and no root was removed.  Indefinite stored stiffness K0<0 creates conservative growth but is not itself an external exchange term.

### SymPy — `S11BB_ENERGY_SOURCES`

    the same signed exchange channels are sources when their computed real quadratic contribution is negative; free coefficients permit this. Indefinite stored moduli generate conservative growth but are not external power sources

## Candidate 55: `WL_UNATTRIBUTED_SINK_TERMS` ↔ `S11BB_UNATTRIBUTED_SINK_TERMS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_UNATTRIBUTED_SINK_TERMS`

    NONE: every external term maps to bulk acoustic transport, interfacial mass conversion including transferred-mass pressure work, or the supplied reciprocal traction.

### SymPy — `S11BB_UNATTRIBUTED_SINK_TERMS`

    none; Q_J^direct=0 and no response kernel was varied

## Candidate 56: `WL_UNATTRIBUTED_EXCHANGE_TERMS` ↔ `S11BB_UNATTRIBUTED_EXCHANGE_TERMS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_UNATTRIBUTED_EXCHANGE_TERMS`

    NONE: Q_J^direct=0; no J*delta_v(deltaW), J*div(delta_v u), differentiated response, or other mechanical face term appears.

### SymPy — `S11BB_UNATTRIBUTED_EXCHANGE_TERMS`

    none; every term maps to pressure transport, material transfer, or the supplied reciprocal traction

## Candidate 57: `WL_PRESSURE_WORK_SIGN_CHECK` ↔ `S11BB_PRESSURE_WORK_SIGN_CHECK`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_PRESSURE_WORK_SIGN_CHECK`

    In the real propagating, impermeable, lx0=0 cut, off-shell pairing of the thickness equation with deltaW_t* gives slab pressure power -Sum_s Re[p_s V_s*]/2, while independent outgoing bulk flux is +Sum_s Re[p_s V_s*]/2.  Their sum/difference check is 0; no on-shell equation or period-average total derivative was used.

### SymPy — `S11BB_PRESSURE_WORK_SIGN_CHECK`

    PASS; off-shell slab pressure contribution paired with partial_t(deltaW) is -sum_s 1/2 Re(delta_p_s V_s*), outgoing bulk contribution has the opposite sign, and their symbolic sum/difference under J=0,Lambda_X0=0 is 0

## Candidate 58: `WL_FULL_TWO_PORT_BALANCE_CHECK` ↔ `S11BB_FULL_TWO_PORT_BALANCE_CHECK`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_FULL_TWO_PORT_BALANCE_CHECK`

    Face by face the off-shell slab pairing is -(1/2)Re[(p_s+lx A_s)V_s*+mu_s J_s*].  Comparison with the independently supplied slab-side exchange has residual 0 at order lx0^0 and residual 0 at order lx0^1; pressure, reciprocal-traction, and transfer channels separately give {0,0,0}.  Amplitudes and {p_s,J_s,A_s} remain algebraically free.

### SymPy — `S11BB_FULL_TWO_PORT_BALANCE_CHECK`

    PASS; face-by-face off-shell slab-minus-supplied-exchange differences are 0 at order Lambda_X0^0 and 0 at order Lambda_X0^1; pressure, mu_s*J and reciprocal-traction channels separately give (0,0,0). This does not test face closure, affinity normalization, or Lambda_X analytic orientation

## Candidate 59: `WL_TRANSVERSE_COUPLING` ↔ `S11BB_TRANSVERSE_COUPLING`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_TRANSVERSE_COUPLING`

    The tested coefficient is the mixed Fourier Hessian coefficient multiplying u_T^* e_W (equivalently the e_W force per transverse displacement amplitude).  B1 contains u only through d=i k.u_L, and every computed O(3) energy cross-invariant contains div(u), so the mixed coefficient is identically 0.  Because it vanishes, a normalization-dependent dimension for this coefficient is undetermined.

### SymPy — `S11BB_TRANSVERSE_COUPLING`

    coefficient defined as the off-diagonal Fourier operator mapping e_W to a transverse u_T equation (and reciprocally u_T to the thickness equation): 0 identically. B1 contains only div u and every scalar energy invariant contains only div u, so projection perpendicular to k annihilates the coupling; because it vanishes, a standalone normalization/dimension is undetermined

## Candidate 60: `WL_TRANSVERSE_DISPERSION` ↔ `S11BB_TRANSVERSE_DISPERSION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_TRANSVERSE_DISPERSION`

    rhoBr*om^2-muR*k^2=0, so om=+/-k*sqrt(muR/rhoBr).  The thickness/interface sector produces no modification on the uniform background.

### SymPy — `S11BB_TRANSVERSE_DISPERSION`

    rho_br0*omega^2=mu_R*k^2, with two transverse polarizations and no thickness/interface modification on the uniform background

## Candidate 61: `WL_TRANSVERSE_DISSIPATION` ↔ `S11BB_TRANSVERSE_DISSIPATION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_TRANSVERSE_DISSIPATION`

    Im(om)=0 for the transverse roots when muR/rhoBr>0 (purely imaginary conservative roots if that ratio is negative).  The coupling and dispersion are independent across the full range of la0,lv0,lx0,om*tauA,om*tauV,om*tauX and the slab-side affinity term.  This uniform-background result does not settle unconditional confinement.

### SymPy — `S11BB_TRANSVERSE_DISSIPATION`

    Im(omega)=0 for real positive mu_R/rho_br0 (or a conservative growing/decaying pair if their ratio is negative); it is independent over the full range of Lambda_A0,Lambda_V0,Lambda_X0,omega*tau_A,omega*tau_V,omega*tau_X and the slab-side affinity. This uniform result does not settle unconditional confinement

## Candidate 62: `WL_HOMOGENEITY_INPLANE_EQUATION` ↔ `S11BB_HOMOGENEITY_INPLANE_EOM`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_HOMOGENEITY_INPLANE_EQUATION`

    True

### SymPy — `S11BB_HOMOGENEITY_INPLANE_EOM`

    PASS

## Candidate 63: `WL_HOMOGENEITY_THICKNESS_EQUATION` ↔ `S11BB_HOMOGENEITY_THICKNESS_EOM`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_HOMOGENEITY_THICKNESS_EQUATION`

    True

### SymPy — `S11BB_HOMOGENEITY_THICKNESS_EOM`

    PASS

## Candidate 64: `WL_HOMOGENEITY_MASS_BALANCE` ↔ `S11BB_HOMOGENEITY_MASS_BALANCE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_HOMOGENEITY_MASS_BALANCE`

    True

### SymPy — `S11BB_HOMOGENEITY_MASS_BALANCE`

    PASS

## Candidate 65: `WL_HOMOGENEITY_AFFINITY` ↔ `S11BB_HOMOGENEITY_AFFINITY`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_HOMOGENEITY_AFFINITY`

    True

### SymPy — `S11BB_HOMOGENEITY_AFFINITY`

    PASS

## Candidate 66: `WL_HOMOGENEITY_CLOSURE` ↔ `S11BB_HOMOGENEITY_CLOSURE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_HOMOGENEITY_CLOSURE`

    True

### SymPy — `S11BB_HOMOGENEITY_CLOSURE`

    PASS

## Candidate 67: `WL_HOMOGENEITY_FACE_RESPONSE` ↔ `S11BB_HOMOGENEITY_FACE_RESPONSE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_HOMOGENEITY_FACE_RESPONSE`

    True

### SymPy — `S11BB_HOMOGENEITY_FACE_RESPONSE`

    PASS

## Candidate 68: `WL_HOMOGENEITY_TWO_PORT_POWER_IDENTITY` ↔ `S11BB_HOMOGENEITY_TWO_PORT_POWER_IDENTITY`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_HOMOGENEITY_TWO_PORT_POWER_IDENTITY`

    True

### SymPy — `S11BB_HOMOGENEITY_TWO_PORT_POWER_IDENTITY`

    PASS

## Candidate 69: `WL_HOMOGENEITY_DISPERSION_DETERMINANT` ↔ `S11BB_HOMOGENEITY_DISPERSION_DETERMINANT`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_HOMOGENEITY_DISPERSION_DETERMINANT`

    True

### SymPy — `S11BB_HOMOGENEITY_DISPERSION_DETERMINANT`

    PASS

## Candidate 70: `WL_HOMOGENEITY_ABLATION_DEMO` ↔ `S11BB_HOMOGENEITY_ABLATION_DEMO`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_HOMOGENEITY_ABLATION_DEMO`

    Deliberately replacing mu_s=muTheta/rhoBr by muTheta in A=mu_s-p/rhom gives equal-dimension test=False (failure detected); restoring division by rhoBr gives True.

### SymPy — `S11BB_HOMOGENEITY_ABLATION_DEMO`

    PASS (failure was detectable): deliberately replacing rho_m by rho_br0 in delta_p/rho gave dimensions [L^1 T^-2 M^0] instead of affinity [L^2 T^-2 M^0], so the checker returned FAIL; restoring rho_m returned PASS

## Candidate 71: `WL_CONTROL_NO_THICKNESS` ↔ `S11BB_CONTROL_NO_THICKNESS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONTROL_NO_THICKNESS`

    Set e_W=V=0 and delete the thickness row/column, but retain J=J_mu*mu_s and p=P_mu*mu_s.  Recomputed K_L=-(piD+piTheta*(-massD/massT))=-aTheta + kD - ((-aTheta + b3 + gTheta*k^2)*(I*om*rhoBr - (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))))/((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))); dispersion=(R_d B_theta-R_theta B_d)=-(k^2*(-aTheta + b3 + gTheta*k^2)*((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))) + (k^2*(aTheta - kD) + om^2*rhoBr)*((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))).  tauA remains in the surviving state-driven transfer.  tauV is irrelevant because its whole V-driven contribution vanishes, and tauX is irrelevant to B4/B5 because the thickness coordinate on which reciprocal traction acts has been removed.

### SymPy — `S11BB_CONTROL_NO_THICKNESS`

    e_W=V=0 gives K_eff=(K_L-D_theta)-(D_theta-b_theta)*R2_d/R2_theta and D_A=R1_d*R2_theta-R1_theta*R2_d; J=J_mu*mu_theta and delta_p=P_mu*mu_theta remain. Thickness response is removed; tau_A remains in the slab-state-driven transfer, while tau_V is irrelevant because the entire V-driven channel vanishes, and tau_X becomes mechanically irrelevant

## Candidate 72: `WL_CONTROL_A_ATTRIBUTION` ↔ `S11BB_CONTROL_A_ATTRIBUTION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONTROL_A_ATTRIBUTION`

    CONFOUNDED: changes can arise jointly from removing deltaW/e_W, face-motion drive, lv*V, pressure work pV, and reciprocal work lx*A*V.  The remaining mu_s-driven J and p forbid attribution to thickness alone; this control cannot separate those simultaneously removed channels.

### SymPy — `S11BB_CONTROL_A_ATTRIBUTION`

    CONFOUNDED: changes can arise jointly from removing the thickness field, face-motion drive, Lambda_V V, pressure work, and reciprocal-traction work; they cannot be attributed to thickness alone. The slab-state-driven transfer channel is retained

## Candidate 73: `WL_CONTROL_NO_GRADIENT_STIFFNESS` ↔ `S11BB_CONTROL_NO_GRADIENT_STIFFNESS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONTROL_NO_GRADIENT_STIFFNESS`

    kapW=0 replaces mE=kW*w0^2+kapW*w0^4*k^2 by kW*w0^2.  Recomputed chiW=w0^2/(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr), and B5 is det(M/.kapW->0)=0: k^2*(-aE + gThetaE*k^2 + cc*w0)*(((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr) - ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)) - k^2*(-aTheta + b3 + gTheta*k^2)*((aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) - ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr)) + (k^2*(aTheta - kD) + om^2*rhoBr)*((-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) - ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + (((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0*(gThetaE*k^2 + cc*w0))/rhoBr)).  C and all three relaxation times remain relevant.

### SymPy — `S11BB_CONTROL_NO_GRADIENT_STIFFNESS`

    set kappa_W=0: chi_W becomes W0/(R3_e|kappa_W=0), D_B=det(R1,R2,R3)|kappa_W=0; only the k^2 thickness stiffness contribution moves; all tau_I remain independent

## Candidate 74: `WL_CONTROL_IMPERMEABLE` ↔ `S11BB_CONTROL_IMPERMEABLE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONTROL_IMPERMEABLE`

    la0=lv0=0 gives J=0,p=Z V,A=mu_s-Z V/rhom while lx remains symbolic; B1 becomes theta+e_W+d=0 for om!=0.  Recomputed dispersion is k^2*(-aE + gThetaE*k^2 + cc*w0)*((-I)*om*rhoBr*(aE - aTheta + (aTheta*lx0*w0)/(rhoBr*(1 - I*om*tauX))) + I*om*rhoBr*(-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*lx0*w0)/(rhoBr*(1 - I*om*tauX)))) - k^2*(-aTheta + b3 + gTheta*k^2)*((-I)*om*rhoBr*(aE - aTheta + (aTheta*lx0*w0)/(rhoBr*(1 - I*om*tauX))) + I*om*rhoBr*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom)/qout - (lx0*om)/(qout*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (lx0*w0*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauX)))) + (k^2*(aTheta - kD) + om^2*rhoBr)*((-I)*om*rhoBr*(-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*lx0*w0)/(rhoBr*(1 - I*om*tauX))) + I*om*rhoBr*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom)/qout - (lx0*om)/(qout*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4 + (lx0*w0*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauX)))).  tauA,tauV are irrelevant; tauX remains if lx0!=0.

### SymPy — `S11BB_CONTROL_IMPERMEABLE`

    Lambda_A0=Lambda_V0=0, Lambda_X symbolic: D_C=det(R1,R2,R3)|Lambda_A0=Lambda_V0=0; mass transfer is cut and tau_A,tau_V become irrelevant, but propagating radiation resistance P_V=Z survives and tau_X remains

## Candidate 75: `WL_CONTROL_NO_CROSS_TERM` ↔ `S11BB_CONTROL_NO_CROSS_TERM`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONTROL_NO_CROSS_TERM`

    cc=0 sends mCross to gThetaE*k^2 and K0 to b3+kW*w0^2.  Recomputed B4 is K_L/.cc->0=-aTheta + kD - ((-aE + gThetaE*k^2)*((I*om*rhoBr - (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr) + ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)))/(-(((-I)*om*rhoBr + (2*gThetaE*k^2*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)) + ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) + (gThetaE*k^2*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4)) - ((-aTheta + b3 + gTheta*k^2)*(((-I)*om*rhoBr + (2*gThetaE*k^2*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))*(aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr) + (I*om*rhoBr - (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) + (gThetaE*k^2*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4)))/(-(((-I)*om*rhoBr + (2*gThetaE*k^2*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)) + ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) + (gThetaE*k^2*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4)); recomputed B5 is k^2*(-aE + gThetaE*k^2)*(((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr) - ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr)) - k^2*(-aTheta + b3 + gTheta*k^2)*(((-I)*om*rhoBr + (2*gThetaE*k^2*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))*(aE - aTheta + (aTheta*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr) - ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) + (gThetaE*k^2*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4)) + (k^2*(aTheta - kD) + om^2*rhoBr)*(((-I)*om*rhoBr + (2*gThetaE*k^2*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + ((b3 + gTheta*k^2)*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr) - ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) + (gThetaE*k^2*((la0*om*rhom^2)/(qout*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + (lx0*rhom^2)/((rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0)/rhoBr + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4)).  The independent gradient cross gThetaE is not cut.

### SymPy — `S11BB_CONTROL_NO_CROSS_TERM`

    C=0: replace b_theta_e by A_theta_e*k^2 in K_eff and D; D_D=det(R1,R2,R3)|C=0; all interface channels/times remain

## Candidate 76: `WL_CONTROL_NO_MU_COUPLING` ↔ `S11BB_CONTROL_NO_MU_COUPLING`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONTROL_NO_MU_COUPLING`

    Set the slab-side term in A to zero while retaining A=-p/rhom, la,lv and bulk coupling.  Then p=P_V V,J=J_V V; the recomputed determinant is k^2*(-aE + gThetaE*k^2 + cc*w0)*((-I)*(aE - aTheta)*om*rhoBr + I*om*rhoBr*(-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0)) - k^2*(-aTheta + b3 + gTheta*k^2)*((aE - aTheta)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + I*om*rhoBr*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4)) + (k^2*(aTheta - kD) + om^2*rhoBr)*((-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0)*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + I*om*rhoBr*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - (I/2)*om*((om*rhom^2*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) - (lx0*om*rhom*(rhom + lv0/(1 - I*om*tauV)))/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))*(1 - I*om*tauX)))*w0^2 + k^2*kapW*w0^4)).  Intrinsic interfacial admissibility for independent (A,V) is unchanged, but the slab port with independent (V,mu_s) contains the free cross power mu_s*J(V) and is PSD for all mu_s only if J_V=0 in addition to Re(F_V)>=0.  This is not the impermeable cut; tauA and tauV remain.

### SymPy — `S11BB_CONTROL_NO_MU_COUPLING`

    eta=0 in A=eta*mu_theta/rho_br0-delta_p/rho_m gives P_mu=J_mu=A_mu=0; D_E=det(R1,R2,R3)|eta=0. Port H_E is recomputed by the same H formula with those zero coefficients; interfacial coefficient admissibility in independent (A,V) is unchanged

## Candidate 77: `WL_CONTROL_NO_RECIPROCAL_TRACTION` ↔ `S11BB_CONTROL_NO_RECIPROCAL_TRACTION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONTROL_NO_RECIPROCAL_TRACTION`

    lx0=0 gives F_V=P_V,F_mu=P_mu.  Recomputed B3 chiW=w0^2/(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - ((I/2)*om^2*rhom^2*(rhom + lv0/(1 - I*om*tauV))*w0^2)/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + k^2*kapW*w0^4 + (la0*om*rhom^2*w0*(gThetaE*k^2 + cc*w0))/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))); recomputed B4 K_L=-aTheta + kD - ((-aE + gThetaE*k^2 + cc*w0)*((I*om*rhoBr - (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(aE - aTheta + (aTheta*la0*om*rhom^2*w0)/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) + ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*la0*om*rhom^2*w0)/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))))/(-((-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*la0*om*rhom^2*w0)/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))) + ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - ((I/2)*om^2*rhom^2*(rhom + lv0/(1 - I*om*tauV))*w0^2)/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + k^2*kapW*w0^4 + (la0*om*rhom^2*w0*(gThetaE*k^2 + cc*w0))/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))) - ((-aTheta + b3 + gTheta*k^2)*((aE - aTheta + (aTheta*la0*om*rhom^2*w0)/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) + (I*om*rhoBr - (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - ((I/2)*om^2*rhom^2*(rhom + lv0/(1 - I*om*tauV))*w0^2)/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + k^2*kapW*w0^4 + (la0*om*rhom^2*w0*(gThetaE*k^2 + cc*w0))/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))))/(-((-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*la0*om*rhom^2*w0)/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))) + ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - ((I/2)*om^2*rhom^2*(rhom + lv0/(1 - I*om*tauV))*w0^2)/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + k^2*kapW*w0^4 + (la0*om*rhom^2*w0*(gThetaE*k^2 + cc*w0))/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))); recomputed B5 determinant=k^2*(-aE + gThetaE*k^2 + cc*w0)*(((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(aE - aTheta + (aTheta*la0*om*rhom^2*w0)/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) - ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*la0*om*rhom^2*w0)/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))) - k^2*(-aTheta + b3 + gTheta*k^2)*((aE - aTheta + (aTheta*la0*om*rhom^2*w0)/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) - ((-I)*om*rhoBr + (2*aTheta*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - ((I/2)*om^2*rhom^2*(rhom + lv0/(1 - I*om*tauV))*w0^2)/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + k^2*kapW*w0^4 + (la0*om*rhom^2*w0*(gThetaE*k^2 + cc*w0))/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))) + (k^2*(aTheta - kD) + om^2*rhoBr)*((-b3 - gTheta*k^2 + gThetaE*k^2 + cc*w0 + ((b3 + gTheta*k^2)*la0*om*rhom^2*w0)/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*((-I)*om*rhoBr - (I*om*rhom*(-((la0*om*rhom)/(qout*(1 - I*om*tauA))) + (lv0*rhom)/(1 - I*om*tauV))*w0)/(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))) + (2*la0*rhom^2*(gThetaE*k^2 + cc*w0))/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA))))) - ((-I)*om*rhoBr + (2*(b3 + gTheta*k^2)*la0*rhom^2)/(rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))*(-(gThetaE*k^2) - cc*w0 + kW*w0^2 - muW*om^2*w0^2 - ((I/2)*om^2*rhom^2*(rhom + lv0/(1 - I*om*tauV))*w0^2)/(qout*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))) + k^2*kapW*w0^4 + (la0*om*rhom^2*w0*(gThetaE*k^2 + cc*w0))/(qout*rhoBr*(1 - I*om*tauA)*(rhom^2 + (la0*om*rhom)/(qout*(1 - I*om*tauA)))))).  Mechanical-operator difference from the B2 symbolic lx0=0 operator has entrywise zero tests={{True, True, True}, {True, True, True}, {True, True, True}}; the two-port identity difference after the same substitution is 0.  tauX is irrelevant, tauA and tauV remain.  Interfacial admissibility reduces to la0>=0 and lv(om)=0 for all real om, hence lv0=0; conditional reciprocity also forces lv0=0.

### SymPy — `S11BB_CONTROL_NO_RECIPROCAL_TRACTION`

    Lambda_X0=0: D_F=det(R1,R2,R3)|Lambda_X0=0 and D_full-D_F=det(R1,R2,DeltaR3); mechanical-operator and power-identity substitution residuals are both 0. Port H_F uses T=p. Interfacial admissibility reduces to Lambda_A0>=0 and Lambda_V0=0; conditional reciprocity likewise forces Lambda_V=0. tau_X is irrelevant

## Candidate 78: `WL_CONTROLS_ON_TRANSVERSE` ↔ `S11BB_CONTROLS_ON_TRANSVERSE`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_CONTROLS_ON_TRANSVERSE`

    Recomputation finds no movement of the transverse coupling (0), dispersion rhoBr*om^2-muR*k^2, or interface dependence under A (field absent but the pre-cut mixed coefficient was zero), B (kapW=0), C (impermeable), or D (cc=0).  The reason discovered from B1 and the full basis is structural: every scalar coupling contains div(u), whose transverse projection is exactly zero; none of A-D creates an allowed transverse-scalar invariant.

### SymPy — `S11BB_CONTROLS_ON_TRANSVERSE`

    recomputed A-D: coupling remains identically 0 and rho_br0*omega^2=mu_R*k^2 in every control. Nothing moves because transverse projection kills div u and every scalar/thickness bilinear; this is a computed uniform-background structural result, not an unconditional-confinement claim

## Candidate 79: `WL_VALIDITY_CONDITIONS` ↔ `S11BB_VALIDITY_CONDITIONS`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_VALIDITY_CONDITIONS`

    Uniform background; |theta|<<1, |e_W|<<1, |k.u|<<1, |k deltaW|<<1; linear face/bulk amplitudes and nonrelativistic velocities; wavelengths resolved by the continuum; homogeneous plane-wave sector only.  Bulk rest-frame linearization additionally requires |v0*qout/om|<<1, using complex modulus as the relative-amplitude norm away from om=0.  The supplied Debye kernels impose no small-om*tau expansion: tauA,tauV,tauX remain independent; low-frequency approximations would separately require |om*tauA|,|om*tauV|,|om*tauX|<<1, while using their full rational forms requires only that their modeled single-relaxation response remain valid.

### SymPy — `S11BB_VALIDITY_CONDITIONS`

    uniform background; |theta|<<1, |e_W|<<1, |grad u|<<1, |deltaW|/W0<<1, small face/bulk velocities, linear constitutive amplitudes, and |v0|*|q_out|/|omega|<<1 for the bulk rest-frame expansion. The scalar complex modulus is the norm used for complex omega,q_out. The supplied Debye kernels retain independent omega*tau_A, omega*tau_V, omega*tau_X and are algebraically defined for all omega off their poles; no microscopic bandwidth criterion was supplied

## Candidate 80: `WL_VALIDITY_FAILURE_REGION` ↔ `S11BB_VALIDITY_FAILURE_REGION`
**Status: REQUIRES REVIEW BEFORE INSERTION**

### Mathematica — `WL_VALIDITY_FAILURE_REGION`

    The background-flow approximation fails where |v0| |qout(om,k)| >= |om|; at k=0 its measure is |v0|/cs, while for fixed nonzero k it fails near om=0 because |qout|~|k|.  At complex om the modulus is meaningful as a local complex-amplitude error norm, but it does not by itself guarantee global spatial boundedness; at om=0 the ratio is undefined and must be assessed before division.  Linearization also fails at large amplitudes/gradients, at face/feedback poles where linear response diverges, and for non-uniform backgrounds where a single k sector is not closed.

### SymPy — `S11BB_VALIDITY_FAILURE_REGION`

    background-flow expansion fails where |v0|*|q_out(omega,k)| is not small compared with |omega|, notably omega->0 at fixed k!=0 and sufficiently far from the sound cone when |q_out/omega| is large; on the cone it reduces to |v0|/c_s0. For complex values the modulus comparison is meaningful as a scalar norm except at omega=0; crossing a dragged cut requires sheet tracking. Kernel-model failure outside a microscopic response bandwidth is NOT_ESTABLISHED because no such bandwidth was supplied
