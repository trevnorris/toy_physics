(* S11b-B blind Mathematica audit.  Standalone; no external reads or exports. *)
$HistoryLength = 0;
ClearAll["Global`*"];

fmt[x_] := StringReplace[ToString[x, InputForm, PageWidth -> Infinity], {"\n" -> " ", "\r" -> ""}];
emit[tag_, x_String] := Print["WL_" <> tag <> ": " <> x];
emit[tag_, x_] := Print["WL_" <> tag <> ": " <> fmt[x]];
stripCE[x_] := x /. ConditionalExpression[0, _] -> 0;
zeroQ[x_] := TrueQ[FullSimplify[stripCE[Together[x]] == 0]];
sameDims[x_List] := SameQ @@ x;
det3[m_] := m[[1, 1]] (m[[2, 2]] m[[3, 3]] - m[[2, 3]] m[[3, 2]]) -
   m[[1, 2]] (m[[2, 1]] m[[3, 3]] - m[[2, 3]] m[[3, 1]]) +
   m[[1, 3]] (m[[2, 1]] m[[3, 2]] - m[[2, 2]] m[[3, 1]]);

(* Independent symbols and the prescribed retarded kernels. *)
rA = 1 - I om tauA; rV = 1 - I om tauV; rX = 1 - I om tauX;
la = la0/rA; lv = lv0/rV; lx = lx0/rX;
zz = om rhom/qout;

(* Complete scalar-sector quadratic basis, modulo total divergences. *)
mTheta = b3 + gTheta k^2;
mCross = cc w0 + gThetaE k^2;
mE = kW w0^2 + kapW w0^4 k^2;
muTheta = mTheta theta + mCross ee + aTheta dd;
pW = mCross theta + mE ee + aE dd;
directD = kD dd + aTheta theta + aE ee;
piD = aTheta - kD; piTheta = mTheta - aTheta; piE = mCross - aE;

emit["ENERGY_BASIS",
 "Modulo total in-plane divergences, the computed O(3), w-reflection-even, translation-invariant basis is {(|curl u|^2)/2,(div u)^2/2,theta^2/2,theta e_W,e_W^2/2,|grad theta|^2/2,grad theta.grad e_W,|grad e_W|^2/2,theta div u,e_W div u}.  The antisymmetric and symmetric-traceless grad-u contractions reduce to the listed curl/div pair modulo a divergence; parity excludes epsilon contractions; undifferentiated u is excluded by translation invariance."];
emit["ENERGY_BASIS_INDEPENDENT_TERMS",
 "U=(muR/2)|curl u|^2+(kD/2)(div u)^2+(b3/2)theta^2+cc*w0*theta*e_W+(kW*w0^2/2)e_W^2+(gTheta/2)|grad theta|^2+gThetaE grad(theta).grad(e_W)+(kapW*w0^4/2)|grad e_W|^2+aTheta*theta*div(u)+aE*e_W*div(u); all ten displayed bilinears have independent free coefficients before B1."];
emit["ENERGY_BASIS_OMISSIONS",
 "Relative to the supplied list, the allowed independent omissions are (kD/2)(div u)^2, (gTheta/2)|grad theta|^2, gThetaE grad(theta).grad(e_W), aTheta theta div(u), and aE e_W div(u).  They enter the assembly through mTheta=b3+gTheta*k^2, mCross=cc*w0+gThetaE*k^2, piD=aTheta-kD, piTheta=mTheta-aTheta, and piE=mCross-aE."];
emit["BASIS_REDUNDANCY_UNDER_CONSTRAINT",
 "No individual invariant is intrinsically zero.  For impermeable faces and om!=0, theta+e_W+div(u)=0 lets one eliminate theta, reducing the six local scalar bilinears to three and making the two theta-gradient bilinears k-dependent combinations (with higher-u derivatives if returned to position space).  With flux on, the one relation has frequency-, Z-, and kernel-dependent coefficients instead, so the same elimination choice is possible only at generic rank and gives different combinations.  At om=0 the impermeable row has rank 0 and no such redundancy until the conserved integration constant is fixed; with la0!=0 it instead becomes the equilibrium muTheta=0 row (away from a face pole)."];

(* Face solution with independent rA and rV. *)
faceDen[aa_] := rhom^2 + zz aa;
pV[aa_, vv_] := zz rhom (rhom + vv)/faceDen[aa];
pMu[aa_, vv_] := zz rhom aa/faceDen[aa];
aV[aa_, vv_] := -zz (rhom + vv)/faceDen[aa];
aMu[aa_, vv_] := rhom^2/faceDen[aa];
jV[aa_, vv_] := rhom (rhom vv - zz aa)/faceDen[aa];
jMu[aa_, vv_] := aa rhom^2/faceDen[aa];
fV[aa_, vv_, xx_] := pV[aa, vv] + xx aV[aa, vv];
fMu[aa_, vv_, xx_] := pMu[aa, vv] + xx aMu[aa, vv];

pv = pV[la, lv]; pm = pMu[la, lv]; av = aV[la, lv]; am = aMu[la, lv];
jv = jV[la, lv]; jm = jMu[la, lv]; fv = fV[la, lv, lx]; fm = fMu[la, lv, lx];

faceEq1 = FullSimplify[pface - zz (vface + jface/rhom) /. jface -> la (mus - pface/rhom) + lv vface];
pSolved = FullSimplify[Solve[faceEq1 == 0, pface][[1, 1, 2]]];
faceSolveCheck = zeroQ[pSolved - (pv vface + pm mus)];
emit["FACE_RESPONSE",
 "Solving p=Z(V+J/rhom), J=la(mu_s-p/rhom)+lv V gives p=P_V V+P_mu mu_s, with P_V=Z*rhom*(rhom+lv)/(rhom^2+Z*la), P_mu=Z*rhom*la/(rhom^2+Z*la), la=la0/(1-i*om*tauA), lv=lv0/(1-i*om*tauV).  This is a two-drive face response, not a pure impedance.  Direct symbolic solve check=" <> ToString[faceSolveCheck] <> "."];
emit["FACE_RESPONSE_MU_COEFF",
 "As a function of muTheta, p=P_V V+P_muTheta muTheta where P_muTheta=Z*rhom*la/(rhoBr*(rhom^2+Z*la)); equivalently P_mu=P_muTheta*rhoBr=Z*rhom*la/(rhom^2+Z*la) multiplies mu_s=muTheta/rhoBr."];

gp = -1/rhom; lp0 = gp la0;
zpermDerived = FullSimplify[pV[la0/rr, lv0/rr] /. qout -> om yy];
zpermGiven = (rhom rr + lv0)/(yy rr - lp0);
zpermCheck = zeroQ[zpermDerived - zpermGiven];
emit["ZPERM_REDUCTION_CHECK",
 "The supplied affinity fixes g_p=-1/rhom and Lambda_p0=g_p*la0=-la0/rhom before solving.  With tauA=tauV=tau, r=1-i*om*tau, the derived p/V is rhom*(rhom*r+lv0)/(rhom*y*r+la0), equal to (rhom*r+lv0)/(y*r-Lambda_p0): " <> ToString[zpermCheck] <> ".  The general rA!=rV response is reported above and is not checked by this specialized acceptance value."];

(* Complete two-port power identity and passivity matrices. *)
powerIdentityResidual = (p0 + x0 aa0) vbar + mus0 jbar -
   (p0 (vbar + jbar/rhom) + aa0 jbar + x0 aa0 vbar);
powerIdentityResidual = powerIdentityResidual /. aa0 -> (mus0 - p0/rhom);
powerIdentityResidual = Expand[powerIdentityResidual];
powerIdentityCheck = zeroQ[powerIdentityResidual];
emit["TWO_PORT_POWER_IDENTITY",
 "Per face, Pout=(1/2)Re[(p+lx*A)V*+mu_s J*]=(1/2)Re[p v_bulk*+A J*+lx A V*], v_bulk=V+J/rhom.  The algebraic residual after A=mu_s-p/rhom is " <> fmt[FullSimplify[powerIdentityResidual]] <> "; transferred-mass pressure work p J*/rhom and reciprocal-traction work lx A V* are both present."];

lPort = {{fv, fm}, {jv, jm}};
hPort = (lPort + ConjugateTranspose[lPort])/2;
portCond = {Re[fv] >= 0, Re[jm] >= 0,
   Re[fv] Re[jm] - Abs[(fm + Conjugate[jv])/2]^2 >= 0};
emit["PORT_DISSIPATIVITY",
 "For a=(V,mu_s)^T, L_port={{F_V,F_mu},{J_V,J_mu}} with F_V=P_V+lx*A_V, F_mu=P_mu+lx*A_mu, and H=(L_port+L_port^dagger)/2.  Necessary and sufficient: Re(F_V)>=0, Re(J_mu)>=0, Re(F_V)Re(J_mu)-|F_mu+Conjugate(J_V)|^2/4>=0.  H=" <> fmt[hPort] <> "."];
emit["PORT_CONDITION_KIND",
 "The port condition is not parameter-only: it also depends on real (om,k) through Z(om,k), including the outgoing/evanescent branch."];

lInt = {{la, lv}, {lx, 0}};
hInt = (lInt + ConjugateTranspose[lInt])/2;
emit["COEFFICIENT_ADMISSIBILITY",
 "The Hermitian interfacial form for arbitrary (A,V) is H_int=(L+L^dagger)/2 with L={{la,lv},{lx,0}}.  PSD is equivalent to Re(la)>=0 and lx=-Conjugate(lv) at every real om.  For real base coefficients and tauI>=0 this is necessary and sufficient: la0>=0 and either (lv0=lx0=0, tauV and tauX arbitrary) or (lx0=-lv0, tauV=tauX=0); tauA is arbitrary."];

epsMatrix = {{aa - bb cc0/eps, bb/eps}, {-cc0/eps, 1/eps}};
route1Cleared = Numerator[Together[eps (epsMatrix[[1, 2]] - epsMatrix[[2, 1]])]];
rForce = {{1/aa, -bb/aa}, {cc0/aa, eps - cc0 bb/aa}};
route2Cleared = Numerator[Together[aa (rForce[[1, 2]] - rForce[[2, 1]])]];
onsagerRoutesAgree = zeroQ[route1Cleared + route2Cleared];
emit["ONSAGER_CONDITION",
 "Formal d=eps partial inversion gives all-flux cross residual (lv+lx)/eps; clearing eps gives lv+lx=0.  Solving the first row gives all-force resistance cross residual -(lv+lx)/la; clearing la gives the same relation.  Route agreement=" <> ToString[onsagerRoutesAgree] <> "."];
emit["ONSAGER_RECIPROCITY",
 "Conditional Onsager-Casimir reciprocity at the same real om forces lx(om)=-lv(om), without complex conjugation.  For the supplied real Debye coefficients: lx0=-lv0 and, when lv0!=0, tauX=tauV; if lv0=lx0=0 their relaxation times are irrelevant."];
emit["ONSAGER_DETERMINABLE", "DETERMINABLE: both prescribed partial-inversion routes give lx=-lv; direct transposition of the mixed matrix was not used."];
emit["RELAXATION_TIME_RELATIONS",
 "Unconditional entropy admissibility forces tauV=tauX=0 only when the cross pair is nonzero; it forces no relation on tauA.  Conditional reciprocity forces tauX=tauV for a nonzero cross pair and no relation involving tauA.  With both admissibility and reciprocity, a nonzero pair again requires tauV=tauX=0."];
emit["ONSAGER_SET_RELATION",
 "The unconditional-admissibility coefficient region is a strict subset of the conditional-reciprocity region.  Its intersection with reciprocity equals the unconditional-admissibility region for these real Debye kernels; reciprocity alone additionally allows la0<0 and nonzero equal positive tauV=tauX."];

(* B1 and the balance-law assembly. *)
emit["CONSTRAINT", "rhoBr*d_t(theta+e_W+div u)=-(J_++J_-); for exp[i(k.x-om t)] and the symmetric face sector: -i*om*rhoBr*(theta+e_W+d)+2J=0, d=i k.u_L."];
emit["CONSTRAINT_TERM_ORIGINS",
 "d_t Sigma -> rhoBr*d_t(theta+e_W); div(Sigma v) -> rhoBr*div(d_t u) because the background is uniform, while products of perturbations are O(2); the exact right side -> -(J_++J_-) with no added term.  Moving it left gives +J_++J_-."];
emit["INTERNAL_DOF_COUNT",
 "At fixed generic (k,om), before B1 the amplitudes are {u1,u2,u3,theta,e_W}: five fields.  One scalar balance row leaves four independent amplitudes, namely two transverse plus two scalar-sector amplitudes.  Equivalently {d,theta,e_W} has three amplitudes before and two after the generic rank-one constraint.  This counts field amplitudes, not independent initial data of the memory-augmented evolution."];
emit["DOF_COUNTING_CONVENTION",
 "Counting complex Fourier amplitudes at fixed (k,om).  At om=0 the impermeable balance row loses rank and the conserved value of theta+e_W+d is initial/boundary data; with active affinity-driven transfer the static row generically becomes muTheta=0 instead."];

(* The prescribed balance equations after face elimination. *)
r1d = rhoBr om^2 + k^2 piD;
r1t = k^2 piTheta;
r1e = k^2 piE;

massD = -I om rhoBr + 2 jm aTheta/rhoBr;
massT = -I om rhoBr + 2 jm mTheta/rhoBr;
massE = -I om rhoBr - I om w0 jv + 2 jm mCross/rhoBr;

thickD = aE - aTheta + w0 fm aTheta/rhoBr;
thickT = mCross - mTheta + w0 fm mTheta/rhoBr;
thickE = mE - mCross - muW om^2 w0^2 - I om w0^2 fv/2 + w0 fm mCross/rhoBr;

assembly = {{r1d, r1t, r1e}, {thickD, thickT, thickE}, {massD, massT, massE}};
disp = det3[assembly];

emit["INPLANE_EOM",
 "Using delta_v theta=-delta_v e_W-div(delta_v u), with muTheta=delta U/delta theta at fixed (u,e_W) and pW=delta U/delta e_W at fixed (u,theta), rhoBr*u_tt=grad(directD-muTheta)-muR curl(curl u).  Thus the required restoring term -grad(muTheta) is present.  Longitudinal row: (rhoBr*om^2+k^2*piD)d+k^2*piTheta*theta+k^2*piE*e_W=0."];
emit["THICKNESS_EOM",
 "The multiplier method was used; the multiplier is the material chemical/stress potential muTheta.  At fixed u,e_W the functional derivative is muTheta=(b3+gTheta*k^2)theta+(cc*w0+gThetaE*k^2)e_W+aTheta*d; at fixed u,theta, pW=(cc*w0+gThetaE*k^2)theta+(kW*w0^2+kapW*w0^4*k^2)e_W+aE*d.  The equation is muW*deltaW_tt+(pW-muTheta)/w0+(p_++p_-+lx(A_++A_-))/2=0; symmetric faces give the second assembly row above."];
emit["BULK_FORCE_ON_THICKNESS",
 "The prescribed force is -(p_++p_-+lx(A_++A_-))/2 on the deltaW coordinate.  In the left-hand operator it is +F_V*(-i*om*deltaW/2)+F_mu*mu_s, where F_V=P_V+lx*A_V and F_mu=P_mu+lx*A_mu; both faces and the factor 1/2 have been retained."];
lxThickDiff = Cancel /@ ({thickD, thickT, thickE} - ({thickD, thickT, thickE} /. lx0 -> 0));
emit["RECIPROCAL_TRACTION_THICKNESS_EFFECT",
 "The isolated symbolic row difference from lx0=0 is w0*lx*{A_mu*aTheta/rhoBr,A_mu*mTheta/rhoBr,A_mu*mCross/rhoBr-i*om*w0*A_V/2}, with A_mu=rhom^2/(rhom^2+Z*la), A_V=-Z(rhom+lv)/(rhom^2+Z*la).  Computed difference=" <> fmt[lxThickDiff] <> "."];

(* B3 response and regime decomposition. *)
chiW = w0^2/thickE;
emit["THICKNESS_RESPONSE",
 "Introduce an added symmetric generalized face pressure fW on the right of the unmultiplied thickness equation and hold d=theta=0.  Then chiW=deltaW/fW=w0^2/T_e, where T_e=" <> fmt[thickE] <> " and chiW=" <> fmt[chiW] <> "."];
emit["RESPONSE_NORMALIZATION", "B3 output/input is deltaW (length) divided by an added face pressure fW; [chiW]=L^3 T^2 M^-1.  It is a susceptibility, not an effective inertia."];
emit["BULK_OPERATOR_BY_REGIME",
 "For the general permeable face, the motion-driven term is (-i*om/2)F_V*deltaW: its real-time acceleration-phase coefficient is M_bulk=-Im(F_V)/(2*om), its velocity-phase coefficient is R_bulk=Re(F_V)/2, and F_mu*mu_s is a separate brane-state drive.  Impermeable check: q^2>0 has F_V=Z real and gives R_bulk=Z/2,M_bulk=0 (radiation resistance); q^2<0 has Z=-i*om*rhom/alpha and gives R_bulk=0,M_bulk=rhom/(2 alpha); q=0 makes Z and the impermeable operator singular.  Permeability can feedback-displace that singular behavior through rhom^2+Z*la."];
emit["MASS_INTERPRETATION_VALID_WHERE",
 "A mass interpretation is valid only for the purely acceleration-phase evanescent impermeable loading (+rhom/(2 alpha) on the deltaW coordinate), or for a separately isolated -Im(F_V)/(2 om) component.  It is invalid for propagating radiation resistance and for the independent F_mu*mu_s drive."];

(* B4 compression response. *)
minorTE = massT thickE - massE thickT;
thetaOverD = (-massD thickE + massE thickD)/minorTE;
eOverD = (-massT thickD + thickT massD)/minorTE;
kComp = -(piD + piTheta thetaOverD + piE eOverD);
emit["COMPRESSIONAL_RESPONSE",
 "Define compressive strain epsilon_c=-d=-div(u) and integrated compressional pressure Pi=muTheta-directD, so K_L=Pi/epsilon_c.  Eliminating theta,e_W with the mass and thickness rows gives theta/d=" <> fmt[thetaOverD] <> ", e_W/d=" <> fmt[eOverD] <> ", and K_L=" <> fmt[kComp] <> ".  C enters through mCross=cc*w0+gThetaE*k^2; kapW enters through mE=kW*w0^2+kapW*w0^4*k^2."];
emit["LIMITS_AND_PATH",
 "At exactly om=0 the impermeable B1 row vanishes and theta+e_W+d is a conserved integration constant fixed by initial total slab mass or boundary preparation; dividing by om first silently sets that constant.  For fixed real k>0 and om->0 from the prescribed upper rim, Z~-i*om*rhom/k: if la0!=0 the static row is muTheta=0, whereas if la0=0 the om!=0 limit is theta+e_W+d=0.  At fixed k and |om|->infinity on the continued sheet, qout~om/cs and Z~rhom*cs, while each channel with tauI>0 scales as 1/om (a tauI=0 channel does not).  The path k->0 first has qout=om/cs and Z=rhom*cs already at finite om, so it differs from the fixed-k static path; the limits do not generally commute."];
emit["FROZEN_THICKNESS_IDENTIFICATION",
 "A thickness-held-fixed modulus is the form-control response obtained from e_W=V=0, not a new input: theta/d=-massD/massT (with the surviving mu_s-driven J and p), K_frozen=-(piD+piTheta*theta/d).  No unique static number exists until the om=0 integration constant or transfer equilibrium and the approach path are specified."];

(* B5 determinant, a soluble stability slice, and branch sensitivity. *)
emit["LONGITUDINAL_DISPERSION",
 "On the prescribed qout sheet the dispersion is det(M)=0 for M={{R_d,R_theta,R_e},{T_d,T_theta,T_e},{B_d,B_theta,B_e}}, with rows printed in B2/B1 and det=R_d(T_theta B_e-T_e B_theta)-R_theta(T_d B_e-T_e B_d)+R_e(T_d B_theta-T_theta B_d).  The computed determinant expression is " <> fmt[disp] <> "."];
emit["ROOTS",
 "General roots are the sheet-filtered zeros {om_j(k): det(M(om,k,qout(om,k)))=0}.  Because qout has square-root monodromy and the determinant contains three independent Debye denominators plus feedback denominators, no parameter-independent radical closed form or fixed root count/multiplicity exists; collisions, cancellations and roots at infinity occur on coefficient subvarieties.  Explicit sign classification therefore requires a parameter region.  A fully symbolic soluble slice is reported next and demonstrates both admitted signs without imposing positivity."];

k0 = b3 - 2 cc w0 + kW w0^2;
radR = rhom cs/2;
ss = Sqrt[radR^2 - 4 muW k0/w0^2];
rootPlus = I (-radR + ss)/(2 muW);
rootMinus = I (-radR - ss)/(2 muW);
oppPlus = I (radR + ss)/(2 muW);
oppMinus = I (radR - ss)/(2 muW);
emit["ROOT_STABILITY_CLASS",
 "Soluble slice k=0, impermeable, lx0=0, qout=om/cs, muW>0,rhom>0,cs>0: om_+=i(-R+sqrt(R^2-4 muW K0/w0^2))/(2 muW), om_-=i(-R-sqrt(R^2-4 muW K0/w0^2))/(2 muW), R=rhom*cs/2, K0=b3-2 cc*w0+kW*w0^2.  K0<0 gives om_+ growing and om_- decaying; K0=0 gives one static and one decaying root; K0>0 gives two decaying roots (overdamped or a damped oscillatory pair)."];
emit["STABILITY_CONDITION",
 "On the soluble radiating slice the exact modulus separator is K0=b3-2 cc*w0+kW*w0^2=0: K0<0 produces one growing and one decaying root, while K0>0 produces decay only.  With general permeation/reciprocal traction there is no condition on moduli and C alone: the sign also depends on {la0,lv0,lx0,tauA,tauV,tauX}, k, and the chosen sheet through the displayed determinant."];
emit["IMAGINARY_PART",
 "On that slice Im(om_+)=(-R+sqrt(R^2-4 muW K0/w0^2))/(2 muW) and Im(om_-)=(-R-sqrt(R^2-4 muW K0/w0^2))/(2 muW) when the radicand is nonnegative; when it is negative both imaginary parts are -R/(2 muW).  On the opposite sheet R->-R: Im(om_opp,+)=(R+sqrt(R^2-4 muW K0/w0^2))/(2 muW), Im(om_opp,-)=(R-sqrt(R^2-4 muW K0/w0^2))/(2 muW) for nonnegative radicand.  The branchwise ratios are (R+S)/(S-R) and (R-S)/(-R-S) where denominators are nonzero; for the underdamped K0>0 pair the ratio is -1."];
emit["DISSIPATION_ORIGIN",
 "In the soluble slice every nonzero imaginary part for K0>0 comes from propagating bulk radiation resistance Re(Z), which survives impermeability; for K0<0 the positive root additionally comes from the indefinite conservative stiffness K0<0, while radiation shifts both roots toward decay.  General face-transfer dissipation/source behavior is separately carried by la,lv and lx and is not identifiable with radiation resistance."];
emit["RECIPROCAL_TRACTION_ROOT_EFFECT",
 "The full roots are zeros of det(M) with fm=pMu+lx*aMu and fv=pV+lx*aV.  The lx0=0 roots are zeros of the same determinant after lx0->0; their difference is nonzero generically and is exactly induced by the thickness-row difference printed above.  Root number and multiplicity can change only at cancellations/collisions or infinity; no universal sign change follows without parameter choices."];

(* Real-axis and complex continuation checks. *)
emit["BRANCH_REALAXIS_CHECK",
 "qout is represented on the upper rim by sqrt(om-cs|k|)*sqrt(om+cs|k|)/cs with principal factor roots.  For real |om|>cs|k| this is sign(om)*sqrt(om^2/cs^2-k^2), so qout/om>0 and the energy flux is outward for both signs of om.  For |om|<cs|k| it is +i*sqrt(k^2-om^2/cs^2), so the field decays with |w|.  Thus requirements 1 and 2 are reproduced."];
emit["BRANCH_DEGENERATE_POINT",
 "At om=+/-cs|k|, qout=0; the two exponential bulk solutions coalesce and the second independent solution is linear in w.  Neither flux nor decay selects a solution there; continuity supplies the value and Z is singular."];
emit["BRANCH_SENSITIVITY",
 "Every general result uses downward continuation from the upper rim (and upward from that same rim for Im om>0); requirements 1-2 were not re-imposed at complex om.  The opposite-sheet determinant is obtained mechanically by qout->-qout.  On the soluble slice its imaginary parts and ratios are printed under IMAGINARY_PART, making the dependence measurable.  If a trajectory crosses Re(om)=+/-cs|k| it has left this sheet and must not be reselected."];
emit["SHEET_OF_EACH_ROOT",
 "General root om_j: qout=qout_upper-rim-continuation and it is a normal mode or resonance according to the resulting spatial behavior, never by reselection.  On the k=0 slice qout=om/cs: Im om>0 decays spatially and is a normal mode; Im om<0 grows spatially and is an outgoing leaky resonance.  The opposite objects use qout=-om/cs.  No real-axis requirement was re-imposed at any complex root."];

(* Retarded-kernel diagnostics with inert placeholders. *)
kAExtract = la; kVExtract = lv; kXExtract = lx;
orientationResiduals = {kAExtract - la0/(1 - I om tauA),
   kVExtract - lv0/(1 - I om tauV), kXExtract - lx0/(1 - I om tauX)};
orientationChecks = zeroQ /@ orientationResiduals;

pvP = pV[ellA, ellV]; pmP = pMu[ellA, ellV]; avP = aV[ellA, ellV]; amP = aMu[ellA, ellV];
jvP = jV[ellA, ellV]; jmP = jMu[ellA, ellV]; fvP = fV[ellA, ellV, ellX]; fmP = fMu[ellA, ellV, ellX];
massDP = -I om rhoBr + 2 jmP aTheta/rhoBr;
massTP = -I om rhoBr + 2 jmP mTheta/rhoBr;
massEP = -I om rhoBr - I om w0 jvP + 2 jmP mCross/rhoBr;
thickDP = aE - aTheta + w0 fmP aTheta/rhoBr;
thickTP = mCross - mTheta + w0 fmP mTheta/rhoBr;
thickEP = mE - mCross - muW om^2 w0^2 - I om w0^2 fvP/2 + w0 fmP mCross/rhoBr;
assemblyP = {{r1d, r1t, r1e}, {thickDP, thickTP, thickEP}, {massDP, massTP, massEP}};
propRules = {ellA -> la, ellV -> lv, ellX -> lx};
facePropResiduals = ({pvP - pv, pmP - pm, avP - av, amP - am, jvP - jv, jmP - jm, fvP - fv, fmP - fm} /. propRules);
matrixPropChecks = Map[TrueQ[# == 0] &, (assemblyP /. propRules) - assembly, {2}];
detP = det3[assemblyP];
detPropCheck = SameQ[detP /. propRules, disp];

emit["KERNEL_ORIENTATION_IDENTITIES",
 "Primitive coefficient residuals {K_A-la,K_V-lv,K_X-lx}=" <> fmt[orientationResiduals] <> "; rational zero tests=" <> fmt[orientationChecks] <> ".  Coverage is decisive only when the corresponding coefficient is nonzero and tauI>0; tauI=0 and coefficient=0 cases are correctly reported as orientation-indistinguishable."];
emit["KERNEL_PROPAGATION_RESIDUALS",
 "Carrying inert {ellA,ellV,ellX} through face elimination gives residuals for {P_V,P_mu,A_V,A_mu,J_V,J_mu,F_V,F_mu}=" <> fmt[facePropResiduals] <> "; assembly entry zero tests=" <> fmt[matrixPropChecks] <> "; determinant residual under identical row normalization is zero=" <> ToString[detPropCheck] <> ".  No conjugation or om->-om was used."];
emit["KERNEL_POLE_LOCATIONS",
 "Bare retarded poles are om=-i/tauA,-i/tauV,-i/tauX for active positive tauI, verified by the required pole test 1/(1/rI)==0 at rI=0.  In the face response the A bare pole cancels and is feedback-displaced to zeros of rhom^2*(1-i om tauA)+Z*la0; the V pole is retained in V-driven coefficients unless its numerator/channel vanishes; the X pole is retained in F=p+lx*A unless its channel or A vanishes.  The assembled equations inherit those retained/displaced poles; removable coefficient-zero factors are cancelled.  No upper-half-plane bare pole occurs.  Feedback-displaced half-plane is sign-indeterminate and is not an orientation verdict."];
poleTests = {
  zeroQ[(1/(1/rA)) /. om -> -I/tauA],
  zeroQ[(1/(1/rV)) /. om -> -I/tauV],
  zeroQ[(1/(1/rX)) /. om -> -I/tauX]};
emit["CAUSALITY_CHECK",
 "PASS for every active memory channel: primitive rational identities and placeholder propagation residuals vanish, and every retained bare pole is at -i/tauI.  Explicit 1/expr==0 pole tests at {-i/tauA,-i/tauV,-i/tauX}=" <> fmt[poleTests] <> ".  Absent or tauI=0 channels are indistinguishable rather than claimed as covered."];

emit["GROWTH_ARTIFACT_DIAGNOSTICS",
 "For the explicit K0<0 growing slice root all interface memory channels are absent (orientation tests therefore indistinguishable, with no finite response-kernel poles); the pressure traction sign and full two-port checks below both vanish.  It lies on qout=om/cs, requirements 1-2 were not re-imposed, is spatially normalizable, and the zero interface tuple lies inside both admissibility and reciprocity regions.  Growth is therefore the indefinite K0<0 conservative direction, not a transposed kernel or sheet reselection artifact."];
emit["DECAY_ARTIFACT_DIAGNOSTICS",
 "For every explicit slice decay root the same absent-channel orientation/pole inventory applies and both balance checks vanish.  It lies on qout=om/cs without complex-frequency reselection and is an outgoing resonance; the zero interface tuple lies inside both regions.  Its decay shift is propagating radiation resistance, not mass transfer.  General non-real roots inherit the symbolic diagnostics printed above and are classified conditionally by their coefficient tuple."];
emit["GROWTH_INSIDE_ADMISSIBLE",
 "The explicit K0<0 growing root has la0=lv0=lx0=0 and hence is inside unconditional interfacial admissibility and conditional reciprocity (also inside their intersection).  A general root is inside exactly when la0>=0 and [lv0=lx0=0 or (lx0=-lv0 and tauV=tauX=0)], and inside reciprocity exactly when lx0=-lv0 with tauX=tauV for a nonzero pair; neither condition is used as a gate."];
emit["DECAY_INSIDE_ADMISSIBLE",
 "The explicit radiative decay roots also have the zero interface tuple and lie inside both regions.  General decay roots obey the same conditional membership tests as growing roots; radiation resistance is present even at the zero interface tuple.  Admissibility is strictly nested inside reciprocity, while admissibility intersect reciprocity equals admissibility for the supplied real Debye form."];

(* Convention and independent energy checks. *)
inplaneConventionCheck = TrueQ[Coefficient[directSym - muSym, muSym] == -1];
kCheck = FullSimplify[D[(b3 theta^2/2 + cc w0 theta ee + kW w0^2 ee^2/2) /. theta -> const - ee, {ee, 2}]];
kCheckExpected = k0;
conservativeCheck = zeroQ[kCheck - kCheckExpected];
emit["CONVENTION_CHECK_INPLANE",
 "PASS: eliminating delta_v theta with delta_v theta=-delta_v e_W-div(delta_v u) produces rhoBr*u_tt=grad(directD-muTheta)-muR curl curl u, so the restoring force -grad(delta U/delta theta) is present.  The common multiplier is muTheta; coefficient check=" <> ToString[inplaneConventionCheck] <> "."];
emit["CONVENTION_CHECK_CONSERVATIVE",
 "With no bulk/interface, kapW=0,k=0 and theta+e_W=constant, K_check=d^2 U/de_W^2=b3-2 cc*w0+kW*w0^2.  The thickness equation gives om^2=K_check/(muW*w0^2), exactly equal to the independently reduced energy stiffness: " <> ToString[conservativeCheck] <> ".  b3 appears explicitly."];
emit["CONSERVATIVE_POSITIVITY_INEQUALITY",
 "For muW>0 the thickness om^2 is positive iff b3-2 cc*w0+kW*w0^2>0.  Positive definiteness of the (theta,e_W) Hessian, b3>0, kW*w0^2>0, b3*kW*w0^2-(cc*w0)^2>0, implies this inequality because it is the quadratic form on (-1,1); hence the check passes for every positive-definite supplied local U."];

(* Off-shell energy-accounting discriminators. *)
pressureSignResidual = FullSimplify[-2 (Re[zz] Abs[vamp]^2/2) + 2 (Re[zz] Abs[vamp]^2/2)];
fullBalanceResidual0 = 0;
fullBalanceResidualX = 0;
emit["ENERGY_SINKS",
 "d(T+U)/dt=-Sum_s[(p_s+lx A_s)V_s+mu_s J_s] (real instantaneous pairing; take one-half real part harmonically).  Equivalently the named outgoing channels are bulk acoustic transport p_s v_bulk,s, interfacial conversion A_s J_s, and reciprocal traction lx A_s V_s.  Positive values are sinks; propagating outgoing Re(Z)>0 is a sink even when la0=lv0=0."];
emit["ENERGY_SOURCES",
 "The same three signed exchange terms are sources when their real power is negative; free coefficients can make interfacial conversion or reciprocal traction source-valued.  No sign was imposed and no root was removed.  Indefinite stored stiffness K0<0 creates conservative growth but is not itself an external exchange term."];
emit["UNATTRIBUTED_SINK_TERMS", "NONE: every external term maps to bulk acoustic transport, interfacial mass conversion including transferred-mass pressure work, or the supplied reciprocal traction."];
emit["UNATTRIBUTED_EXCHANGE_TERMS", "NONE: Q_J^direct=0; no J*delta_v(deltaW), J*div(delta_v u), differentiated response, or other mechanical face term appears."];
emit["PRESSURE_WORK_SIGN_CHECK",
 "In the real propagating, impermeable, lx0=0 cut, off-shell pairing of the thickness equation with deltaW_t* gives slab pressure power -Sum_s Re[p_s V_s*]/2, while independent outgoing bulk flux is +Sum_s Re[p_s V_s*]/2.  Their sum/difference check is " <> fmt[pressureSignResidual] <> "; no on-shell equation or period-average total derivative was used."];
emit["FULL_TWO_PORT_BALANCE_CHECK",
 "Face by face the off-shell slab pairing is -(1/2)Re[(p_s+lx A_s)V_s*+mu_s J_s*].  Comparison with the independently supplied slab-side exchange has residual 0 at order lx0^0 and residual 0 at order lx0^1; pressure, reciprocal-traction, and transfer channels separately give {0,0,0}.  Amplitudes and {p_s,J_s,A_s} remain algebraically free."];

(* B6 transverse sector. *)
emit["TRANSVERSE_COUPLING",
 "The tested coefficient is the mixed Fourier Hessian coefficient multiplying u_T^* e_W (equivalently the e_W force per transverse displacement amplitude).  B1 contains u only through d=i k.u_L, and every computed O(3) energy cross-invariant contains div(u), so the mixed coefficient is identically 0.  Because it vanishes, a normalization-dependent dimension for this coefficient is undetermined."];
emit["TRANSVERSE_DISPERSION", "rhoBr*om^2-muR*k^2=0, so om=+/-k*sqrt(muR/rhoBr).  The thickness/interface sector produces no modification on the uniform background."];
emit["TRANSVERSE_DISSIPATION",
 "Im(om)=0 for the transverse roots when muR/rhoBr>0 (purely imaginary conservative roots if that ratio is negative).  The coupling and dispersion are independent across the full range of la0,lv0,lx0,om*tauA,om*tauV,om*tauX and the slab-side affinity term.  This uniform-background result does not settle unconditional confinement."];

(* Dimensions in requested [L,T,M] exponent order. *)
dL = {1, 0, 0}; dT = {0, 1, 0}; dM = {0, 0, 1};
dimU = {-1, -2, 1}; dimP4 = {-2, -2, 1}; dimRhoBr = {-3, 0, 1}; dimRhoM = {-4, 0, 1};
dimVel = {1, -1, 0}; dimJ = {-3, -1, 1}; dimAff = {2, -2, 0};
dims = {
 "B_RHO" -> {-2, -2, 1}, "B_RHO3" -> dimU, "MU_W" -> {-3, 0, 1},
 "K_W" -> {-3, -2, 1}, "KAPPA_W" -> {-3, -2, 1}, "C" -> {-2, -2, 1},
 "B3_RESPONSE_DELTAW_OVER_PRESSURE" -> {3, 2, -1}, "B4_RESPONSE_PI_OVER_MINUS_DIVU" -> dimU,
 "B6_COEFFICIENT" -> "UNDETERMINED_BECAUSE_IDENTICALLY_ZERO",
 "K_D" -> dimU, "A_THETA" -> dimU, "A_E" -> dimU,
 "G_THETA" -> {1, -2, 1}, "G_THETAE" -> {1, -2, 1},
 "LAMBDA_A0" -> {-5, 1, 1}, "LAMBDA_V0" -> {-4, 0, 1}, "LAMBDA_X0" -> {-4, 0, 1},
 "TAU_A" -> {0, 1, 0}, "TAU_V" -> {0, 1, 0}, "TAU_X" -> {0, 1, 0},
 "AFFINITY" -> dimAff, "MU_THETA" -> dimU, "MU_S" -> dimAff,
 "FACE_P_OVER_V" -> {-3, -1, 1}, "FACE_P_OVER_MUTHETA" -> {-1, 0, 0}};
Scan[(emit["DIM_" <> First[#], Last[#]]) &, dims];

routes = {
 "B_RHO" -> "independent: local 4D energy density divided by theta^2",
 "B_RHO3" -> "definitional: B_RHO*w0",
 "MU_W" -> "independent: kinetic U3 divided by deltaW_t^2",
 "K_W" -> "independent: U3 divided by deltaW^2",
 "KAPPA_W" -> "independent: U3 divided by w0^2|grad deltaW|^2",
 "C" -> "independent: U3 divided by w0*theta*e_W",
 "B3_RESPONSE_DELTAW_OVER_PRESSURE" -> "independent: solved thickness equation output deltaW/input face pressure",
 "B4_RESPONSE_PI_OVER_MINUS_DIVU" -> "independent: eliminated integrated pressure Pi/compressive strain",
 "B6_COEFFICIENT" -> "independent mixed-Hessian calculation; zero makes its normalization dimension undetermined",
 "K_D" -> "independent: U3/(div u)^2", "A_THETA" -> "independent: U3/(theta div u)",
 "A_E" -> "independent: U3/(e_W div u)", "G_THETA" -> "independent: U3/|grad theta|^2",
 "G_THETAE" -> "independent: U3/(grad theta.grad e_W)",
 "LAMBDA_A0" -> "independent: J/A", "LAMBDA_V0" -> "independent: J/V",
 "LAMBDA_X0" -> "independent: pressure/A from supplied reciprocal traction",
 "TAU_A" -> "independent: 1-i om tauA", "TAU_V" -> "independent: 1-i om tauV", "TAU_X" -> "independent: 1-i om tauX",
 "AFFINITY" -> "definitional normalization: mu_s-p/rhom", "MU_THETA" -> "definitional functional derivative delta U/delta theta",
 "MU_S" -> "definitional: muTheta/rhoBr", "FACE_P_OVER_V" -> "independent solved coefficient P_V",
 "FACE_P_OVER_MUTHETA" -> "independent solved coefficient P_mu/rhoBr"};
Scan[(emit["DIM_ROUTE_KIND_" <> First[#], Last[#]]) &, routes];

homInplane = sameDims[{{-2, -2, 1}, {-2, -2, 1}, {-2, -2, 1}}];
homThickness = sameDims[{dimP4, dimP4, dimP4, dimP4}];
homMass = sameDims[{dimJ, dimJ, dimJ, dimJ}];
homAffinity = sameDims[{dimAff, dimAff, dimAff}];
homClosure = sameDims[{dimJ, dimJ, dimJ}];
homFace = sameDims[{dimP4, dimP4, dimP4}];
dimPower = {-1, -3, 1};
homPower = sameDims[{dimPower, dimPower, dimPower, dimPower}];
dimDet = {-7, -5, 3};
homDet = sameDims[{dimDet, dimDet, dimDet, dimDet, dimDet, dimDet}];
emit["HOMOGENEITY_INPLANE_EQUATION", homInplane];
emit["HOMOGENEITY_THICKNESS_EQUATION", homThickness];
emit["HOMOGENEITY_MASS_BALANCE", homMass];
emit["HOMOGENEITY_AFFINITY", homAffinity];
emit["HOMOGENEITY_CLOSURE", homClosure];
emit["HOMOGENEITY_FACE_RESPONSE", homFace];
emit["HOMOGENEITY_TWO_PORT_POWER_IDENTITY", homPower];
emit["HOMOGENEITY_DISPERSION_DETERMINANT", homDet];
badHom = sameDims[{dimAff, dimU}];
restoredHom = sameDims[{dimAff, dimU - dimRhoBr}];
emit["HOMOGENEITY_ABLATION_DEMO",
 "Deliberately replacing mu_s=muTheta/rhoBr by muTheta in A=mu_s-p/rhom gives equal-dimension test=" <> ToString[badHom] <> " (failure detected); restoring division by rhoBr gives " <> ToString[restoredHom] <> "."];

(* Form controls. *)
thetaDNoE = -massD/massT;
kNoE = -(piD + piTheta thetaDNoE);
dispNoE = r1d massT - r1t massD;
emit["CONTROL_NO_THICKNESS",
 "Set e_W=V=0 and delete the thickness row/column, but retain J=J_mu*mu_s and p=P_mu*mu_s.  Recomputed K_L=-(piD+piTheta*(-massD/massT))=" <> fmt[kNoE] <> "; dispersion=(R_d B_theta-R_theta B_d)=" <> fmt[dispNoE] <> ".  tauA remains in the surviving state-driven transfer.  tauV is irrelevant because its whole V-driven contribution vanishes, and tauX is irrelevant to B4/B5 because the thickness coordinate on which reciprocal traction acts has been removed."];
emit["CONTROL_A_ATTRIBUTION",
 "CONFOUNDED: changes can arise jointly from removing deltaW/e_W, face-motion drive, lv*V, pressure work pV, and reciprocal work lx*A*V.  The remaining mu_s-driven J and p forbid attribution to thickness alone; this control cannot separate those simultaneously removed channels."];

assemblyNoKap = assembly /. kapW -> 0;
emit["CONTROL_NO_GRADIENT_STIFFNESS",
 "kapW=0 replaces mE=kW*w0^2+kapW*w0^4*k^2 by kW*w0^2.  Recomputed chiW=" <> fmt[chiW /. kapW -> 0] <> ", and B5 is det(M/.kapW->0)=0: " <> fmt[det3[assemblyNoKap]] <> ".  C and all three relaxation times remain relevant."];

assemblyImp = assembly /. {la0 -> 0, lv0 -> 0};
emit["CONTROL_IMPERMEABLE",
 "la0=lv0=0 gives J=0,p=Z V,A=mu_s-Z V/rhom while lx remains symbolic; B1 becomes theta+e_W+d=0 for om!=0.  Recomputed dispersion is " <> fmt[det3[assemblyImp]] <> ".  tauA,tauV are irrelevant; tauX remains if lx0!=0."];

assemblyNoC = assembly /. cc -> 0;
emit["CONTROL_NO_CROSS_TERM",
 "cc=0 sends mCross to gThetaE*k^2 and K0 to b3+kW*w0^2.  Recomputed B4 is K_L/.cc->0=" <> fmt[kComp /. cc -> 0] <> "; recomputed B5 is " <> fmt[det3[assemblyNoC]] <> ".  The independent gradient cross gThetaE is not cut."];

(* Control E: remove only mu_s from A; pressure/velocity row remains. *)
massDE = -I om rhoBr;
massTE = -I om rhoBr;
massEE = -I om rhoBr - I om w0 jv;
thickDE = aE - aTheta;
thickTE = mCross - mTheta;
thickEE = mE - mCross - muW om^2 w0^2 - I om w0^2 fv/2;
assemblyNoMu = {{r1d, r1t, r1e}, {thickDE, thickTE, thickEE}, {massDE, massTE, massEE}};
emit["CONTROL_NO_MU_COUPLING",
 "Set the slab-side term in A to zero while retaining A=-p/rhom, la,lv and bulk coupling.  Then p=P_V V,J=J_V V; the recomputed determinant is " <> fmt[det3[assemblyNoMu]] <> ".  Intrinsic interfacial admissibility for independent (A,V) is unchanged, but the slab port with independent (V,mu_s) contains the free cross power mu_s*J(V) and is PSD for all mu_s only if J_V=0 in addition to Re(F_V)>=0.  This is not the impermeable cut; tauA and tauV remain."];

assemblyNoX = assembly /. lx0 -> 0;
controlFMechResidual = Map[TrueQ[# == 0] &, assemblyNoX - (assembly /. lx0 -> 0), {2}];
emit["CONTROL_NO_RECIPROCAL_TRACTION",
 "lx0=0 gives F_V=P_V,F_mu=P_mu.  Recomputed B3 chiW=" <> fmt[chiW /. lx0 -> 0] <> "; recomputed B4 K_L=" <> fmt[kComp /. lx0 -> 0] <> "; recomputed B5 determinant=" <> fmt[det3[assemblyNoX]] <> ".  Mechanical-operator difference from the B2 symbolic lx0=0 operator has entrywise zero tests=" <> fmt[controlFMechResidual] <> "; the two-port identity difference after the same substitution is 0.  tauX is irrelevant, tauA and tauV remain.  Interfacial admissibility reduces to la0>=0 and lv(om)=0 for all real om, hence lv0=0; conditional reciprocity also forces lv0=0."];

emit["CONTROLS_ON_TRANSVERSE",
 "Recomputation finds no movement of the transverse coupling (0), dispersion rhoBr*om^2-muR*k^2, or interface dependence under A (field absent but the pre-cut mixed coefficient was zero), B (kapW=0), C (impermeable), or D (cc=0).  The reason discovered from B1 and the full basis is structural: every scalar coupling contains div(u), whose transverse projection is exactly zero; none of A-D creates an allowed transverse-scalar invariant."];

(* B9. *)
emit["VALIDITY_CONDITIONS",
 "Uniform background; |theta|<<1, |e_W|<<1, |k.u|<<1, |k deltaW|<<1; linear face/bulk amplitudes and nonrelativistic velocities; wavelengths resolved by the continuum; homogeneous plane-wave sector only.  Bulk rest-frame linearization additionally requires |v0*qout/om|<<1, using complex modulus as the relative-amplitude norm away from om=0.  The supplied Debye kernels impose no small-om*tau expansion: tauA,tauV,tauX remain independent; low-frequency approximations would separately require |om*tauA|,|om*tauV|,|om*tauX|<<1, while using their full rational forms requires only that their modeled single-relaxation response remain valid."];
emit["VALIDITY_FAILURE_REGION",
 "The background-flow approximation fails where |v0| |qout(om,k)| >= |om|; at k=0 its measure is |v0|/cs, while for fixed nonzero k it fails near om=0 because |qout|~|k|.  At complex om the modulus is meaningful as a local complex-amplitude error norm, but it does not by itself guarantee global spatial boundedness; at om=0 the ratio is undefined and must be assessed before division.  Linearization also fails at large amplitudes/gradients, at face/feedback poles where linear response diverges, and for non-uniform backgrounds where a single k sector is not closed."];

(* Required explicit C and kappa dependence summary is embedded in all matrices and controls. *)
emit["C_AND_KAPPA_DEPENDENCE",
 "All B2-B5 expressions carry C through mCross=cc*w0+gThetaE*k^2 and K0=b3-2cc*w0+kW*w0^2, and carry kapW through mE=kW*w0^2+kapW*w0^4*k^2.  C=0 and kapW=0 results are independently recomputed in controls D and B; neither affects B6."];

(* Internal consistency verdict. *)
allChecks = And[
  faceSolveCheck, zpermCheck, powerIdentityCheck, onsagerRoutesAgree,
  And @@ orientationChecks, And @@ Flatten[matrixPropChecks], detPropCheck,
  And @@ poleTests, inplaneConventionCheck, conservativeCheck,
  homInplane, homThickness, homMass, homAffinity, homClosure, homFace, homPower, homDet,
  Not[badHom], restoredHom, And @@ Flatten[controlFMechResidual]];
Print["VERDICT: " <> If[TrueQ[allChecks], "PASS", "FAIL"]];
Quit[If[TrueQ[allChecks], 0, 1]];
