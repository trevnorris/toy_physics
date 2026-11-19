(* Free exterior connection with flux 2 Pi n, unit-sphere angular averages.
   t=r/R, f=R Sqrt[u/K] Psi, b=r^2 B_r-n/2, m=gamma u/K^2.
   In units 4 Pi K^2/(u R), B3 is g^2 w/2, w=1-1/Lambda. *)
ClearAll["Global`*"];
emit[label_, x_] := Print[label <> " = " <> ToString[x, InputForm]];
avg[q_] := Integrate[q, {z, 0, 1}];
rad[q_, l_] := Integrate[q, {t, 1, l}];

(* Direct volume integral of the gradient of the vacuum potential with
   R=1, epsilon=1 and b(1)=P_l(cos theta). A source-free finite-energy
   competitor minus this extension has zero flux trace and zero inner
   product with its gradient (integration by parts). *)
vacuum[l_] := Module[{potential, radial, transverse},
  potential = LegendreP[l, x]/((l+1) t^(l+1));
  radial = Integrate[t^2 D[potential, t]^2/2, {t, 1, Infinity}, {x,-1,1}];
  transverse = Integrate[(1-x^2) D[potential,x]^2/2,
                         {t,1,Infinity}, {x,-1,1}];
  {radial, transverse, radial+transverse}
];
vacuumModes = vacuum /@ {0,1,2};
emit["vacuum radial,tangential,total costs / (4 Pi gamma epsilon^2/R), l=0,1,2",
     vacuumModes];

(* Axially symmetric section in a real-amplitude regular gauge, A_phi
   includes the native connection. No material trial varies radially. *)
kin[f_, az_] := Together[z (1-z) D[f,z]^2 + f^2 az^2/(4 z (1-z))];
neutralK = avg[kin[z,2(1-z)+2 e z (1-z)] + kin[1-z,2 z+2 e z(1-z)]];
neutralN = avg[z^2+(1-z)^2];
neutralQ = avg[(z^2+(1-z)^2)^2];
fUnit = Sqrt[z] (1+(1-z)/2);
unitK = avg[kin[fUnit,1-z]];
unitN = avg[fUnit^2];
unitQ = avg[fUnit^4];
emit["neutral and unit trial moments (kinetic,norm,quartic)",
     {{neutralK,neutralN,neutralQ},{unitK,unitN,unitQ}}];
w[l_] := 1-1/l;
vol[l_] := rad[t^2,l];
trialGain[l_,g_,intGain_,q_] := intGain^2/(g^2 w[l] vol[l] q);
unitUpper[l_,g_,y_] := y-trialGain[l,g,(l-1)(g unitN-unitK),unitQ];
oldL=5/4; oldG=3/2; oldY=2/9; eps=-2/3;
oldGain = rad[oldG neutralN-(neutralK /. e->eps/t),oldL];
oldDipole = 4 oldY eps^2 vacuumModes[[2,3]];
oldNeutral = oldDipole-trialGain[oldL,oldG,oldGain,neutralQ];
oldUnit = unitUpper[oldL,oldG,oldY];
emit["original equal-gain neutral trial (dipole,integrated gain,exterior,energy)",
     {eps,oldGain,oldDipole,oldNeutral}];
emit["original neutral trial,uniform unit trial,and difference",
     N[{oldNeutral,oldUnit,oldUnit-oldNeutral},16]];

(* On a shell, the integrated chiral-square identities bound the
   material angular gradients below by (signed curvature) times density
   for either choice of chirality, even when the curvature changes sign.
   For opposite chiralities let Z=rho_s-rho_w, S=rho_s+rho_w.
   Charge conservation gives <b>=0: the cross term b Z equals b (Z-<Z>)
   under the integral. At m=1/2, completing the b-square leaves
   t^2/2(<S^2-Z^2>+<Z>^2); the first term is nonnegative. *)
covarianceResidual = Expand[bb (zz-zMean) + bb^2/(2t^2) +
    t^2(zz-zMean)^2/2 - (bb+t^2(zz-zMean))^2/(2t^2)];
densityResidual = Factor[(rhoS+rhoW)^2-(rhoS-rhoW)^2];
emit["magnetic square residual,and remaining local density polynomial",
     {covarianceResidual,densityResidual}];

(* Complete the radial kinetic term with D_t f - p f. The connection
   drops out of Re(conj(f)D_t f). Vanishing p at both ends removes the
   boundary term, even for free traces of f. Set
     v=-(t^2 p)'-t^2 p^2.
   The resulting lower bound contributes v times the shell mean density. *)
aspect=5/3; strong=9/5; weak=1/10; m=1/2;
y = m/(2 strong^2 w[aspect]);
c = (aspect-1)/vol[aspect];
p = Factor[Integrate[-a(1-c s^2),{s,1,t}]/t^2];
vPlus = -D[t^2 p,t]-t^2 p^2;
vMinus = D[t^2 p,t]-t^2 p^2;
residualA = Factor[a-vPlus-t^2(a c+p^2)];
d = (aspect-1)^2 (aspect+2)/(4(1+aspect+aspect^2));
radialChecks = {p /. t->1, p /. t->aspect, residualA,
  Resolve[ForAll[{a,t}, a>=0 && 1<=t<=aspect, -a d <= t p <= a d],Reals]};
emit["aspect,strong,weak,m,y,c,d,and radial p", {aspect,strong,weak,m,y,c,d,p}];
emit["radial endpoint values,linear-coefficient residual,and uniform bound",
     radialChecks];
(* t>=1, |t p|<=a d implies 0<=a c+p²<=a c+a² d². Minimizing
   t² Zmean²/2 -(a-vPlus) Zmean yields a lower bound even for signed
   Zmean. For the tail use Zmean=<S> and ignore positive angular terms.
   Neutral uses opposite p for strong/weak densities. *)
bindingBound[a0_] := Factor[vol[aspect](a0 c+a0^2 d^2)^2/
                            (strong^2 w[aspect])];
aZero=strong-1;
weakLeftover = Factor[2-strong-weak+vPlus+vMinus];
weakLower = 2-strong-weak-2(aZero d)^2;
emit["neutral residual weak coefficient,and positive lower bound",
     {weakLeftover,weakLower}];
neutralLower = -bindingBound[aZero];
tailLower = 4y-bindingBound[strong];
(* For n opposite to the strong native sign, the two angular curvatures
   have uniform parts of magnitude 3/2 and 1/2. Choose p=0 here. *)
oppA = strong-3/2;
oppositeLower = y-(oppA/strong)^2;
oppositeWeak = 2-strong-weak;
alignedUpper = unitUpper[aspect,strong,y];
floors = {neutralLower,oppositeLower,tailLower};
emit["alternate aligned unit upper; floors for neutral,opposite unit,|n|>=2",
     Prepend[floors,alignedUpper]];
emit["alternate floors minus aligned unit upper", floors-alignedUpper];
emit["alternate numerical energies", N[Prepend[floors,alignedUpper],14]];
(* Raising m preserves the fluctuation bound above. Its mean monopole
   contribution is exactly y(m)n², so every interval comparison is affine
   in m. Endpoint tests cover the closed interval. *)
yAt[m0_] := m0/(2 strong^2 w[aspect]);
intervalGaps[m0_] := Module[{shift=yAt[m0]-y},
  floors+{0,shift,4 shift}-(alignedUpper+shift)];
endpoints={m,51/100};
endpointGaps=intervalGaps /@ endpoints;
emit["m interval endpoints,y values,and endpoint margins",
     {endpoints,yAt /@ endpoints,endpointGaps}];
(* C_h=0 control. A constant neutral field is globally defined. Its
   trial uses the stronger channel. Every charged state has |n|>=1 and
   the same crude material bound used in the tail, with a=strong. *)
flatUpper = -trialGain[aspect,strong,(aspect-1)strong,1];
flatGap = y-bindingBound[strong]-flatUpper;
native[h_,s_] := -2h s;
selected[h_,s_] := Sign[native[h,s]];
mirrorRows=Flatten[Table[
  {h,s,selected[h,s],selected[h,-s],selected[-h,-s]},
  {h,{-1,1}},{s,{-1,1}}],1];
energyScale=2 Pi stiffness^2 strong^2 w[aspect]/(quarticU radius);
yPhysical=Factor[(Pi m0 stiffness^2/quarticU/radius)/energyScale];
emit["zero-native neutral upper,all charged lower,and gap at lower m",
     {flatUpper,y-bindingBound[strong],flatGap}];
emit["favored helicity,orientation,winner,fixed-helicity mirror,swapped mirror",mirrorRows];
emit["physical B3, y", {energyScale,yPhysical}];
checks={covarianceResidual==0, densityResidual==4 rhoS rhoW,
        radialChecks=={0,0,0,True}, weakLower>0,oppositeWeak>=0,
        strong unitN>unitK,oldGain>0,oldNeutral<oldUnit,
        Min[endpointGaps]>0,flatGap>0,yPhysical==yAt[m0],
        And@@(#[[3]]==-#[[4]] && #[[3]]==#[[5]] & /@ mirrorRows)};
emit["checks",checks];
If[!TrueQ[And@@checks],Exit[1]];
