(* A selected helicity line is an additional constitutive hypothesis.
   Here 'a' converts the orthonormal W one-form to the Berry connection;
   outward orientation and A=i e^dagger d e fix all signs. *)
ClearAll["Global`*"];
emit[label_,value_] := Print[label," = ",ToString[value,InputForm]];
checks = {};
compare[label_,l_,r_,assumptions_:True] := Module[{res},
  res = FullSimplify[l-r,assumptions];
  emit["check " <> label <> " {left,right,residual}",{l,r,res}];
  AppendTo[checks,AllTrue[Flatten[{res}],#===0&]]
];
realSimp[expr_] := Simplify[ComplexExpand[expr]];
inputDir = "software/em_charge_attribute/reports/";
files = inputDir <> # & /@ {"mouth_count_unwinding.md","trapped_wave_charge.md"};
reportLines = StringSplit[Import[#,"Text"],"\n"]& /@ files;
emit["saved inputs and SHA256", {#,FileHash[#,"SHA256","HexString"]}& /@ files];
emit["saved flare and doublet equations",{
  SelectFirst[reportLines[[1]],StringContainsQ[#,"w=H_*+"]&],
  SelectFirst[reportLines[[2]],StringStartsQ[StringTrim[#],"E=\\omega(I_+"]&]}];

radial={Sin[th] Cos[ph],Sin[th] Sin[ph],Cos[th]};
polar=D[radial,th];
azimuth=Simplify[D[radial,ph]/Sin[th]];
e[h_]:=(polar+I h azimuth)/Sqrt[2];
connection[state_]:=Table[realSimp[I Conjugate[state].D[state,coord]],{coord,{th,ph}}];
curve[conn_]:=Simplify[D[conn[[2]],th]-D[conn[[1]],ph]];
localData[h_]:=Module[{ev=e[h],conn,f,cnum,cn,cs,rn,rs,overlap},
  conn=connection[ev];
  f=curve[conn];
  cnum=Integrate[f,{th,0,Pi},{ph,0,2 Pi}]/(2 Pi);
  (* Remove the polar-frame phase at the corresponding regular pole. *)
  cn=Simplify[conn-{0,h}];
  cs=Simplify[conn-{0,-h}];
  rn=realSimp[Exp[I h ph] ev /. th->0];
  rs=realSimp[Exp[-I h ph] ev /. th->Pi];
  overlap=Integrate[cn[[2]]-cs[[2]],{ph,0,2 Pi}];
  emit["selected {h,e,A_theta/phi,F_theta_phi,Chern number}",{h,ev,conn,f,cnum}];
  emit["two-chart {h,regular north/south states at poles,north/south connections,overlap integral}",
    {h,rn,rs,cn,cs,overlap}];
  compare["north/south Stokes for h="<>ToString[h],overlap/(2 Pi),cnum];
  compare["selected helicity eigenvalue for h="<>ToString[h],
    realSimp[I Cross[radial,ev]],h ev];
  {f,cnum}
];
lineData=localData /@ {-1,1};
(* Non-Abelian Hermitian connection corresponding to d-iA for
   A=i V^dagger dV; curvature dA-i[A,A]. *)
basis=Transpose[{e[-1],e[1]}];
matrixA[coord_]:=Map[realSimp,I ConjugateTranspose[basis].D[basis,coord],{2}];
at=matrixA[th]; ap=matrixA[ph];
fmat=Map[realSimp,D[ap,th]-D[at,ph]-I (at.ap-ap.at),{2}];
emit["complete doublet {A_theta,A_phi,F_theta_phi,trace curvature,trace Chern number}",
  {at,ap,fmat,Tr[fmat],Integrate[Tr[fmat],{th,0,Pi},{ph,0,2 Pi}]/(2 Pi)}];
compare["doublet diagonal versus separate circular lines",Diagonal[fmat],lineData[[All,1]]];

(* Parameter sphere of one closed material section, with t=r/R.
   Pull back the signed area form using the oriented triple product
   det(raw,raw_theta,raw_phi)/|raw|^3. *)
raw=t radial-{0,0,beta};
norm2=Simplify[raw.raw];
triple=Simplify[raw.Cross[D[raw,th],D[raw,ph]]];
area=triple/norm2^(3/2);
nz=(t Cos[th]-beta)/Sqrt[norm2];
interior=t>0 && beta>=0 && 0<th<Pi && norm2>0;
emit["section {projected vector, squared norm, oriented triple, signed area density, n_z}",
  {raw,norm2,triple,area,nz}];
compare["area from normalized polar direction",area,-D[nz,th],interior];
factors=Simplify[#/Sin[th]]& /@ lineData[[All,1]]; (* coeff of target signed area *)
sector[assumptions_]:=Module[{endpoints,degree,cnums},
  endpoints=FullSimplify[{nz /. th->0,nz /. th->Pi},assumptions];
  degree=Simplify[(endpoints[[1]]-endpoints[[2]])/2];
  cnums=degree factors 2;
  emit["regular section {conditions,n_z at north/south,degree,Chern(-,+),exterior flux(-,+)}",
    {assumptions,endpoints,degree,cnums,2 Pi cnums/a}];
  cnums
];
inside=sector[t>0 && 0<=beta<t];
outside=sector[t>0 && beta>t];
emit["flux change with increasing section radius across t=beta, h order {-1,+1}",
  2 Pi (inside-outside)/a];
(* At equality the vector vanishes at the north pole. Compute distinct
   meridian limits without assigning a Chern number to the punctured map. *)
zeroAtNorth=Simplify[raw /. {th->0,beta->t}];
limitNormal=FullSimplify[
  Limit[(raw/Sqrt[norm2] /. beta->t),th->0,Direction->"FromAbove",Assumptions->t>0],
  t>0 && Element[ph,Reals]];
meridians={limitNormal /. ph->0,limitNormal /. ph->Pi/2};
emit["t=beta {north vector,meridional limit,limits at phi=0,Pi/2}",
  {zeroAtNorth,limitNormal,meridians}];
compare["squared norm at regular-section north pole",norm2 /. th->0,(t-beta)^2];

(* Reflect the projected normal with section orientation held fixed.
   The helicity projector depends both on h and on the oriented normal.
   This explicitly distinguishes holding helicity and holding the line. *)
mirrorTriple=Simplify[(-raw).Cross[D[-raw,th],D[-raw,ph]]];
crossMat[n_]:=Transpose[Cross[n,#]& /@ IdentityMatrix[3]];
projector[n_,h_]:=(IdentityMatrix[3]-Outer[Times,n,n]+I h crossMat[n])/2;
doMirror[h_]:=Module[{pr,mir,rev},
  pr=projector[radial,h]; mir=projector[-radial,h]; rev=projector[-radial,-h];
  compare["projector from circular state for h="<>ToString[h],
    Map[realSimp,Outer[Times,e[h],Conjugate[e[h]]],{2}],pr];
  emit["mirror {h,projector difference at fixed h at north,projector difference reversing h}",
    {h,Map[realSimp,(mir-pr) /. th->0,{2}],Simplify[rev-pr]}]
];
Scan[doMirror,{-1,1}];
compare["mirror reverses oriented area numerator",mirrorTriple,-triple];
emit["mirrored Chern(-,+) by sector {t>beta,t<beta}, helicity held",
  -( {inside,outside} )];
emit["mirrored Chern(-,+) by sector {t>beta,t<beta}, helicity reversed",
  (-Reverse[#]&)/@{inside,outside}];

(* Localized relative phase stiffness, only on the ordered surface:
   v=d eta-a W_parallel+A_h. On the complement of phase zeros both
   connections can be taken regular on each small disk, so only d eta
   survives its shrinking positive boundary integral. The exterior
   boundary orientation of the punctured surface is the opposite one. *)
cv={et1,et2}-a {ww1,ww2}+{ah1,ah2};
emit["surface relative term {covariant components in orthonormal frame,energy density}",
  {cv,stiff amp^2 cv.cv/2}];
loop=Integrate[D[m ph,ph],{ph,0,2 Pi}];
regularBound=Limit[2 Pi smallR maxNorm,smallR->0,Direction->"FromAbove"];
smoothCurl=-a surfaceFlux+2 Pi chern;
orientedBoundary=-(loop /. m->totalWind);
fluxRule=First[Solve[smoothCurl==orientedBoundary,surfaceFlux]];
emit["punctured-surface Stokes {positive d-eta loop,bound on shrinking regular loop,bulk dv,oriented boundary,normalized flux}",
  {loop,regularBound,smoothCurl,orientedBoundary,Factor[a surfaceFlux/(2 Pi) /. fluxRule]}];
emit["smooth everywhere-nonvanishing relative phase, normalized flux",
  Factor[a surfaceFlux/(2 Pi) /. fluxRule /. totalWind->0]];
compare["oriented puncture boundary versus dv after integrated constraint",
  smoothCurl /. fluxRule,orientedBoundary];
emit["all consistency checks",checks];
Exit[If[And@@checks,0,1]];
