(* Conditional compact completion of the curl-field exterior. flux is
   the measured D flux, a is a nonzero real fixed conversion factor.
   Kphase=stiff is uniform in the modal calculation, around neutral W=0.
   No temporal phase kinetic term or core boundary law is assumed. *)
ClearAll["Global`*"];
emit[label_,value_] := Print[label," = ",ToString[value,InputForm]];
checks = {};
compare[label_,left_,right_,assumptions_:True] := Module[{residual},
  residual = FullSimplify[left-right,assumptions];
  emit["check " <> label <> " {left,right,residual}",{left,right,residual}];
  AppendTo[checks,AllTrue[Flatten[{residual}],#===0&]]
];

(* The azimuthal component is expressed in the unit-vector basis. *)
wNorth = flux (1-Cos[th])/(4 Pi r Sin[th]);
wSouth = -flux (1+Cos[th])/(4 Pi r Sin[th]);
overlapConnection = Simplify[wNorth-wSouth];
dAngle = Simplify[a r Sin[th] overlapConnection];
angle = Integrate[dAngle,ph]; (* irrelevant integration constant set to zero *)
transition = Exp[I angle];
turnFactor = Simplify[(transition /. ph->ph+2 Pi)/transition];
realPeriodic = Reduce[{Cos[v]==1,Sin[v]==0},v,Reals];
fluxUnit = 2 Pi/a;
emit["{W_N phi,W_S phi,phase angular derivative,phase difference,transition,period factor}",
  {wNorth,wSouth,dAngle,angle,transition,turnFactor}];
emit["allowed real v=a flux from the real and imaginary parts of exp(i v)=1", realPeriodic];
emit["{flux expression for winding n,period factor at winding n}",{n fluxUnit,
  FullSimplify[turnFactor /. flux->n fluxUnit,Element[n,Integers]]}];
compare["covariant relative gradient is the same on the overlap",dAngle/(r Sin[th]),
  a overlapConnection];
compare["unit phase closes at integer winding",turnFactor /. flux->n fluxUnit,1,Element[n,Integers]];

(* On each disk a smooth nonzero complex phase z has the smooth winding
   one-form dz/(2 Pi I z). Full-disk regularity lets us use Stokes on its
   boundary. The same orientation is used for the two boundary windings. *)
windForm[field_,parameter_] := D[field,parameter]/(2 Pi I field);
diskExterior[field_] := Together[D[windForm[field,yy],xx]-D[windForm[field,xx],yy]];
forms = Table[windForm[zn[xx,yy],param],{param,{xx,yy}}];
diskCurls = diskExterior /@ {zn[xx,yy],zs[xx,yy]};
diskIntegrals = Integrate[#, {xx,-1,1},{yy,-Sqrt[1-xx^2],Sqrt[1-xx^2]}]& /@ diskCurls;
ratioWind = windForm[north[ph]/south[ph],ph];
differenceWind = windForm[north[ph],ph]-windForm[south[ph],ph];
emit["nonvanishing disk phase {winding form components,exterior derivatives N/S,disk integrals N/S}",
  {forms,diskCurls,diskIntegrals}];
compare["winding form of ratio on common positive-azimuth equator",
  Together[ratioWind],Together[differenceWind]];
emit["transition winding if both phases extend regularly over their caps",
  diskIntegrals[[1]]-diskIntegrals[[2]]];
emit["transition winding from prescribed exterior flux",
  Integrate[Simplify[windForm[transition,ph]],{ph,0,2 Pi}]];

(* In a neutral uniform region take k along z. W multiplies cos(k z)
   and the phase perturbation multiplies sin(k z). Their relative
   spatial gradient thus lies in the same real Fourier quadrature.
   The common positive Fourier norm is omitted. *)
wv = {wx,wy,wz};
vel = {vx,vy,vz};
kv = {0,0,k};
gradAlpha = kv phaseAmp;
rel = gradAlpha-a wv;
curlAmp = Cross[kv,wv];
potential = (curlAmp.curlAmp)/(2 eps)+stiff rel.rel/2;
kinetic = b vel.vel/2;
stiffness = Table[D[potential,wv[[i]],wv[[j]]],{i,3},{j,3}];
inertia = Table[D[kinetic,vel[[i]],vel[[j]]],{i,3},{j,3}];
mix = Table[D[potential,wv[[i]],phaseAmp],{i,3}];
transverse = (stiffness-omega2 inertia)[[1;;2,1;;2]];
pole = Factor[omega2 /. First[Solve[transverse[[1,1]]==0,omega2]]];
staticRoot = Factor[k2 /. First[Solve[(transverse[[1,1]] /. {omega2->0,k^2->k2})==0,k2]]];
emit["neutral uniform Fourier {kinetic,potential,phase-W stiffness block,W stiffness,W inertia}",
  {kinetic,potential,mix,stiffness,inertia}];
emit["{transverse operator,dispersion omega^2,gap squared,omega^2 with Kphase=0,static k^2 root}",
  {transverse,pole,pole /. k->0,pole /. stiff->0,staticRoot}];
compare["both transverse polarizations solve the dispersion relation",
  transverse /. omega2->pole,ConstantArray[0,{2,2}]];
compare["longitudinal phase adjustment leaves transverse energy",
  D[potential /. phaseAmp->a wz/k,wx],D[potential,wx]];

(* An independent compact order phase living on a material circle of
   radius ringR may wind even when the ambient W is identically zero. *)
ringPhase = Exp[I n ph];
ringWinding = Integrate[Simplify[windForm[ringPhase,ph]],{ph,0,2 Pi}];
ambientW = {0,0,0};
ambientD = Curl[ambientW,{xx,yy,zz}];
normal = {Sin[th] Cos[ph],Sin[th] Sin[ph],Cos[th]};
sphereFlux = Integrate[(ambientD.normal) r^2 Sin[th],{th,0,Pi},{ph,0,2 Pi}];
emit["ring control {ring phase,azimuthal unit gradient of lifted angle,phase winding,ambient W,D,sphere flux}",
  {ringPhase,n/ringR,ringWinding,ambientW,ambientD,sphereFlux}];
compare["ring phase closes at integer winding",Exp[I n 2 Pi],1,Element[n,Integers]];
emit["all consistency checks",checks];
Exit[If[And@@checks,0,1]];
