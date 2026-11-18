(* Work in a fixed flat laboratory patch for the conditional advective example.
   The live material ledger does not assign this momentum law to the drain. *)
xyz = {x,y,z};
rho0 = R0[x,y,z]; drho = R1[t,x,y,z];
v0 = Through[{vx,vy,vz}[x,y,z]]; dv = Through[{ux,uy,uz}[t,x,y,z]];
(* Build the directional derivatives of the mass and momentum fluxes from their products. *)
massFlux0 = rho0 v0;
dMassFlux = drho v0 + rho0 dv;
dMom = drho v0[[1]] + rho0 dv[[1]];
dMomFlux = dMassFlux v0[[1]] + massFlux0 dv[[1]];
(* Eliminate the outward carried mass loss with the base and incremental mass balances. *)
j0 = -Div[massFlux0,xyz];
dj = -D[drho,t] - Div[dMassFlux,xyz];
transport = D[dMom,t] + Div[dMomFlux,xyz] + j0 dv[[1]] + dj v0[[1]];
local = rho0 D[dv[[1]],t];
advection = rho0 (v0.Grad[dv[[1]],xyz] + dv.Grad[v0[[1]],xyz]);
densityTerm = drho v0.Grad[v0[[1]],xyz];
remainder = Expand[transport-local-advection-densityTerm];
Print["Conditional mass-loss elimination residual in sideways momentum: ",remainder];

(* Time-domain check: a unit step of incremental shear strain gives an
   instantaneous stress jump, followed by relaxation at fixed strain. *)
step = DSolveValue[{tau sig'[t]+sig[t]==0,sig[0]==modulus},sig[t],t];
step2 = step /. tau->2 tau;
values = FullSimplify[{step,step2} /. t->tau,Assumptions->tau>0];
atStart = {step,step2} /. t->0;
atEnd = Limit[{step,step2},t->Infinity,Assumptions->tau>0];
Print["Example stress histories at onset, late time, t=tau: ",InputForm[{atStart,atEnd,values}]];
If[remainder=!=0 || atStart=!={modulus,modulus} || atEnd=!={0,0} ||
   !TrueQ[FullSimplify[values[[2]]>values[[1]],Assumptions->modulus>0]],Exit[1]];
