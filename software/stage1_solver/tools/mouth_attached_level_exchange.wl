(* Phase-referenced parcel crossing of a hard-core localized bright receiver.
   Basis e0,e1,g0,g1. Both arms share the free parcel and receiver evolution.
   The phase of the determinant is followed from zero dwell, modulo pi for
   half the argument. Energies and times use hbar=1. *)
ClearAll["Global`*"];
emit[key_, value_] := WriteString["stdout", key, " = ", ToString[value, InputForm], "\n"];
residuals = {};
extendedOnly = Environment["MOUTH_ATTACHMENT_EXTENDED_RAMP"] === "1";
record[value_] := AppendTo[residuals, Max[Abs[Flatten[{value}]]]];
id2 = IdentityMatrix[2];
bright = {1,4};
gap = 6;
state[n_, c_:0] := {{1-n,c},{Conjugate[c],n}};
amplitude[u_, rho_] := Table[
  Sum[u[[2 (a-1)+i,2 (b-1)+j]] rho[[j,i]],{i,2},{j,2}], {a,2},{b,2}];
common[base_, probed_] := Arg[Det[ConjugateTranspose[base].probed]]/2;
population[q_, rho_] := Re[Tr[q.rho]];
occupied = DiagonalMatrix[{0,1,0,1}];
pExcited = DiagonalMatrix[{1,1,0,0}];
vac = state[0];
coupling = DiagonalMatrix[{0,0,0,0}]; coupling[[1,4]]=coupling[[4,1]]=1;
hfree[d_] := gap pExcited + (gap+d) occupied;
arm[g_,d_,t_] := MatrixExp[I t hfree[d]].MatrixExp[-I t (hfree[d]+g coupling)];
survival[g_,d_,t_] := With[{w=Sqrt[g^2+d^2/4]},
  Exp[-I d t/2] (Cos[w t]+I d Sin[w t]/(2 w))];
response[g_,d_,t_,n_,c_] := Module[{u,v,m,a,p,detFormula,rho0,rho1,
     input=0.73, dn,energy,width,op,mirror,angles,tilted,mixture},
  u=N[arm[g,d,t],25]; v=amplitude[u,vac]; m=amplitude[u,state[n,c]];
  a=survival[g,d,t]; p=1-Abs[a]^2;
  detFormula=(n+(1-n) a) ((1-n)+n Conjugate[a])+p Abs[c]^2;
  record[{u[[1,1]]-a,u[[4,4]]-Conjugate[a],Det[m]-detFormula}];
  (* Distinct diagonal probe for energetics; tracing an initially coherent
     receiver against this input also gives the population formula below. *)
  rho0=KroneckerProduct[DiagonalMatrix[{input,1-input}],state[n,c]];
  rho1=u.rho0.ConjugateTranspose[u];
  dn=population[occupied,rho1-rho0]; energy=population[hfree[d],rho1-rho0];
  record[{dn-p (input-n),energy-d dn,population[pExcited,rho1-rho0]+dn}];
  mirror=KroneckerProduct[PauliMatrix[1],id2];
  tilted=amplitude[mirror.u.ConjugateTranspose[mirror],state[n,c]];
  mixture=amplitude[mirror.u.ConjugateTranspose[mirror],vac];
  record[common[v,m]-common[mixture,tilted]];
  angles=Arg[Diagonal[ConjugateTranspose[v].m]];
  {N[common[v,m],14], N[Min[SingularValueList[m]],12],N[p,12],N[dn,12],
   N[energy,12], N[angles,13]}];
archived=522640709968/10^13;
If[!extendedOnly,
emit["abrupt columns","{half determinant phase, worst amplitude, swap probability, receiver population change for x=.73, total free energy change, spin phases}"];
Do[
 emit["rectangular {G,Delta,T,n} "<>ToString[{g,d,t,n},InputForm],
      response[N[g,30],N[d,30],N[t,30],N[n,30],0]],
 {g,{Sqrt[2]}},{d,{-2,0,2,4}},{t,{33/200,33/100}},{n,{0,archived,1/2,1}}];
emit["retained coherence d=2,t=.33,n=.25,c=.2i",
     response[N[Sqrt[2],30],2,33/100,1/4,I/5]];
(* Undressed lead spin and receiver have gaps 6 and 6+d. The first
   nonzero common phase is cubic for a sudden short passage. *)
series=Assuming[{Element[{g,d,n,t},Reals],g>0,d>0},
  FullSimplify[Series[survival[g,d,t],{t,0,3}]]];
emit["series for amplitude of e0 in rectangular zone",series]];

(* A passive spatial coupling envelope g(z)=gmax sin^2(pi z/L)
   sampled along a positive constant-speed trajectory; the free momentum
   supplies level mismatch in nonadiabatic passages. Use exact propagation
   over each short constant cell and check phase refinement. *)
level[g_,d_,k_:1] := Sign[d] (Sqrt[d^2+4 k g^2]-Abs[d])/2;
relativeEnd[g0_,d_,duration_,steps_] := Module[{dt=N[duration/steps],freq,
    env,ham,trans,phase,id=N[id2], result=N[id2]},
  Do[
    env=g0 Sin[Pi (i-.5)/steps]^2;
    freq=Sqrt[env^2+d^2/4];
    ham={{-d/2,env},{env,d/2}};
    trans=Cos[freq dt] id-I Sin[freq dt] ham/freq;
    result=trans.result,
    {i,steps}];
  DiagonalMatrix[Exp[I {-d/2,d/2} duration]].result];
adiabaticSample[g0_,d_,duration_,n_] := Module[{fine,coarse,adiabaticAngle,
   a,phaseError,reference,pred,exact},
  fine=relativeEnd[g0,d,duration,4096];
  coarse=relativeEnd[g0,d,duration,2048];
  adiabaticAngle=duration NIntegrate[level[g0 Sin[Pi x]^2,d],{x,0,1}];
  a=fine[[1,1]]; phaseError=Arg[a Exp[-I adiabaticAngle]];
  pred=(1-n)+n Exp[-I adiabaticAngle];
  reference=Arg[a];
  exact=Arg[(n+(1-n) a) ((1-n)+n Conjugate[a]) Conjugate[a]]/2;
  {duration,Abs[fine[[2,1]]]^2,phaseError,Abs[a-coarse[[1,1]]],
    exact,Arg[pred],Abs[pred]}];
emit["smooth d=2 g=sqrt(2), n=.25; columns",
    "{T, swap probability, e0 phase error vs adiabatic radians, 2048/4096 amplitude difference, measured half phase, ideal phase, ideal fringe amplitude}"];
smooth=Table[adiabaticSample[N[Sqrt[2]],2.,tt,1/4],
  {tt,If[extendedOnly,{128},{8,16,32,64,128}]}];
Scan[emit["smooth",#]&,smooth];
(* Independent evolution of the last smooth pulse checks whether the
   asymptotic threshold is reached. The rectangular segment product uses
   piecewise-constant midpoint coupling and this ODE uses the exact envelope. *)
duration=128;
ode=NDSolveValue[{I x'[t]==Sqrt[2] Sin[Pi t/duration]^2 y[t],
     I y'[t]==Sqrt[2] Sin[Pi t/duration]^2 x[t]+2 y[t],
     x[0]==1,y[0]==0},{x,y},{t,0,duration},WorkingPrecision->28,
     AccuracyGoal->16,PrecisionGoal->16,MaxSteps->Infinity];
transfer=Abs[ode[[2]][duration]]^2;
rkError=Abs[ode[[1]][duration]-relativeEnd[N[Sqrt[2]],2.,duration,4096][[1,1]]];
emit["extended T=128 continuous-envelope ODE {swap probability, amplitude difference from 4096 cells}",
     {transfer,rkError}];
(* A harmonic receiver has a second rung |e,1> <-> |g,2> of strength
   sqrt(2)G. Test with independently sorted block eigenvalues; compare
   the level centroid shift for known n=1 vs 0 on the two completions. *)
ladderTest[g_,d_] := Module[{branches,asIfHard,soft,s=Sign[d],b1,b2},
 branches[k_]:=Sort[Eigenvalues[N[{{0,Sqrt[k] g},{Sqrt[k] g,d}},25]]];
 b1=If[d>0,-First[branches[1]],-Last[branches[1]]];
 b2=If[d>0,-First[branches[2]],-Last[branches[2]]];
 asIfHard={b1,b1}; soft={b1-b2,b1};
 record[{Mean[asIfHard]-level[g,d],Mean[soft]-(2 level[g,d]-level[g,d,2])/2}];
 {N[asIfHard,12],N[soft,12],N[Mean[soft],12],N[Mean[soft]/b1,12]}];
If[!extendedOnly,
emit["centroid comparison columns","{hard-core occupied-minus-empty energy shifts on e,g; oscillator shifts on e,g; oscillator mean; mean/hard-core}"];
Do[emit["spectrum G=1 Delta="<>ToString[d],ladderTest[1,d]],{d,{-12,-4,4,12,40}}];
large=Assuming[{d>0,g>0},FullSimplify[Series[(2 level[g,d]-level[g,d,2])/2,{g,0,4}]]];
emit["oscillator common energy shift series at positive detuning",large];
emit["max independent matrix/energy/phase/spectrum residual",Max[residuals]]];
Exit[If[(extendedOnly || (Max[residuals]<10^-12 &&
    Abs[Last[smooth][[3]]]<Abs[First[smooth][[3]]])) &&
  Max[Last[smooth][[2]],transfer]<10^-6 && rkError<10^-5 &&
  Max[smooth[[All,4]]]<10^-5,0,1]];
