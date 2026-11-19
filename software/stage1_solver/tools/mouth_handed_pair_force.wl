(* Two separate controls: neutral tangent-plane response, then fixed
   discrete states of two well-separated mouths at prescribed centers.
   The material scalar is real and passive (kappa>0). Remove only the
   independent self energies, identical in the preparations compared. *)
ClearAll["Global`*"];
emit[label_, obj_] := Print[label <> " = " <> ToString[obj,InputForm]];

generator = {{0,-1},{1,0}};
reflection = DiagonalMatrix[{1,-1}];
response={{aa,cc-I b},{cc+I b,ee}};
commutator=response.generator-generator.response;
isotropicRules=First[Solve[Thread[Flatten[commutator]==0], {cc,ee}]];
isotropic=response /. isotropicRules;
circular[h_] := {1,I h}/Sqrt[2];
gain[mat_,h_] := Expand[Conjugate[circular[h]].mat.circular[h]];
gains[mat_] := gain[mat,#]& /@ {1,-1};
mirror=reflection.isotropic.Transpose[reflection];
invariant=Solve[Thread[Flatten[mirror-isotropic]==0],b];
example=First[Solve[Thread[gains[isotropic]=={9/5,1/10}],{aa,b}]];
emit["tangent rotation commutator, proper-isotropic response",
     {commutator,isotropic}];
emit["circular gain (+,-), reflected gains, reflection-invariant condition",
     {gains[isotropic],gains[mirror],invariant}];
emit["example mean,handed coefficient,split,normalized contrast,ratio",
     {aa,b,Subtract@@gains[isotropic],b/aa,Divide@@gains[isotropic]} /. example];

(* (-Laplacian+mu^2) G = delta in three dimensions. Obtain the
   decaying solution from a Yukawa ansatz and fix its normalization by
   integrating the PDE over a small sphere. *)
kernel0=amp Exp[-mu d]/d;
ode=Factor[-D[d^2 D[kernel0,d],d]/d^2+mu^2 kernel0];
fluxAtOrigin=Limit[-4 Pi d^2 D[kernel0,d],d->0,Direction->"FromAbove"];
kernel=kernel0 /. First[Solve[fluxAtOrigin==1,amp]];
(* Stationarity makes the field energy equal to half the source pairing,
   so the minimized total energy is minus half the pairing. Both
   source orders contribute once to the following cross term. *)
pairPhi= {lambda chi2 kernel/kappa,lambda chi1 kernel/kappa};
scalar=Factor[-lambda {chi1,chi2}.pairPhi/2];
n1=-s1 chi1; n2=-s2 chi2;
energy=Factor[coul n1 n2/d+scalar];
outward=Factor[-D[energy,d]];
emit["radial Green equation residual,origin flux,kernel",{ode,fluxAtOrigin,kernel}];
emit["integrated-out scalar cross energy,combined cross energy,outward force",
     {scalar,energy,outward}];

(* For equal s, aligned chi: c=coul>0, x=mu d, screened amplitude
   v=lambda² exp(-x)/(4 Pi kappa)>=0. For antialignment both signs
   reverse, so 2U is the aligned-minus-anti energy. *)
alignRules={s1->1,s2->1,chi1->1,chi2->1};
strip[q_] := Expand[q] /. lambda^2 -> 4 Pi kappa v Exp[mu d];
scaledEnergy=Simplify[strip[d (energy /. alignRules)]];
scaledForce=Simplify[strip[d^2 (outward /. alignRules)]] /. mu->x/d // Simplify;
gap=2 scaledEnergy;
identity=Simplify[scaledEnergy-scaledForce];
finite=Reduce[coul>0 && v>=0 && x>0 && gap<0 && scaledForce>0,
              {coul,v,x},Reals];
massless=({gap,scaledForce} /. x->0);
infinite=Reduce[coul>0 && v>=0 && First[massless]<0 && Last[massless]>0,
                {coul,v},Reals];
emit["same orientation aligned: d*(U_aligned-U_anti), d² F, d U_aligned-d² F",
     {gap,scaledForce,identity}];
emit["finite screening strict lower-and-repulsive parameter set",finite];
emit["infinite range difference and force, strict simultaneous parameter set",
     {massless,infinite}];
checks={gains[mirror]==Reverse[gains[isotropic]],(b /. First[invariant])==0,
        ode==0,fluxAtOrigin==4 Pi amp,identity==v x,finite===False,infinite===False};
emit["checks",checks];
If[!TrueQ[And@@checks],Exit[1]];
