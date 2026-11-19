(* Weak held-patch response of a jointly replaced two-face spin pair.
   Units set hbar=1. The input is a field f coupled to S_top,x-S_bottom,x.
   A fresh pair occupies the patch at every event. With probability keep
   it receives the old joint state by SWAP, then acquires an energy phase.
   Otherwise it enters in the instantaneous ground state. *)
ClearAll["Global`*"];
emit[name_, value_] := Print[name, " = ", ToString[value, InputForm]];
spin = Table[PauliMatrix[a]/2,{a,3}];
top = KroneckerProduct[#,IdentityMatrix[2]]& /@ spin;
bottom = KroneckerProduct[IdentityMatrix[2],#]& /@ spin;
staggered = top-bottom;
exchange = 6; alignment = 4; rate = 4; keep = 1/2;
ham = exchange Sum[top[[a]].bottom[[a]],{a,3}]-alignment staggered[[3]];
probe = staggered[[1]];
energies = Eigenvalues[ham];
groundEnergy = Min[energies];
ground = Normalize[First[NullSpace[ham-groundEnergy IdentityMatrix[4]]]];
groundProjector = Outer[Times,ground,Conjugate[ground]];
(* The perturbation tilts the pinning field about y through f/alignment. *)
comm[a_,b_] := a.b-b.a;
groundSlope = -I comm[top[[2]]+bottom[[2]],groundProjector]/alignment;
susceptibility0 = Tr[probe.groundSlope];
emit["J, h0, replacement rate, state handoff probability",
     {exchange,alignment,rate,keep}];
emit["ground energy, staggered polarization, transverse static susceptibility",
     {groundEnergy,Tr[staggered[[3]].groundProjector],susceptibility0}];

phaseUnitary[angle_] := groundProjector+Exp[-I angle](IdentityMatrix[4]-groundProjector);
phaseSlope[angle_] := (1-Exp[-I angle]) groundSlope;
handoff[state_,angle_] := phaseUnitary[angle].state.ConjugateTranspose[phaseUnitary[angle]];
handoffSlope[angle_] := phaseSlope[angle].groundProjector.ConjugateTranspose[phaseUnitary[angle]] +
     phaseUnitary[angle].groundProjector.ConjugateTranspose[phaseSlope[angle]];
event[state_,angles_] := keep Mean[handoff[state,#]& /@ angles]+(1-keep)Tr[state] groundProjector;
eventSlope[angles_] := keep Mean[handoffSlope /@ angles]+(1-keep) groundSlope;
liouville[state_,angles_] := -I comm[ham,state]+rate (event[state,angles]-state);
basis = Flatten[Table[Normal[SparseArray[{{a,b}->1},{4,4}]],{a,4},{b,4}],1];
matrix[op_] := Transpose[Flatten[op[#]]& /@ basis];
susceptibility[angles_] := Module[{drive, response, generator},
     generator = matrix[Function[state,liouville[state,angles]]];
     drive = I comm[probe,groundProjector]+rate eventSlope[angles];
     response = LinearSolve[s IdentityMatrix[16]-generator, Flatten[drive]];
     Factor[Tr[probe.ArrayReshape[response,{4,4}]]/susceptibility0]];
(* Two meters read the real orientation and its conjugate transition
   quadrature. Normalize both against the same small pre-event tilt. *)
gap = (probe.ground).((ham-groundEnergy IdentityMatrix[4]).(probe.ground)) /
       ((probe.ground).(probe.ground));
quadrature = -I comm[ham,probe]/gap;
anglesByCase = {{Pi/2},{Pi/2,-Pi/2}};
readings = Table[With[{after=event[groundSlope,angles]},
   {Tr[probe.after],Tr[quadrature.after]}/susceptibility0],{angles,anglesByCase}];
responseFunctions = susceptibility /@ anglesByCase;
emit["energy gap and single-event (orientation, quadrature) for fixed and random phase",
     {gap,readings}];
emit["held susceptibilities divided by the static value, s=-i omega",responseFunctions];
lagTimes = (-D[#,s]/.s->0)& /@ responseFunctions;
emit["low-frequency phase slope for fixed and random phase",lagTimes];
emit["complex response at omega=4, divided by static",responseFunctions /. s->-4 I];

(* Fixed-field energy never grows in a jump for any incoming state.
   This checks the full observable identity on all sixteen matrix units. *)
energyResiduals = Table[Tr[ham.(event[b,angles]-b)] -
     (1-keep) (groundEnergy Tr[b]-Tr[ham.b]),{angles,anglesByCase},{b,basis}];
zeroMatrix = ConstantArray[0,{4,4}];
stationarySlope = Table[liouville[groundSlope,angles]+I comm[probe,groundProjector]+
     rate eventSlope[angles],{angles,anglesByCase}];
absorption = FullSimplify[ComplexExpand[Im[# /. s->-I w]],
             Assumptions->w>0]& /@ responseFunctions;
emit["imaginary part for positive omega",absorption];
checks = {And@@(#===zeroMatrix& /@ stationarySlope), Flatten[energyResiduals]===ConstantArray[0,32],
    And@@(TrueQ[FullSimplify[#>0,Assumptions->w>0]]& /@ absorption),
    (responseFunctions/.s->0)==={1,1}, First[readings][[1]]===Last[readings][[1]],
    First[lagTimes]=!=Last[lagTimes]};
emit["stationarity, energy, absorbed work, static, equal directional transfer, distinct lag",checks];
If[!TrueQ[And@@checks],Exit[1]];
