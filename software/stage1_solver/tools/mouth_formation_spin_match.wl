(* Spin deposited per solid angle, directions in one ambient frame.
   signMap=+1 matches the same radial direction and -1 the antipode.
   This is endpoint angular-impulse accounting for the candidate textures. *)
ClearAll["Global`*"];
emit[label_,data_] := Print[label <> " = " <> ToString[data,InputForm]];
radial={Sin[theta]Cos[phi],Sin[theta]Sin[phi],Cos[theta]};
solid[q_] := Integrate[q Sin[theta], {theta,0,Pi}, {phi,0,2Pi}];
charge[s_,chi_] := -s chi;
spin[s_,chi_,magnitude_,direction_] := magnitude s chi direction;
first=spin[s,chi1,l1,radial];
second=spin[-s,chi2,l2,radial];
netSpin={solid[first],solid[second]};
emit["whole-sphere spin on each mirrored mouth",netSpin];

(* When both momenta are torques integrated over a birth, an equal and
   opposite impulse at every registered pair of points requires J1+J2=0.
   The squared mismatch keeps directional cancellations on the sphere
   from hiding residual impulses. *)
residual[signMap_]:= first+spin[-s,chi2,l2,signMap radial];
rawMismatch=Factor[solid[residual[signMap].residual[signMap]]/(4Pi)];
qPair={charge[s,chi1],charge[-s,chi2]};
emit["two charges,total charge,pointwise mismatch squared (sphere mean)",
     {qPair,Total[qPair],rawMismatch}];
(* Fix s=chi1=1 as a convention. Simultaneously flipping either
   orientation or both chi merely reverses charge signs. *)
rows=Flatten[Table[Module[{r={s->1,chi1->1,chi2->chi,signMap->p},norm},
  norm=rawMismatch /. r;
  {p,chi,qPair /. r,Factor[norm],
   Reduce[l1>=0 && l2>=0 && norm==0,{l1,l2},Reals]}],
  {p,{1,-1}},{chi,{1,-1}}],1];
emit["point registration,chi2,charges,mean square,nonnegative magnitudes with exact match",
     rows];
matchingCharge=Table[{p,Resolve[ForAll[{l1,l2},l1>0 && l2>0,
  And@@Flatten[Table[
    Implies[(rawMismatch/.{signMap->p,s->si,chi1->a,chi2->b})==0,
            b==p a && charge[-si,b]==-p charge[si,a] && l1==l2],
    {si,{-1,1}},{a,{-1,1}},{b,{-1,1}}]]],Reals]}, {p,{-1,1}}];
emit["positive-magnitude local-match implication for both orientations",matchingCharge];

(* Independent births: list sign assignments, with no probability law.
   Two completed same-direction-matched births labeled A,B each have
   two possible signs on the member with orientation s=+1. *)
birthPair[chi_,p_]:={charge[1,chi],charge[-1,p chi]};
independent=Flatten[Table[
  {a,b,birthPair[a,1],birthPair[b,1],
   Times@@{First[birthPair[a,1]],First[birthPair[b,1]]}},
  {a,{-1,1}},{b,{-1,1}}],1];
compensated={birthPair[1,-1],birthPair[-1,-1]};
emit["independent matched births: chiA,chiB,pairA,pairB,same-orientation product",
     independent];
emit["two antipodally matched births with zero four-mouth net charge",
     {compensated,Total[Flatten[compensated]]}];
checks={netSpin==ConstantArray[0,{2,3}],And@@matchingCharge[[All,2]],
        Sort[DeleteDuplicates[independent[[All,5]]]]=={-1,1},
        Total[Flatten[compensated]]==0};
emit["checks",checks];
If[!TrueQ[And@@checks],Exit[1]];
