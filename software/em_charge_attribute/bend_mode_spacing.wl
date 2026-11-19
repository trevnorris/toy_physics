(* Work with measured squared frequencies; eliminate the common gap and
   independently compute the pressure of a gapped mode from its field current. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[label_, value_] := Print[label, " = ", ToString[value,InputForm]];
check[label_, a_, b_] := Module[{res=Factor[a-b]},
  emit[label,{a,b,res}]; If[res=!=0,Exit[1]]];
inv = First[Solve[{w0sq == gapSq + c^2 kSq,
                    w1sq == gapSq + 9 c^2 kSq},{gapSq,kSq}]];
emit["{gap squared,length squared} from measured frequencies",
     Factor[{gapSq,Pi^2/(4 kSq)}/.inv]];
(* Average over the fast phase first. ux is the displacement's axial
   derivative amplitude and ut is its time derivative amplitude. *)
shape = Sin[Pi x/(2 ell)];
ut = omega shape;
ux = D[shape,x];
dens = rho (ut^2+c^2 ux^2+gap^2 shape^2)/4;
current = rho (ut^2+c^2 ux^2-gap^2 shape^2)/4;
disp = omega^2 -> gap^2 + c^2 Pi^2/(4 ell^2);
e = Integrate[dens,{x,0,ell},Assumptions->ell>0] /.disp;
ends = Factor[({-(current/.x->0),current/.x->ell}/e) /. disp];
emit["endpoint loads per energy at fixed material",ends];
check["internal reaction", Total[ends],0];
check["fraction of gapless end pressure",Factor[ell Last[ends]],
  (c Pi/(2 ell))^2/(gap^2+(c Pi/(2 ell))^2)];
Exit[0];
