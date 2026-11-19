(* Split infinitesimal axial dilation work by energy channel, at fixed action. *)
ClearAll["Global`*"];
$Assumptions = r>0 && len>0 && act>0 && vr>0 && vz>0 && gap>=0;
check[name_, expr_, expected_] := Module[{v = FullSimplify[expr]},
  Print[name, " = ", ToString[v,InputForm]];
  If[!TrueQ[FullSimplify[v==expected]],Exit[1]]
];
freq[r_, l_] := Sqrt[gap^2+vr^2/r^2+vz^2/l^2];
ampEnergy = act freq[r,len];
(* The pressure includes the volume energy of depleted feed and the bulk
   pressure term in the cylindrical model. Cap tension is independent of len. *)
channels[x_] := {act freq[r,len Exp[x]], 4 Pi sig r^2 len Exp[x],
  8 Pi sig r^3/3, 4 Pi p r^3 len Exp[x]/3, 16 Pi bend len Exp[x]/9};
slopes = (D[channels[x],x]/len)/.x->0;
push = -First[slopes];
check["fraction of E/length carried as wave thrust",
  push len/ampEnergy, (vz/len/freq[r,len])^2];
balanced = thrust /. First[Solve[-thrust+Total[Rest[slopes]]==0,thrust]];
check["equilibrium axial work contributions", Rest[slopes],
  {4 Pi sig r^2,0,4 Pi p r^3/3,16 Pi bend/9}];
check["balanced thrust increment at doubled radius",
  (balanced/.r->2 r)-balanced,12 Pi sig r^2+28 Pi p r^3/3];
(* The optional curvature energy stays constant under purely radial
   dilation at fixed length. Test whether it alone selects a free radius. *)
curvatureOnly = act freq[r Exp[y],len]+16 Pi bend len/9;
check["outward radial drive in curvature-only example",
  (-D[curvatureOnly,y]/r)/.y->0,act vr^2/(r^3 freq[r,len])];
