(* Differentiate the gradient energy with respect to a common virtual tilt.
   A tilted order level can also change width and can slip relative to the
   half-order crossing. Use disjoint planar face integrals (wide-slab limit). *)
ClearAll["Global`*"];
$HistoryLength = 0;
$Assumptions = delta > 0 && kap > 0 && Element[{s, h, g, probe}, Reals];
check[label_, actual_, expected_] := Module[{value = Factor[actual]},
  Print[label <> " = " <> ToString[value, InputForm]];
  If[FullSimplify[value == expected] =!= True, Exit[1]]];

(* x=tanh(w/(2 delta)) at one face. Translating its level sets with
   slope v(x) gives radial derivative -v(x) dchi/dw at linear order. *)
grad = (1-x^2)/(4 delta);
jac = 2 delta/(1-x^2);
energy[v_] := Integrate[Cancel[kap grad^2 v^2 jac/2], {x,-1,1}];
lowerSlope = s + probe - h/2 + g*x^2;
upperSlope = s + probe + h/2 + g*x^2;
e = energy[lowerSlope]+energy[upperSlope];
readout = D[e,probe] /. probe->0;
check["translation stiffness from area energy", D[readout,s], kap/(3 delta)];
check["width gradient contribution to common traction", D[readout,h], 0];
check["mean half-order height slope", (lowerSlope+upperSlope)/2 /. {x->0,probe->0}, s];
check["shape / translation coefficient from virtual work",
  D[readout,g]/D[readout,s], 1/5];
check["two-face traction from virtual work", readout, kap (s+g/5)/(3 delta)];
Exit[0];
