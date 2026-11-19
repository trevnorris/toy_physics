(* Use the signed areas of radial sections and the continuity defect to test a
   gauge change. Work in units of the charge assigned to the upper mouth. *)
ClearAll["Global`*"];
loss[r_] := 30 (r - 2)^2 (3 - r)^2;
ordered[r_] = Integrate[loss[u], {u, r, 3}];
dir[c_] := (r c - b)/Sqrt[r^2 + b^2 - 2 r b c];
degree = (dir[1] - dir[-1])/2;
sections = {FullSimplify[degree, 0 < r < b], FullSimplify[degree, r > b > 0]};
Print["degrees inside and outside the moving point = ", InputForm[sections]];
index = Last[sections] - First[sections];
gauge[b_] := (b-2) (3-b);
cases = {{"layer alone", 1-ordered[r], 0}, {"completed", 1-ordered[r], 1},
         {"localized primitive", ordered[r] (1-ordered[r]), 1}};
Do[
  {label, p, withPoint} = item;
  regular = First[sections] Integrate[D[p,r], {r,2,b}] +
            Last[sections] Integrate[D[p,r], {r,b,3}];
  point = withPoint index (p /. r -> b);
  defect = D[regular + point, b];
  pairing = -Integrate[gauge[b] defect, {b,2,3}];
  halves = Factor[{regular, point} /. b -> 5/2];
  edges = Factor[{regular+point /. b -> 2, regular+point /. b -> 3}];
  Print[label, " gauge pairing, midpoint components, limiting charges = ",
        InputForm[{pairing, halves, edges}]];
  If[!TrueQ[Simplify[pairing == Integrate[gauge'[b] (regular+point), {b,2,3}]]] ||
     (withPoint == 1 && !TrueQ[Simplify[defect == 0]]), Exit[1]],
  {item, cases}];
