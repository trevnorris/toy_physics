(* Outward momentum flux through a sphere about one target in the steady,
   constant-density inviscid potential-flow reference. k = Q_source/(4 Pi).
   vt is the signed outward target radial speed. Remote sink is at +d. *)
ClearAll["Global`*"];
cs[l_] := (d^2+a^2-l^2)/(2 a d);
bn = k (d cs[l]-a)/l^3;
bx = k (d-a cs[l])/l^3;
uniform = k/d^2;
(* The area element after substituting distance to the remote source. *)
sphere[term_] := Module[{primitive = Integrate[Expand[2 Pi a l term/d],l]},
  Factor[(primitive /. l -> d+a) - (primitive /. l -> d-a)]];
pairPressure = sphere[-rho vt bn cs[l]];
pairCarriage = sphere[rho vt (bx+bn cs[l])];
uniformPressure = Integrate[-2 Pi a^2 rho vt uniform t^2,{t,-1,1}];
uniformCarriage = Integrate[2 Pi a^2 rho vt uniform (1+t^2),{t,-1,1}];
Print["outward cross pressure, pair and stream = ",ToString[{pairPressure,uniformPressure},InputForm]];
Print["outward cross carriage, pair and stream = ",ToString[{pairCarriage,uniformCarriage},InputForm]];
If[Together[pairPressure-uniformPressure] =!= 0 ||
   Together[pairCarriage-uniformCarriage] =!= 0, Exit[1]];

(* Pressure alone from the undisturbed companion has a higher-order gradient.
   Evaluate its leading moment by expanding the angular integrand at a/d=0. *)
pressureDensity = -Pi a^2 rho k^2 t/(d^4 (1-2 h t+h^2)^2);
leadingPressure = FullSimplify[(a/d) Integrate[
  SeriesCoefficient[pressureDensity,{h,0,1}],{t,-1,1}]];
(* Total source-only momentum flux vanishes identically by steady Euler
   balance in the source-free ball. Check divergence directly in Cartesian coordinates. *)
vars = {x,y,z};
vel = D[k/Sqrt[(x-d)^2+y^2+z^2],{vars}];
stressRow = rho vel[[1]] vel - rho (vel.vel) {1,0,0}/2;
divRow = FullSimplify[Sum[D[stressRow[[i]],vars[[i]]],{i,3}],
  Assumptions -> Element[{x,y,z,d},Reals] && (x-d)^2+y^2+z^2>0];
Print["leading source-only outward pressure and total-flux divergence = ",
  ToString[{leadingPressure, divRow},InputForm]];
If[Together[leadingPressure + 8 Pi rho a^3 k^2/(3 d^5)] =!= 0 || divRow =!= 0, Exit[1]];
