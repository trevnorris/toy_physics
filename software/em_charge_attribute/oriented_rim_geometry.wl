(* Independent check via the surface swept out by the unit normal in S^3. *)
ClearAll["Global`*"];
check[name_, actual_, expected_] := Module[{value = FullSimplify[actual]},
  Print[name, " = ", InputForm[value]];
  If[!TrueQ[FullSimplify[value == expected]], Exit[1]]
];

normal = {Sin[t] Sin[b] Cos[f], Sin[t] Sin[b] Sin[f], Sin[t] Cos[b], Cos[t]};
tangents = (D[normal, #] & /@ {t, b, f}) /. f -> 0;
(* Rotation about the third horizontal axis removes the azimuth dependence. *)
gram = FullSimplify[tangents . Transpose[tangents]];
area = FullSimplify[Sqrt[Det[gram]], 0 < b < Pi && 0 < t < Pi/2];
check["normal-space area element", area, Sin[t]^2 Sin[b]];
cap[angle_] := 2 Pi Integrate[area, {b, 0, Pi}, {t, 0, angle}];
(* Define positive signed turning for the rising upper face. Reflection of
   the graph reverses its three horizontal normal components at fixed +w. *)
reflection = Det[-IdentityMatrix[3]];
check["upper and lower vertical lips", {cap[Pi/2], reflection cap[Pi/2]}, {Pi^2, -Pi^2}];
check["45 degree lip as fraction", cap[Pi/4]/cap[Pi/2], 1/2 - 1/Pi];
check["cap traversed from flat back to flat", cap[angle] - cap[angle], 0];

$Assumptions = a > 0 && c > 0;
(* With an independently adjustable height profile, the strength can vary. *)
profile = c ArcSin[a/r]/a;
check["far coefficients of mirror profiles", Limit[r {profile, -profile}, r -> Infinity], {c, -c}];
check["vertical inclinations of mirror profiles",
  Limit[ArcTan /@ {D[profile, r], -D[profile, r]}, r -> a, Direction -> "FromAbove"], {-Pi/2, Pi/2}];

(* Parametrize by height w instead: radius R(w), slope v=dR/dw.
   The area density has no explicit w; at the lip R=a, v=0. *)
lagrangian = rad^2 Sqrt[1 + v^2];
beltrami = FullSimplify[lagrangian - v D[lagrangian, v]];
check["area first integral at the waist", beltrami /. {rad -> a, v -> 0}, a^2];
mag = Sqrt[(rad/a)^4 - 1];
check["conserved area first integral off the waist",
  FullSimplify[beltrami /. v -> mag, rad > a > 0], a^2];
far = Limit[rad^2/mag, rad -> Infinity];
check["upper and lower far bend coefficients", {far, -far}, {a^2, -a^2}];
check["factor for twice the mouth radius", (far /. a -> 2 a)/far, 4];
