(* Recover the field from intact-face pressures and the spherical Green function. *)
ClearAll["Global`*"];
$HistoryLength = 0;
$Assumptions = r >= a > 0 && k > 0 && stiff > 0 && rho > 0 && rhof > 0 &&
  rate > 0 && Element[{aa, bb, eps, transfer}, Reals];
check[label_, actual_, expected_] := Module[{value = FullSimplify[actual]},
  Print[label <> " = " <> ToString[value, InputForm]];
  If[FullSimplify[value == expected] =!= True, Exit[1]]];

(* Feed pressure on the two reflecting faces, oriented toward increasing z.
   The receiver is on the upper face. Uniform ambient stresses are subtracted. *)
feed = (aa + bb Exp[-k r] Cos[k z])/r;
pf = -rhof (D[feed,r]^2 + D[feed,z]^2)/2;
feedShell = FullSimplify[4 Pi r^2 ((pf /. z->Pi/k) - (pf /. z->0))];
bulk = -rate/(2 Pi^2 rho (r^2+w^2)) -
  eps rate a^2 (4 w^2/(r^2+w^2)-1)/(4 Pi^2 rho (r^2+w^2)^2);
pb = -rho (D[bulk,r]^2 + D[bulk,w]^2)/2;
bulkShell = FullSimplify[-4 Pi r^2 (pb /. w->0)];
outside[shell_] := Integrate[shell /. r->u, {u,r,Infinity},
  Assumptions->$Assumptions];
feedReturn = outside[feedShell];
bulkReturn = outside[bulkShell];
check["feed pressure load outside r", feedReturn, 8 Pi rhof aa bb Exp[-k r]/r];
check["bulk pressure load outside r", bulkReturn,
  2 rate^2/(Pi^3 rho) (1/(3 r^3)-2 eps a^2/(5 r^5)+eps^2 a^4/(7 r^7))];

(* Compact recoil balances the integrated distributed pressure for each
   exterior with vanishing outer carrier flux and no separate recipient.
   For r>=a, interior shells contribute 1/r, exterior shells 1/u. *)
return = feedReturn + bulkReturn;
compact = -(return /. r->a);
enclosed = compact + (return /. r->a) - return;
bulkHeight = FullSimplify[(-bulkReturn/r + outside[bulkShell/r])/(4 Pi stiff)];
check["bulk midpoint displacement", bulkHeight, -rate^2/(2 Pi^4 rho stiff) *
  (1/(12 r^4)-eps a^2/(15 r^6)+eps^2 a^4/(56 r^8))];
(* The feed Green integral is -2 rhof aa bb/stiff times Integrate[Exp[-k u]/u^3].
   Its derivative and the zero at infinity fix the integration constant. *)
feedHeight = -2 rhof aa bb ExpIntegralE[3,k r]/(stiff r^2);
recur = ExpIntegralE[3,k r] -> (Exp[-k r]-k r ExpIntegralE[2,k r])/2;
check["combined radial force balance",
  Expand[-4 Pi stiff r^2 D[feedHeight+bulkHeight,r] - enclosed] /. recur, 0];
check["height at infinity", Limit[feedHeight+bulkHeight,r->Infinity], 0];
check["1/r coefficient of total displacement", Limit[r (feedHeight+bulkHeight),r->Infinity], 0];
check["r^4 times total displacement at infinity",
  Limit[r^4 (feedHeight+bulkHeight),r->Infinity], -rate^2/(24 Pi^4 rho stiff)];
check["1/r coefficient with separate transfer",
  Limit[r (feedHeight+bulkHeight-transfer/(4 Pi stiff r)),r->Infinity],
  -transfer/(4 Pi stiff)];
Exit[0];
