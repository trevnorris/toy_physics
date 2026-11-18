(* Vary the harmonic wall action; then use stationarity of the normalized
   displacement to differentiate the quadratic boundary response. *)
ClearAll["Global`*"];
check[name_, left_, right_] := Module[{res = FullSimplify[left - right, $Assumptions]},
  Print[name, " | ", ToString[left, InputForm], " | ", ToString[right, InputForm],
    " | residual ", ToString[res, InputForm]];
  If[res =!= 0, Exit[1]]
];
angular = Integrate[Cos[th]^2 Sin[th], {th, 0, Pi}, {az, 0, 2 Pi}];
check["angular norm", angular, 4 Pi/3];
lag = angular (mu[w] v^2 - ten[w] d^2 - k1[w] x^2)/2;
(* Lower endpoint contributes minus the gradient derivative to delta S. *)
fieldOnRim = -D[lag, d] /. d -> q'[0] /. w -> 0;
check["field force on rim", fieldOnRim, angular ten[0] q'[0]];

pot = angular (ten[w] D[p[w, lam], w]^2 + (k1[w] + lam mu[w]) p[w, lam]^2);
source = -D[ten[w] D[p[w, lam], w], w] + (k1[w] + lam mu[w]) p[w, lam];
(* Derivative of the quadratic form before eliminating the interior and endpoint terms. *)
deriv = Expand[D[pot, lam] - 2 angular D[p[w, lam], lam] source -
  D[2 angular ten[w] D[p[w, lam], w] D[p[w, lam], lam], w]];
check["stationary response derivative density", deriv, angular mu[w] p[w, lam]^2];
forced[y_] := -D[ten[w] D[y, w], w] + (k1[w] + lam mu[w]) y;
workA = Expand[D[ten[w] q[w] p0'[w], w] + q[w] forced[p0[w]]];
workB = Expand[D[ten[w] p0[w] q'[w], w] + p0[w] forced[q[w]]];
check["interchange of source and probe in weak work", workA, workB];

(* Solve two remote/mouth trace equations in a free-ended uniform control. *)
material = len > 0 && stiff > 0 && dens > 0;
$Assumptions = s > 0 && material;
wave = s Sqrt[dens/stiff];
field[w_] := h1 Exp[wave w] + h2 Exp[-wave w];
coeffs = First@Solve[{field[0] == 1, field'[len] == 0}, {h1, h2}];
transfer = FullSimplify[-angular stiff field'[0] /. coeffs, $Assumptions];
check["remote force", FullSimplify[stiff field'[len] /. coeffs, $Assumptions], 0];
check["drive force per displacement", transfer, angular s Sqrt[stiff dens] Tanh[wave len]];
check["slow force per acceleration", Limit[transfer/s^2, s -> 0, Assumptions -> material],
  angular dens len];
check["fast force per acceleration", Limit[transfer/s^2, s -> Infinity, Assumptions -> material], 0];
