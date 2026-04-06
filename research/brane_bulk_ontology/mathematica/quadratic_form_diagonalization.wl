(* ---------------------------------------------------------------------- *)
(* Script: quadratic_form_diagonalization.wl                              *)
(* Minimal algebraic check of the positive-definite two-mode toy model.   *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

assumePositive = {AT > 0, AL > 0};

alpha2 = 3/4;
alpha = Sqrt[alpha2];
uT = Symbol["uT"];
uVector = {uT, alpha*uT};
matM = DiagonalMatrix[{AT, AL}];

energy = Expand[(1/2) * (uVector . matM . uVector)];
energyTarget = (1/8) (4*AT + 3*AL) uT^2;
evals = Eigenvalues[matM];
positiveCheck = Assuming[assumePositive, Simplify[energyTarget/uT^2 > 0]];

Print["--- Two-Mode Input ---"];
Print["alpha^2 = ", alpha2];
Print["u_L = alpha u_T = ", alpha // TraditionalForm, " u_T"];
Print[""];

Print["--- Effective Single-Mode Energy ---"];
Print["Derived E(uT) ="];
Print[energy // TraditionalForm];
Print["Target calibrated form ="];
Print[energyTarget // TraditionalForm];
Print["Difference (derived - target) ="];
Print[Simplify[energy - energyTarget] // TraditionalForm];
Print[""];

Print["--- Eigenvalues of the Underlying Quadratic Form ---"];
Print[evals // TraditionalForm];
Print["Positive for AT > 0 and AL > 0:"];
Print[positiveCheck];

(* 
Output:
--- Two-Mode Input ---
alpha^2 = 3/4
u_L = alpha u_T = TraditionalForm[Sqrt[3]/2] u_T

--- Effective Single-Mode Energy ---
Derived E(uT) =
TraditionalForm[(3*AL*uT^2)/8 + (AT*uT^2)/2]
Target calibrated form =
TraditionalForm[((3*AL + 4*AT)*uT^2)/8]
Difference (derived - target) =
TraditionalForm[0]

--- Eigenvalues of the Underlying Quadratic Form ---
TraditionalForm[{AL, AT}]
Positive for AT > 0 and AL > 0:
True
*)
