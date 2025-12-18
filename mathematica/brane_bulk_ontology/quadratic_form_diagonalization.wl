(* ---------------------------------------------------------------------- *)
(* Script: QuadraticForm_Diagonalization.wl *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

(* 1. Parameters *)
(* We use AT and AL instead of A_T, A_L to avoid underscores *)
assumePositive = {AT > 0, AL > 0};

alpha2 = 3/4;       (* Rational number *)
alpha = Sqrt[alpha2]; (* This will automatically handle the Imaginary unit I *)

(* 2. Basis Vector *)
(* Python: u = [u_T, alpha * u_T] *)
(* Mathematica uses Lists {} to represent vectors *)
uT = Symbol["uT"];
uVector = {uT, alpha * uT};

(* 3. Euclidean kernel M *)
(* Python: M = sp.diag(A_T, A_L) *)
matM = DiagonalMatrix[{AT, AL}];

(* 4. Calculate Energy E *)
(* Python: E = (u.T * M * u)[0] *)
(* In Mathematica, the Dot product (.) automatically handles the *)
(* vector-matrix-vector multiplication: Row . Matrix . Column *)
energy = uVector . matM . uVector;

energySimplified = Simplify[energy];

Print["--- Energy E(uT) ---"];
Print["Unsimplified:"];
Print[energy // TraditionalForm];
Print["Simplified:"];
Print[energySimplified // TraditionalForm];
Print[""];

(* 5. Diagonalize (Eigenvalues) *)
evals = Eigenvalues[matM];

Print["--- Eigenvalues of M ---"];
Print[evals // TraditionalForm];

(*"
Output:

--- Energy E(uT) ---
Unsimplified:
TraditionalForm[(3*AL*uT^2)/4 + AT*uT^2]
Simplified:
TraditionalForm[((3*AL + 4*AT)*uT^2)/4]

--- Eigenvalues of M ---
TraditionalForm[{AL, AT}]
"*)
