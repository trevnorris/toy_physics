(* ---------------------------------------------------------------------- *)
(* Script: brane_mode_resonance.wl *)
(* Calculation of fundamental 4D mode frequencies and aspect ratios *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

(* 1. Bessel Root Calculation *)
(* Python: sp.nsolve(J0, 2.4) *)
(* Mathematica: BesselJZero[0, 1] gives the exact symbolic first zero of J0 *)
x01 = BesselJZero[0, 1];
x01_numeric = N[x01]; (* Compute numeric value *)

Print["--- 1. Bessel Root ---"];
Print["First J0 root x_01 (Exact): ", x01];
Print["First J0 root x_01 (Numeric): ", x01_numeric];
Print[""];

(* ---------------------------------------------------------------------- *)
(* 2. Wavenumbers and Mode Frequency *)
(* ---------------------------------------------------------------------- *)
(* Define wavenumbers *)
k_r = x01 / a;
k_w = Pi / L;

(* Calculate Omega *)
(* Note: cs is the sound speed *)
omega = cs * Sqrt[k_r^2 + k_w^2];

Print["--- 2. Wavenumbers & Frequency ---"];
Print["Radial wavenumber k_r: ", k_r // TraditionalForm];
Print["Axial wavenumber k_w: ", k_w // TraditionalForm];
Print["Mode frequency omega: "];
Print[omega // TraditionalForm];
Print[""];

(* ---------------------------------------------------------------------- *)
(* 3. Fundamental 4D mode h(r, w, t) *)
(* ---------------------------------------------------------------------- *)
(* Python: h = A * sp.besselj(0, k_r * r) * sp.sin(k_w * w) * ... *)
(* Mathematica: I is the imaginary unit *)

h = A * BesselJ[0, k_r * r] * Sin[k_w * w] * Exp[-I * omega * t];

Print["--- 3. Fundamental 4D Mode h(r, w, t) ---"];
Print[h // TraditionalForm];
Print[""];

(* ---------------------------------------------------------------------- *)
(* 4. Aspect Ratio L/a check *)
(* From Paper IV result: L/a = sqrt(2)*pi/x01 *)
(* ---------------------------------------------------------------------- *)
LOverAExact = Sqrt[2] * Pi / x01;
LOverANumeric = N[LOverAExact];

Print["--- 4. Geometric Aspect Ratio ---"];
Print["Predicted L/a (Exact): "];
Print[LOverAExact // TraditionalForm];
Print["Predicted L/a (Numeric): ", LOverANumeric];

(*"
Output:

--- 1. Bessel Root ---
First J0 root x_01 (Exact): BesselJZero[0, 1]
First J0 root x_01 (Numeric): x01_numeric

--- 2. Wavenumbers & Frequency ---
Radial wavenumber k_r: TraditionalForm[k_r]
Axial wavenumber k_w: TraditionalForm[k_w]
Mode frequency omega:
TraditionalForm[cs*Sqrt[(k_r)^2 + (k_w)^2]]

--- 3. Fundamental 4D Mode h(r, w, t) ---
TraditionalForm[(A*BesselJ[0, r*(k_r)]*Sin[w*(k_w)])/E^(I*cs*t*Sqrt[(k_r)^2 + (k_w)^2])]

--- 4. Geometric Aspect Ratio ---
Predicted L/a (Exact):
TraditionalForm[(Sqrt[2]*Pi)/BesselJZero[0, 1]]
Predicted L/a (Numeric): 1.8474865771201279
"*)
