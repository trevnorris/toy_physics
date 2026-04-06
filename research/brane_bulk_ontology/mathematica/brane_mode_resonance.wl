(* ---------------------------------------------------------------------- *)
(* Script: brane_mode_resonance.wl                                        *)
(* Basic checks for the fundamental 4D throat mode and aspect ratio.      *)
(* ---------------------------------------------------------------------- *)
ClearAll["Global`*"]

x01 = BesselJZero[0, 1];
x01Numeric = N[x01];

kr = x01/a;
kw = Pi/L;
omega11Squared = cs^2 (kr^2 + kw^2);
omega11 = Symbol["omega11"];
h11 = BesselJ[0, kr*r] * Sin[kw*w] * Cos[omega11*t];
LOverAExact = Sqrt[2] * Pi / x01;
LOverANumeric = N[LOverAExact];

Print["--- 1. Bessel Root ---"];
Print["First J0 root x_01 (Exact): ", x01];
Print["First J0 root x_01 (Numeric): ", x01Numeric];
Print[""];

Print["--- 2. Fundamental Mode Data ---"];
Print["Radial wavenumber k_r ="];
Print[kr // TraditionalForm];
Print["Axial wavenumber k_w ="];
Print[kw // TraditionalForm];
Print["Dispersion relation omega_11^2 ="];
Print[omega11Squared // TraditionalForm];
Print[""];

Print["--- 3. Fundamental Throat Mode h_11(t, r, w) ---"];
Print[h11 // TraditionalForm];
Print[""];

Print["--- 4. Preferred Aspect Ratio ---"];
Print["Predicted L/a (Exact) ="];
Print[LOverAExact // TraditionalForm];
Print["Predicted L/a (Numeric): ", LOverANumeric];

(* 
Output:
--- 1. Bessel Root ---
First J0 root x_01 (Exact): BesselJZero[0, 1]
First J0 root x_01 (Numeric): 2.404825557695773

--- 2. Fundamental Mode Data ---
Radial wavenumber k_r =
TraditionalForm[BesselJZero[0, 1]/a]
Axial wavenumber k_w =
TraditionalForm[Pi/L]
Dispersion relation omega_11^2 =
TraditionalForm[cs^2*(Pi^2/L^2 + BesselJZero[0, 1]^2/a^2)]

--- 3. Fundamental Throat Mode h_11(t, r, w) ---
TraditionalForm[BesselJ[0, (r*BesselJZero[0, 1])/a]*Cos[omega11*t]*Sin[(Pi*w)/L]]

--- 4. Preferred Aspect Ratio ---
Predicted L/a (Exact) =
TraditionalForm[(Sqrt[2]*Pi)/BesselJZero[0, 1]]
Predicted L/a (Numeric): 1.8474865771201279
*)
