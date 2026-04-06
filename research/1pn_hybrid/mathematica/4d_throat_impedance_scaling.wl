(* ---------------------------------------------------------------------- *)
(* Script: 4d_throat_impedance_scaling.wl                                  *)
(* Purpose: Verify the appendix-level throat-impedance bookkeeping used in *)
(* the paper's near-linear mass-radius discussion.                         *)
(*                                                                        *)
(* Audit note: this is a phenomenological scaling model anchored to the    *)
(* preferred aspect ratio from Paper V. It checks that finite impedance    *)
(* can keep k_eff close to 1 while preserving Z_th ~ a^-2.                *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

lambdaStar = (Sqrt[2] * Pi)/BesselJZero[0, 1];
sigmaFS = FullSimplify[1/(1 + lambdaStar^2)];
a0 = 1;
aa = Unique["a"];

Zth[a_] := c/a^2;
Mmodel[a_] := a * (1 + sigmaFS/(a^2 + a0^2));
kEffExpr = FullSimplify[aa * D[Log[Mmodel[aa]], aa]];
kEffDisplay = kEffExpr /. aa -> a;
kEff[a_] := kEffExpr /. aa -> a;

sampleRadii = {1, 2, 4, 8};

Print["--- Throat Impedance Scaling ---"];
Print["lambda_* = ", N[lambdaStar, 16]];
Print["sigma_fs = ", N[sigmaFS, 16]];
Print["Z_th(a) = ", Zth[a]];
Print["k_eff(a) = ", kEffDisplay];
Print["k_eff(1) = ", N[kEff[1], 16]];
Print["k_eff(2) = ", N[kEff[2], 16]];
Print["k_eff(4) = ", N[kEff[4], 16]];
Print["k_eff(8) = ", N[kEff[8], 16]];
Print["Limit[k_eff, a->Infinity] = ", Limit[kEff[a], a -> Infinity]];

(*"
Output:

--- Throat Impedance Scaling ---
lambda_* = 1.84748657712012805104268508521531981984`16.
sigma_fs = 0.22659260685243721683646187392498156403`16.
Z_th(a) = c/a^2
k_eff(a) = (2*(1 + a^2)^2*Pi^2 + (2 + a^2 + a^4)*BesselJZero[0, 1]^2)/(2*(1 + a^2)^2*Pi^2 + (2 + 3*a^2 + a^4)*BesselJZero[0, 1]^2)
k_eff(1) = 0.89823346841488270087027110485748949739`16.
k_eff(2) = 0.9306339333797371292256273632632720372`16.
k_eff(4) = 0.97524018419127971973858043005711623294`16.
k_eff(8) = 0.99315903045582128266328230336510598539`16.
Limit[k_eff, a->Infinity] = 1
"*)
