(* ::Package:: *)

(*
========================================================
Paper VIII (referee add-on v4): axial gauge + brane-source consistency
========================================================

This version fixes the v3 reduction-rule bug that prevented the zero-mode
replacement rules from applying cleanly (Table iterator / rule-flattening issue).

What this add-on establishes (referee-walkable):
  (1) Axial gauge reachability: a (formal) gauge parameter chi can be constructed
      such that Aw + d_w chi = 0, hence Aw=0 is reachable.
  (2) Residual gauge freedom in axial gauge is w-independent chi0(t,x,y,z), i.e.
      ordinary 3+1 gauge freedom on the brane.
  (3) For a brane-localized source J^mu(x) DiracDelta[w] with Jw=0,
      5D current conservation d_M J^M = 0 reduces to 3+1 continuity.
  (4) Under the photon zero-mode truncation
        Aw = 0,   d_w A_mu = 0,
      the gauge-invariant 5D Maxwell operator d_M(Z F^{MN}) - mu0 J^N reduces to
        Op^nu(w) = Z(w) * (d_mu f^{mu nu}) - mu0 * Jb^nu * DiracDelta[w]
      and integrating over w gives the standard 3+1 Maxwell equations
        d_mu f^{mu nu} = mu0eff * Jb^nu,
      with mu0eff = mu0 / Zint and Zint = ∫dw Z(w).

Notes:
  - This add-on uses the *gauge-invariant* Maxwell equations (no gauge-fixing term).
  - No global symbol names contain underscores; underscores appear only in WL patterns.
  - Any paper-like notation appears only in printed strings.

Conventions:
  - Coordinates: (t,x,y,z,w)
  - Metric signature: (-,+,+,+,+)
  - A_M are covariant components (A0,Ax,Ay,Az,Aw).
  - J^M are contravariant components (J0,Jx,Jy,Jz,Jw).
========================================================
*)

Print["\n========================================"]; 
Print["Paper VIII add-on v4: axial gauge + brane-source consistency"]; 
Print["========================================\n"]; 

ClearAll["Global`*"];

(* ---------- Assumptions + simplification helpers ---------- *)
$Assumptions = {Element[{mu0, lambdaConf}, Reals], mu0 > 0, lambdaConf > 0};
Simp[expr_] := FullSimplify[Together[expr], Assumptions -> $Assumptions];
ZeroQAssume[expr_] := TrueQ[FullSimplify[expr == 0, Assumptions -> $Assumptions]];

(* ---------- Coordinates + metric ---------- *)
coords5 = {t, x, y, z, w};
dim5 = Length[coords5];
eta5 = DiagonalMatrix[{-1, 1, 1, 1, 1}];
etaU5 = eta5;

Print["--- Objects ---"]; 
Print[HoldForm[coords5], " -> ", coords5];
Print[HoldForm[eta5], " -> ", eta5];

(* Localization profile (Gaussian) *)
Zprof[w_] := Exp[-w^2/lambdaConf^2];
Print[HoldForm[Zprof[w]], " -> ", Zprof[w]];

Zint = Simp[Integrate[Zprof[w], {w, -Infinity, Infinity}]];
mu0eff = Simp[mu0/Zint];

Print["\nZint = Integrate[Z(w) dw] -> ", Zint];
Print["mu0eff = mu0/Zint -> ", mu0eff];

(* ---------- Derivative operators ---------- *)
dDown5[i_][expr_] := D[expr, coords5[[i]]];

(* ---------- Field components (covariant) ---------- *)
ADown5 = {
  A0[t, x, y, z, w],
  Ax[t, x, y, z, w],
  Ay[t, x, y, z, w],
  Az[t, x, y, z, w],
  Aw[t, x, y, z, w]
};

(* ---------- (1) Axial gauge reachability ---------- *)

Print["\n--- (1) Axial gauge reachability: choose chi so that Aw + d_w chi = 0 ---\n"]; 

chiFun = chi[t, x, y, z, w];

(* Constructive choice: chiAx = -∫_0^w Aw(...,wp) dwp *)
chiAx = -Integrate[Aw[t, x, y, z, wp], {wp, 0, w}];
AwAxial = Simp[Aw[t, x, y, z, w] + D[chiAx, w]];

Print["Constructive choice (string):  chiAx = -∫_0^w Aw(t,x,y,z,wp) dwp"]; 
Print["Check Aw' = Aw + d_w chiAx  (should be 0):"]; 
Print[AwAxial];
If[ZeroQAssume[AwAxial],
  Print["OK: Axial gauge Aw=0 is reachable by a gauge transformation (formal check)."],
  Print["WARNING: Aw' did not simplify to 0 (symbolic integration/derivative issue)."]
];

(* Residual gauge freedom: d_w chi = 0 => chi = chi0(t,x,y,z). *)
Print["\nResidual gauge freedom in axial gauge: d_w chi = 0  =>  chi is w-independent."]; 
chi0expr = chi0[t, x, y, z];
residualCheck = Simp[D[chi0expr, w]];
Print["Check: d_w chi0(t,x,y,z) (should be 0): ", residualCheck];
If[ZeroQAssume[residualCheck],
  Print["OK: Residual gauge parameter chi0(t,x,y,z) is the usual 3+1 gauge freedom on the brane."],
  Print["WARNING: Unexpected nonzero residual."]
];

(* ---------- (2) Brane-localized current conservation ---------- *)

Print["\n--- (2) Brane-localized current: d_M J^M reduces to 3+1 continuity ---\n"]; 

(* Brane source ansatz: J^mu(x,w) = Jb^mu(x) DiracDelta[w], Jw=0 *)
J0[t_, x_, y_, z_, w_] := J0b[t, x, y, z] DiracDelta[w];
Jx[t_, x_, y_, z_, w_] := Jxb[t, x, y, z] DiracDelta[w];
Jy[t_, x_, y_, z_, w_] := Jyb[t, x, y, z] DiracDelta[w];
Jz[t_, x_, y_, z_, w_] := Jzb[t, x, y, z] DiracDelta[w];
Jw[t_, x_, y_, z_, w_] := 0;

JUp5 = {J0[t, x, y, z, w], Jx[t, x, y, z, w], Jy[t, x, y, z, w], Jz[t, x, y, z, w], Jw[t, x, y, z, w]};

(* Flat divergence (since we are working in flat 5D here): d_M J^M *)
divJ5 = Simp[
  D[J0[t, x, y, z, w], t] +
  D[Jx[t, x, y, z, w], x] +
  D[Jy[t, x, y, z, w], y] +
  D[Jz[t, x, y, z, w], z] +
  D[Jw[t, x, y, z, w], w]
];

Print["d_M J^M (5D) = "]; 
Print[divJ5];

divJ4 = Simp[Integrate[divJ5, {w, -Infinity, Infinity}]];
Print["Integrate over w: ∫dw d_M J^M = "]; 
Print[divJ4];

Print["So if d_M J^M = 0 distributionally, then the induced 3+1 current is conserved:"]; 
Print["  d_t J0b + d_x Jxb + d_y Jyb + d_z Jzb = 0."]; 

(* ---------- (3) Zero-mode reduction of gauge-invariant 5D Maxwell EOM ---------- *)

Print["\n--- (3) Zero-mode reduction of gauge-invariant Maxwell EOM: d_M(Z F^{MN}) = mu0 J^N ---\n"]; 

(* Define F_{MN} and F^{MN} from covariant A_M *)
Fdown5[m_, n_] := dDown5[m][ADown5[[n]]] - dDown5[n][ADown5[[m]]];
Fup5[m_, n_] := Simp[Sum[etaU5[[m, a]] etaU5[[n, b]] Fdown5[a, b], {a, 1, dim5}, {b, 1, dim5}]];

(* Maxwell operator: d_M(Z F^{MN}) - mu0 J^N *)
MaxwellOp5[n_] := Simp[Sum[dDown5[m][Zprof[w] * Fup5[m, n]], {m, 1, dim5}] - mu0 * JUp5[[n]]];

(* Robust zero-mode reduction rules (no Table destructuring):
   - Aw = 0 and all derivatives of Aw are 0
   - A_mu(t,x,y,z,w) -> A_mu^b(t,x,y,z)
   - Any derivative with at least one w-derivative is 0
   - Pure (t,x,y,z) derivatives map to derivatives of brane fields
*)
pairsMu = {{A0, A0b}, {Ax, Axb}, {Ay, Ayb}, {Az, Azb}};

awRules = {
  Aw[tt_, xx_, yy_, zz_, ww_] :> 0,
  Derivative[__][Aw][tt_, xx_, yy_, zz_, ww_] :> 0
};

muRules = Flatten@Table[
  With[{Acomp = pair[[1]], AcompB = pair[[2]]},
    {
      (* field itself *)
      Acomp[tt_, xx_, yy_, zz_, ww_] :> AcompB[tt, xx, yy, zz],

      (* any derivative with at least one w-derivative is set to 0 *)
      Derivative[dt_, dx_, dy_, dz_, dw_][Acomp][tt_, xx_, yy_, zz_, ww_] :> 0 /; dw > 0,

      (* pure brane derivatives: drop w and map to brane field *)
      Derivative[dt_, dx_, dy_, dz_, 0][Acomp][tt_, xx_, yy_, zz_, ww_] :> Derivative[dt, dx, dy, dz][AcompB][tt, xx, yy, zz]
    }
  ],
  {pair, pairsMu}
];

zeroModeRules = Join[awRules, muRules];

nuNames = {"t", "x", "y", "z", "w"};

maxOp5zm = Table[
  Simp[MaxwellOp5[n] /. zeroModeRules],
  {n, 1, 5}
];

Print["Maxwell operator under zero-mode + brane source (distributional in w):"]; 
Do[
  Print["  N=", nuNames[[n]], " -> ", maxOp5zm[[n]]],
  {n, 1, 5}
];

(* w-component should vanish under Aw=0 and d_w A_mu=0 *)
Print["\nCheck w-component under the ansatz (should be 0):"]; 
Print["  N=w -> ", maxOp5zm[[5]]];
If[ZeroQAssume[maxOp5zm[[5]]],
  Print["OK: w-component vanishes under axial gauge + zero-mode (no off-brane flux, Jw=0)."],
  Print["WARNING: w-component did not simplify to 0 (unexpected)."]
];

(* Integrate over w to obtain the effective 3+1 equations *)
maxOp4int = Table[
  Simp[Integrate[maxOp5zm[[n]], {w, -Infinity, Infinity}]],
  {n, 1, 4}
];

Print["\nIntegrate over w (N = t,x,y,z):"]; 
Do[
  Print["  ∫dw Op(N=", nuNames[[n]], ") -> ", maxOp4int[[n]]],
  {n, 1, 4}
];

(* ---------- Explicit 3+1 operator check: maxOp4int == Zint*d_mu f^{mu nu} - mu0*Jb^nu ---------- *)

coords4 = {t, x, y, z};
dim4 = Length[coords4];
eta4 = DiagonalMatrix[{-1, 1, 1, 1}];
etaU4 = eta4;

dDown4[i_][expr_] := D[expr, coords4[[i]]];

ADown4 = {A0b[t, x, y, z], Axb[t, x, y, z], Ayb[t, x, y, z], Azb[t, x, y, z]};

JbUp4 = {J0b[t, x, y, z], Jxb[t, x, y, z], Jyb[t, x, y, z], Jzb[t, x, y, z]};

fdown4[m_, n_] := dDown4[m][ADown4[[n]]] - dDown4[n][ADown4[[m]]];
fup4[m_, n_] := Simp[Sum[etaU4[[m, a]] etaU4[[n, b]] fdown4[a, b], {a, 1, dim4}, {b, 1, dim4}]];

braneDivF[nu_] := Simp[Sum[dDown4[mu][fup4[mu, nu]], {mu, 1, dim4}]];

expectedInt[nu_] := Simp[Zint * braneDivF[nu] - mu0 * JbUp4[[nu]]];

Print["\nCheck integrated equations vs explicit 3+1 form:"]; 
Do[
  Module[{res = Simp[maxOp4int[[nu]] - expectedInt[nu]]},
    Print["  nu=", nuNames[[nu]], "  residual (should be 0): ", res];
    If[ZeroQAssume[res],
      Null,
      Print["  WARNING: residual not zero (convention mismatch / simplification issue)."]
    ]
  ],
  {nu, 1, 4}
];

Print["\nInterpretation:"]; 
Print["  Integrated EOM:  Zint * (d_mu f^{mu nu}) - mu0 * Jb^nu = 0."]; 
Print["  Divide by Zint:  d_mu f^{mu nu} = (mu0/Zint) Jb^nu = mu0eff * Jb^nu."]; 

Print["\n========== End add-on v4 =========="]; 

(*"
Output:

========================================
Paper VIII add-on v4: axial gauge + brane-source consistency
========================================

--- Objects ---
HoldForm[coords5] -> {t, x, y, z, w}
HoldForm[eta5] -> {{-1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}}
HoldForm[Zprof[w]] -> E^(-(w^2/lambdaConf^2))

Zint = Integrate[Z(w) dw] -> lambdaConf*Sqrt[Pi]
mu0eff = mu0/Zint -> mu0/(lambdaConf*Sqrt[Pi])

--- (1) Axial gauge reachability: choose chi so that Aw + d_w chi = 0 ---

Constructive choice (string):  chiAx = -∫_0^w Aw(t,x,y,z,wp) dwp
Check Aw' = Aw + d_w chiAx  (should be 0):
0
OK: Axial gauge Aw=0 is reachable by a gauge transformation (formal check).

Residual gauge freedom in axial gauge: d_w chi = 0  =>  chi is w-independent.
Check: d_w chi0(t,x,y,z) (should be 0): 0
OK: Residual gauge parameter chi0(t,x,y,z) is the usual 3+1 gauge freedom on the brane.

--- (2) Brane-localized current: d_M J^M reduces to 3+1 continuity ---

d_M J^M (5D) =
DiracDelta[w]*Derivative[0, 0, 0, 1][Jzb][t, x, y, z] + DiracDelta[w]*Derivative[0, 0, 1, 0][Jyb][t, x, y, z] + DiracDelta[w]*Derivative[0, 1, 0, 0][Jxb][t, x, y, z] + DiracDelta[w]*Derivative[1, 0, 0, 0][J0b][t, x, y, z]
Integrate over w: ∫dw d_M J^M =
Derivative[0, 0, 0, 1][Jzb][t, x, y, z] + Derivative[0, 0, 1, 0][Jyb][t, x, y, z] + Derivative[0, 1, 0, 0][Jxb][t, x, y, z] + Derivative[1, 0, 0, 0][J0b][t, x, y, z]
So if d_M J^M = 0 distributionally, then the induced 3+1 current is conserved:
  d_t J0b + d_x Jxb + d_y Jyb + d_z Jzb = 0.

--- (3) Zero-mode reduction of gauge-invariant Maxwell EOM: d_M(Z F^{MN}) = mu0 J^N ---

Maxwell operator under zero-mode + brane source (distributional in w):
  N=t -> -((mu0*DiracDelta[w]*J0b[t, x, y, z] + Derivative[0, 0, 0, 2][A0b][t, x, y, z] + Derivative[0, 0, 2, 0][A0b][t, x, y, z] + Derivative[0, 2, 0, 0][A0b][t, x, y, z] - Derivative[1, 0, 0, 1][Azb][t, x, y, z] - Derivative[1, 0, 1, 0][Ayb][t, x, y, z] - Derivative[1, 1, 0, 0][Axb][t, x, y, z])/E^(w^2/lambdaConf^2))
  N=x -> -((mu0*DiracDelta[w]*Jxb[t, x, y, z] - Derivative[0, 0, 0, 2][Axb][t, x, y, z] - Derivative[0, 0, 2, 0][Axb][t, x, y, z] + Derivative[0, 1, 0, 1][Azb][t, x, y, z] + Derivative[0, 1, 1, 0][Ayb][t, x, y, z] - Derivative[1, 1, 0, 0][A0b][t, x, y, z] + Derivative[2, 0, 0, 0][Axb][t, x, y, z])/E^(w^2/lambdaConf^2))
  N=y -> -((mu0*DiracDelta[w]*Jyb[t, x, y, z] - Derivative[0, 0, 0, 2][Ayb][t, x, y, z] + Derivative[0, 0, 1, 1][Azb][t, x, y, z] + Derivative[0, 1, 1, 0][Axb][t, x, y, z] - Derivative[0, 2, 0, 0][Ayb][t, x, y, z] - Derivative[1, 0, 1, 0][A0b][t, x, y, z] + Derivative[2, 0, 0, 0][Ayb][t, x, y, z])/E^(w^2/lambdaConf^2))
  N=z -> -((mu0*DiracDelta[w]*Jzb[t, x, y, z] + Derivative[0, 0, 1, 1][Ayb][t, x, y, z] - Derivative[0, 0, 2, 0][Azb][t, x, y, z] + Derivative[0, 1, 0, 1][Axb][t, x, y, z] - Derivative[0, 2, 0, 0][Azb][t, x, y, z] - Derivative[1, 0, 0, 1][A0b][t, x, y, z] + Derivative[2, 0, 0, 0][Azb][t, x, y, z])/E^(w^2/lambdaConf^2))
  N=w -> 0

Check w-component under the ansatz (should be 0):
  N=w -> 0
OK: w-component vanishes under axial gauge + zero-mode (no off-brane flux, Jw=0).

Integrate over w (N = t,x,y,z):
  ∫dw Op(N=t) -> -(mu0*J0b[t, x, y, z]) + lambdaConf*Sqrt[Pi]*(-Derivative[0, 0, 0, 2][A0b][t, x, y, z] - Derivative[0, 0, 2, 0][A0b][t, x, y, z] - Derivative[0, 2, 0, 0][A0b][t, x, y, z] + Derivative[1, 0, 0, 1][Azb][t, x, y, z] + Derivative[1, 0, 1, 0][Ayb][t, x, y, z] + Derivative[1, 1, 0, 0][Axb][t, x, y, z])
  ∫dw Op(N=x) -> -(mu0*Jxb[t, x, y, z]) + lambdaConf*Sqrt[Pi]*(Derivative[0, 0, 0, 2][Axb][t, x, y, z] + Derivative[0, 0, 2, 0][Axb][t, x, y, z] - Derivative[0, 1, 0, 1][Azb][t, x, y, z] - Derivative[0, 1, 1, 0][Ayb][t, x, y, z] + Derivative[1, 1, 0, 0][A0b][t, x, y, z] - Derivative[2, 0, 0, 0][Axb][t, x, y, z])
  ∫dw Op(N=y) -> -(mu0*Jyb[t, x, y, z]) + lambdaConf*Sqrt[Pi]*(Derivative[0, 0, 0, 2][Ayb][t, x, y, z] - Derivative[0, 0, 1, 1][Azb][t, x, y, z] - Derivative[0, 1, 1, 0][Axb][t, x, y, z] + Derivative[0, 2, 0, 0][Ayb][t, x, y, z] + Derivative[1, 0, 1, 0][A0b][t, x, y, z] - Derivative[2, 0, 0, 0][Ayb][t, x, y, z])
  ∫dw Op(N=z) -> -(mu0*Jzb[t, x, y, z]) + lambdaConf*Sqrt[Pi]*(-Derivative[0, 0, 1, 1][Ayb][t, x, y, z] + Derivative[0, 0, 2, 0][Azb][t, x, y, z] - Derivative[0, 1, 0, 1][Axb][t, x, y, z] + Derivative[0, 2, 0, 0][Azb][t, x, y, z] + Derivative[1, 0, 0, 1][A0b][t, x, y, z] - Derivative[2, 0, 0, 0][Azb][t, x, y, z])

Check integrated equations vs explicit 3+1 form:
  nu=t  residual (should be 0): 0
  nu=x  residual (should be 0): 0
  nu=y  residual (should be 0): 0
  nu=z  residual (should be 0): 0

Interpretation:
  Integrated EOM:  Zint * (d_mu f^{mu nu}) - mu0 * Jb^nu = 0.
  Divide by Zint:  d_mu f^{mu nu} = (mu0/Zint) Jb^nu = mu0eff * Jb^nu.

========== End add-on v4 ==========
"*)
