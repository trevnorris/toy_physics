(* ::Package:: *)
(*
  paper7_brane_momentum_appendix_min.wl

  Minimal Mathematica harness for Paper 7 Appendix:
  - weighted brane projection
  - projected continuity + leakage source
  - exact longitudinal identity (Helmholtz split)
  - projected brane momentum equation (structure) + Reynolds stress + momentum-leakage source

  Design goals:
  - keep everything HOLDING (no evaluation / no Integrate::idiv / no Div::sclr warnings)
  - print compact, paper-ready equations

  Usage:
    wolframscript -file paper7_brane_momentum_appendix_min.wl
*)

ClearAll["Global`*"];

(* ----------------------------- *)
(* Small pretty-printer helpers  *)
(* ----------------------------- *)

Section[title_String] := (
  Print["\n=============================================================================="];
  Print[title];
  Print["------------------------------------------------------------------------------"];
);

EqPrint[name_String, expr_] := (
  Print[name];
  Print[TraditionalForm[expr]];
  Print[""];
);

EqBlock[name_String, eqns_List] := (
  Print[name];
  Scan[(Print[TraditionalForm[#]]) &, eqns];
  Print[""];
);

(* ----------------------------- *)
(* Symbols / placeholders        *)
(* ----------------------------- *)

coords3 = {x, y, z};

(* Brane projection weight *)
W = Symbol["W"]; (* used as W[w] *)

(* Bulk / brane observables (formal placeholders) *)
rho      = Symbol["rho"];          (* bulk density rho(x,y,z,w,t) *)
Jw       = Symbol["Jw"];           (* bulk w-current component *)
Jxyz     = {Jx, Jy, Jz};            (* bulk xyz current components *)

rhobrane = Symbol["rhobrane"];   (* brane density *)
Jbrane   = Symbol["Jbrane"];     (* brane xyz current (vector) *)
vbrane   = Symbol["vbrane"];     (* brane velocity (vector) *)

Sleak  = Symbol["Sleak"];        (* brane continuity source from leakage *)

(* Helmholtz pieces (brane) *)
Phi = Symbol["Phi"];               (* longitudinal potential *)
vT  = Symbol["vT"];                (* transverse component (vector) *)

(* Momentum reduction placeholders *)
Pibrane = Symbol["Pibrane"];     (* projected momentum flux tensor Pibrane[i,a] *)
Rstress  = Symbol["R"];            (* Reynolds / stress correction R[i,a] *)
Smom    = Symbol["Smom"];        (* momentum leakage source Smom[i] *)

(* Force-side placeholders (keep generic) *)
Vconf = Symbol["Vconf"];
h     = Symbol["h"];
Q     = Symbol["Q"];
Ee     = Symbol["Ee"];
B     = Symbol["B"];

m = Symbol["m"]; q = Symbol["q"];

(* ----------------------------- *)
(* Appendix: continuity + Poisson hook *)
(* ----------------------------- *)

Section["APPENDIX: Brane continuity and exact longitudinal identity"];

EqBlock["Projection (definitions)", {
  HoldForm[rhobrane == Integrate[W[w]*rho, {w, -Infinity, Infinity}]],
  HoldForm[Jbrane == Integrate[W[w]*Jxyz, {w, -Infinity, Infinity}]],
  HoldForm[Jxyz == {Jx, Jy, Jz}],
  HoldForm[vbrane == Jbrane/rhobrane]
}];

EqBlock["Projected continuity + leakage source (identity)", {
  HoldForm[D[rhobrane, t] + Div[Jbrane, coords3] == Sleak],
  HoldForm[Sleak == -Limit[W[w]*Jw, w -> Infinity] + Limit[W[w]*Jw, w -> -Infinity] + Integrate[D[W[w], w]*Jw, {w, -Infinity, Infinity}]]
}];

EqBlock["Helmholtz decomposition (brane)", {
  HoldForm[vbrane == Grad[Phi, coords3] + vT],
  HoldForm[Div[vT, coords3] == 0]
}];

EqBlock["Exact longitudinal identity (from continuity + Helmholtz)", {
  (* Start from continuity in terms of vbrane = Jbrane/rhobrane: *)
  HoldForm[D[rhobrane, t] + Div[rhobrane*vbrane, coords3] == Sleak],

  (* Expand Div[rho*(Grad Phi + vT)] using the standard product rule identity.
     The last term drops if Div[vT]=0. *)
  HoldForm[rhobrane*Laplacian[Phi, coords3] == Sleak - D[rhobrane, t] - Grad[rhobrane, coords3].(Grad[Phi, coords3] + vT) - rhobrane*Div[vT, coords3]],
  HoldForm[rhobrane*Laplacian[Phi, coords3] == Sleak - D[rhobrane, t] - Grad[rhobrane, coords3].(Grad[Phi, coords3] + vT)]
}];

(* ----------------------------- *)
(* Appendix: brane momentum reduction (structure) *)
(* ----------------------------- *)

Section["APPENDIX: Projected brane momentum equation (structure + correction terms)"];

EqBlock["Bulk momentum (conservative form; i in {x,y,z}, j in {x,y,z,w})", {
  HoldForm[m*(D[rho*v[i], t] + Subscript[del, j][rho*v[i]*v[j]]) == q*rho*(Ee[i] + v[j]*B[i, j]) - rho*Subscript[del, i][Vconf + h[rho] + Q]]
}];

EqBlock["Projected brane momentum density and flux tensor", {
  HoldForm[Jbrane[i] == Integrate[W[w]*rho*v[i], {w, -Infinity, Infinity}]],
  HoldForm[Pibrane[i, a] == Integrate[W[w]*rho*v[i]*v[a], {w, -Infinity, Infinity}]],
  HoldForm[Rstress[i, a] == Pibrane[i, a] - rhobrane*vbrane[i]*vbrane[a]]
}];

EqBlock["Momentum leakage source (from the w-flux term)", {
  HoldForm[Smom[i] == -Limit[W[w]*rho*v[i]*v[w], w -> Infinity] + Limit[W[w]*rho*v[i]*v[w], w -> -Infinity] + Integrate[D[W[w], w]*rho*v[i]*v[w], {w, -Infinity, Infinity}]]
}];

EqBlock["Projected brane momentum equation (i in {x,y,z}; a in {x,y,z})", {
  HoldForm[m*(D[Jbrane[i], t] + Subscript[del, a][Pibrane[i, a]]) == Integrate[W[w]*(q*rho*(Ee[i] + v[j]*B[i, j]) - rho*Subscript[del, i][Vconf + h[rho] + Q]), {w, -Infinity, Infinity}] + m*Smom[i]],
  HoldForm[Pibrane[i, a] == rhobrane*vbrane[i]*vbrane[a] + Rstress[i, a]]
}];

EqPrint["(Optional alternate form) Reynolds/stress correction as a variance in w",
  HoldForm[Rstress[i, a] == Integrate[W[w]*rho*(v[i] - vbrane[i])*(v[a] - vbrane[a]), {w, -Infinity, Infinity}]]
];

Print["END"];

(*"
Output:

==============================================================================
APPENDIX: Brane continuity and exact longitudinal identity
------------------------------------------------------------------------------
Projection (definitions)
TraditionalForm[HoldForm[rhobrane == Integrate[W[w]*rho, {w, -Infinity, Infinity}]]]
TraditionalForm[HoldForm[Jbrane == Integrate[W[w]*Jxyz, {w, -Infinity, Infinity}]]]
TraditionalForm[HoldForm[Jxyz == {Jx, Jy, Jz}]]
TraditionalForm[HoldForm[vbrane == Jbrane/rhobrane]]

Projected continuity + leakage source (identity)
TraditionalForm[HoldForm[D[rhobrane, t] + Div[Jbrane, coords3] == Sleak]]
TraditionalForm[HoldForm[Sleak == -Limit[W[w]*Jw, w -> Infinity] + Limit[W[w]*Jw, w -> -Infinity] + Integrate[D[W[w], w]*Jw, {w, -Infinity, Infinity}]]]

Helmholtz decomposition (brane)
TraditionalForm[HoldForm[vbrane == Grad[Phi, coords3] + vT]]
TraditionalForm[HoldForm[Div[vT, coords3] == 0]]

Exact longitudinal identity (from continuity + Helmholtz)
TraditionalForm[HoldForm[D[rhobrane, t] + Div[rhobrane*vbrane, coords3] == Sleak]]
TraditionalForm[HoldForm[rhobrane*Laplacian[Phi, coords3] == Sleak - D[rhobrane, t] - Grad[rhobrane, coords3] . (Grad[Phi, coords3] + vT) - rhobrane*Div[vT, coords3]]]
TraditionalForm[HoldForm[rhobrane*Laplacian[Phi, coords3] == Sleak - D[rhobrane, t] - Grad[rhobrane, coords3] . (Grad[Phi, coords3] + vT)]]


==============================================================================
APPENDIX: Projected brane momentum equation (structure + correction terms)
------------------------------------------------------------------------------
Bulk momentum (conservative form; i in {x,y,z}, j in {x,y,z,w})
TraditionalForm[HoldForm[m*(D[rho*v[i], t] + Subscript[del, j][rho*v[i]*v[j]]) == q*rho*(Ee[i] + v[j]*B[i, j]) - rho*Subscript[del, i][Vconf + h[rho] + Q]]]

Projected brane momentum density and flux tensor
TraditionalForm[HoldForm[Jbrane[i] == Integrate[W[w]*rho*v[i], {w, -Infinity, Infinity}]]]
TraditionalForm[HoldForm[Pibrane[i, a] == Integrate[W[w]*rho*v[i]*v[a], {w, -Infinity, Infinity}]]]
TraditionalForm[HoldForm[Rstress[i, a] == Pibrane[i, a] - rhobrane*vbrane[i]*vbrane[a]]]

Momentum leakage source (from the w-flux term)
TraditionalForm[HoldForm[Smom[i] == -Limit[W[w]*rho*v[i]*v[w], w -> Infinity] + Limit[W[w]*rho*v[i]*v[w], w -> -Infinity] + Integrate[D[W[w], w]*rho*v[i]*v[w], {w, -Infinity, Infinity}]]]

Projected brane momentum equation (i in {x,y,z}; a in {x,y,z})
TraditionalForm[HoldForm[m*(D[Jbrane[i], t] + Subscript[del, a][Pibrane[i, a]]) == Integrate[W[w]*(q*rho*(Ee[i] + v[j]*B[i, j]) - rho*Subscript[del, i][Vconf + h[rho] + Q]), {w, -Infinity, Infinity}] + m*Smom[i]]]
TraditionalForm[HoldForm[Pibrane[i, a] == rhobrane*vbrane[i]*vbrane[a] + Rstress[i, a]]]

(Optional alternate form) Reynolds/stress correction as a variance in w
TraditionalForm[HoldForm[Rstress[i, a] == Integrate[W[w]*rho*(v[i] - vbrane[i])*(v[a] - vbrane[a]), {w, -Infinity, Infinity}]]]

END
"*)
