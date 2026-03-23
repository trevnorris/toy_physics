(* ::Package:: *)

(*
Paper 7 Hard-Mode 4D ("Stiff Water") -- Master Derivation Harness (Paper-Polished)

Goals (unchanged physics):
  1) Encode the frozen model stack (GNLS + stiff EOS + projection/leakage + geometry forces + localized Maxwell).
  2) Print a referee-clean "paper form" bundle first (compact, notation-consistent).
  3) Optionally print an appendix with expanded component details.

This file is a presentation polish of 4d_master_derivations.wl. The algebra and model
choices are unchanged; the printouts are rewritten to avoid giant expansions and
notation drift.

Toggles:
  - $UseGauge: include minimal coupling in the internal definitions.
  - $PrintAppendix: if True, prints expanded component forms after the paper bundle.
  - $SimplifyMode: "Off" | "Light" | "Full".
*)

(* ----------------------------- *)
(* 0) Global setup + utilities   *)
(* ----------------------------- *)

(* Preserve toggle overrides if set before loading this file.
   You can set either the $-prefixed variables (recommended) or the non-$ aliases:
      $PrintAppendix / PrintAppendix
      $UseGauge      / UseGauge
      $SimplifyMode  / SimplifyMode
*)
Module[{sm, pa, ug},
  sm = Which[
    ValueQ[$SimplifyMode], $SimplifyMode,
    ValueQ[SimplifyMode], SimplifyMode,
    True, "Light"
  ];
  pa = Which[
    ValueQ[$PrintAppendix], $PrintAppendix,
    ValueQ[PrintAppendix], PrintAppendix,
    True, False
  ];
  ug = Which[
    ValueQ[$UseGauge], $UseGauge,
    ValueQ[UseGauge], UseGauge,
    True, True
  ];
  ClearAll["Global`*"];
  $HistoryLength = 1;
  $SimplifyMode  = sm;
  $PrintAppendix = TrueQ[pa];
  $UseGauge      = TrueQ[ug];
];

$PrintAppendix = True;
$UseGauge = True;
$SimplifyMode = "Light";

Print["Toggles: $UseGauge=", $UseGauge, "  |  $PrintAppendix=", $PrintAppendix, "  |  $SimplifyMode=", $SimplifyMode];

(* Pretty section printer *)
Section[title_String] := (
  Print["\n", Style[StringRepeat["=", 78], 12, GrayLevel[0.5]]];
  Print[Style[title, 16, Bold]];
  Print[Style[StringRepeat["-", 78], 12, GrayLevel[0.5]]];
);

EqPrint[name_String, expr_] := (
  Print[Style[name, 13, Bold]];
  Print[TraditionalForm[expr]];
  Print[""];
);

EqBlock[name_String, eqns_List] := (
  Print[Style[name, 13, Bold]];
  Scan[(Print[TraditionalForm[#]]) &, eqns];
  Print[""];
);

SmartSimplify[expr_] := Which[
  $SimplifyMode === "Off", expr,
  $SimplifyMode === "Light", Simplify[expr],
  True, FullSimplify[expr]
];

(* ---- TeX export helpers ---- *)
TeXString[expr_] := ToString[TeXForm[expr]];

ExportTeX[file_String, expr_] := Export[file, TeXString[expr], "String"];

(* Export a list of equations to a minimal LaTeX align block.

   IMPORTANT: Do NOT strip HoldForm before TeX conversion.
   If HoldForm is removed, Mathematica will try to evaluate objects like:
     - Div[ vT, {x,y,z} ]     (vT is symbolic)
     - Integrate[ J^N, dw ]   (J is unspecified)
     - Geometry-force integrals with symbolic rho
   which produces Div::sclr / Integrate::idiv warnings and can change the
   expressions. For export we want *purely formal* printing.
*)
SetAttributes[ExportPaperTeX, HoldAll];
ExportPaperTeX[file_String, eqns_List, replacements_: {}] := Module[
  {wrapHold, held, lines, body},
  wrapHold[e_] := If[Head[e] === HoldForm, e, HoldForm[e]];
  held = wrapHold /@ eqns;
  lines = TeXString /@ held;
  body = StringRiffle[("  " <> # <> "\\\\") & /@ lines, "\n"];
  body = StringReplace[body, replacements];
  Export[file,
    StringRiffle[{"\\begin{align}", body, "\\end{align}"}, "\n"],
    "String"
  ]
];

(* Default (optional) string replacements for nicer LaTeX symbols/macros. *)
DefaultTeXReplacements = {
  "\\text{psib}" -> "\\bar{\\psi}",
  "\\text{psi}"  -> "\\psi",
  "\\text{rho}"  -> "\\rho",
  "psib" -> "\\bar{\\psi}",
  "psi"  -> "\\psi",
  "rho"  -> "\\rho",
  "Vconf" -> "V_{\\rm conf}",
  "L_psi" -> "\\mathcal{L}_\\psi",
  "L_EM"  -> "\\mathcal{L}_{\\rm EM}"
};


(* ----------------------------- *)
(* 1) Coordinates, parameters    *)
(* ----------------------------- *)

Section["PAPER FORM: Conventions, coordinates, frozen EOS"]; 

(* Independent variables *)
t = Symbol["t"]; x = Symbol["x"]; y = Symbol["y"]; z = Symbol["z"]; w = Symbol["w"];
coords4 = {x, y, z, w};
coords5 = {t, x, y, z, w};

(* Parameters/constants *)
hbar = Symbol["hbar"]; hPlanck = Symbol["h"]; m = Symbol["m"]; etaQ = Symbol["etaQ"]; eStar = Symbol["eStar"]; qstar = Symbol["qstar"]; eEff = Symbol["eEff"]; qeff = Symbol["qeff"]; K = Symbol["K"];
mu0  = Symbol["mu0"]; xi = Symbol["xi"]; rhoBg = Symbol["rhoBg"];

(* Geometry knobs *)
a = Symbol["a"]; L = Symbol["L"];

(* Confinement/gate knobs (Family-1 baseline) *)
OmegaOut = Symbol["OmegaOut"]; OmegaIn = Symbol["OmegaIn"];
deltaPar = Symbol["deltaPar"]; epsW = Symbol["epsW"];

(* Localization length for Z(w) *)
lambdaConf = Symbol["lambdaConf"];

(* Geometry energy knobs *)
Pvac = Symbol["Pvac"]; sigma = Symbol["sigma"]; kappaB = Symbol["kappaB"];

(* Index conventions (paper):
   - spatial i,j in {x,y,z,w}
   - spacetime M,N in {0,x,y,z,w} with 0 == t
   - metric signature eta = diag(-,+,+,+,+)
*)
eta = DiagonalMatrix[{-1, 1, 1, 1, 1}];

EqBlock["Index + metric conventions", {
  HoldForm[coords4 == {x, y, z, w}],
  HoldForm[coords5 == {t, x, y, z, w}],
  HoldForm[eta == DiagonalMatrix[{-1, 1, 1, 1, 1}]]
}];

EqBlock["Charge ontology (paper)", {
  HoldForm[qstar == etaQ*eStar],
  HoldForm[etaQ^2 == 1],
  HoldForm[eStar > 0]
}];

(* Frozen EOS ladder: P = K rho^5, U = K/4 rho^5, h = 5K/4 rho^4 *)
Uof[r_] := (K/4) r^5;
hof[r_] := D[Uof[r], r];
Pof[r_] := r*hof[r] - Uof[r];
cs2of[r_] := D[Pof[r], r];

EqBlock["EOS ladder", {
  HoldForm[P[r]] -> Pof[r],
  HoldForm[U[r]] -> Uof[r],
  HoldForm[h[r]] -> hof[r],
  HoldForm[cs2[r]] -> cs2of[r]
}];

EqBlock["EOS checks", {
  HoldForm[P[r] == K r^5] -> SmartSimplify[Pof[r] == K r^5],
  HoldForm[h[r] == (5 K/4) r^4] -> SmartSimplify[hof[r] == (5 K/4) r^4],
  HoldForm[cs2[r] == 5 K r^4] -> SmartSimplify[cs2of[r] == 5 K r^4]
}];

(* ----------------------------- *)
(* 2) Fields + covariant pieces  *)
(* ----------------------------- *)

Section["PAPER FORM: Matter field, minimal coupling, currents"]; 

(* Fields (bulk) *)
psiF  = psi[t, x, y, z, w];
psibF = psib[t, x, y, z, w];

(* Gauge field components A_M (M=0..4) *)
A0F = A0[t, x, y, z, w];
AxF = Ax[t, x, y, z, w];
AyF = Ay[t, x, y, z, w];
AzF = Az[t, x, y, z, w];
AwF = Aw[t, x, y, z, w];

Avec = {AxF, AyF, AzF, AwF};

rhoF = psibF*psiF;

(* Gauge-covariant derivatives (internal definitions used for checks) *)
Dtcov[expr_, useGauge_:$UseGauge] := D[expr, t] + If[useGauge, I*qstar/hbar*A0F*expr, 0];
DtcovBar[expr_, useGauge_:$UseGauge] := D[expr, t] - If[useGauge, I*qstar/hbar*A0F*expr, 0];

Di[expr_, var_, Acomp_, useGauge_:$UseGauge] := D[expr, var] - If[useGauge, I*qstar/hbar*Acomp*expr, 0];
DiBar[expr_, var_, Acomp_, useGauge_:$UseGauge] := D[expr, var] + If[useGauge, I*qstar/hbar*Acomp*expr, 0];

(* Covariant Laplacian D_i D_i (i in {x,y,z,w}) acting on a scalar *)
D2[expr_, useGauge_:$UseGauge] := Module[{vars = coords4, Acomps = Avec},
  Sum[
    D[Di[expr, vars[[i]], Acomps[[i]], useGauge], vars[[i]]] -
      If[useGauge, I*qstar/hbar*Acomps[[i]]*Di[expr, vars[[i]], Acomps[[i]], useGauge], 0]
  , {i, 1, 4}]
];

(* Gauge-covariant mass/number current components *)
Jcomp[var_, Acomp_, useGauge_:$UseGauge] := (hbar/(2 I m))*(
  psibF*Di[psiF, var, Acomp, useGauge] - psiF*DiBar[psibF, var, Acomp, useGauge]
);

Jx = Jcomp[x, AxF]; Jy = Jcomp[y, AyF]; Jz = Jcomp[z, AzF]; Jw = Jcomp[w, AwF];

(* Paper definitions (compact, notation-consistent) *)
EqBlock["Minimal coupling (paper definitions)", {
  HoldForm[qstar == etaQ*eStar],
  HoldForm[Subscript[D, t][psi] == Subscript[del, t][psi] + (I qstar/hbar) Subscript[A, 0] psi],
  HoldForm[Subscript[D, t][psib] == Subscript[del, t][psib] - (I qstar/hbar) Subscript[A, 0] psib],
  HoldForm[Subscript[D, i][psi] == Subscript[del, i][psi] - (I qstar/hbar) Subscript[A, i] psi],
  HoldForm[Subscript[D, i][psib] == Subscript[del, i][psib] + (I qstar/hbar) Subscript[A, i] psib],
  HoldForm[rho == psib psi]
}];

EqBlock["Current + continuity (paper definitions)", {
  HoldForm[Jmass[i] == (hbar/m) Im[psib Subscript[D, i][psi]]],
  HoldForm[Jch[i] == qstar Jmass[i]],
  HoldForm[Qtot == qstar],
  HoldForm[Subscript[del, t][rho] + Subscript[del, i][rho v[i]] == 0]
}];

(* Optional: if gauge is disabled, print the explicit specialization used for "no-gauge" baseline runs.
   (Paper equations remain gauge-covariant; this block is only a readability aid.) *)
If[!$UseGauge,
  EqBlock["No-gauge specialization (baseline; set A_M = 0)", {
    HoldForm[Subscript[A, 0] == 0],
    HoldForm[Subscript[A, i] == 0],
    HoldForm[Subscript[D, t] == Subscript[del, t]],
    HoldForm[Subscript[D, i] == Subscript[del, i]],
    HoldForm[v[i] == (hbar/m) Subscript[del, i][theta]],
    HoldForm[Omega[i, j] == 0]
  }];
];

(* For appendix / explicit work, keep the fully explicit internal component definitions: *)
If[$PrintAppendix,
  EqPrint["(Appendix) Explicit density and current components", {
    HoldForm[rho] -> rhoF,
    HoldForm[Jx] -> Jx,
    HoldForm[Jy] -> Jy,
    HoldForm[Jz] -> Jz,
    HoldForm[Jw] -> Jw
  }];
];

(* ----------------------------- *)
(* 3) Frozen confinement V_conf  *)
(* ----------------------------- *)

Section["PAPER FORM: Confinement potential V_conf (Family-1 baseline)"]; 

Sstep[u_] := (1/2) (1 + Tanh[u]);
SmoothAbs[u_] := Sqrt[u^2 + epsW^2];

R3 = Sqrt[x^2 + y^2 + z^2];
Gr = Exp[-R3^4/a^4];
Gw = 1 - Sstep[(SmoothAbs[w] - L/2)/deltaPar];
Ggate = Gr*Gw;

OmegaW2 = OmegaOut^2 - (OmegaOut^2 - OmegaIn^2)*Ggate;
Vbrane = (1/2)*m*OmegaW2*w^2;

(* Optional wall/cap term:
   - Baseline harness: set to 0 (so the internal algebra matches the prior harness).
   - To include explicit walls/endcaps, replace VwallcapFn with a smooth, differentiable barrier.
*)
VwallcapFn[x_, y_, z_, w_, a_, L_] := 0;
Vwallcap = VwallcapFn[x, y, z, w, a, L];
Vconf = Vbrane + Vwallcap;

EqBlock["Family-1 modulated brane trap (with optional wall/cap add-on)", {
  HoldForm[Vconf] -> HoldForm[Vbrane + Vwallcap],
  HoldForm[Vbrane] -> Vbrane,
  HoldForm[Vwallcap] -> Vwallcap
}];

EqBlock["Required geometry derivatives (symbolic)", {
  HoldForm[dVda] -> D[Vconf, a],
  HoldForm[dVdL] -> D[Vconf, L]
}];

(* ----------------------------- *)
(* 4) Matter action -> GNLS      *)
(* ----------------------------- *)

Section["PAPER FORM: Matter Lagrangian density and GNLS (stiff EOS)"]; 

(* Full matter Lagrangian density (internal expression used for derivations) *)
LpsiInternal = (I*hbar/2) (psibF*Dtcov[psiF] - psiF*DtcovBar[psibF])
  - (hbar^2/(2 m)) Sum[
      DiBar[psibF, coords4[[k]], Avec[[k]]]*Di[psiF, coords4[[k]], Avec[[k]]],
      {k, 1, 4}
    ]
  - Vconf*rhoF
  - (K/4) rhoF^5;

(* Paper-visible compact L_psi (time-symplectic term minus Hamiltonian density) *)
LpsiPaper = HoldForm[
  L_psi ==
    (I*hbar/2) (psib Subscript[D, t][psi] - psi Subscript[D, t][psib])
    - (hbar^2/(2 m)) ( Subscript[D, x][psib] Subscript[D, x][psi]
                     + Subscript[D, y][psib] Subscript[D, y][psi]
                     + Subscript[D, z][psib] Subscript[D, z][psi]
                     + Subscript[D, w][psib] Subscript[D, w][psi] )
    - Vconf rho
    - (K/4) rho^5
];

EqPrint["Matter Lagrangian density L_psi (paper form)", LpsiPaper];

(* Paper GNLS (stiff water): i hbar D_t psi = [-(hbar^2/2m) D_i D_i + Vconf + (5K/4) rho^4] psi *)
GNLSPaper = HoldForm[
  I*hbar*Subscript[D, t][psi] ==
    (-(hbar^2/(2 m)) ( Subscript[D, x][Subscript[D, x][psi]]
                     + Subscript[D, y][Subscript[D, y][psi]]
                     + Subscript[D, z][Subscript[D, z][psi]]
                     + Subscript[D, w][Subscript[D, w][psi]] )
     + Vconf + (5 K/4) rho^4) psi
];
EqPrint["GNLS (paper form)", GNLSPaper];

(* Optional: show the internal Lagrangian in appendix (can be large). *)
If[$PrintAppendix,
  EqPrint["(Appendix) L_psi internal (expanded covariant derivative form)", HoldForm[LpsiInternal] -> LpsiInternal];
];

(* ----------------------------- *)
(* 5) Continuity identity        *)
(* ----------------------------- *)

Section["PAPER FORM: Continuity identity"]; 

ContinuityPaper = HoldForm[Subscript[del, t][rho] + Subscript[del, i][rho v[i]] == 0];
EqPrint["Continuity (paper form)", ContinuityPaper];

(* Internal target form in explicit components (appendix) *)
If[$PrintAppendix,
  EqPrint["(Appendix) 4D continuity (explicit components)",
    HoldForm[D[rhoF, t] + D[Jx, x] + D[Jy, y] + D[Jz, z] + D[Jw, w] == 0]
  ];
];

(* ----------------------------- *)
(* 6) Madelung + Euler           *)
(* ----------------------------- *)

Section["PAPER FORM: Madelung variables, Q, Euler equation, vorticity identity"]; 

(* Hydrodynamic variables for substitution checks (internal) *)
rhoMField   = rhoM[t, x, y, z, w];
thetaMField = thetaM[t, x, y, z, w];

psiM  = Sqrt[rhoMField] * Exp[I*thetaMField];
psibM = Sqrt[rhoMField] * Exp[-I*thetaMField];

vComp[var_, Acomp_] := (hbar/m) (D[thetaMField, var] - (qstar/hbar) Acomp);
vx = vComp[x, AxF]; vy = vComp[y, AyF]; vz = vComp[z, AzF]; vw = vComp[w, AwF];

EqBlock["Madelung definitions (paper)", {
  HoldForm[psi == Sqrt[rho] Exp[I theta]],
  HoldForm[v[i] == (hbar/m) (Subscript[del, i][theta] - (qstar/hbar) Subscript[A, i])],
  HoldForm[h[rho] == (5 K/4) rho^4]
}];

(* Quantum potential: keep canonical compact form in the paper bundle *)
Qcompact = HoldForm[Q == -(hbar^2/(2 m)) Laplacian[Sqrt[rho], {x, y, z, w}]/Sqrt[rho]];
EqPrint["Quantum potential Q (canonical compact form)", Qcompact];

(* Euler-like equation: compact index form (no component expansion) *)
EulerPaper = HoldForm[
  m (Subscript[del, t] + v[j] Subscript[del, j])[v[i]] ==
    qstar (E[i] + v[j] B[i, j]) - Subscript[del, i][Vconf + h[rho] + Q]
];
EqPrint["Euler-like equation (compact index form)", EulerPaper];

EqBlock["EM-like fields (paper definitions)", {
  HoldForm[E[i] == -Subscript[del, t][Subscript[A, i]] - Subscript[del, i][Subscript[A, 0]]],
  HoldForm[B[i, j] == Subscript[del, i][Subscript[A, j]] - Subscript[del, j][Subscript[A, i]]]
}];

VortIdentityPaper = HoldForm[
  Omega[i, j] == Subscript[del, i][v[j]] - Subscript[del, j][v[i]] == -(qstar/m) B[i, j]
];
EqPrint["Vorticity <-> gauge identity (paper form)", VortIdentityPaper];

EqBlock["Fluxoid quantization (paper form)", {
  HoldForm[LoopInt[Subscript[del, i][theta] - (qstar/hbar) Subscript[A, i], Cloop] == 2*Pi*nWind],
  HoldForm[GammaC == nWind*hPlanck/m - (qstar/m)*PhiC],
  HoldForm[PhiC == 0 -> GammaC == nWind*hPlanck/m]
}];

(* Appendix: explicit vector/component form (kept off by default) *)
If[$PrintAppendix,
  (* Madelung identity: J_i = rho v_i. For referee readability we print the residuals
     J_i - rho v_i (these should simplify to 0). *)
  JMadelungResiduals[] := Module[{JxM, JyM, JzM, JwM},
    JxM = SmartSimplify[Jcomp[x, AxF] /. {psiF -> psiM, psibF -> psibM}];
    JyM = SmartSimplify[Jcomp[y, AyF] /. {psiF -> psiM, psibF -> psibM}];
    JzM = SmartSimplify[Jcomp[z, AzF] /. {psiF -> psiM, psibF -> psibM}];
    JwM = SmartSimplify[Jcomp[w, AwF] /. {psiF -> psiM, psibF -> psibM}];
    {
      HoldForm[Jx - rhoMField*vx] -> SmartSimplify[JxM - rhoMField*vx],
      HoldForm[Jy - rhoMField*vy] -> SmartSimplify[JyM - rhoMField*vy],
      HoldForm[Jz - rhoMField*vz] -> SmartSimplify[JzM - rhoMField*vz],
      HoldForm[Jw - rhoMField*vw] -> SmartSimplify[JwM - rhoMField*vw]
    }
  ];
  EqBlock["(Appendix) Madelung identity residuals (should be 0)", JMadelungResiduals[]];

  (* Expanded Q (component form) *)
  Qexpanded = SmartSimplify[
    -(hbar^2/(2 m)) * (D[Sqrt[rhoMField], x, x] + D[Sqrt[rhoMField], y, y] + D[Sqrt[rhoMField], z, z] + D[Sqrt[rhoMField], w, w]) / Sqrt[rhoMField]
  ];
  EqPrint["(Appendix) Q expanded in components (rhoM)", HoldForm[Q] -> Qexpanded];
];

(* ----------------------------- *)
(* 7) EM sector: compact + checks*)
(* ----------------------------- *)

Section["PAPER FORM: EM sector (localized Maxwell + gauge fixing + jellium)"]; 

(* Internal Maxwell objects (used if you call MaxwellEOM[] in the appendix) *)
A5 = {A0F, AxF, AyF, AzF, AwF};

(* Localization profile (frozen family): Z(w) = exp(-w^2/lambdaConf^2) *)
Zfac = Exp[-w^2/lambdaConf^2];

Fmn = Table[
  D[A5[[nuIdx]], coords5[[muIdx]]] - D[A5[[muIdx]], coords5[[nuIdx]]],
  {muIdx, 1, 5}, {nuIdx, 1, 5}
];
Fup = eta . Fmn . eta;

(* Contraction: (1/2) sum_{MN} F_{MN} F^{MN} so that -(Z/(2 mu0))*Fcontract == -(Z/(4 mu0)) F_{MN}F^{MN} *)
Fcontract = (1/2) Sum[Fmn[[muIdx, nuIdx]]*Fup[[muIdx, nuIdx]], {muIdx, 1, 5}, {nuIdx, 1, 5}];

divA = Sum[eta[[muIdx, muIdx]]*D[A5[[muIdx]], coords5[[muIdx]]], {muIdx, 1, 5}];

(* Charge bookkeeping: matter current carries the fixed topological sign, explicit A.J terms are background only. *)
J0psi = qstar*rhoF;
JpsiVec = qstar*{Jx, Jy, Jz, Jw};
Jpsi5 = Join[{J0psi}, JpsiVec];

J0ext = -qstar*rhoBg;
Jext5 = Join[{J0ext}, {0, 0, 0, 0}];
Jtot5 = Jpsi5 + Jext5;
J0tot = Jtot5[[1]];

LemInternal = -(Zfac/(2 mu0)) * Fcontract - (1/(2 xi mu0)) * divA^2 - Sum[Jext5[[muIdx]]*A5[[muIdx]], {muIdx, 1, 5}];

(* Paper L_EM (compact invariant form) *)
M = Symbol["M"]; Nn = Symbol["Nn"]; MNn = Symbol["MNn"]; mu = Symbol["mu"]; muNn = Symbol["muNn"]; (* symbolic indices for printing only; avoid built-in N[] *)
LemPaper = HoldForm[
  L_EM ==
    -(Z[w]/(4 mu0)) (Subscript[F, MNn] * Superscript[F, MNn])
    - (1/(2 xi mu0)) (divA)^2
    - Subscript[A, M] * Superscript[Jext, M]
];

EqBlock["Charge bookkeeping + jellium (paper)", {
  HoldForm[qstar == etaQ*eStar],
  HoldForm[Superscript[Jpsi, 0] == qstar*rho],
  HoldForm[Superscript[Jext, 0] == -qstar*rhoBg],
  HoldForm[Superscript[Jtot, 0] == qstar*(rho - rhoBg)],
  HoldForm[Z[w] == Exp[-w^2/lambdaConf^2]],
  HoldForm[Subscript[F, MNn] == Subscript[del, M][Subscript[A, Nn]] - Subscript[del, Nn][Subscript[A, M]]],
  HoldForm[divA == Subscript[del, M][Superscript[A, M]]]
}];

EqPrint["L_EM (paper form; compact invariant expression)", LemPaper];

(* Maxwell EOM (paper form) *)
MaxwellPaper = HoldForm[
  Subscript[del, M][ Z[w] * Superscript[F, MNn] ] + (1/xi) * Superscript[del, Nn][divA] == mu0 * (Superscript[Jpsi, Nn] + Superscript[Jext, Nn])
];
EqPrint["Maxwell EOM (paper form)", MaxwellPaper];

(* Appendix: component equations from the internal Lagrangian, if desired *)
If[$PrintAppendix,
  EqPrint["(Appendix) L_EM internal (expanded)", HoldForm[LemInternal] -> LemInternal];

  MaxwellPaperEOMComponents = Table[
    With[{nn = n},
      HoldForm[
        Sum[D[Zfac*Fup[[muIdx, nn]], coords5[[muIdx]]], {muIdx, 1, 5}] + (1/xi)*eta[[nn, nn]]*D[divA, coords5[[nn]]] == mu0*Jtot5[[nn]]
      ]
    ],
    {n, 1, 5}
  ];
  EqBlock["(Appendix) Maxwell EOM components (compact sums)", MaxwellPaperEOMComponents];
];



(* ----------------------------- *)
(* 7b) EM localization reduction *)
(* ----------------------------- *)

Section["PAPER FORM: EM localization reduction (controlled zero-mode)"];

(* Normalization of the frozen localization profile Z(w) = exp(-w^2/lambdaConf^2). *)
Zint = SmartSimplify[
  Integrate[Exp[-w^2/lambdaConf^2], {w, -Infinity, Infinity}, Assumptions -> (lambdaConf > 0)]
];

EqBlock["Localization normalization", {
  HoldForm[Zint == Integrate[Z[w], {w, -Infinity, Infinity}]],
  HoldForm[Zint] -> Zint
}];

(* Controlled brane-dominant / zero-w-mode assumptions used for the paper reduction:
   - A_w = 0
   - A_\[Nu] has negligible w-dependence on the support of Z (so F_{w\[Nu]} \approx 0)
   - Work in Lorenz gauge on-shell (divA = 0), so the gauge-fixing term is inactive.
*)
EqBlock["Zero-mode / brane-dominant assumptions (paper)", {
  HoldForm[Subscript[A, w] == 0],
  HoldForm[Subscript[del, w][Subscript[A, Nn]] == 0],
  HoldForm[divA == 0]
}];

(* Integrated Maxwell equation on the brane (paper form).

   Starting point (Lorenz gauge):
     \partial_M ( Z F^{M N} ) = mu0 J^N.

   For N restricted to brane directions {0,x,y,z} and with F^{w N} = 0,
     Z(w) \partial_\[Mu] F^{\[Mu] N} = mu0 J^N.

   Integrating over w yields the effective brane equation
     \partial_\[Mu] F^{\[Mu] N} = mu0eff * J_eff^N,
   where
     mu0eff = mu0 / Zint,
     J_eff^N = \int_{-\infty}^{\infty} J^N dw.
*)
EqBlock["Effective brane Maxwell equation (controlled reduction)", {
  HoldForm[Jeff[Nn] == Integrate[Superscript[Jtot, Nn], {w, -Infinity, Infinity}]],
  HoldForm[mu0eff == mu0/Zint],
  HoldForm[Subscript[del, mu][Superscript[F, muNn]] == mu0eff*Jeff[Nn]]
}];

EqBlock["Canonical normalization / effective charge (paper)", {
  HoldForm[Acanon[mu] == Sqrt[Zint]*Subscript[A, mu]],
  HoldForm[qeff == qstar/Sqrt[Zint]],
  HoldForm[eEff == eStar/Sqrt[Zint]],
  HoldForm[qeff == etaQ*eEff],
  HoldForm[eEff == eStar/Sqrt[lambdaConf*Sqrt[Pi]]]
}];

(* Action-level reduction (referee-proof): show the same scaling directly at the Lagrangian-density level.

   In Lorenz gauge on-shell (divA=0), the gauge-fixing term drops out and we can reduce the physical kinetic
   term cleanly.

   Key invariant split:
      F_{MN}F^{MN} = F_{mu nu}F^{mu nu} + 2 F_{w mu}F^{w mu}.
   Under the controlled zero-mode ansatz (A_w=0, \partial_w A_\mu = 0), we have F_{w\mu}=0.

   Integrating the kinetic prefactor over w then gives a single prefactor Zint = \int Z(w) dw.
*)
EqBlock["Action-level localization reduction (paper)", {
  HoldForm[divA == 0],
  HoldForm[Subscript[F, MNn] Superscript[F, MNn] == Subscript[F, munu] Superscript[F, munu] + 2 Subscript[F, wmu] Superscript[F, wmu]],
  HoldForm[Subscript[F, wmu] == 0],
  HoldForm[JextEff[mu] == Integrate[Superscript[Jext, mu], {w, -Infinity, Infinity}]],
  HoldForm[L_EM4 == -(Zint/(4*mu0))*(Subscript[F, munu] Superscript[F, munu]) - Subscript[A, mu] JextEff[mu]],
  HoldForm[mu0eff == mu0/Zint]
}];

(* Electrostatic specialization (paper hook):
   A_i = 0, A_0 = Phi(x,y,z), \partial_t=0, \partial_w Phi=0  =>  3D Poisson-type equation.
*)
EqPrint["Electrostatic brane limit (hook)", HoldForm[Laplacian[PhiEM, {x, y, z}] == -mu0eff*Jeff[0]]];

If[$PrintAppendix,
  Section["APPENDIX: EM localization reduction sanity check (electrostatic zero-mode)"];

  (* Define a brane potential Phi(t,x,y,z) and substitute into the internal Maxwell residual. *)
  PhiEMF = PhiEM[t, x, y, z];

  subElectro = {A0F -> PhiEMF, AxF -> 0, AyF -> 0, AzF -> 0, AwF -> 0};

  electroResidual = SmartSimplify[
    (
      Sum[D[Zfac*Fup[[muIdx, 1]], coords5[[muIdx]]], {muIdx, 1, 5}] +
      (1/xi)*eta[[1, 1]]*D[divA, coords5[[1]]] -
      mu0*J0tot
    ) /. subElectro /. {D[PhiEMF, t] -> 0, D[PhiEMF, t, t] -> 0}
  ];

  EqPrint["Time-component Maxwell residual under electrostatic zero-mode", HoldForm[resid0] -> electroResidual];

  (* w-integrated form for Phi independent of w (PhiF has no w argument, so this is automatic). *)
  EqPrint["w-integrated form (PhiEM independent of w)",
    HoldForm[Zint*Laplacian[PhiEM, {x, y, z}] == -mu0*Integrate[J0tot, {w, -Infinity, Infinity}]]
  ];
];

(* ----------------------------- *)
(* 8) Projection + leakage       *)
(* ----------------------------- *)

Section["PAPER FORM: Brane projection and leakage source term"]; 

(* Keep W(w) abstract in the paper bundle; show HO ground-state W0 in appendix. *)
rhoBraneDef = HoldForm[rho_brane == Integrate[W[w] * rho, {w, -Infinity, Infinity}]];
JBraneDef = HoldForm[J_brane == Integrate[W[w] * J_xyz, {w, -Infinity, Infinity}]];
EqBlock["Weighted projection (paper)", {
  rhoBraneDef,
  HoldForm[J_xyz == {Jx, Jy, Jz}],
  JBraneDef,
  HoldForm[Integrate[W[w], {w, -Infinity, Infinity}] == 1]
}];

(* Projected continuity identity with leakage/boundary term *)
LeakageSourcePaper = HoldForm[
  S_leak == -Limit[W[w] * Jw, w -> Infinity] + Limit[W[w] * Jw, w -> -Infinity] + Integrate[D[W[w], w] * Jw, {w, -Infinity, Infinity}]
];
ProjContinuityPaper = HoldForm[Subscript[del, t][rho_brane] + Div[J_brane, {x, y, z}] == S_leak];

EqBlock["Projected continuity (paper identity)", {ProjContinuityPaper, LeakageSourcePaper}];

(* Sanity simplification: for fast-decaying fields/basis, boundary term vanishes. *)
EqBlock["Fast-decay simplification (paper sanity check)", {
  HoldForm[Limit[W[w] * Jw, w -> Infinity] == 0],
  HoldForm[Limit[W[w] * Jw, w -> -Infinity] == 0],
  HoldForm[S_leak == Integrate[D[W[w], w] * Jw, {w, -Infinity, Infinity}]]
}];

(* Appendix: show that HO ground-state weight is normalized and consistent with the frozen harmonic trap. *)
If[$PrintAppendix,
  ellW = Symbol["ellW"]; (* HO length scale *)
  chi0[ww_] := (1/(Pi^(1/4)*Sqrt[ellW]))*Exp[-ww^2/(2 ellW^2)];
  W0[ww_] := chi0[ww]^2;

  EqBlock["(Appendix) Harmonic trap ground-state weight W0(w)", {
    HoldForm[W[w] == W0[w]],
    HoldForm[W0[w]] -> W0[w],
    HoldForm[Integrate[W0[w], {w, -Infinity, Infinity}]] -> SmartSimplify[Integrate[W0[w], {w, -Infinity, Infinity}, Assumptions -> (ellW > 0)]]
  }];
];

(* ----------------------------- *)
(* 9) Geometry force ledger      *)
(* ----------------------------- *)

Section["PAPER FORM: Geometry force ledger and wall law"]; 

(* Explicit (paper-friendly) tube measures from the master sheet *)
Vgeom[a_, L_] := (4 Pi/3) a^3 L;
Ageom[a_, L_] := 4 Pi a^2 L + (8 Pi/3) a^3;

EgeomPaper = HoldForm[E_geom[a, L] == Pvac Vgeom[a, L] + sigma Ageom[a, L] + kappaB*a^2/L];
EqPrint["Geometry energy (paper form)", EgeomPaper];

GeometryDerivsPaper = {
  HoldForm[D[kappaB*a^2/L, a] == 2*kappaB*a/L],
  HoldForm[D[kappaB*a^2/L, L] == -(kappaB*a^2/L^2)]
};
EqBlock["Geometry derivatives (paper)", GeometryDerivsPaper];

ForceLedgerPaper = {
  HoldForm[F_a[t] == -D[H_tot, a]],
  HoldForm[F_L[t] == -D[H_tot, L]],
  HoldForm[F_a == F_a_psi + F_a_EM + F_a_geom + F_a_other],
  HoldForm[F_L == F_L_psi + F_L_EM + F_L_geom + F_L_other]
};
EqBlock["Force ledger identity (ASCII-safe)", ForceLedgerPaper];

FluidForceKernelPaper = {
  HoldForm[F_a_psi == -Integrate[rho * D[Vconf, a], {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity}, {w, -Infinity, Infinity}]],
  HoldForm[F_L_psi == -Integrate[rho * D[Vconf, L], {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity}, {w, -Infinity, Infinity}]]
};
EqBlock["Matter contribution (Hellmann-Feynman form)", FluidForceKernelPaper];

WallLawPaper = {
  HoldForm[Ma a''[t] + Ca a'[t] == F_a[t]],
  HoldForm[ML L''[t] + CL L'[t] == F_L[t]]
};
EqBlock["Dynamic wall law (paper form)", WallLawPaper];

(* ----------------------------- *)
(* 10) Optional notes for Poisson*)
(* ----------------------------- *)



(* ----------------------------- *)
(* 9b) Brane kinematics / sectors*)
(* ----------------------------- *)

Section["PAPER FORM: Brane sector decomposition (Helmholtz + Poisson hook)"];

EqPrint["Brane velocity (operational)", HoldForm[v_brane == J_brane/rho_brane]];

EqBlock["Brane Helmholtz decomposition (symbolic)", {
  HoldForm[v_brane == Grad[Phi, {x, y, z}] + vT],
  HoldForm[Div[vT, {x, y, z}] == 0]
}];

(* Exact identity: longitudinal potential obeys a Poisson-like equation sourced by leakage and compressibility.

   Start from projected continuity:
      del_t rho_brane + Div(J_brane) = S_leak,
   and J_brane = rho_brane v_brane with v_brane = Grad(Phi)+vT, Div(vT)=0.

   Expanding Div(rho_brane Grad(Phi)) yields:
      rho_brane Laplacian(Phi) = S_leak - del_t rho_brane - Grad(rho_brane)·(Grad(Phi)+vT).
*)
EqBlock["Exact longitudinal identity (from projected continuity)", {
  HoldForm[J_brane == rho_brane*v_brane],
  HoldForm[Subscript[del, t][rho_brane] + Div[rho_brane*v_brane, {x, y, z}] == S_leak],
  HoldForm[rho_brane*Laplacian[Phi, {x, y, z}] == S_leak - Subscript[del, t][rho_brane] - Grad[rho_brane, {x, y, z}] . (Grad[Phi, {x, y, z}] + vT)]
}];

EqPrint["Poisson candidate (regime statement)",
  HoldForm[Laplacian[Phi, {x, y, z}] == 4*Pi*Geff*rho3D]
];

Section["PAPER NOTE: Poisson sector is a regime/limit (not forced)"]; 

Print["Regime statement (paper text):"]; 
Print[
  "In a quasi-static, longitudinal-dominant brane regime with small leakage and subdominant Q, " <>
  "the longitudinal brane potential (from a Helmholtz split) is the natural Poisson candidate." 
];
Print[""]; 


(* Collect the main paper-form equations for convenient TeX export. *)
PaperEqns = Flatten[{
  (* Core matter sector *)
  LpsiPaper,
  GNLSPaper,
  ContinuityPaper,
  HoldForm[psi == Sqrt[rho] Exp[I theta]],
  HoldForm[v[i] == (hbar/m) (Subscript[del, i][theta] - (qstar/hbar) Subscript[A, i])],
  Qcompact,
  EulerPaper,
  VortIdentityPaper,
  HoldForm[LoopInt[Subscript[del, i][theta] - (qstar/hbar) Subscript[A, i], Cloop] == 2*Pi*nWind],
  HoldForm[GammaC == nWind*hPlanck/m - (qstar/m)*PhiC],

  (* EM sector *)
  LemPaper,
  HoldForm[qstar == etaQ*eStar],
  HoldForm[Superscript[Jtot, 0] == qstar (rho - rhoBg)],
  HoldForm[Z[w] == Exp[-w^2/lambdaConf^2]],
  HoldForm[Subscript[F, MNn] == Subscript[del, M][Subscript[A, Nn]] - Subscript[del, Nn][Subscript[A, M]]],
  HoldForm[divA == Subscript[del, M][Superscript[A, M]]],
  MaxwellPaper,

  (* EM localization reduction (controlled) *)
  HoldForm[Zint == Integrate[Z[w], {w, -Infinity, Infinity}]],
  HoldForm[mu0eff == mu0/Zint],
  HoldForm[Jeff[Nn] == Integrate[Superscript[Jtot, Nn], {w, -Infinity, Infinity}]],
  HoldForm[Subscript[del, mu][Superscript[F, muNn]] == mu0eff*Jeff[Nn]],
  HoldForm[qeff == qstar/Sqrt[Zint]],
  HoldForm[eEff == eStar/Sqrt[Zint]],
  HoldForm[JextEff[mu] == Integrate[Superscript[Jext, mu], {w, -Infinity, Infinity}]],
  HoldForm[L_EM4 == -(Zint/(4*mu0))*(Subscript[F, munu] Superscript[F, munu]) - Subscript[A, mu] JextEff[mu]],

  (* Brane sector decomposition / Poisson hook *)
  HoldForm[v_brane == J_brane/rho_brane],
  HoldForm[v_brane == Grad[Phi, {x, y, z}] + vT],
  HoldForm[Div[vT, {x, y, z}] == 0],
  HoldForm[rho_brane*Laplacian[Phi, {x, y, z}] == S_leak - Subscript[del, t][rho_brane] - Grad[rho_brane, {x, y, z}] . (Grad[Phi, {x, y, z}] + vT)],
  HoldForm[Laplacian[Phi, {x, y, z}] == 4*Pi*Geff*rho3D],

  (* Projection/leakage *)
  rhoBraneDef,
  HoldForm[J_xyz == {Jx, Jy, Jz}],
  JBraneDef,
  ProjContinuityPaper,
  LeakageSourcePaper,

  (* Geometry ledger *)
  EgeomPaper,
  ForceLedgerPaper,
  FluidForceKernelPaper,
  WallLawPaper
}];


(* ----------------------------- *)
(* 11) Extra appendices (optional)*)
(* ----------------------------- *)

If[$PrintAppendix,
  Section["APPENDIX: Minimal w-mode truncation and Gaussian scaling (optional)"];

  (* Two-mode HO truncation: psi ~ psi0 chi0 + eps psi1 chi1 *)
  eps = Symbol["eps"];
  psi0F  = psi0[t, x, y, z];
  psib0F = psib0[t, x, y, z];
  psi1F  = psi1[t, x, y, z];
  psib1F = psib1[t, x, y, z];

  chi1[ww_] := (Sqrt[2]/(Pi^(1/4)*Sqrt[ellW]))*(ww/ellW)*Exp[-ww^2/(2 ellW^2)];

  psiAns  = psi0F*chi0[w] + eps*psi1F*chi1[w];
  psibAns = psib0F*chi0[w] + eps*psib1F*chi1[w];
  rhoAns = SmartSimplify[Expand[psiAns*psibAns]];

  rhoBraneAns = SmartSimplify[Integrate[W0[w]*rhoAns, {w, -Infinity, Infinity}, Assumptions -> (ellW > 0)]];

  EqBlock["Two-mode ansatz: rho and projected rho_brane", {
    HoldForm[psi] -> psiAns,
    HoldForm[rho] -> rhoAns,
    HoldForm[rho_brane] -> rhoBraneAns
  }];

  crossIntegral = SmartSimplify[Integrate[W0[w]*chi0[w]*chi1[w], {w, -Infinity, Infinity}, Assumptions -> (ellW > 0)]];
  EqPrint["Parity check (should be 0): Integrate[W0 chi0 chi1]", HoldForm[cross] -> crossIntegral];

  (* Leakage kernel structure (no gauge) for the 2-mode ansatz: J_w is mode-mixing. *)
  JwAnsNoGauge = SmartSimplify[(hbar/(2 I m))*(psibAns*D[psiAns, w] - psiAns*D[psibAns, w])];
  EqPrint["Two-mode ansatz: J_w structure (no gauge)", HoldForm[Jw] -> JwAnsNoGauge];

  (* (Optional) Projected leakage source for this symmetric two-mode ansatz.
     With W0 even and Jw even (here Jw ~ chi0^2), the integrand W0'(w) Jw is odd and integrates to 0.
     This is a useful sanity check: in the minimal symmetric truncation, the weighted brane continuity closes.
  *)
  SleakAns = SmartSimplify[Integrate[D[W0[w], w]*JwAnsNoGauge, {w, -Infinity, Infinity}, Assumptions -> (ellW > 0)]];
  EqPrint["(Appendix) Two-mode projected leakage source (symmetric W0)", HoldForm[S_leak] -> SleakAns];

  (* Symmetry-breaking example (useful for Paper 7 text):
       - In the symmetric baseline (even W0, even Jw), S_leak = ∫ W0'(w) Jw(w) dw vanishes by parity.
       - If the effective brane weight is shifted off-center, W(w)=W0(w-w0), then parity is broken and
         S_leak turns on at O(w0).

     This is a clean demonstration that leakage is not “forced”; it is a controllable sector.
  *)
  w0 = Symbol["w0"]; (* brane offset in w; can be positive or negative *)
  Wshift[ww_] := W0[ww - w0];
  SleakShift = SmartSimplify[
    Integrate[D[Wshift[w], w]*JwAnsNoGauge, {w, -Infinity, Infinity},
      Assumptions -> (ellW > 0 && Element[w0, Reals])
    ]
  ];
  EqBlock["(Appendix) Leakage turns on under parity breaking (shifted weight)", {
    HoldForm[W[w] == W0[w - w0]],
    HoldForm[S_leak_shift] -> SleakShift,
    HoldForm[S_leak_shift_series] -> Normal[Series[SleakShift, {w0, 0, 1}]]
  }];

  (* (Optional) Controlled separable-trap projection of the stiff-EOS GNLS onto the brane-dominant mode.

     Controlled simplification (used only for this appendix derivation):
       V_conf(x,y,z,w) ~ V3(t,x,y,z) + (hbar^2/(2 m ellW^4)) w^2,
     i.e. a separable HO trap in w with ground-state width ellW.

     This yields an effective 3D GNLS for psi0 with an explicit stiff-EOS overlap constant.
  *)
  V3F = V3[t, x, y, z];
  VwRef = (hbar^2*w^2)/(2*m*ellW^4);

  C10 = SmartSimplify[Integrate[chi0[w]^10, {w, -Infinity, Infinity}, Assumptions -> (ellW > 0)]];
  E0w = SmartSimplify[Integrate[chi0[w]*(-hbar^2/(2*m)*D[chi0[w], {w, 2}] + VwRef*chi0[w]), {w, -Infinity, Infinity}, Assumptions -> (ellW > 0)]];
  E1w = SmartSimplify[Integrate[chi1[w]*(-hbar^2/(2*m)*D[chi1[w], {w, 2}] + VwRef*chi1[w]), {w, -Infinity, Infinity}, Assumptions -> (ellW > 0)]];

  EqBlock["(Appendix) Separable HO-w trap constants", {
    HoldForm[C10 == Integrate[chi0[w]^10, {w, -Infinity, Infinity}]],
    HoldForm[C10] -> C10,
    HoldForm[E0w] -> E0w,
    HoldForm[E1w] -> E1w
  }];

  rho03 = psi0F*psib0F;
  GNLS0proj = HoldForm[
    I*hbar*D[psi0F, t] == (-(hbar^2/(2*m))*Laplacian[psi0F, {x, y, z}] + V3F + E0w + (5*K/4)*C10*rho03^4)*psi0F
  ];
  EqPrint["(Appendix) Projected brane-mode GNLS (psi0; separable HO w-trap)", GNLS0proj];

  (* Gaussian scaling appendix (kept as in the prior harness) *)
  rhoGaussian = (Symbol["M"]/ (Pi^2*a^3*L)) * Exp[-(x^2 + y^2 + z^2)/a^2] * Exp[-w^2/L^2];

  QGaussian = SmartSimplify[-(hbar^2/(2 m))*(D[Sqrt[rhoGaussian], x, x] + D[Sqrt[rhoGaussian], y, y] + D[Sqrt[rhoGaussian], z, z] + D[Sqrt[rhoGaussian], w, w]) / Sqrt[rhoGaussian]];
  EqPrint["Gaussian quantum potential Q(x,y,z,w)", HoldForm[Qgauss] -> QGaussian];

  UInternalGaussian = SmartSimplify[(K/4) * Integrate[rhoGaussian^5, {x, -Infinity, Infinity}, {y, -Infinity, Infinity}, {z, -Infinity, Infinity}, {w, -Infinity, Infinity}, Assumptions -> (a > 0 && L > 0)]];
  EqPrint["Gaussian internal energy (stiff EOS)", HoldForm[Uint] -> UInternalGaussian];

  EqBlock["Restoring forces from Uint(a,L)", {
    HoldForm[F_a_from_U] -> SmartSimplify[-D[UInternalGaussian, a]],
    HoldForm[F_L_from_U] -> SmartSimplify[-D[UInternalGaussian, L]]
  }];
];

Section["END: Paper-polished master derivation harness loaded"]; 

(*
NOTE on LaTeX export:
  - To export the main paper equations, collect them into a list (e.g. PaperEqns)
    and call ExportPaperTeX["paper_eqns.tex", PaperEqns, replacements].

Example replacements (optional):
  replacements = {
    "\\text{psib}" -> "\\bar{\\psi}",
    "psib" -> "\\bar{\\psi}",
    "psi" -> "\\psi",
    "rho" -> "\\rho"
  };
*)
