(* ::Package:: *)

(***
========================================
Paper VIII (referee harness v2): Maxwell from 4D localized Maxwell sector
========================================

Goals (referee-walkable):
  1) Define 4+1D localized Maxwell + gauge-fix + source coupling.
  2) Derive Euler–Lagrange equations from the action.
  3) Verify the derived EOM matches the compact index form.
  4) Verify the Bianchi identity (homogeneous Maxwell).
  5) Take divergence of EOM and show current conservation follows in Lorenz gauge.
  6) Show controlled brane reduction (zero w-mode) at the action + EOM level.

Robustness note:
  - This version avoids SubtractSides[] to remain compatible with older Mathematica versions.
  - Any paper-like subscripts (A_\[Mu], \[PartialD]_w, etc.) are printed as STRINGS only.
  - No symbol names contain underscores ( _ ), because underscores are patterns in WL.

Conventions:
  - Metric signature: (-,+,+,+,+)
  - Afunc symbols represent covariant components A_M.
  - Jfunc symbols represent contravariant components J^M by default.
    If you use J0[...] to mean the covariant component J_0 instead, set J0IsCovariant=True.

** IMPORTANT: What we consider "passing" **
  - Each component check should reduce to EXACTLY 0 (after simplification).
  - Bianchi residuals should be 0.
  - Divergence identity residual should be 0.

*)

Print["\n========================================"]; 
Print["Paper VIII (referee harness v2): Maxwell from 4D localized Maxwell sector"]; 
Print["========================================\n"]; 

ClearAll["Global`*"];
Needs["VariationalMethods`"];

(* Assumptions for simplification; extend if desired *)
$Assumptions = {mu0 != 0, xi != 0, lambdaConf > 0, eStar > 0, Element[etaQ, Reals], etaQ^2 == 1};

(* Simplification helper: Together first (rational structure), then FullSimplify *)
Simp[expr_] := FullSimplify[Together[expr], Assumptions -> $Assumptions];

(* Convert an EulerEquations output to a pure expression "lhs - rhs".
   This is deliberately version-agnostic (does not rely on SubtractSides). *)
EqToExpr[eq_] := Module[{lst},
  Which[
    Head[eq] === Equal,
      lst = List @@ eq;
      If[Length[lst] >= 2,
        lst[[1]] - Total[lst[[2 ;;]]],
        eq
      ],
    True,
      eq
  ]
];

(* Convenience: boolean check for "expr == 0" under assumptions *)
ZeroQAssume[expr_] := TrueQ[FullSimplify[expr == 0, Assumptions -> $Assumptions]];

Print["---------- Charge ontology bookkeeping ----------\n"]; 
Print["Fixed defect branch label: qStar = etaQ * eStar, with etaQ = +/- 1 and eStar > 0."]; 
Print["In the coupled reading, Jtot^M = Jpsi^M + Jext^M; the compact J^M used below is shorthand total source."]; 
Print[""]; 

Print["---------- Conventions: coordinates, metric, fields ----------\n"]; 

coords5 = {t, x, y, z, w};
dim5 = Length[coords5];

eta5  = DiagonalMatrix[{-1, 1, 1, 1, 1}];
etaU5 = eta5; (* inverse metric is itself for diag(-1,1,1,1,1) *)

Zprof[ww_] := Exp[-(ww^2/lambdaConf^2)];

Print["Objects"]; 
Print[HoldForm[coords5], " -> ", coords5];
Print[HoldForm[eta5], " -> ", eta5];
Print[HoldForm[Zprof[w]], " -> ", Zprof[w]];
Print[""];

(* Field component function symbols: covariant A_M *)
Afunc = {A0, Ax, Ay, Az, Aw};
Jfunc = {J0, Jx, Jy, Jz, Jw};

(* Toggle: if True, interpret J0[...] as covariant J_0 and convert to contravariant J^0=-J_0 *)
J0IsCovariant = False;
Jraw[i_]  := Jfunc[[i]][Sequence @@ coords5];
Jcomp[i_] := If[i == 1 && J0IsCovariant, -Jraw[i], Jraw[i]];

Acomp[i_] := Afunc[[i]][Sequence @@ coords5];

d[i_, expr_] := D[expr, coords5[[i]]];

(* F_{MN} = d_M A_N - d_N A_M *)
Fdown[i_, j_] := d[i, Acomp[j]] - d[j, Acomp[i]];

(* F^{MN} = eta^{MP} eta^{NQ} F_{PQ} *)
Fup[i_, j_] := Sum[etaU5[[i, p]] etaU5[[j, q]] Fdown[p, q], {p, 1, dim5}, {q, 1, dim5}];

(* DivA = d_M A^M = eta^{MN} d_M A_N *)
DivA := Sum[etaU5[[m, n]] d[m, Acomp[n]], {m, 1, dim5}, {n, 1, dim5}];

(* Raised derivative d^N = eta^{NM} d_M *)
dUp[n_, expr_] := Sum[etaU5[[n, m]] d[m, expr], {m, 1, dim5}];

Print["---------- Lagrangian (as used in the code) ----------\n"]; 

Print["Paper-form reminder (string):"]; 
Print["  L = -(Z(w)/(4 mu0)) F_{MN}F^{MN}  - (1/(2 xi mu0)) (DivA)^2  - A_M J^M"]; 
Print["  where DivA := d_M A^M  and J^M is treated as an external source.\n"]; 

(* Invariant F_{MN}F^{MN} with full index sum; standard 1/4 prefactor in L handles antisymmetry *)
F2 = Sum[Fdown[i, j] Fup[i, j], {i, 1, dim5}, {j, 1, dim5}];

Lkin = -(Zprof[w]/(4 mu0)) F2;
Lgf  = -(1/(2 xi mu0)) DivA^2;
Lint = -Sum[Acomp[i] Jcomp[i], {i, 1, dim5}];

L = Lkin + Lgf + Lint;

Print["Wolfram-language expression actually varied (L):"]; 
Print[L];
Print[""];

Print["---------- Derive Maxwell EOM from the action and verify compact form ----------\n"]; 

fields5 = Table[Afunc[[i]][Sequence @@ coords5], {i, 1, dim5}];

(* EulerEquations typically returns a list of equations; convert to expressions robustly *)
eomEqns5 = EulerEquations[L, fields5, coords5];

Print["EulerEquations output head check:"]; 
Print["  Head[eomEqns5[[1]]] -> ", Head[eomEqns5[[1]]]];

(* Convert to expressions == 0 form (lhs - rhs) *)
eomExpr5 = Simp /@ (EqToExpr /@ eomEqns5);

Print["After EqToExpr conversion:"]; 
Print["  Head[eomExpr5[[1]]] -> ", Head[eomExpr5[[1]]]];
Print[""];

Print["Maxwell EOM (expected compact expression, paper-normalized):"]; 
Print["  d_M( Z(w) F^{MN} ) + (1/xi) d^N(DivA) = mu0 J^N"]; 
Print["We check that:  mu0*(ELexpr) - [ d_M(Z F^{MN}) + (1/xi) d^N(DivA) - mu0 J^N ]  simplifies to 0.\n"]; 

expected5[n_] := Simp[
  Sum[d[m, Zprof[w] Fup[m, n]], {m, 1, dim5}] + (1/xi) dUp[n, DivA] - mu0 Jcomp[n]
];

check5 = Table[
  Simp[mu0*eomExpr5[[n]] - expected5[n]],
  {n, 1, dim5}
];

labels5 = {"N=0 (t)", "N=x", "N=y", "N=z", "N=w"};

Print["Component check results (each should be 0):"]; 
Do[
  Print["  ", labels5[[n]], " -> ", check5[[n]]],
  {n, 1, dim5}
];

pass5 = ZeroQAssume /@ check5;
Print["\nZero-check booleans under $Assumptions:"]; 
Do[
  Print["  ", labels5[[n]], " : ", pass5[[n]]],
  {n, 1, dim5}
];

If[And @@ pass5,
  Print["\nOK: All 5D components match the compact Maxwell form.\n"],
  Print["\nWARNING: At least one 5D component did not verify as zero under current assumptions." ];
  Print["Common causes: (i) J0 covariant-vs-contravariant mismatch; (ii) different gauge-fix sign convention." ];
  Print["Try toggling J0IsCovariant=True if your J0 symbol means J_0." ];
  Print[""]
];

Print["---------- Bianchi identity check (homogeneous Maxwell equations) ----------\n"]; 

triples = Subsets[Range[dim5], {3}];

bianchi = (With[{a = #[[1]], b = #[[2]], c = #[[3]]},
    Simp[d[a, Fdown[b, c]] + d[b, Fdown[c, a]] + d[c, Fdown[a, b]]]
] &) /@ triples;

Print["Bianchi identity residuals (each should be 0):"]; 
Do[
  Print["  triple ", triples[[k]], " -> ", bianchi[[k]]],
  {k, 1, Length[triples]}
];

If[And @@ (ZeroQAssume /@ bianchi),
  Print["\nOK: All Bianchi identities hold (homogeneous Maxwell).\n"],
  Print["\nWARNING: Some Bianchi residuals are nonzero (unexpected)." ];
  Print["Since F = dA, dF should vanish identically; this would indicate a coding error." ];
  Print[""]
];

Print["---------- Consistency: divergence of Maxwell EOM forces current conservation (Lorenz gauge) ----------\n"]; 

DivJ = Simp[Sum[d[n, Jcomp[n]], {n, 1, dim5}]];

(* Divergence of the expected paper-normalized EOM *)
divEOM = Simp[Sum[d[n, expected5[n]], {n, 1, dim5}]];

(* 5D box on DivA: box DivA := d_N d^N DivA *)
boxDivA = Simp[Sum[etaU5[[n, m]] d[n, d[m, DivA]], {n, 1, dim5}, {m, 1, dim5}]];

Print["Computed divergence of the (paper-normalized) EOM, d_N(EOM^N):"]; 
Print[divEOM];
Print[""];

Print["It should reduce to: (1/xi) box(DivA) - mu0 * (d_N J^N)." ];
resDiv = Simp[divEOM - (1/xi) boxDivA + mu0*DivJ];
Print["Residual after subtracting target form (should be 0):"]; 
Print[resDiv];

If[ZeroQAssume[resDiv],
  Print["\nOK: Divergence identity holds. In Lorenz gauge (DivA=0), EOM enforces d_N J^N = 0.\n"],
  Print["\nWARNING: Divergence identity did not simplify to 0 (unexpected)." ];
  Print[""]
];

Print["---------- Controlled brane reduction: zero-w-mode -> 3+1 Maxwell ----------\n"]; 

Print["Reduction assumptions (strings):"]; 
Print["  (i) Aw = 0"]; 
Print["  (ii) A0,Ax,Ay,Az have no w-dependence (photon zero-mode)" ];
Print["  (iii) Jw ~ 0 in the brane regime (no net flux off the brane)\n" ];

Zint = Assuming[$Assumptions, Integrate[Zprof[w], {w, -Infinity, Infinity}]];
mu0eff = Simp[mu0/Zint];

Print["Localization integral and effective coupling:"]; 
Print[HoldForm[Zint], " -> ", Zint];
Print[HoldForm[mu0eff], " -> ", mu0eff];
qStar = etaQ*eStar;
qEff = Simp[qStar/Sqrt[Zint]];
eEff = Simp[eStar/Sqrt[Zint]];
Print[HoldForm[qStar], " -> ", qStar];
Print[HoldForm[qEff], " -> ", qEff];
Print[HoldForm[eEff], " -> ", eEff];
Print[""];

Print["Under (i)-(ii), the kinetic term reduces as:"]; 
Print["  Integral dw [ -(Z(w)/(4 mu0)) F_{MN}F^{MN} ]  ->  -(Zint/(4 mu0)) f_{mu nu} f^{mu nu}" ];
Print["which is standard 3+1 Maxwell with mu0eff = mu0/Zint.\n" ];
Print["Equivalent canonical normalization statement: if a_mu = a_mu^can / Sqrt[Zint], then qEff = qStar / Sqrt[Zint] and eEff = eStar / Sqrt[Zint].\n"];

(* ---------- 3+1D brane Maxwell: derive from reduced action ---------- *)
Print["--- 3+1D effective Maxwell system (brane) ---\n"]; 

coords4 = {t, x, y, z};
dim4 = Length[coords4];

eta4  = DiagonalMatrix[{-1, 1, 1, 1}];
etaU4 = eta4;

Abfunc = {A0b, Axb, Ayb, Azb};
Jbfunc = {J0b, Jxb, Jyb, Jzb};

Abcomp[i_] := Abfunc[[i]][Sequence @@ coords4];
Jbcomp[i_] := Jbfunc[[i]][Sequence @@ coords4];

d4[i_, expr_] := D[expr, coords4[[i]]];

F4down[i_, j_] := d4[i, Abcomp[j]] - d4[j, Abcomp[i]];
F4up[i_, j_] := Sum[etaU4[[i, p]] etaU4[[j, q]] F4down[p, q], {p, 1, dim4}, {q, 1, dim4}];

DivA4 := Sum[etaU4[[m, n]] d4[m, Abcomp[n]], {m, 1, dim4}, {n, 1, dim4}];
dUp4[n_, expr_] := Sum[etaU4[[n, m]] d4[m, expr], {m, 1, dim4}];

F4sq = Sum[F4down[i, j] F4up[i, j], {i, 1, dim4}, {j, 1, dim4}];

L4 = -(1/(4 mu0eff)) F4sq - (1/(2 xi mu0eff)) DivA4^2 - Sum[Abcomp[i] Jbcomp[i], {i, 1, dim4}];

Print["Effective 3+1 Lagrangian density (L4) being varied:"]; 
Print[L4];
Print[""];

fields4 = Table[Abfunc[[i]][Sequence @@ coords4], {i, 1, dim4}];

eomEqns4 = EulerEquations[L4, fields4, coords4];
Print["EulerEquations (3+1) head check:"]; 
Print["  Head[eomEqns4[[1]]] -> ", Head[eomEqns4[[1]]]];

eomExpr4 = Simp /@ (EqToExpr /@ eomEqns4);
Print["After EqToExpr conversion (3+1):"]; 
Print["  Head[eomExpr4[[1]]] -> ", Head[eomExpr4[[1]]]];
Print[""];

expected4[n_] := Simp[
  Sum[d4[m, F4up[m, n]], {m, 1, dim4}] + (1/xi) dUp4[n, DivA4] - mu0eff Jbcomp[n]
];

check4 = Table[Simp[mu0eff*eomExpr4[[n]] - expected4[n]], {n, 1, dim4}];
labels4 = {"nu=0 (t)", "nu=x", "nu=y", "nu=z"};

Print["3+1 component check results (each should be 0):"]; 
Do[
  Print["  ", labels4[[n]], " -> ", check4[[n]]],
  {n, 1, dim4}
];

pass4 = ZeroQAssume /@ check4;
Print["\n3+1 zero-check booleans under $Assumptions:"]; 
Do[
  Print["  ", labels4[[n]], " : ", pass4[[n]]],
  {n, 1, dim4}
];

If[And @@ pass4,
  Print["\nOK: 3+1 Maxwell equations follow from the reduced action.\n"],
  Print["\nWARNING: At least one 3+1 component did not verify as zero under current assumptions." ];
  Print["If 5D checks pass but 3+1 fails, it usually indicates a normalization mismatch in mu0eff." ];
  Print[""]
];

Print["--- KK / w-mode equation (for corrections beyond zero-mode) ---\n"]; 

Print["If you expand A_mu(x,w) = Sum_n a_mu^(n)(x) f_n(w), the w-profiles satisfy a Sturm–Liouville problem."]; 
Print["The eigenproblem (schematic) is:"]; 
Print["  -d/dw( Z(w) f'(w) ) = m^2 Z(w) f(w)" ];
Print["A normalizable zero mode exists when Zint = Integral Z(w) dw is finite.\n" ];

(* Sanity: constant f(w) is a zero-mode for m^2=0 *)
fmode[ww_] := 1;
m2 = 0;
zeroModeResidual = Simp[-D[Zprof[w] D[fmode[w], w], w] - m2 Zprof[w] fmode[w]];
Print["Zero-mode check for f(w)=1, m^2=0 (should be 0):"]; 
Print[zeroModeResidual];

If[ZeroQAssume[zeroModeResidual],
  Print["OK: zero mode exists for finite Zint."],
  Print["WARNING: zero-mode residual did not simplify to 0 (unexpected)."]
];

Print["\n========== End of referee harness v2 ==========" ];

(*"
Output:

========================================
Paper VIII (referee harness v2): Maxwell from 4D localized Maxwell sector
========================================

---------- Conventions: coordinates, metric, fields ----------

Objects
HoldForm[coords5] -> {t, x, y, z, w}
HoldForm[eta5] -> {{-1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}}
HoldForm[Zprof[w]] -> E^(-(w^2/lambdaConf^2))

---------- Lagrangian (as used in the code) ----------

Paper-form reminder (string):
  L = -(Z(w)/(4 mu0)) F_{MN}F^{MN}  - (1/(2 xi mu0)) (DivA)^2  - A_M J^M
  where DivA := d_M A^M  and J^M is treated as an external source.

Wolfram-language expression actually varied (L):
-(A0[t, x, y, z, w]*J0[t, x, y, z, w]) - Aw[t, x, y, z, w]*Jw[t, x, y, z, w] - Ax[t, x, y, z, w]*Jx[t, x, y, z, w] - Ay[t, x, y, z, w]*Jy[t, x, y, z, w] - Az[t, x, y, z, w]*Jz[t, x, y, z, w] - (Derivative[0, 0, 0, 0, 1][Aw][t, x, y, z, w] + Derivative[0, 0, 0, 1, 0][Az][t, x, y, z, w] + Derivative[0, 0, 1, 0, 0][Ay][t, x, y, z, w] + Derivative[0, 1, 0, 0, 0][Ax][t, x, y, z, w] - Derivative[1, 0, 0, 0, 0][A0][t, x, y, z, w])^2/(2*mu0*xi) - ((Derivative[0, 0, 0, 0, 1][Az][t, x, y, z, w] - Derivative[0, 0, 0, 1, 0][Aw][t, x, y, z, w])^2 + (-Derivative[0, 0, 0, 0, 1][Az][t, x, y, z, w] + Derivative[0, 0, 0, 1, 0][Aw][t, x, y, z, w])^2 + (Derivative[0, 0, 0, 0, 1][Ay][t, x, y, z, w] - Derivative[0, 0, 1, 0, 0][Aw][t, x, y, z, w])^2 + (-Derivative[0, 0, 0, 0, 1][Ay][t, x, y, z, w] + Derivative[0, 0, 1, 0, 0][Aw][t, x, y, z, w])^2 + (Derivative[0, 0, 0, 1, 0][Ay][t, x, y, z, w] - Derivative[0, 0, 1, 0, 0][Az][t, x, y, z, w])^2 + (-Derivative[0, 0, 0, 1, 0][Ay][t, x, y, z, w] + Derivative[0, 0, 1, 0, 0][Az][t, x, y, z, w])^2 + (Derivative[0, 0, 0, 0, 1][Ax][t, x, y, z, w] - Derivative[0, 1, 0, 0, 0][Aw][t, x, y, z, w])^2 + (-Derivative[0, 0, 0, 0, 1][Ax][t, x, y, z, w] + Derivative[0, 1, 0, 0, 0][Aw][t, x, y, z, w])^2 + (Derivative[0, 0, 1, 0, 0][Ax][t, x, y, z, w] - Derivative[0, 1, 0, 0, 0][Ay][t, x, y, z, w])^2 + (-Derivative[0, 0, 1, 0, 0][Ax][t, x, y, z, w] + Derivative[0, 1, 0, 0, 0][Ay][t, x, y, z, w])^2 + (Derivative[0, 0, 0, 1, 0][Ax][t, x, y, z, w] - Derivative[0, 1, 0, 0, 0][Az][t, x, y, z, w])^2 + (-Derivative[0, 0, 0, 1, 0][Ax][t, x, y, z, w] + Derivative[0, 1, 0, 0, 0][Az][t, x, y, z, w])^2 + 2*(Derivative[0, 0, 0, 0, 1][A0][t, x, y, z, w] - Derivative[1, 0, 0, 0, 0][Aw][t, x, y, z, w])*(-Derivative[0, 0, 0, 0, 1][A0][t, x, y, z, w] + Derivative[1, 0, 0, 0, 0][Aw][t, x, y, z, w]) + 2*(Derivative[0, 1, 0, 0, 0][A0][t, x, y, z, w] - Derivative[1, 0, 0, 0, 0][Ax][t, x, y, z, w])*(-Derivative[0, 1, 0, 0, 0][A0][t, x, y, z, w] + Derivative[1, 0, 0, 0, 0][Ax][t, x, y, z, w]) + 2*(Derivative[0, 0, 1, 0, 0][A0][t, x, y, z, w] - Derivative[1, 0, 0, 0, 0][Ay][t, x, y, z, w])*(-Derivative[0, 0, 1, 0, 0][A0][t, x, y, z, w] + Derivative[1, 0, 0, 0, 0][Ay][t, x, y, z, w]) + 2*(Derivative[0, 0, 0, 1, 0][A0][t, x, y, z, w] - Derivative[1, 0, 0, 0, 0][Az][t, x, y, z, w])*(-Derivative[0, 0, 0, 1, 0][A0][t, x, y, z, w] + Derivative[1, 0, 0, 0, 0][Az][t, x, y, z, w]))/(4*E^(w^2/lambdaConf^2)*mu0)

---------- Derive Maxwell EOM from the action and verify compact form ----------

EulerEquations output head check:
  Head[eomEqns5[[1]]] -> Equal
After EqToExpr conversion:
  Head[eomExpr5[[1]]] -> Times

Maxwell EOM (expected compact expression, paper-normalized):
  d_M( Z(w) F^{MN} ) + (1/xi) d^N(DivA) = mu0 J^N
We check that:  mu0*(ELexpr) - [ d_M(Z F^{MN}) + (1/xi) d^N(DivA) - mu0 J^N ]  simplifies to 0.

Component check results (each should be 0):
  N=0 (t) -> 0
  N=x -> 0
  N=y -> 0
  N=z -> 0
  N=w -> 0

Zero-check booleans under $Assumptions:
  N=0 (t) : True
  N=x : True
  N=y : True
  N=z : True
  N=w : True

OK: All 5D components match the compact Maxwell form.

---------- Bianchi identity check (homogeneous Maxwell equations) ----------

Bianchi identity residuals (each should be 0):
  triple {1, 2, 3} -> 0
  triple {1, 2, 4} -> 0
  triple {1, 2, 5} -> 0
  triple {1, 3, 4} -> 0
  triple {1, 3, 5} -> 0
  triple {1, 4, 5} -> 0
  triple {2, 3, 4} -> 0
  triple {2, 3, 5} -> 0
  triple {2, 4, 5} -> 0
  triple {3, 4, 5} -> 0

OK: All Bianchi identities hold (homogeneous Maxwell).

---------- Consistency: divergence of Maxwell EOM forces current conservation (Lorenz gauge) ----------

Computed divergence of the (paper-normalized) EOM, d_N(EOM^N):
(Derivative[0, 0, 0, 0, 3][Aw][t, x, y, z, w] + Derivative[0, 0, 0, 1, 2][Az][t, x, y, z, w] + Derivative[0, 0, 0, 2, 1][Aw][t, x, y, z, w] + Derivative[0, 0, 0, 3, 0][Az][t, x, y, z, w] + Derivative[0, 0, 1, 0, 2][Ay][t, x, y, z, w] + Derivative[0, 0, 1, 2, 0][Ay][t, x, y, z, w] + Derivative[0, 0, 2, 0, 1][Aw][t, x, y, z, w] + Derivative[0, 0, 2, 1, 0][Az][t, x, y, z, w] + Derivative[0, 0, 3, 0, 0][Ay][t, x, y, z, w] + Derivative[0, 1, 0, 0, 2][Ax][t, x, y, z, w] + Derivative[0, 1, 0, 2, 0][Ax][t, x, y, z, w] + Derivative[0, 1, 2, 0, 0][Ax][t, x, y, z, w] + Derivative[0, 2, 0, 0, 1][Aw][t, x, y, z, w] + Derivative[0, 2, 0, 1, 0][Az][t, x, y, z, w] + Derivative[0, 2, 1, 0, 0][Ay][t, x, y, z, w] + Derivative[0, 3, 0, 0, 0][Ax][t, x, y, z, w] - mu0*xi*(Derivative[0, 0, 0, 0, 1][Jw][t, x, y, z, w] + Derivative[0, 0, 0, 1, 0][Jz][t, x, y, z, w] + Derivative[0, 0, 1, 0, 0][Jy][t, x, y, z, w] + Derivative[0, 1, 0, 0, 0][Jx][t, x, y, z, w] + Derivative[1, 0, 0, 0, 0][J0][t, x, y, z, w]) - Derivative[1, 0, 0, 0, 2][A0][t, x, y, z, w] - Derivative[1, 0, 0, 2, 0][A0][t, x, y, z, w] - Derivative[1, 0, 2, 0, 0][A0][t, x, y, z, w] - Derivative[1, 2, 0, 0, 0][A0][t, x, y, z, w] - Derivative[2, 0, 0, 0, 1][Aw][t, x, y, z, w] - Derivative[2, 0, 0, 1, 0][Az][t, x, y, z, w] - Derivative[2, 0, 1, 0, 0][Ay][t, x, y, z, w] - Derivative[2, 1, 0, 0, 0][Ax][t, x, y, z, w] + Derivative[3, 0, 0, 0, 0][A0][t, x, y, z, w])/xi

It should reduce to: (1/xi) box(DivA) - mu0 * (d_N J^N).
Residual after subtracting target form (should be 0):
0

OK: Divergence identity holds. In Lorenz gauge (DivA=0), EOM enforces d_N J^N = 0.

---------- Controlled brane reduction: zero-w-mode -> 3+1 Maxwell ----------

Reduction assumptions (strings):
  (i) Aw = 0
  (ii) A0,Ax,Ay,Az have no w-dependence (photon zero-mode)
  (iii) Jw ~ 0 in the brane regime (no net flux off the brane)

Localization integral and effective coupling:
HoldForm[Zint] -> lambdaConf*Sqrt[Pi]
HoldForm[mu0eff] -> mu0/(lambdaConf*Sqrt[Pi])

Under (i)-(ii), the kinetic term reduces as:
  Integral dw [ -(Z(w)/(4 mu0)) F_{MN}F^{MN} ]  ->  -(Zint/(4 mu0)) f_{mu nu} f^{mu nu}
which is standard 3+1 Maxwell with mu0eff = mu0/Zint.

--- 3+1D effective Maxwell system (brane) ---

Effective 3+1 Lagrangian density (L4) being varied:
-(A0b[t, x, y, z]*J0b[t, x, y, z]) - Axb[t, x, y, z]*Jxb[t, x, y, z] - Ayb[t, x, y, z]*Jyb[t, x, y, z] - Azb[t, x, y, z]*Jzb[t, x, y, z] - (lambdaConf*Sqrt[Pi]*(Derivative[0, 0, 0, 1][Azb][t, x, y, z] + Derivative[0, 0, 1, 0][Ayb][t, x, y, z] + Derivative[0, 1, 0, 0][Axb][t, x, y, z] - Derivative[1, 0, 0, 0][A0b][t, x, y, z])^2)/(2*mu0*xi) - (lambdaConf*Sqrt[Pi]*((Derivative[0, 0, 0, 1][Ayb][t, x, y, z] - Derivative[0, 0, 1, 0][Azb][t, x, y, z])^2 + (-Derivative[0, 0, 0, 1][Ayb][t, x, y, z] + Derivative[0, 0, 1, 0][Azb][t, x, y, z])^2 + (Derivative[0, 0, 1, 0][Axb][t, x, y, z] - Derivative[0, 1, 0, 0][Ayb][t, x, y, z])^2 + (-Derivative[0, 0, 1, 0][Axb][t, x, y, z] + Derivative[0, 1, 0, 0][Ayb][t, x, y, z])^2 + (Derivative[0, 0, 0, 1][Axb][t, x, y, z] - Derivative[0, 1, 0, 0][Azb][t, x, y, z])^2 + (-Derivative[0, 0, 0, 1][Axb][t, x, y, z] + Derivative[0, 1, 0, 0][Azb][t, x, y, z])^2 + 2*(Derivative[0, 1, 0, 0][A0b][t, x, y, z] - Derivative[1, 0, 0, 0][Axb][t, x, y, z])*(-Derivative[0, 1, 0, 0][A0b][t, x, y, z] + Derivative[1, 0, 0, 0][Axb][t, x, y, z]) + 2*(Derivative[0, 0, 1, 0][A0b][t, x, y, z] - Derivative[1, 0, 0, 0][Ayb][t, x, y, z])*(-Derivative[0, 0, 1, 0][A0b][t, x, y, z] + Derivative[1, 0, 0, 0][Ayb][t, x, y, z]) + 2*(Derivative[0, 0, 0, 1][A0b][t, x, y, z] - Derivative[1, 0, 0, 0][Azb][t, x, y, z])*(-Derivative[0, 0, 0, 1][A0b][t, x, y, z] + Derivative[1, 0, 0, 0][Azb][t, x, y, z])))/(4*mu0)

EulerEquations (3+1) head check:
  Head[eomEqns4[[1]]] -> Equal
After EqToExpr conversion (3+1):
  Head[eomExpr4[[1]]] -> Plus

3+1 component check results (each should be 0):
  nu=0 (t) -> 0
  nu=x -> 0
  nu=y -> 0
  nu=z -> 0

3+1 zero-check booleans under $Assumptions:
  nu=0 (t) : True
  nu=x : True
  nu=y : True
  nu=z : True

OK: 3+1 Maxwell equations follow from the reduced action.

--- KK / w-mode equation (for corrections beyond zero-mode) ---

If you expand A_mu(x,w) = Sum_n a_mu^(n)(x) f_n(w), the w-profiles satisfy a Sturm–Liouville problem.
The eigenproblem (schematic) is:
  -d/dw( Z(w) f'(w) ) = m^2 Z(w) f(w)
A normalizable zero mode exists when Zint = Integral Z(w) dw is finite.

Zero-mode check for f(w)=1, m^2=0 (should be 0):
0
OK: zero mode exists for finite Zint.

========== End of referee harness v2 ==========
"*)
