
(* 4D -> 3+1 Maxwell -> two-fluid -> MHD reduction (symbolic)
   Referee-proof Wolfram Language script: v11
   - avoids protected symbols
   - avoids Solve pitfalls (uses manual + LinearSolve cross-checks)
   - uses explicit div/curl operators to prevent CAS false negatives
   - includes Maxwell=>continuity and Poynting theorem checks that close the loop
   - treats species charges as fixed plasma inputs; no circulation-based charge law is used
   - records the equivalent reduced source-charge packaging qEff = qStar/Sqrt[Zint]
   - treats suppression of (Aw, Jw, F_{mu w}) as a controlled far-field brane reduction
   - now also checks Gaussian mode-tower identities, projected topology-EMF bookkeeping,
     open-system leakage bookkeeping, and geometry-force formulas used in the paper
*)

Print["--- 4D -> 3+1 Maxwell -> two-fluid -> MHD reduction (symbolic) ---"];

(* =========
   0) Setup
   ========= *)

ClearAll[
  mu0Eff, mu0eff, mu0, eps0, c, Zint,
  eC, eEff, qIon, qElectron, qStar, qEff, etaQ, eStar,
  nDen, mi, me,
  cross3, dot3, outer3, grad3, div3, curl3, divTensor3,
  checkScalar, checkVector, checkListZero, checkTrue, checkEqScalar, checkEqVector,
  vVec, JVec, viVec, veVec, BVec, gpEVec, EVec,
  x, y, z, t
];

mu0Eff = mu0/Zint;
mu0eff = mu0Eff;
eEff = eC;
qIon = +eEff;
qElectron = -eEff;
qStar = etaQ eStar;
qEff = qStar/Sqrt[Zint];

Print["Effective brane coupling:  mu0Eff = mu0/Zint = ", mu0Eff];
Print["Equivalent canonical charge packaging:  qEff = qStar/Sqrt[Zint] with qStar = etaQ eStar and eStar > 0."];
Print["Fixed species charge labels:  qIon = +eEff, qElectron = -eEff, with eEff the fixed brane-observed elementary charge magnitude."];
Print["Controlled reduction note: suppressing (Aw, Jw, F_{mu w}) is a far-field brane reduction assumption used to recover standard MHD; it does not deny microscopic mixed-core structure."];
Print["Charge-ontology note: this harness treats species charges as fixed inputs and does not use any circulation-based charge law."];

(* =========
   1) Helper operators (explicit, not built-ins)
   ========= *)

cross3[a_List, b_List] := {
  a[[2]] b[[3]] - a[[3]] b[[2]],
  a[[3]] b[[1]] - a[[1]] b[[3]],
  a[[1]] b[[2]] - a[[2]] b[[1]]
};

dot3[a_List, b_List] := a[[1]] b[[1]] + a[[2]] b[[2]] + a[[3]] b[[3]];

outer3[a_List, b_List] := Table[a[[i]] b[[j]], {i, 1, 3}, {j, 1, 3}];

grad3[f_, {xx_, yy_, zz_}] := {D[f, xx], D[f, yy], D[f, zz]};

div3[F_List, {xx_, yy_, zz_}] := D[F[[1]], xx] + D[F[[2]], yy] + D[F[[3]], zz];

curl3[F_List, {xx_, yy_, zz_}] := {
  D[F[[3]], yy] - D[F[[2]], zz],
  D[F[[1]], zz] - D[F[[3]], xx],
  D[F[[2]], xx] - D[F[[1]], yy]
};

divTensor3[T_?MatrixQ, {xx_, yy_, zz_}] := {
  D[T[[1, 1]], xx] + D[T[[1, 2]], yy] + D[T[[1, 3]], zz],
  D[T[[2, 1]], xx] + D[T[[2, 2]], yy] + D[T[[2, 3]], zz],
  D[T[[3, 1]], xx] + D[T[[3, 2]], yy] + D[T[[3, 3]], zz]
};

Print["(Helpers defined: cross3, dot3, outer3, grad3, div3, curl3, divTensor3)"];

(* =========
   2) PASS/FAIL harness
   ========= *)

checkScalar[name_String, expr_, assum_: True] := Module[{res},
  res = FullSimplify[expr, assum];
  If[res === 0,
    Print["PASS: ", name],
    Print["FAIL: ", name, "\n  residual -> ", res]
  ];
];

checkVector[name_String, expr_List, assum_: True] := Module[{res},
  res = FullSimplify[expr, assum];
  If[res === {0, 0, 0},
    Print["PASS: ", name],
    Print["FAIL: ", name, "\n  residual -> ", res]
  ];
];

checkListZero[name_String, exprs_List, assum_: True] := Module[{res},
  res = FullSimplify[exprs, assum];
  If[VectorQ[res, (# === 0) &],
    Print["PASS: ", name],
    Print["FAIL: ", name, "\n  residual -> ", res]
  ];
];

checkTrue[name_String, cond_, assum_: True] := Module[{res},
  res = FullSimplify[cond, assum];
  If[TrueQ[res],
    Print["PASS: ", name],
    Print["FAIL: ", name, "\n  condition -> ", res]
  ];
];

checkEqScalar[name_String, lhs_, rhs_, assum_: True] := checkScalar[name, lhs - rhs, assum];
checkEqVector[name_String, lhs_List, rhs_List, assum_: True] := checkVector[name, lhs - rhs, assum];

(* =========
   3) Two-fluid definitions -> solve vi, ve in terms of v, J
      Charge labels are fixed inputs: qIon = +eEff, qElectron = -eEff.
   ========= *)

vVec = {vX, vY, vZ};
JVec = {JX, JY, JZ};

Print["Mass density rho = (mi+me) nDen = ", (mi + me) nDen];
Print["Current definition:  J = eEff nDen (vi - ve), with eEff carried symbolically as eC below."];
Print["COM definition:      v = (mi vi + me ve)/(mi+me)"];

(* Manual solution (vector form) *)
viVec = vVec + (me/(eC nDen (mi + me))) JVec;
veVec = vVec - (mi/(eC nDen (mi + me))) JVec;

Print["Solved velocities (manual):"];
Print["  vi = ", viVec];
Print["  ve = ", veVec];

(* LinearSolve cross-check per component *)
With[
  {M = {{eC nDen, -eC nDen}, {mi, me}}},
  Module[{solX, solY, solZ, viLS, veLS},
    solX = LinearSolve[M, {JX, (mi + me) vX}];
    solY = LinearSolve[M, {JY, (mi + me) vY}];
    solZ = LinearSolve[M, {JZ, (mi + me) vZ}];
    viLS = {solX[[1]], solY[[1]], solZ[[1]]};
    veLS = {solX[[2]], solY[[2]], solZ[[2]]};

    checkEqVector["Velocity cross-check: vi(LinearSolve) == vi(manual)", viLS, viVec,
      eC != 0 && nDen != 0 && mi != 0 && me != 0 && (mi + me) != 0
    ];
    checkEqVector["Velocity cross-check: ve(LinearSolve) == ve(manual)", veLS, veVec,
      eC != 0 && nDen != 0 && mi != 0 && me != 0 && (mi + me) != 0
    ];
  ]
];

(* Consistency checks *)
checkEqVector[
  "Consistency: J == e n (vi - ve)",
  JVec,
  eC nDen (viVec - veVec),
  eC != 0 && nDen != 0 && (mi + me) != 0
];

checkEqVector[
  "Consistency: v == (mi vi + me ve)/(mi+me)",
  vVec,
  (mi viVec + me veVec)/(mi + me),
  (mi + me) != 0
];

(* me -> 0 limit *)
Print["me -> 0 limit:"];
Print["  vi -> ", Simplify[viVec /. me -> 0, eC != 0 && nDen != 0 && mi != 0]];
Print["  ve -> ", Simplify[veVec /. me -> 0, eC != 0 && nDen != 0 && mi != 0]];

checkEqVector["me->0: vi == v", Simplify[viVec /. me -> 0], vVec, eC != 0 && nDen != 0 && mi != 0];
checkEqVector[
  "me->0: ve == v - J/(e n)",
  Simplify[veVec /. me -> 0],
  vVec - (1/(eC nDen)) JVec,
  eC != 0 && nDen != 0 && mi != 0
];

(* =========
   4) Ohm law from electron force balance (pressure + E,B only)
   ========= *)

Print["--- Ohm law derivation + checks ---"];

BVec = {BX, BY, BZ};
gpEVec = {gpEX, gpEY, gpEZ};

(* Electron balance (inertialess form): 0 = -e n (E + ve×B) - ∇p_e  *)
EVec = -cross3[veVec, BVec] - (1/(eC nDen)) gpEVec;

Print["E (exact, in v,J): ", EVec];

checkEqVector[
  "Ohm (exact): E + v×B == (mi/(mi+me)) (J×B)/(e n) - (∇p_e)/(e n)",
  FullSimplify[EVec + cross3[vVec, BVec], eC != 0 && nDen != 0 && (mi + me) != 0],
  FullSimplify[(mi/(mi + me)) (cross3[JVec, BVec]/(eC nDen)) - gpEVec/(eC nDen),
    eC != 0 && nDen != 0 && (mi + me) != 0],
  eC != 0 && nDen != 0 && (mi + me) != 0
];

checkEqVector[
  "Ohm (me->0): E + v×B == (J×B)/(e n) - (∇p_e)/(e n)",
  FullSimplify[(EVec + cross3[vVec, BVec]) /. me -> 0, eC != 0 && nDen != 0 && mi != 0],
  FullSimplify[(cross3[JVec, BVec]/(eC nDen) - gpEVec/(eC nDen)), eC != 0 && nDen != 0 && mi != 0],
  eC != 0 && nDen != 0 && mi != 0
];

checkEqVector[
  "Ideal limit: (me->0, J=0, ∇p_e=0) => E == -v×B",
  FullSimplify[(EVec /. {me -> 0, JX -> 0, JY -> 0, JZ -> 0, gpEX -> 0, gpEY -> 0, gpEZ -> 0}),
    eC != 0 && nDen != 0 && mi != 0],
  -cross3[vVec, BVec],
  eC != 0 && nDen != 0 && mi != 0
];

(* =========
   5) Vector calculus checks (explicit operators)
   ========= *)

Print["--- Higher-value checks (vector calculus, explicit operators) ---"];

(* Induction from Faraday + ideal Ohm *)
Module[{coords, vFun, bFun, Eideal, induction},
  coords = {x, y, z};
  vFun = {v1[x, y, z, t], v2[x, y, z, t], v3f[x, y, z, t]};
  bFun = {b1[x, y, z, t], b2[x, y, z, t], b3f[x, y, z, t]};
  Eideal = -cross3[vFun, bFun];
  induction = D[bFun, t] + curl3[Eideal, coords] - (D[bFun, t] - curl3[cross3[vFun, bFun], coords]);
  checkVector["Induction: Faraday + (E=-v×B) => ∂t B = curl(v×B)", FullSimplify[induction]];
];

(* div(curl(F)) == 0 *)
Module[{coords, Ffun},
  coords = {x, y, z};
  Ffun = {g1[x, y, z, t], g2[x, y, z, t], g3[x, y, z, t]};
  checkScalar["Identity: div(curl(F)) == 0 (for smooth fields)", div3[curl3[Ffun, coords], coords]];
];

(* div(curl(v×B)) == 0 *)
Module[{coords, vFun, bFun},
  coords = {x, y, z};
  vFun = {v1[x, y, z, t], v2[x, y, z, t], v3f[x, y, z, t]};
  bFun = {b1[x, y, z, t], b2[x, y, z, t], b3f[x, y, z, t]};
  checkScalar["Divergence: div(curl(v×B)) == 0 (hence ∂t div(B)=0)", div3[curl3[cross3[vFun, bFun], coords], coords]];
];

(* div(∂t E) == ∂t div(E) *)
Module[{coords, Efun, lhs, rhs},
  coords = {x, y, z};
  Efun = {e1[x, y, z, t], e2[x, y, z, t], e3[x, y, z, t]};
  lhs = div3[D[Efun, t], coords];
  rhs = D[div3[Efun, coords], t];
  checkScalar["Identity: div(∂t E) == ∂t div(E)", lhs - rhs];
];

(* =========
   6) Even-higher-value: induction split (ideal + Hall + Biermann) + Biermann identity
   ========= *)

Print["--- Even-higher-value check: induction including Hall + Biermann terms ---"];

Module[
  {coords, vFun, bFun, nFun, pEfun, JFun, Eohm, hall, bier, rhsInduction, bierIdLHS, bierIdRHS},
  coords = {x, y, z};
  vFun = {v1[x, y, z, t], v2[x, y, z, t], v3f[x, y, z, t]};
  bFun = {b1[x, y, z, t], b2[x, y, z, t], b3f[x, y, z, t]};
  JFun = {j1[x, y, z, t], j2[x, y, z, t], j3f[x, y, z, t]};
  nFun = n0[x, y, z, t];
  pEfun = pe[x, y, z, t];

  hall = (1/(eC nFun)) cross3[JFun, bFun];
  bier = -(1/(eC nFun)) {D[pEfun, x], D[pEfun, y], D[pEfun, z]};

  (* Non-ideal Ohm (no resistivity, no electron inertia): E = -v×B + Hall + Biermann.
     Species charge labels remain fixed plasma inputs throughout. *)
  Eohm = -cross3[vFun, bFun] + hall + bier;

  (* Faraday => ∂t B = -curl(E) *)
  rhsInduction = -curl3[Eohm, coords];

  checkVector[
    "Non-ideal induction splits cleanly into ideal + Hall + (∇p_e) terms",
    FullSimplify[
      rhsInduction - (curl3[cross3[vFun, bFun], coords] - curl3[hall, coords] - curl3[bier, coords])
    ]
  ];

  (* Biermann identity: curl((∇p)/(e n)) = (1/e) ∇(1/n)×∇p *)
  bierIdLHS = curl3[(1/(eC nFun)) {D[pEfun, x], D[pEfun, y], D[pEfun, z]}, coords];
  bierIdRHS = (1/eC) cross3[{D[1/nFun, x], D[1/nFun, y], D[1/nFun, z]}, {D[pEfun, x], D[pEfun, y], D[pEfun, z]}];

  checkVector[
    "Biermann identity: curl((∇p_e)/(e n)) = (1/e) ∇(1/n)×∇p_e",
    FullSimplify[bierIdLHS - bierIdRHS]
  ];
];

(* =========
   7) Very strong: conservative (flux) form of the induction equation
   ========= *)

Print["--- Very-strong check: conservative (flux) form of induction ---"];

Module[{coords, vFun, bFun, fluxT, id},
  coords = {x, y, z};
  vFun = {v1[x, y, z, t], v2[x, y, z, t], v3f[x, y, z, t]};
  bFun = {b1[x, y, z, t], b2[x, y, z, t], b3f[x, y, z, t]};

  (* Standard MHD magnetic flux tensor: F = B⊗v - v⊗B with components F_ij = B_i v_j - v_i B_j *)
  fluxT = outer3[bFun, vFun] - outer3[vFun, bFun];

  (* Identity: curl(v×B) + div(F) == 0 *)
  id = curl3[cross3[vFun, bFun], coords] + divTensor3[fluxT, coords];
  checkVector["Flux-form identity: curl(v×B) + div(B⊗v - v⊗B) == 0", FullSimplify[id]];

  (* Equivalence of PDE residuals: (∂tB + divF) - (∂tB - curl(v×B)) == divF + curl(v×B) == 0 *)
  checkVector[
    "Flux-form induction equivalence: (∂tB = curl(v×B)) <=> (∂tB + divF = 0)",
    FullSimplify[
      (D[bFun, t] + divTensor3[fluxT, coords]) - (D[bFun, t] - curl3[cross3[vFun, bFun], coords])
    ]
  ];
];

(* =========
   8) Momentum structure checks (two-fluid sum)
   ========= *)

Print["--- Higher-value checks (momentum structure) ---"];

Module[
  {viF, veF, Bf, Efield, gpI, gpE, aI, aE, nF, ionEq, elecEq, sumEq, JdefF},
  viF = {vi1, vi2, vi3};
  veF = {ve1, ve2, ve3};
  Bf = {bX, bY, bZ};
  Efield = {eX, eY, eZ};
  gpI = {gpiX, gpiY, gpiZ};
  gpE = {gpeX, gpeY, gpeZ};
  aI = {aIX, aIY, aIZ};
  aE = {aEX, aEY, aEZ};
  nF = nDen;

  (* Two-fluid momentum (schematic, collisionless, no viscosity):
     mi n aI = + e n (E + vi×B) - ∇p_i
     me n aE = - e n (E + ve×B) - ∇p_e
  *)
  ionEq = mi nF aI - ( eC nF (Efield + cross3[viF, Bf]) - gpI );
  elecEq = me nF aE - ( -eC nF (Efield + cross3[veF, Bf]) - gpE );

  sumEq = FullSimplify[ionEq + elecEq];

  (* Strong check: the sum has no dependence on Efield *)
  checkVector[
    "Momentum add identity: (ionEq + elecEq) independent of E",
    FullSimplify[sumEq - (sumEq /. {Efield -> {0, 0, 0}})]
  ];

  checkVector["Cross bilinearity: vi×B - ve×B == (vi-ve)×B", FullSimplify[cross3[viF, Bf] - cross3[veF, Bf] - cross3[viF - veF, Bf]]];

  JdefF = eC nF (viF - veF);

  checkVector["Current rewrite: e n (vi-ve)×B == J×B", FullSimplify[eC nF cross3[viF - veF, Bf] - cross3[JdefF, Bf]]];

  (* Momentum sum is not an identity; it follows from imposing ionEq==0 and elecEq==0.
     We therefore (i) show the sum residual is exactly the single-fluid residual,
     and (ii) show that substituting the two-fluid equations makes it vanish. *)
  Module[{momRes, ionRules, elecRules},
    momRes = FullSimplify[mi nF aI + me nF aE - (cross3[JdefF, Bf] - (gpI + gpE))];

    checkVector[
      "Momentum sum residual equals (ionEq+elecEq)",
      FullSimplify[momRes - sumEq]
    ];

    ionRules = Thread[aI -> (1/(mi nF)) ( eC nF (Efield + cross3[viF, Bf]) - gpI )];
    elecRules = Thread[aE -> (1/(me nF)) ( -eC nF (Efield + cross3[veF, Bf]) - gpE )];

    checkVector[
      "Momentum sum: applying ionEq=0 and elecEq=0 gives mi n aI + me n aE = J×B - ∇(p_i+p_e)",
      FullSimplify[momRes /. Join[ionRules, elecRules]]
    ];
  ];
];

(* =========
   9) Maxwell => continuity (charge conservation)
   ========= *)

Print["--- Even-higher-value checks (Maxwell => continuity) ---"];

Module[
  {coords, EEM, BEM, JEM, rhoEM,
   divEexpr, divJexpr, divEdotexpr, dtDivEexpr,
   ampereVec, ampereDiv, contRes,
   ruleDivJFromAmp, ruleCommute, ruleGaussTime, contFromMaxwell,
   gaussProp},
  coords = {x, y, z};

  EEM = {ExM[x, y, z, t], EyM[x, y, z, t], EzM[x, y, z, t]};
  BEM = {BxM[x, y, z, t], ByM[x, y, z, t], BzM[x, y, z, t]};
  JEM = {JxM[x, y, z, t], JyM[x, y, z, t], JzM[x, y, z, t]};
  rhoEM = rhoC[x, y, z, t];

  divEexpr = div3[EEM, coords];
  divJexpr = div3[JEM, coords];
  divEdotexpr = div3[D[EEM, t], coords];
  dtDivEexpr = D[divEexpr, t];

  checkScalar["Identity: div(curl B) == 0", div3[curl3[BEM, coords], coords]];
  checkScalar["Identity: div(∂t E) == ∂t div(E)", divEdotexpr - dtDivEexpr];

  (* Ampere-Maxwell (SI): curl B = mu0 J + mu0 eps0 ∂t E -> residual *)
  ampereVec = curl3[BEM, coords] - mu0 JEM - mu0 eps0 D[EEM, t];
  ampereDiv = FullSimplify[div3[ampereVec, coords]];

  (* Divergence(Ampere) reduces to -μ0 divJ - μ0 eps0 div(∂t E) *)
  checkScalar[
    "Divergence(Ampere) reduces to -μ0 div(J) - μ0 eps0 div(∂t E)",
    FullSimplify[ampereDiv - (-mu0 divJexpr - mu0 eps0 divEdotexpr)]
  ];

  (* Continuity residual *)
  contRes = D[rhoEM, t] + divJexpr;

  (* Use divergence(Ampere)=0 -> divJ = -eps0 div(∂t E) and Gauss: divE = rho/eps0 *)
  ruleDivJFromAmp = (divJexpr -> -eps0 divEdotexpr);
  ruleCommute = (divEdotexpr -> dtDivEexpr);
  ruleGaussTime = (dtDivEexpr -> D[rhoEM/eps0, t]);

  contFromMaxwell = FullSimplify[contRes /. ruleDivJFromAmp /. ruleCommute /. ruleGaussTime];

  checkScalar["Maxwell => continuity: ∂t rho + div J == 0", contFromMaxwell];

  (* Gauss constraint propagation (uses only divergence(Ampere) + commutation) *)
  gaussProp = FullSimplify[
    (D[divEexpr - rhoEM/eps0, t] + contRes/eps0) /. ruleDivJFromAmp /. ruleCommute
  ];

  checkScalar[
    "Gauss propagation: ∂t(divE - rho/eps0) + (∂t rho + divJ)/eps0 == 0 (via div(Ampere))",
    gaussProp
  ];
];

(* =========
   10) Poynting theorem (Maxwell energy conservation)
   ========= *)

Print["--- Even-higher-value check: Poynting theorem (Maxwell energy conservation) ---"];

Module[
  {coords, EEM, BEM, JEM,
   curlEexpr, curlBexpr, dEdt, dBdt,
   uDot, Sdiv, poyRes,
   faradayRes, ampereRes, poyFactor, rulesMaxwell},
  coords = {x, y, z};

  EEM = {ExM[x, y, z, t], EyM[x, y, z, t], EzM[x, y, z, t]};
  BEM = {BxM[x, y, z, t], ByM[x, y, z, t], BzM[x, y, z, t]};
  JEM = {JxM[x, y, z, t], JyM[x, y, z, t], JzM[x, y, z, t]};

  curlEexpr = curl3[EEM, coords];
  curlBexpr = curl3[BEM, coords];
  dEdt = D[EEM, t];
  dBdt = D[BEM, t];

  (* Identity: div(E×B) = B·curl(E) - E·curl(B) *)
  checkScalar[
    "Identity: div(E×B) = B·curl(E) - E·curl(B)",
    FullSimplify[div3[cross3[EEM, BEM], coords] - (dot3[BEM, curlEexpr] - dot3[EEM, curlBexpr])]
  ];

  (* Check: ∂t( (1/2)E^2 ) = E·∂tE and same for B *)
  checkScalar[
    "Identity: ∂t( (1/2)E^2 ) = E·∂tE",
    FullSimplify[D[(1/2) dot3[EEM, EEM], t] - dot3[EEM, dEdt]]
  ];
  checkScalar[
    "Identity: ∂t( (1/2)B^2 ) = B·∂tB",
    FullSimplify[D[(1/2) dot3[BEM, BEM], t] - dot3[BEM, dBdt]]
  ];

  (* div(S) with S=(1/mu0)E×B computed via the identity above *)
  Sdiv = (1/mu0) (dot3[BEM, curlEexpr] - dot3[EEM, curlBexpr]);

  (* ∂t u = eps0 E·∂tE + (1/mu0) B·∂tB *)
  uDot = eps0 dot3[EEM, dEdt] + (1/mu0) dot3[BEM, dBdt];

  (* Poynting residual: ∂t u + div S + J·E *)
  poyRes = FullSimplify[uDot + Sdiv + dot3[JEM, EEM]];

  (* Define Maxwell residual vectors (to avoid fragile substitutions after expansion) *)
  faradayRes = curlEexpr + dBdt;                          (* Faraday: curlE + ∂tB = 0 *)
  ampereRes  = curlBexpr - mu0 JEM - mu0 eps0 dEdt;       (* Ampere:  curlB - μ0J - μ0ε0∂tE = 0 *)

  (* Factorization identity:
     ∂t u + div S + J·E == (1/mu0)( B·FaradayRes - E·AmpereRes )
     so it vanishes when Maxwell residuals vanish.
  *)
  poyFactor = (1/mu0) (dot3[BEM, faradayRes] - dot3[EEM, ampereRes]);

  checkScalar[
    "Poynting factorization: (∂t u + div S + J·E) - (1/mu0)(B·FaradayRes - E·AmpereRes) == 0",
    FullSimplify[poyRes - poyFactor]
  ];

  rulesMaxwell = Join[Thread[faradayRes -> {0, 0, 0}], Thread[ampereRes -> {0, 0, 0}]];

  checkScalar[
    "Poynting theorem: under Maxwell, ∂t u + div S + J·E = 0",
    FullSimplify[poyFactor /. rulesMaxwell]
  ];
];

(* =========
   11) Energetic sanity checks: Hall term neutrality + ideal transfer + Ohm projection
   ========= *)

Print["--- Energetic sanity checks (Hall term, ideal transfer, Ohm projection) ---"];

Module[{Jv, Bv, vv, gpEv, n0s, EohmMe0, eta, EohmRes, workMe0, workRes},
  Jv = {JX, JY, JZ};
  Bv = {BX, BY, BZ};
  vv = {vX, vY, vZ};
  gpEv = {gpEX, gpEY, gpEZ};
  n0s = nDen;

  checkScalar["Hall workless identity: J·(J×B) = 0", FullSimplify[dot3[Jv, cross3[Jv, Bv]]]];

  checkScalar["Ideal transfer identity: J·(-v×B) = v·(J×B)", FullSimplify[dot3[Jv, -cross3[vv, Bv]] - dot3[vv, cross3[Jv, Bv]]]];

  (* Explicitly show Hall term drops out of J·E for generalized Ohm (me->0) *)
  EohmMe0 = -cross3[vv, Bv] + (1/(eC n0s)) cross3[Jv, Bv] - (1/(eC n0s)) gpEv;

  workMe0 = FullSimplify[dot3[Jv, EohmMe0] - (dot3[vv, cross3[Jv, Bv]] - (1/(eC n0s)) dot3[Jv, gpEv])];
  checkScalar["Ohm work projection (me->0): J·E = v·(J×B) - (J·∇p_e)/(e n)", workMe0, eC != 0 && n0s != 0];

  (* Optional resistive extension sanity: add eta J and check J·(eta J)=eta |J|^2 *)
  eta = etaOhm;
  EohmRes = EohmMe0 + eta Jv;
  workRes = FullSimplify[dot3[Jv, EohmRes] - (dot3[Jv, EohmMe0] + eta dot3[Jv, Jv])];
  checkScalar["Resistive work identity: J·(eta J) = eta |J|^2", workRes];
];



(* =========
   11.5) Very-strong check: resistive dissipation (Joule heating) decomposition
   ========= *)

Print["--- Very-strong check: resistive dissipation (Joule heating) ---"];

Module[{Jv, Bv, vv, eta, Eideal, Eres, rhs},
  Jv = {JX, JY, JZ};
  Bv = {BX, BY, BZ};
  vv = {vX, vY, vZ};
  eta = etaOhm;

  Eideal = -cross3[vv, Bv];
  Eres = Eideal + eta Jv;

  (* Local power transfer: -J·E = -v·(J×B) - eta |J|^2 (ideal exchange + dissipation) *)
  rhs = -dot3[vv, cross3[Jv, Bv]] - eta dot3[Jv, Jv];

  checkScalar[
    "Resistive decomposition: -J·E = -v·(J×B) - eta |J|^2",
    FullSimplify[-dot3[Jv, Eres] - rhs]
  ];
];

(* =========
   11.6) Very-strong check: magnetic helicity (ideal MHD)
   ========= *)

Print["--- Very-strong check: magnetic helicity (ideal MHD) ---"];

Module[{coords, Avec, phi, Bvec, Evec, helicityRes, vvec},
  coords = {x, y, z};

  (* Vector potential and gauge scalar (arbitrary smooth fields) *)
  Avec = {AxH[x, y, z, t], AyH[x, y, z, t], AzH[x, y, z, t]};
  phi = phiH[x, y, z, t];

  (* Definitions *)
  Bvec = curl3[Avec, coords];
  Evec = -D[Avec, t] - grad3[phi, coords];

  (* Helicity density identity:
     ∂t(A·B) + div( φ B + E×A ) = -2 E·B
     Here B = curl A => div B = 0 automatically.
  *)
  helicityRes =
    FullSimplify[
      D[dot3[Avec, Bvec], t] + div3[phi Bvec + cross3[Evec, Avec], coords] + 2 dot3[Evec, Bvec]
    ];

  checkScalar[
    "Helicity density identity: ∂t(A·B)+div(φ B + E×A) = -2 E·B",
    helicityRes
  ];

  (* Ideal MHD implies E·B = 0 via E = -v×B *)
  vvec = {v1H[x, y, z, t], v2H[x, y, z, t], v3H[x, y, z, t]};
  checkScalar[
    "Ideal Ohm implies E·B=0: (-v×B)·B = 0",
    FullSimplify[dot3[-cross3[vvec, Bvec], Bvec]]
  ];
];

(* =========
   11.7) Exact checks: localized Maxwell reduction + Gaussian mode tower
   ========= *)

Print["--- Exact checks: localized Maxwell reduction + Gaussian mode tower ---"];

Module[{wLoc, lamLoc, zGauss, zIntGauss, modeResiduals, normResiduals, evenCouplingResiduals},
  zGauss[ww_] := Exp[-ww^2/lamLoc^2];
  zIntGauss = Assuming[lamLoc > 0, FullSimplify[Integrate[zGauss[wLoc], {wLoc, -Infinity, Infinity}]]];

  checkEqScalar["Gaussian localization integral: Zint = λ Sqrt[π]", zIntGauss, lamLoc Sqrt[Pi], lamLoc > 0];
  checkEqScalar["Effective Maxwell coupling under Gaussian localization: μ0eff = μ0/(λ Sqrt[π])", mu0/zIntGauss, mu0/(lamLoc Sqrt[Pi]), lamLoc > 0];
  checkEqScalar["Canonical charge packaging under Gaussian localization: qEff = qStar/Sqrt[λ Sqrt[π]]", qStar/Sqrt[zIntGauss], qStar/Sqrt[lamLoc Sqrt[Pi]], lamLoc > 0];

  modeResiduals = Table[
    FullSimplify[
      -D[zGauss[wLoc] D[HermiteH[nMode, wLoc/lamLoc], wLoc], wLoc] -
      ((2 nMode)/(lamLoc^2)) zGauss[wLoc] HermiteH[nMode, wLoc/lamLoc],
      lamLoc > 0
    ],
    {nMode, 0, 6}
  ];
  checkListZero["Gaussian soft-wall Sturm-Liouville modes: Hermite profiles satisfy the weighted eigenproblem for n=0..6", modeResiduals, lamLoc > 0];

  normResiduals = Table[
    FullSimplify[
      Integrate[
        zGauss[wLoc] HermiteH[nMode, wLoc/lamLoc]^2,
        {wLoc, -Infinity, Infinity},
        Assumptions -> lamLoc > 0
      ] - lamLoc Sqrt[Pi] 2^nMode nMode!,
      lamLoc > 0
    ],
    {nMode, 0, 6}
  ];
  checkListZero["Gaussian mode norms: ||H_n(w/λ)||_Z^2 = λ Sqrt[π] 2^n n! for n=0..6", normResiduals, lamLoc > 0];

  checkListZero["Odd-mode brane coupling selection rule: H_{2m+1}(0)=0 for m=0..4", Table[HermiteH[2 m + 1, 0], {m, 0, 4}]];

  evenCouplingResiduals = Table[
    FullSimplify[
      (HermiteH[2 m, 0]^2)/(lamLoc Sqrt[Pi] 2^(2 m) (2 m)!) -
      (1/(lamLoc Sqrt[Pi])) Binomial[2 m, m]/4^m,
      lamLoc > 0
    ],
    {m, 0, 4}
  ];
  checkListZero["Even-mode brane couplings: c_{2m} = (1/(λ Sqrt[π])) Binomial[2m,m]/4^m for m=0..4", evenCouplingResiduals, lamLoc > 0];
  checkEqScalar["Zero-mode coupling matches controlled reduction: c0 = 1/Zint", 1/(lamLoc Sqrt[Pi]), 1/zIntGauss, lamLoc > 0];
];

(* =========
   11.8) Exact checks: projected generalized Ohm bookkeeping
   ========= *)

Print["--- Exact checks: projected generalized Ohm bookkeeping ---"];

Module[
  {
    barE, barv, barB, barVe, barJ, barGradPe, barExtra, overVcrossB,
    pressureProj, resistProj, inertiaProj, ewTopo, ecov, ep, er, ei, etopo,
    hallExact, hallClosure, etaProj, iBar, exactProjected, splitClosed, splitResidual, finalResidual
  },
  barE = {Eb1, Eb2, Eb3};
  barv = {vb1, vb2, vb3};
  barB = {Bb1, Bb2, Bb3};
  barVe = {veb1, veb2, veb3};
  barJ = {Jb1, Jb2, Jb3};
  barGradPe = {gpb1, gpb2, gpb3};
  barExtra = {xpb1, xpb2, xpb3};
  overVcrossB = {ovcb1, ovcb2, ovcb3};
  pressureProj = {praw1, praw2, praw3};
  resistProj = {rraw1, rraw2, rraw3};
  inertiaProj = {iraw1, iraw2, iraw3};
  ewTopo = {ew1, ew2, ew3};
  etaProj = etaOhm;
  iBar = {ibar1, ibar2, ibar3};

  exactProjected = -overVcrossB + ewTopo + pressureProj + resistProj + inertiaProj + barExtra;
  hallExact = cross3[barv - barVe, barB];

  splitResidual = FullSimplify[
    (exactProjected + cross3[barv, barB]) -
    (hallExact + (cross3[barVe, barB] - overVcrossB) + ewTopo + pressureProj + resistProj + inertiaProj + barExtra)
  ];
  checkVector["Projected exact Ohm split: add-subtract vbar×Bbar exposes Hall + covariance + projected remainder", splitResidual];

  ecov = cross3[barVe, barB] - overVcrossB;
  ep = pressureProj + barGradPe/(eC nDen);
  er = resistProj - etaProj barJ;
  ei = inertiaProj + (me/eC) iBar;
  etopo = ecov + ewTopo + ep + er + ei;
  hallClosure = cross3[barJ, barB]/(eC nDen);
  splitClosed = hallExact + ecov + ewTopo + pressureProj + resistProj + inertiaProj + barExtra;

  finalResidual = FullSimplify[
    (splitClosed /. Thread[(barv - barVe) -> barJ/(eC nDen)]) -
    (hallClosure - barGradPe/(eC nDen) + etaProj barJ - (me/eC) iBar + etopo + barExtra),
    eC != 0 && nDen != 0
  ];
  checkVector["Projected generalized Ohm law closes with E_topo bookkeeping", finalResidual, eC != 0 && nDen != 0];

  checkVector[
    "Covariance EMF vanishes in the zero-mode reduction: overline(ve×B) = vbar_e×Bbar",
    FullSimplify[ecov /. Thread[overVcrossB -> cross3[barVe, barB]]]
  ];
];

(* =========
   11.9) Exact checks: open-system leakage + subscale helicity bookkeeping
   ========= *)

Print["--- Exact checks: open-system leakage + subscale helicity bookkeeping ---"];

Module[{wLoc, wKernel, sW, lhsLeak, rhsLeak, hProj, hRes, fProj, fRes, ebProj, ebRes, hSub, fSub, xLoc},
  wKernel[ww_] := Exp[-ww^2];
  sW[ww_] := (1 + ww) Exp[-ww^2/2];

  lhsLeak = FullSimplify[-Integrate[wKernel[wLoc] D[sW[wLoc], wLoc], {wLoc, -Infinity, Infinity}]];
  rhsLeak = FullSimplify[
    -(Limit[wKernel[wLoc] sW[wLoc], wLoc -> Infinity] - Limit[wKernel[wLoc] sW[wLoc], wLoc -> -Infinity]) +
    Integrate[D[wKernel[wLoc], wLoc] sW[wLoc], {wLoc, -Infinity, Infinity}]
  ];
  checkEqScalar["Projected leakage identity: -∫W ∂w S = -[WS] + ∫W' S", lhsLeak, rhsLeak];

  hProj = hProjF[xLoc, t];
  hRes = hResF[xLoc, t];
  fProj = fProjF[xLoc, t];
  fRes = fResF[xLoc, t];
  ebProj = ebProjF[xLoc, t];
  ebRes = ebResF[xLoc, t];
  hSub = hProj - hRes;
  fSub = fProj - fRes;

  checkScalar[
    "Subscale helicity bookkeeping: subtracting projected and resolved budgets yields the covariance source",
    FullSimplify[
      (D[hProj, t] + D[fProj, xLoc] + 2 ebProj) -
      (D[hRes, t] + D[fRes, xLoc] + 2 ebRes) -
      (D[hSub, t] + D[fSub, xLoc] + 2 (ebProj - ebRes))
    ]
  ];
];

(* =========
   11.10) Exact checks: geometry-force ledger formulas
   ========= *)

Print["--- Exact checks: geometry-force ledger formulas ---"];

Module[{wLoc, lamLoc, zGauss, kappa, alpha, aVar, lVar, pVac, sigmaWall, vTube, aTube, eGeom},
  zGauss[ww_] := Exp[-ww^2/lamLoc^2];

  checkEqScalar["Gaussian width derivative: ∂λ Z = (2 w^2 / λ^3) Z", D[zGauss[wLoc], lamLoc], (2 wLoc^2/lamLoc^3) zGauss[wLoc], lamLoc > 0];
  checkTrue["Gaussian-width EM force density is compressive: -(∂λ Z) <= 0 for λ>0", -D[zGauss[wLoc], lamLoc] <= 0, lamLoc > 0 && Element[wLoc, Reals]];

  checkEqScalar["Ratio-penalty force on a: F_a^{ratio} = κ α (L - α a)", -D[(kappa/2) (lVar - alpha aVar)^2, aVar], kappa alpha (lVar - alpha aVar)];
  checkEqScalar["Ratio-penalty force on L: F_L^{ratio} = -κ (L - α a)", -D[(kappa/2) (lVar - alpha aVar)^2, lVar], -kappa (lVar - alpha aVar)];

  vTube = (4 Pi/3) aVar^3 lVar;
  aTube = 4 Pi aVar^2 lVar + (8 Pi/3) aVar^3;
  eGeom = pVac vTube + sigmaWall aTube;

  checkEqScalar[
    "Tube geometry force on a from E_geom(a,L)",
    -D[eGeom, aVar],
    -(4 Pi pVac aVar^2 lVar + sigmaWall (8 Pi aVar lVar + 8 Pi aVar^2))
  ];
  checkEqScalar[
    "Tube geometry force on L from E_geom(a,L)",
    -D[eGeom, lVar],
    -((4 Pi/3) pVac aVar^3 + 4 Pi sigmaWall aVar^2)
  ];
];

(* =========
   12) Paper-friendly display (no math used for checks)
   ========= *)

Print["--- Paper-friendly identities (display only) ---"];
Print["Ideal induction equation (display):"];
Print[HoldForm[D[BDisplay, t] == curl3[cross3[vDisplay, BDisplay], {x, y, z}]]];

Print["Done."];

(*"
Output:

--- 4D -> 3+1 Maxwell -> two-fluid -> MHD reduction (symbolic) ---
Effective brane coupling:  mu0Eff = mu0/Zint = mu0/Zint
Equivalent canonical charge packaging:  qEff = qStar/Sqrt[Zint] with qStar = etaQ eStar and eStar > 0.
Fixed species charge labels:  qIon = +eEff, qElectron = -eEff, with eEff the fixed brane-observed elementary charge magnitude.
Controlled reduction note: suppressing (Aw, Jw, F_{mu w}) is a far-field brane reduction assumption used to recover standard MHD; it does not deny microscopic mixed-core structure.
Charge-ontology note: this harness treats species charges as fixed inputs and does not use any circulation-based charge law.
(Helpers defined: cross3, dot3, outer3, grad3, div3, curl3, divTensor3)
Mass density rho = (mi+me) nDen = (me + mi)*nDen
Current definition:  J = eEff nDen (vi - ve), with eEff carried symbolically as eC below.
COM definition:      v = (mi vi + me ve)/(mi+me)
Solved velocities (manual):
  vi = {(JX*me)/(eC*(me + mi)*nDen) + vX, (JY*me)/(eC*(me + mi)*nDen) + vY, (JZ*me)/(eC*(me + mi)*nDen) + vZ}
  ve = {-((JX*mi)/(eC*(me + mi)*nDen)) + vX, -((JY*mi)/(eC*(me + mi)*nDen)) + vY, -((JZ*mi)/(eC*(me + mi)*nDen)) + vZ}
PASS: Velocity cross-check: vi(LinearSolve) == vi(manual)
PASS: Velocity cross-check: ve(LinearSolve) == ve(manual)
PASS: Consistency: J == e n (vi - ve)
PASS: Consistency: v == (mi vi + me ve)/(mi+me)
me -> 0 limit:
  vi -> {vX, vY, vZ}
  ve -> {-(JX/(eC*nDen)) + vX, -(JY/(eC*nDen)) + vY, -(JZ/(eC*nDen)) + vZ}
PASS: me->0: vi == v
PASS: me->0: ve == v - J/(e n)
--- Ohm law derivation + checks ---
E (exact, in v,J): {-(gpEX/(eC*nDen)) - BZ*(-((JY*mi)/(eC*(me + mi)*nDen)) + vY) + BY*(-((JZ*mi)/(eC*(me + mi)*nDen)) + vZ), -(gpEY/(eC*nDen)) + BZ*(-((JX*mi)/(eC*(me + mi)*nDen)) + vX) - BX*(-((JZ*mi)/(eC*(me + mi)*nDen)) + vZ), -(gpEZ/(eC*nDen)) - BY*(-((JX*mi)/(eC*(me + mi)*nDen)) + vX) + BX*(-((JY*mi)/(eC*(me + mi)*nDen)) + vY)}
PASS: Ohm (exact): E + v×B == (mi/(mi+me)) (J×B)/(e n) - (∇p_e)/(e n)
PASS: Ohm (me->0): E + v×B == (J×B)/(e n) - (∇p_e)/(e n)
PASS: Ideal limit: (me->0, J=0, ∇p_e=0) => E == -v×B
--- Higher-value checks (vector calculus, explicit operators) ---
PASS: Induction: Faraday + (E=-v×B) => ∂t B = curl(v×B)
PASS: Identity: div(curl(F)) == 0 (for smooth fields)
PASS: Divergence: div(curl(v×B)) == 0 (hence ∂t div(B)=0)
PASS: Identity: div(∂t E) == ∂t div(E)
--- Even-higher-value check: induction including Hall + Biermann terms ---
PASS: Non-ideal induction splits cleanly into ideal + Hall + (∇p_e) terms
PASS: Biermann identity: curl((∇p_e)/(e n)) = (1/e) ∇(1/n)×∇p_e
--- Very-strong check: conservative (flux) form of induction ---
PASS: Flux-form identity: curl(v×B) + div(B⊗v - v⊗B) == 0
PASS: Flux-form induction equivalence: (∂tB = curl(v×B)) <=> (∂tB + divF = 0)
--- Higher-value checks (momentum structure) ---
PASS: Momentum add identity: (ionEq + elecEq) independent of E
PASS: Cross bilinearity: vi×B - ve×B == (vi-ve)×B
PASS: Current rewrite: e n (vi-ve)×B == J×B
PASS: Momentum sum residual equals (ionEq+elecEq)
PASS: Momentum sum: applying ionEq=0 and elecEq=0 gives mi n aI + me n aE = J×B - ∇(p_i+p_e)
--- Even-higher-value checks (Maxwell => continuity) ---
PASS: Identity: div(curl B) == 0
PASS: Identity: div(∂t E) == ∂t div(E)
PASS: Divergence(Ampere) reduces to -μ0 div(J) - μ0 eps0 div(∂t E)
PASS: Maxwell => continuity: ∂t rho + div J == 0
PASS: Gauss propagation: ∂t(divE - rho/eps0) + (∂t rho + divJ)/eps0 == 0 (via div(Ampere))
--- Even-higher-value check: Poynting theorem (Maxwell energy conservation) ---
PASS: Identity: div(E×B) = B·curl(E) - E·curl(B)
PASS: Identity: ∂t( (1/2)E^2 ) = E·∂tE
PASS: Identity: ∂t( (1/2)B^2 ) = B·∂tB
PASS: Poynting factorization: (∂t u + div S + J·E) - (1/mu0)(B·FaradayRes - E·AmpereRes) == 0
PASS: Poynting theorem: under Maxwell, ∂t u + div S + J·E = 0
--- Energetic sanity checks (Hall term, ideal transfer, Ohm projection) ---
PASS: Hall workless identity: J·(J×B) = 0
PASS: Ideal transfer identity: J·(-v×B) = v·(J×B)
PASS: Ohm work projection (me->0): J·E = v·(J×B) - (J·∇p_e)/(e n)
PASS: Resistive work identity: J·(eta J) = eta |J|^2
--- Very-strong check: resistive dissipation (Joule heating) ---
PASS: Resistive decomposition: -J·E = -v·(J×B) - eta |J|^2
--- Very-strong check: magnetic helicity (ideal MHD) ---
PASS: Helicity density identity: ∂t(A·B)+div(φ B + E×A) = -2 E·B
PASS: Ideal Ohm implies E·B=0: (-v×B)·B = 0
--- Exact checks: localized Maxwell reduction + Gaussian mode tower ---
PASS: Gaussian localization integral: Zint = λ Sqrt[π]
PASS: Effective Maxwell coupling under Gaussian localization: μ0eff = μ0/(λ Sqrt[π])
PASS: Canonical charge packaging under Gaussian localization: qEff = qStar/Sqrt[λ Sqrt[π]]
PASS: Gaussian soft-wall Sturm-Liouville modes: Hermite profiles satisfy the weighted eigenproblem for n=0..6
PASS: Gaussian mode norms: ||H_n(w/λ)||_Z^2 = λ Sqrt[π] 2^n n! for n=0..6
PASS: Odd-mode brane coupling selection rule: H_{2m+1}(0)=0 for m=0..4
PASS: Even-mode brane couplings: c_{2m} = (1/(λ Sqrt[π])) Binomial[2m,m]/4^m for m=0..4
PASS: Zero-mode coupling matches controlled reduction: c0 = 1/Zint
--- Exact checks: projected generalized Ohm bookkeeping ---
PASS: Projected exact Ohm split: add-subtract vbar×Bbar exposes Hall + covariance + projected remainder
PASS: Projected generalized Ohm law closes with E_topo bookkeeping
PASS: Covariance EMF vanishes in the zero-mode reduction: overline(ve×B) = vbar_e×Bbar
--- Exact checks: open-system leakage + subscale helicity bookkeeping ---
PASS: Projected leakage identity: -∫W ∂w S = -[WS] + ∫W' S
PASS: Subscale helicity bookkeeping: subtracting projected and resolved budgets yields the covariance source
--- Exact checks: geometry-force ledger formulas ---
PASS: Gaussian width derivative: ∂λ Z = (2 w^2 / λ^3) Z
PASS: Gaussian-width EM force density is compressive: -(∂λ Z) <= 0 for λ>0
PASS: Ratio-penalty force on a: F_a^{ratio} = κ α (L - α a)
PASS: Ratio-penalty force on L: F_L^{ratio} = -κ (L - α a)
PASS: Tube geometry force on a from E_geom(a,L)
PASS: Tube geometry force on L from E_geom(a,L)
--- Paper-friendly identities (display only) ---
Ideal induction equation (display):
HoldForm[D[BDisplay, t] == curl3[cross3[vDisplay, BDisplay], {x, y, z}]]
Done.
"*)
