ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

Print["STAGE 004 PROJECTED MAXWELL BUNDLE INDEX MATHEMATICA AUDIT"];

Clear[w, lam, mu0];
scaleAssumptions = lam > 0;
couplingAssumptions = lam > 0 && mu0 > 0;

(* M1: integration by parts with decaying test functions.
   Combine W*f' + W'*f into one integrand so Mathematica recognizes it as
   d/dw(W*f) and evaluates to the boundary term, which vanishes for the
   decaying Gaussian profile. *)
decayingWindow = Exp[-w^2/lam^2];
oddProbe = w*Exp[-w^2/lam^2];
m1Integrand = decayingWindow * D[oddProbe, w] + D[decayingWindow, w] * oddProbe;
m1Left = Integrate[m1Integrand, {w, -Infinity, Infinity}, Assumptions -> scaleAssumptions];
m1Right = 0;
m1Residual = FullSimplify[m1Left - m1Right, Assumptions -> scaleAssumptions];
Print["M1 residual = ", fmt[m1Residual]];
If[FullSimplify[m1Left - m1Right, Assumptions -> scaleAssumptions] =!= 0,
  Print["FAIL: M1 density-level integration by parts"]; Exit[1]
];
Print["PASS: M1 density-level integration by parts"];

(* M2: cyclic Bianchi identity for F = dA implies Maxwell-Faraday by
   Schwarz symmetry of mixed partials. Non-tautological: a sign error in
   the E,B<->F map produces a nonzero residual. *)
Clear[t, x, y, z, AA];
spaceTimeAssumptions = Element[{t, x, y, z}, Reals];
coordList = {t, x, y, z};
potentialList = {AA[0][t, x, y, z], AA[1][t, x, y, z],
                  AA[2][t, x, y, z], AA[3][t, x, y, z]};

(* fieldStrength evaluated explicitly per index pair, avoiding Module
   pre-evaluation pitfalls where Part-indexing on a pattern parameter
   inside a Do/Module body fires Part::pkspec1 with unbound gensym locals. *)
fStr[i_Integer, j_Integer] := D[potentialList[[j + 1]], coordList[[i + 1]]] -
                              D[potentialList[[i + 1]], coordList[[j + 1]]];

(* Precompute the six field-strength components needed below.  Using
   concrete integer literals here forces full evaluation before any
   Module/Do scope opens, so the cyclic Bianchi sums are built from
   already-evaluated D[AA[k], coord] expressions. *)
F01 = fStr[0, 1]; F02 = fStr[0, 2]; F03 = fStr[0, 3];
F12 = fStr[1, 2]; F13 = fStr[1, 3]; F23 = fStr[2, 3];
F10 = -F01; F20 = -F02; F30 = -F03;
F21 = -F12; F31 = -F13; F32 = -F23;

(* M2: cyclic Bianchi identity per triple, with the three field-strength
   components substituted in literally so no pattern-matching happens
   inside the loop body. *)
bianchiChecks = {
  {{0, 2, 3}, D[F23, t] + D[F30, y] + D[F02, z]},
  {{0, 3, 1}, D[F31, t] + D[F10, z] + D[F03, x]},
  {{0, 1, 2}, D[F12, t] + D[F20, x] + D[F01, y]}
};
Do[Module[{triple, cyc, residual},
    {triple, cyc} = entry;
    residual = FullSimplify[cyc, Assumptions -> spaceTimeAssumptions];
    Print["M2 cyclic Bianchi ", triple, " residual = ", fmt[residual]];
    If[residual =!= 0,
      Print["FAIL: M2 cyclic Bianchi ", triple]; Exit[1]
    ];
    Print["PASS: M2 cyclic Bianchi ", triple];
  ],
  {entry, bianchiChecks}
];

(* Maxwell-Faraday checks: each component is the standard cyclic identity
   ∂_t B_k + ε_{klm} ∂_l E_m = 0 with E_i = -F_{0i} and B_1=F_{23}, B_2=F_{31}, B_3=F_{12}.
   The original Module/Do/Part-indexing pattern lost half its terms to
   pre-evaluation; precomputing each component as an immediate expression
   from already-evaluated F-components avoids the trap. *)
mf1 = D[F23, t] + D[-F03, y] - D[-F02, z];
mf2 = D[F31, t] + D[-F01, z] - D[-F03, x];
mf3 = D[F12, t] + D[-F02, x] - D[-F01, y];

maxwellFaradayChecks = {{1, mf1}, {2, mf2}, {3, mf3}};
Do[Module[{compIdx, expr, residual},
    {compIdx, expr} = entry;
    residual = FullSimplify[expr, Assumptions -> spaceTimeAssumptions];
    Print["M2 Maxwell-Faraday component ", compIdx, " residual = ", fmt[residual]];
    If[residual =!= 0,
      Print["FAIL: M2 Maxwell-Faraday component ", compIdx]; Exit[1]
    ];
    Print["PASS: M2 Maxwell-Faraday component ", compIdx];
  ],
  {entry, maxwellFaradayChecks}
];

(* M3: Gaussian normalization. *)
localizedProfile = Exp[-w^2/lam^2];
profileArea =
  Integrate[localizedProfile, {w, -Infinity, Infinity},
    Assumptions -> scaleAssumptions];
m3Target = Sqrt[Pi]*lam;
m3Residual = FullSimplify[profileArea - m3Target, Assumptions -> scaleAssumptions];
Print["M3 residual = ", fmt[m3Residual]];
If[FullSimplify[profileArea - m3Target, Assumptions -> scaleAssumptions] =!= 0,
  Print["FAIL: M3 Gaussian normalization"]; Exit[1]
];
Print["PASS: M3 Gaussian normalization"];

(* M4: Gaussian squared norm. *)
profileSelfMass =
  Integrate[localizedProfile^2, {w, -Infinity, Infinity},
    Assumptions -> scaleAssumptions];
m4Target = Sqrt[2*Pi]*lam/2;
m4Residual = FullSimplify[profileSelfMass - m4Target, Assumptions -> scaleAssumptions];
Print["M4 residual = ", fmt[m4Residual]];
If[FullSimplify[profileSelfMass - m4Target, Assumptions -> scaleAssumptions] =!= 0,
  Print["FAIL: M4 Gaussian squared norm"]; Exit[1]
];
Print["PASS: M4 Gaussian squared norm"];

(* M5: matched-kernel overlap. *)
matchedWeight = localizedProfile/profileArea;
overlapValue =
  Integrate[matchedWeight*localizedProfile, {w, -Infinity, Infinity},
    Assumptions -> scaleAssumptions];
m5Target = Sqrt[2]/2;
m5Residual = FullSimplify[overlapValue - m5Target, Assumptions -> scaleAssumptions];
Print["M5 residual = ", fmt[m5Residual]];
If[FullSimplify[overlapValue - m5Target, Assumptions -> scaleAssumptions] =!= 0,
  Print["FAIL: M5 matched-kernel overlap"]; Exit[1]
];
Print["PASS: M5 matched-kernel overlap"];

(* M6: delta-source projection compared with reduced coupling. *)
pointSourceCoupling = mu0*(matchedWeight /. w -> 0)/overlapValue;
volumeReducedCoupling = mu0/profileArea;
m6Target = Sqrt[2];
m6Residual =
  FullSimplify[pointSourceCoupling/volumeReducedCoupling - m6Target,
    Assumptions -> couplingAssumptions];
Print["M6 residual = ", fmt[m6Residual]];
If[FullSimplify[pointSourceCoupling/volumeReducedCoupling - m6Target,
    Assumptions -> couplingAssumptions] =!= 0,
  Print["FAIL: M6 delta-source projection/reduction ratio"]; Exit[1]
];
Print["PASS: M6 delta-source projection/reduction ratio"];

Print["STATUS: PASS"];
Exit[0];
