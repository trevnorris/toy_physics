(*
========================================================
Paper VIII (referee add-on v2): KK spectrum, Coulomb law, and propagator
for the localized Maxwell sector with Gaussian profile Z(w)

Changes vs v1:
  - Pass Assumptions directly into Integrate where needed.
  - Avoid symbolic-parameter Integrate for norm[n] (which can mis-diagnose
    convergence for non-integer n). Instead:
        (i) state the standard closed form normClosed[n],
        (ii) spot-check normClosed[n] by direct integration for many integer n.
  - Verify the even-mode closed form c_{2k} by spot-checks (k=0..K),
    and keep the analytic expression used downstream.
========================================================
*)

Print["\n========================================"];
Print["Paper VIII add-on v2: KK spectrum + Coulomb + propagator (Gaussian Z)"];
Print["========================================\n"];

ClearAll["Global`*"];

(* ---------- Formatting helpers ---------- *)
PrintHeader[title_String] := Print[
  "\n" <> StringJoin[ConstantArray["=", 60]] <> "\n" <>
  title <> "\n" <>
  StringJoin[ConstantArray["-", 60]]
];

PrintEq[label_String, expr_] := (Print[""]; Print[label <> ":"]; Print[OutputForm[expr]]);

zeroQ[expr_, asm_] := TrueQ[FullSimplify[expr == 0, Assumptions -> asm]];

(* ---------- Assumptions ---------- *)
asmBase = {
  Element[{mu0, lambdaConf}, Reals],
  mu0 > 0, lambdaConf > 0
};

(* ---------- Localization profile ---------- *)
Zprof[ww_] := Exp[-ww^2/lambdaConf^2];

PrintHeader["Profile and basic integrals"];
Zint = FullSimplify[
  Integrate[Zprof[w], {w, -Infinity, Infinity}, Assumptions -> asmBase],
  Assumptions -> asmBase
];
mu0eff = FullSimplify[mu0/Zint, Assumptions -> asmBase];

PrintEq["Zint = Integrate[Z(w) dw]", Zint];
PrintEq["mu0eff = mu0/Zint", mu0eff];

(* ---------- Sturm-Liouville eigenproblem (Gaussian Z) ----------
   KK separation gives w-profiles f_n(w) satisfying:
      -d/dw( Z(w) f'(w) ) = m_n^2 Z(w) f(w)
*)

PrintHeader["Sturm-Liouville eigenproblem (Gaussian Z)"];

(* Operator L[f] := -(1/Z) d/dw( Z f' ) *)
Lop[f_] := FullSimplify[
  -(1/Zprof[w]) * D[Zprof[w] * D[f, w], w],
  Assumptions -> asmBase
];

(* Candidate eigenfunctions and eigenvalues *)
m2[n_] := FullSimplify[(2*n)/lambdaConf^2, Assumptions -> asmBase];
fn[n_Integer?NonNegative][ww_] := HermiteH[n, ww/lambdaConf];

(* Verify eigen-equation for symbolic integer n is tricky; do a clean spot-check. *)
Print["Eigen-residual spot-checks (should all be 0):"];
Table[
  With[{nn = nTest},
    res = FullSimplify[Lop[fn[nn][w]] - m2[nn]*fn[nn][w], Assumptions -> asmBase];
    Print["  n=", nn, " -> ", OutputForm[res]];
  ],
  {nTest, 0, 8}
];

(* ---------- Norms and brane couplings ---------- *)
PrintHeader["Weighted norms and brane couplings"];

(*
For integer n >= 0 (physicists' Hermite polynomials),
the standard orthogonality gives:

  ∫_{-∞}^{∞} exp(-x^2) H_n(x)^2 dx = 2^n n! sqrt(pi).

With x = w/lambdaConf and dw = lambdaConf dx, and Z(w)=exp(-(w/lambdaConf)^2),
we get:

  normClosed[n] := ∫ dw Z(w) f_n(w)^2 = lambdaConf * 2^n n! sqrt(pi).
*)

normClosed[n_Integer?NonNegative] := lambdaConf*Sqrt[Pi]*2^n*Factorial[n];

PrintEq["normClosed[n] (standard Hermite orthogonality result)", normClosed[n]];

(* Spot-check normClosed by direct integration (no symbolic-parameter pitfalls). *)
Print["\nSpot-check norms via direct integration (integer n):"];
Table[
  Module[{nn = nTest, integ, diff},
    integ = Quiet@FullSimplify[
      Integrate[Zprof[w] * (fn[nn][w])^2, {w, -Infinity, Infinity},
        Assumptions -> asmBase
      ],
      Assumptions -> asmBase
    ];
    diff = FullSimplify[integ - normClosed[nn], Assumptions -> asmBase];
    Print[
      "  n=", nn,
      "  integral=", OutputForm[integ],
      "  diff=", OutputForm[diff]
    ];
    If[!zeroQ[diff, asmBase],
      Print["  WARNING: norm mismatch at n=", nn, " (unexpected)."]
    ];
  ],
  {nTest, 0, 10}
];

Print["\nOK if all diffs above are 0: normClosed[n] is verified for these n."];

(* Brane coupling for a delta-localized brane source at w=0 is proportional to f_n(0). *)
f0[n_Integer?NonNegative] := fn[n][0];

PrintEq["f_n(0) = HermiteH[n,0] (odd n vanish)", f0[n]];

(* Spectral weight for brane-to-brane exchange: c_n := f_n(0)^2 / normClosed[n] *)
cCoef[n_Integer?NonNegative] := FullSimplify[(f0[n]^2)/normClosed[n], Assumptions -> asmBase];

PrintEq["cCoef[n] := f_n(0)^2 / normClosed[n]", cCoef[n]];

PrintEq["Check cCoef[0] - 1/Zint (should be 0)", FullSimplify[cCoef[0] - 1/Zint, Assumptions -> asmBase]];

(* Show first few coefficients explicitly *)
Print["\nFirst few (n, m_n, c_n) values:"];
Table[
  Print[
    "n=", nn,
    "   m_n=", OutputForm[FullSimplify[Sqrt[m2[nn]], Assumptions -> asmBase]],
    "   c_n=", OutputForm[FullSimplify[cCoef[nn], Assumptions -> asmBase]]
  ],
  {nn, 0, 10}
];

(* Closed form for even modes n = 2k:
     H_{2k}(0) = (-1)^k (2k)!/k!
   => c_{2k} = Binomial[2k,k]/(2^(2k) lambdaConf sqrt(pi)),  and c_{odd}=0.
*)
cEvenClosed[k_Integer?NonNegative] := Binomial[2*k, k]/(2^(2*k) * lambdaConf * Sqrt[Pi]);

PrintEq["Closed form: c_{2k} = Binomial[2k,k]/(2^(2k) lambdaConf Sqrt[Pi])", cEvenClosed[k]];

Print["\nSpot-check c_{2k} closed form (k=0..8):"];
Table[
  Module[{kk = kTest, lhs, rhs, diff},
    lhs = cCoef[2*kk];
    rhs = cEvenClosed[kk];
    diff = FullSimplify[lhs - rhs, Assumptions -> asmBase];
    Print[
      "  k=", kk,
      "  cCoef[", 2*kk, "]=", OutputForm[lhs],
      "  closed=", OutputForm[rhs],
      "  diff=", OutputForm[diff]
    ];
    If[!zeroQ[diff, asmBase],
      Print["  WARNING: even-mode coefficient mismatch at k=", kk, " (unexpected)."]
    ];
  ],
  {kTest, 0, 8}
];

(* ---------- Brane Coulomb law + Yukawa corrections (static exchange picture) ---------- *)
PrintHeader["Brane potential: Coulomb + Yukawa tower"];

(*
For a static point charge on the brane (w=0), the brane-observed potential has
the spectral representation:

  A0(r) = (mu0*q)/(4 Pi r) * Sum_{n>=0} c_n * Exp[-m_n r]

Zero mode (n=0) gives:
  A0_0(r) = (mu0*q)/(4 Pi r) * (1/Zint) = (mu0eff*q)/(4 Pi r)

Even-n modes correct this by Yukawa terms with range ~ 1/m_n.
Odd-n modes do not couple for a symmetric brane source because f_{odd}(0)=0.
*)

ClearAll[r, q];
asmR = Join[asmBase, {Element[{r, q}, Reals], r > 0}];

A0zero[r_] := FullSimplify[(mu0*q)/(4*Pi*r) * cCoef[0], Assumptions -> asmR];
PrintEq["A0_zero(r) = (mu0 q)/(4 Pi r) * c0 = (mu0eff q)/(4 Pi r)", A0zero[r]];

(* Correction ratio Delta(r) := (A0 - A0_zero)/A0_zero, truncated at kmax even modes *)
corrRatioTrunc[r_, kmax_Integer?NonNegative] := FullSimplify[
  Sum[
    (cCoef[2*kk]/cCoef[0]) * Exp[-Sqrt[m2[2*kk]]*r],
    {kk, 1, kmax}
  ],
  Assumptions -> asmR
];

Print["\nCorrection ratio Delta(r) := (A0 - A0_zero)/A0_zero  (finite truncations):"];
Table[
  Print["kmax=", km, "  Delta(r) = ", OutputForm[corrRatioTrunc[r, km]]],
  {km, 1, 4}
];

Print["\nLeading correction (k=1 => n=2) has coefficient 1/2 and mass m=2/lambdaConf:"];
leadCoeff = FullSimplify[(cCoef[2]/cCoef[0]), Assumptions -> asmR];
leadMass = FullSimplify[Sqrt[m2[2]], Assumptions -> asmR];
PrintEq["(c2/c0)", leadCoeff];
PrintEq["m2", leadMass];
PrintEq["Leading Delta(r) approx", FullSimplify[leadCoeff*Exp[-leadMass*r], Assumptions -> asmR]];

(* ---------- Brane-to-brane propagator (momentum space; Lorentz scalar) ---------- *)
PrintHeader["Brane-to-brane propagator depends only on k^2"];

(*
Mode-expanded Euclidean propagator (schematic):

  D_eff(k2) = mu0 * Sum_{n>=0} c_n / (k2 + m_n^2)

Each term depends on k2 and a scalar mass m_n^2. So the Maxwell sector induces
no direction-dependent anisotropy in the brane sector: Deff is a function of k2 only.

Any preferred-frame effect would have to enter through the source/matter sector or
time-dependent backgrounds, not from this localized Maxwell kernel itself.
*)

ClearAll[k2];
asmK = Join[asmBase, {Element[k2, Reals], k2 >= 0}];

DeffTrunc[k2_, nmax_Integer?NonNegative] := FullSimplify[
  mu0 * Sum[cCoef[nn]/(k2 + m2[nn]), {nn, 0, nmax}],
  Assumptions -> asmK
];

Print["Truncated Euclidean propagator Deff(k2) examples:"];
Table[
  Print["nmax=", nm, "  Deff(k2) = ", OutputForm[DeffTrunc[k2, nm]]],
  {nm, 0, 8}
];

Print["\nNote: Deff depends on k2 only (a Lorentz scalar in the brane sector)."];

Print["\n========== End KK/Coulomb/propagator add-on v2 ==========\n"];

(*"
Output:

========================================
Paper VIII add-on v2: KK spectrum + Coulomb + propagator (Gaussian Z)
========================================


============================================================
Profile and basic integrals
------------------------------------------------------------

Zint = Integrate[Z(w) dw]:
OutputForm[lambdaConf*Sqrt[Pi]]

mu0eff = mu0/Zint:
OutputForm[mu0/(lambdaConf*Sqrt[Pi])]

============================================================
Sturm-Liouville eigenproblem (Gaussian Z)
------------------------------------------------------------
Eigen-residual spot-checks (should all be 0):
  n=0 -> OutputForm[0]
  n=1 -> OutputForm[0]
  n=2 -> OutputForm[0]
  n=3 -> OutputForm[0]
  n=4 -> OutputForm[0]
  n=5 -> OutputForm[0]
  n=6 -> OutputForm[0]
  n=7 -> OutputForm[0]
  n=8 -> OutputForm[0]

============================================================
Weighted norms and brane couplings
------------------------------------------------------------

normClosed[n] (standard Hermite orthogonality result):
OutputForm[normClosed[n]]

Spot-check norms via direct integration (integer n):
  n=0  integral=OutputForm[lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]
  n=1  integral=OutputForm[2*lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]
  n=2  integral=OutputForm[8*lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]
  n=3  integral=OutputForm[48*lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]
  n=4  integral=OutputForm[384*lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]
  n=5  integral=OutputForm[3840*lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]
  n=6  integral=OutputForm[46080*lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]
  n=7  integral=OutputForm[645120*lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]
  n=8  integral=OutputForm[10321920*lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]
  n=9  integral=OutputForm[185794560*lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]
  n=10  integral=OutputForm[3715891200*lambdaConf*Sqrt[Pi]]  diff=OutputForm[0]

OK if all diffs above are 0: normClosed[n] is verified for these n.

f_n(0) = HermiteH[n,0] (odd n vanish):
OutputForm[f0[n]]

cCoef[n] := f_n(0)^2 / normClosed[n]:
OutputForm[cCoef[n]]

Check cCoef[0] - 1/Zint (should be 0):
OutputForm[0]

First few (n, m_n, c_n) values:
n=0   m_n=OutputForm[0]   c_n=OutputForm[1/(lambdaConf*Sqrt[Pi])]
n=1   m_n=OutputForm[Sqrt[2]/lambdaConf]   c_n=OutputForm[0]
n=2   m_n=OutputForm[2/lambdaConf]   c_n=OutputForm[1/(2*lambdaConf*Sqrt[Pi])]
n=3   m_n=OutputForm[Sqrt[6]/lambdaConf]   c_n=OutputForm[0]
n=4   m_n=OutputForm[(2*Sqrt[2])/lambdaConf]   c_n=OutputForm[3/(8*lambdaConf*Sqrt[Pi])]
n=5   m_n=OutputForm[Sqrt[10]/lambdaConf]   c_n=OutputForm[0]
n=6   m_n=OutputForm[(2*Sqrt[3])/lambdaConf]   c_n=OutputForm[5/(16*lambdaConf*Sqrt[Pi])]
n=7   m_n=OutputForm[Sqrt[14]/lambdaConf]   c_n=OutputForm[0]
n=8   m_n=OutputForm[4/lambdaConf]   c_n=OutputForm[35/(128*lambdaConf*Sqrt[Pi])]
n=9   m_n=OutputForm[(3*Sqrt[2])/lambdaConf]   c_n=OutputForm[0]
n=10   m_n=OutputForm[(2*Sqrt[5])/lambdaConf]   c_n=OutputForm[63/(256*lambdaConf*Sqrt[Pi])]

Closed form: c_{2k} = Binomial[2k,k]/(2^(2k) lambdaConf Sqrt[Pi]):
OutputForm[cEvenClosed[k]]

Spot-check c_{2k} closed form (k=0..8):
  k=0  cCoef[0]=OutputForm[1/(lambdaConf*Sqrt[Pi])]  closed=OutputForm[1/(lambdaConf*Sqrt[Pi])]  diff=OutputForm[0]
  k=1  cCoef[2]=OutputForm[1/(2*lambdaConf*Sqrt[Pi])]  closed=OutputForm[1/(2*lambdaConf*Sqrt[Pi])]  diff=OutputForm[0]
  k=2  cCoef[4]=OutputForm[3/(8*lambdaConf*Sqrt[Pi])]  closed=OutputForm[3/(8*lambdaConf*Sqrt[Pi])]  diff=OutputForm[0]
  k=3  cCoef[6]=OutputForm[5/(16*lambdaConf*Sqrt[Pi])]  closed=OutputForm[5/(16*lambdaConf*Sqrt[Pi])]  diff=OutputForm[0]
  k=4  cCoef[8]=OutputForm[35/(128*lambdaConf*Sqrt[Pi])]  closed=OutputForm[35/(128*lambdaConf*Sqrt[Pi])]  diff=OutputForm[0]
  k=5  cCoef[10]=OutputForm[63/(256*lambdaConf*Sqrt[Pi])]  closed=OutputForm[63/(256*lambdaConf*Sqrt[Pi])]  diff=OutputForm[0]
  k=6  cCoef[12]=OutputForm[231/(1024*lambdaConf*Sqrt[Pi])]  closed=OutputForm[231/(1024*lambdaConf*Sqrt[Pi])]  diff=OutputForm[0]
  k=7  cCoef[14]=OutputForm[429/(2048*lambdaConf*Sqrt[Pi])]  closed=OutputForm[429/(2048*lambdaConf*Sqrt[Pi])]  diff=OutputForm[0]
  k=8  cCoef[16]=OutputForm[6435/(32768*lambdaConf*Sqrt[Pi])]  closed=OutputForm[6435/(32768*lambdaConf*Sqrt[Pi])]  diff=OutputForm[0]

============================================================
Brane potential: Coulomb + Yukawa tower
------------------------------------------------------------

A0_zero(r) = (mu0 q)/(4 Pi r) * c0 = (mu0eff q)/(4 Pi r):
OutputForm[(mu0*q)/(4*lambdaConf*Pi^(3/2)*r)]

Correction ratio Delta(r) := (A0 - A0_zero)/A0_zero  (finite truncations):
kmax=1  Delta(r) = OutputForm[1/(2*E^((2*r)/lambdaConf))]
kmax=2  Delta(r) = OutputForm[1/(2*E^((2*r)/lambdaConf)) + 3/(8*E^((2*Sqrt[2]*r)/lambdaConf))]
kmax=3  Delta(r) = OutputForm[(8/E^((2*r)/lambdaConf) + 6/E^((2*Sqrt[2]*r)/lambdaConf) + 5/E^((2*Sqrt[3]*r)/lambdaConf))/16]
kmax=4  Delta(r) = OutputForm[(35/E^((4*r)/lambdaConf) + 64/E^((2*r)/lambdaConf) + 48/E^((2*Sqrt[2]*r)/lambdaConf) + 40/E^((2*Sqrt[3]*r)/lambdaConf))/128]

Leading correction (k=1 => n=2) has coefficient 1/2 and mass m=2/lambdaConf:

(c2/c0):
OutputForm[1/2]

m2:
OutputForm[2/lambdaConf]

Leading Delta(r) approx:
OutputForm[1/(2*E^((2*r)/lambdaConf))]

============================================================
Brane-to-brane propagator depends only on k^2
------------------------------------------------------------
Truncated Euclidean propagator Deff(k2) examples:
nmax=0  Deff(k2) = OutputForm[mu0/(k2*lambdaConf*Sqrt[Pi])]
nmax=1  Deff(k2) = OutputForm[mu0/(k2*lambdaConf*Sqrt[Pi])]
nmax=2  Deff(k2) = OutputForm[((8 + 3*k2*lambdaConf^2)*mu0)/(2*k2*lambdaConf*(4 + k2*lambdaConf^2)*Sqrt[Pi])]
nmax=3  Deff(k2) = OutputForm[((8 + 3*k2*lambdaConf^2)*mu0)/(2*k2*lambdaConf*(4 + k2*lambdaConf^2)*Sqrt[Pi])]
nmax=4  Deff(k2) = OutputForm[((256 + 5*k2*lambdaConf^2*(28 + 3*k2*lambdaConf^2))*mu0)/(8*k2*lambdaConf*(4 + k2*lambdaConf^2)*(8 + k2*lambdaConf^2)*Sqrt[Pi])]
nmax=5  Deff(k2) = OutputForm[((256 + 5*k2*lambdaConf^2*(28 + 3*k2*lambdaConf^2))*mu0)/(8*k2*lambdaConf*(4 + k2*lambdaConf^2)*(8 + k2*lambdaConf^2)*Sqrt[Pi])]
nmax=6  Deff(k2) = OutputForm[((6144 + 7*k2*lambdaConf^2*(576 + 5*k2*lambdaConf^2*(20 + k2*lambdaConf^2)))*mu0)/(16*k2*lambdaConf*(4 + k2*lambdaConf^2)*(8 + k2*lambdaConf^2)*(12 + k2*lambdaConf^2)*Sqrt[Pi])]
nmax=7  Deff(k2) = OutputForm[((6144 + 7*k2*lambdaConf^2*(576 + 5*k2*lambdaConf^2*(20 + k2*lambdaConf^2)))*mu0)/(16*k2*lambdaConf*(4 + k2*lambdaConf^2)*(8 + k2*lambdaConf^2)*(12 + k2*lambdaConf^2)*Sqrt[Pi])]
nmax=8  Deff(k2) = OutputForm[(3*(262144 + k2*lambdaConf^2*(192896 + 7*k2*lambdaConf^2*(6096 + 5*k2*lambdaConf^2*(104 + 3*k2*lambdaConf^2))))*mu0)/(128*k2*lambdaConf*(4 + k2*lambdaConf^2)*(8 + k2*lambdaConf^2)*(12 + k2*lambdaConf^2)*(16 + k2*lambdaConf^2)*Sqrt[Pi])]

Note: Deff depends on k2 only (a Lorentz scalar in the brane sector).

========== End KK/Coulomb/propagator add-on v2 ==========
"*)
