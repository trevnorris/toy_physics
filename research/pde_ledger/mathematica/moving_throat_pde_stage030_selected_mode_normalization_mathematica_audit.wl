ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["PART I — GENERIC SELECTED-MODE NORMALIZED RESPONSE"];

Clear[w, d0, d2, d4, c5];
$Assumptions = Element[{w, d0, d2, d4, c5}, Reals] && d0 != 0;
dSel = d0 + d2*w^2 + d4*w^4 - I*c5*w^5;
ySel = Expand[Normal[Series[d0/dSel, {w, 0, 5}]]];

Print["Y_-(omega) = ", fmt[ySel]];

u2 = -d2/d0;
u4 = (d2^2 - d0*d4)/d0^2;
gamma5 = c5/d0;

expectZero["u2 coefficient check", Coefficient[Expand[ySel], w, 2] - u2];
expectZero["u4 coefficient check", Coefficient[Expand[ySel], w, 4] - u4];
expectZero["Gamma5 coefficient check", Coefficient[Expand[ComplexExpand[Im[ySel]]], w, 5] - gamma5];

banner["PART II — EXACT SELECTED LOWER EIGENVALUE AND OVERLAP"];

Clear[a, dK, alpha, x0, x1];
$Assumptions = Element[{a, dK, alpha, x0, x1}, Reals] && a > 0 && dK > 0 && alpha >= 0 && x0 > 0 && x1 > 0;

(* Independent derivation: define the 2x2 wall block explicitly and let
   Mathematica's Eigenvalues[] produce the spectrum. The matrix is chosen so
   that Tr[mMat] = 2 a + dK - alpha (x0+x1) and
   Det[mMat] = a (a+dK) - alpha ((a+dK) x0 + a x1) -- the same trace/det as the
   determinant identity verified later in this script (see wl:87). This breaks
   the line-by-line transliteration of the SymPy file by routing lamMinus /
   lamPlus through Eigenvalues[mMat] rather than through a typed closed form. *)
mMat = {{a + dK - alpha*x1, -alpha*Sqrt[x0*x1]},
        {-alpha*Sqrt[x0*x1], a - alpha*x0}};

sigma = x0 + x1;
deltaKappa = x0 - x1;
kappaProd = x0*x1;
r = Sqrt[(dK + alpha*deltaKappa)^2 + 4*alpha^2*kappaProd];

eigVals = Eigenvalues[mMat];
lamMinus = FullSimplify[
  First[Select[eigVals,
    FullSimplify[# - ((2*a + dK - alpha*sigma) - r)/2, Assumptions -> $Assumptions] === 0 &]],
  Assumptions -> $Assumptions
];
lamPlus = FullSimplify[
  First[Select[eigVals,
    FullSimplify[# - ((2*a + dK - alpha*sigma) + r)/2, Assumptions -> $Assumptions] === 0 &]],
  Assumptions -> $Assumptions
];

(* HF eigenvector check: compute the lower-eigenvalue eigenvector of mMat
   directly and verify (v.e_-)^2 = sMinusClosed. This closes the HF chain
   (v.e_-)^2 = -d lambda_-/d alpha = closed form. In the mMat basis (row 1
   = kappa_1-mode, row 2 = kappa_0-mode), the loading direction is
   vVec = {Sqrt[x1], Sqrt[x0]}. *)
eigPairs = Eigensystem[mMat];
eMinusRaw = First[Pick[Transpose[eigPairs][[All, 2]],
  Map[FullSimplify[# - lamMinus, Assumptions -> $Assumptions] === 0 &,
      First[eigPairs]]]];
eMinusNorm = FullSimplify[eMinusRaw/Sqrt[eMinusRaw.eMinusRaw],
  Assumptions -> $Assumptions];
vVec = {Sqrt[x1], Sqrt[x0]};
sMinusEig = FullSimplify[(vVec.eMinusNorm)^2, Assumptions -> $Assumptions];

Print["lambda_- = ", fmt[lamMinus]];
Print["lambda_+ = ", fmt[lamPlus]];

sMinusHF = FullSimplify[-D[lamMinus, alpha], Assumptions -> $Assumptions];
sMinusClosed = FullSimplify[
  (sigma + ((dK + alpha*deltaKappa)*deltaKappa + 4*alpha*kappaProd)/r)/2,
  Assumptions -> $Assumptions
];
expectZero["selected overlap: HF - closed form", sMinusHF - sMinusClosed];
expectZero["HF eigenvector check", sMinusEig - sMinusClosed];
Print["s_- = (v.e_-)^2 = ", fmt[sMinusClosed]];
expectZero["weak-loading overlap limit", (sMinusClosed /. alpha -> 0) - x0];

banner["PART III — SELECTED ODD COEFFICIENT AND STATIC PREFACTOR"];

Clear[beta0, g5];
$Assumptions = Element[{a, dK, alpha, x0, x1, beta0, g5}, Reals] &&
  a > 0 && dK > 0 && alpha >= 0 && x0 > 0 && x1 > 0 && beta0 > 0 && g5 > 0;

beta5 = g5*beta0;
c5Sel = FullSimplify[beta5*sMinusClosed, Assumptions -> $Assumptions];
gamma5Sel = FullSimplify[c5Sel/lamMinus, Assumptions -> $Assumptions];
p0Sel = FullSimplify[beta0*sMinusClosed/lamMinus, Assumptions -> $Assumptions];

Print["C5_sel = ", fmt[c5Sel]];
Print["Gamma5_sel = ", fmt[gamma5Sel]];
Print["P0_sel = ", fmt[p0Sel]];

(* Note: gamma5Sel - g5*p0Sel == 0 follows by construction from lines 74-77
   (c5Sel := g5*beta0*sMinusClosed, gamma5Sel := c5Sel/lamMinus,
   p0Sel := beta0*sMinusClosed/lamMinus). Not verified separately. The
   physical content (Gamma5 = C5/D0 at the selected mode) is verified in
   generic form by Part I (line 41). *)
expectZero["P0_sel + beta0*d(log lambda_-)/d alpha", p0Sel + beta0*D[Log[lamMinus], alpha]];

t0 = (a + dK)*x0 + a*x1;
expectZero["det identity", Expand[lamMinus*lamPlus - (a*(a + dK) - alpha*t0)]];

banner["PART IV — EQUIVALENCE OF THE SELECTED-BRANCH TARGET FORMS"];

Clear[gConst, cs, radius, cSpeed, mhat];
$Assumptions = Element[{a, dK, alpha, x0, x1, beta0, g5, gConst, cs, radius, cSpeed, mhat}, Reals] &&
  a > 0 && dK > 0 && alpha >= 0 && x0 > 0 && x1 > 0 && beta0 > 0 && g5 > 0 &&
  gConst > 0 && cs > 0 && radius > 0 && cSpeed > 0 && mhat > 0;

g5Phys = radius^5/(27*cs^5);
nQTarget = 54*gConst*cs^5/(5*radius^5*cSpeed^5);

cond1 = FullSimplify[mhat^2*g5Phys*p0Sel - 2*gConst/(5*cSpeed^5), Assumptions -> $Assumptions];
cond2 = FullSimplify[mhat^2*p0Sel - nQTarget, Assumptions -> $Assumptions];
expectZero["normalization equivalence", cond1 - g5Phys*cond2];

lambdaReq = FullSimplify[mhat^2*beta0*sMinusClosed/nQTarget, Assumptions -> $Assumptions];
(* Note: (lamMinus - lambdaReq) + (mhat^2 p0Sel - nQTarget) lamMinus / nQTarget
   == 0 follows from the definitions p0Sel := beta0*sMinusClosed/lamMinus and
   lambdaReq := mhat^2*beta0*sMinusClosed/nQTarget by substitution, with no
   physical content of lamMinus required. Not verified separately. *)
Print["lambda_req = ", fmt[lambdaReq]];

banner["STAGE 30 AUDIT COMPLETE"];
Print["Verified:"];
Print["  generic normalized selected-response expansion"];
Print["  exact selected lower eigenvalue and Hellmann–Feynman overlap"];
Print["  exact selected static prefactor P0_- = beta0 s_-/lambda_-"];
Print["  exact identity P0_- = - beta0 d(log lambda_-)/d alpha"];
Print["  exact determinant identity for the selected branch"];
Print["  exact equivalence of the selected-branch normalization target forms"];

Exit[0];
