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

expectTrue[name_String, cond_] := Module[{res},
  res = FullSimplify[cond, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res], pass[name], fail[name, res]];
];

banner["STAGE 088 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE"];

(* Independent Mathematica derivation: start from the paper Y_Q^cons input
   directly and use Limit + subtraction to extract (c0, c1). This is a
   different algebraic path than the SymPy script (which works in the
   rho_alpha-parameterized form first); both engines arrive at the same
   numerical answer. *)

Clear[rhoAlpha, omega, omegaQ, u, cMix];
$Assumptions =
  Element[{rhoAlpha, omega, omegaQ, u, cMix}, Reals] &&
  rhoAlpha > 0 && omegaQ > 0 && cMix > 0 && -1 < u < 1;

(* Paper Input: Y_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2). *)
yQpaper = 3/4 + (1/4)/(1 - omega^2/omegaQ^2);
yQpaperU = yQpaper /. omega^2/omegaQ^2 -> u;

(* Independent extraction: c1 from pole residue at u = 1, then c0 by
   subtraction at u = 0. This is a single non-tautological probe of the
   paper form. *)
c1Paper = FullSimplify[Limit[(1 - u)*yQpaperU, u -> 1]];
c0Paper = FullSimplify[(yQpaperU /. u -> 0) - c1Paper];

Print["Y_Q^cons (paper)   = ", fmt[yQpaper]];
Print["c0 (subtract-pole) = ", fmt[c0Paper]];
Print["c1 (pole residue)  = ", fmt[c1Paper]];

expectZero["c0_paper - 3/4", c0Paper - 3/4];
expectZero["c1_paper - 1/4", c1Paper - 1/4];
expectZero["c0_paper + c1_paper - 1", c0Paper + c1Paper - 1];

(* Loading-ratio extraction from coefficients: rho_alpha = 1/c0, zeta = c1/c0. *)
rhoMin = FullSimplify[1/c0Paper];
zetaMin = FullSimplify[c1Paper/c0Paper];

Print["rho_alpha (= 1/c0) = ", fmt[rhoMin]];
Print["zeta_req (= c1/c0) = ", fmt[zetaMin]];

expectZero["rho_min - 4/3", rhoMin - 4/3];
expectZero["zeta_min - 1/3", zetaMin - 1/3];

(* Reconstruct the contact-plus-pole form from extracted (c0, c1) and confirm
   it matches the paper precursor. This is the actual coefficient-matching
   claim of the stage. *)
yRhoFromCoeffs = c0Paper + c1Paper/(1 - omega^2/omegaQ^2);
expectZero["paper form - reconstruction from extracted (c0, c1)",
           yQpaper - yRhoFromCoeffs];

(* Also confirm the general rho-parameterized form rebuilds yQpaper after
   the substitution rhoAlpha -> rhoMin. *)
yRhoParam = 1/rhoAlpha + ((rhoAlpha - 1)/rhoAlpha)/(1 - omega^2/omegaQ^2);
expectZero["rho-parameterized form (rhoAlpha -> rhoMin) - paper form",
           (yRhoParam /. rhoAlpha -> rhoMin) - yQpaper];

(* Stage-085 product identity: Pi_tr = rho_alpha * C_mix (verified upstream
   in the stage 085 Mathematica audit files). Substitute rho_min. *)
piFromRho = FullSimplify[rhoMin*cMix];
Print["Pi_tr (= rho_min * C_mix) = ", fmt[piFromRho]];
expectZero["Pi_tr_from_rho - (4/3) C_mix", piFromRho - (4/3)*cMix];
expectTrue["1 < rho_min < 2 (symmetric-lowest-twin regime)", 1 < rhoMin < 2];

Print[""];
Print["Stage 088 Mathematica audit passed."];

Exit[0];
