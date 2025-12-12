(* ============================================================
   verify_w_zero_mode.wl

   Purpose:
     Analyze the transverse (bulk) eigenproblem in w:

       -χ''(w) = λ χ(w),   on w ∈ [0, Lw]
       with λ = k_w^2

     and determine whether a true zero mode exists (k_w = 0).

   Why this matters (ties to your prior script):
     If photon propagation along the brane has dispersion
       ω = c_s(r) Sqrt[k^2 + k_w^2],
     then k_w ≠ 0 generally introduces (k_w/k)^2 corrections
     (chromaticity / cutoff-like behavior).
     A genuine zero mode (k_w = 0) gives a clean, gapless branch.

   This script:
     1) Prints eigenvalue spectra for common BCs (analytic).
     2) Flags whether k_w = 0 is allowed.
     3) Includes a simple Robin-BC zero-mode test (k_w=0 case treated properly).

   ============================================================ *)

ClearAll["Global`*"];

(* ---- User knobs ---- *)
LwVal = 1;      (* choose a convenient normalization; change as desired *)
mMax  = 6;      (* how many modes to print *)

(* ---- Helper: pretty printer for a mode list ---- *)
printModes[name_, kwList_] := Module[{},
  Print["\n--- ", name, " ---"];
  Print["k_w values (first ", Length[kwList], "): ", kwList // N];
  Print["λ = k_w^2 (first ", Length[kwList], "): ", (kwList^2) // N];
  Print["Zero mode present?  ", If[MemberQ[kwList, 0], "YES (k_w=0)", "NO"]];
];

(* ============================================================
   1) Analytic spectra for standard BC choices
   ============================================================ *)

(* Dirichlet-Dirichlet: χ(0)=0, χ(L)=0  =>  k_w = m π/L, m=1,2,... (NO zero mode) *)
kwDirichlet[L_] := Table[m Pi/L, {m, 1, mMax}];

(* Neumann-Neumann: χ'(0)=0, χ'(L)=0  =>  k_w = m π/L, m=0,1,2,... (YES zero mode) *)
kwNeumann[L_] := Table[m Pi/L, {m, 0, mMax - 1}];

(* Periodic: χ(0)=χ(L), χ'(0)=χ'(L)  =>  k_w = 2π m/L, m=0,1,2,... (YES zero mode) *)
kwPeriodic[L_] := Table[2 Pi m/L, {m, 0, mMax - 1}];

(* Mixed Dirichlet-Neumann: χ(0)=0, χ'(L)=0 => k_w=(m+1/2)π/L (NO zero mode) *)
kwMixedDN[L_] := Table[(m + 1/2) Pi/L, {m, 0, mMax - 1}];

(* Mixed Neumann-Dirichlet: χ'(0)=0, χ(L)=0 => k_w=(m+1/2)π/L (NO zero mode) *)
kwMixedND[L_] := Table[(m + 1/2) Pi/L, {m, 0, mMax - 1}];

Print["\n--- Setup ---"];
Print["Eigenproblem:  -χ''(w) = λ χ(w),  λ = k_w^2,  w ∈ [0, Lw]"];
Print["Using Lw = ", LwVal, ", printing up to mMax = ", mMax, " modes."];

printModes["Dirichlet-Dirichlet  (χ(0)=0, χ(L)=0)", kwDirichlet[LwVal]];
printModes["Neumann-Neumann      (χ'(0)=0, χ'(L)=0)", kwNeumann[LwVal]];
printModes["Periodic             (χ(0)=χ(L), χ'(0)=χ'(L))", kwPeriodic[LwVal]];
printModes["Mixed D-N            (χ(0)=0, χ'(L)=0)", kwMixedDN[LwVal]];
printModes["Mixed N-D            (χ'(0)=0, χ(L)=0)", kwMixedND[LwVal]];

Print["\nKey takeaway so far:"];
Print["  • Dirichlet-Dirichlet and mixed BCs forbid k_w=0 (no true zero mode)."];
Print["  • Neumann-Neumann and Periodic BCs allow k_w=0 (true zero mode exists)."];

(* ============================================================
   2) Robin BC: χ'(0)=α χ(0),  χ'(L)=β χ(L)
   Zero-mode analysis must treat k=0 separately.
   ============================================================ *)

Print["\n======================================================="];
Print[" Robin BC zero-mode test"];
Print["======================================================="];

(* For k_w = 0, equation becomes χ''=0 => χ(w)=A + B w.
   Apply Robin BC:
     χ'(0)=B = α χ(0)=α A  => B=α A
     χ'(L)=B = β χ(L)=β (A + B L)=β (A + α A L)=β A (1 + α L)

   For nontrivial A:
     α A = β A (1 + α L)  => α = β (1 + α L)
   Special cases:
     - Neumann is α=β=0  => YES zero mode
     - Many other Robin choices => NO zero mode
*)

robinZeroModeCondition[L_, α_, β_] := Simplify[α == β (1 + α L)];

(* A practical boolean check: allow exact or numeric *)
robinHasZeroMode[L_, α_, β_] := Module[{cond},
  cond = robinZeroModeCondition[L, α, β];
  Which[
    cond === True, True,
    cond === False, False,
    True, False (* if symbolic undecided, default to "not guaranteed" *)
  ]
];

Print["Robin BC: χ'(0)=α χ(0), χ'(L)=β χ(L)"];
Print["Zero-mode condition (k_w=0):  α = β (1 + α L)"];

(* Demonstrate a few representative parameter choices *)
robinExamples = {
  {"Neumann (α=0, β=0)", 0, 0},
  {"Symmetric Robin (α=1, β=1)", 1, 1},
  {"One-sided Robin (α=0, β=1)", 0, 1},
  {"Dirichlet-like limit (α→∞, β→∞)", Infinity, Infinity}
};

Do[
  Module[{label = robinExamples[[i, 1]],
          α = robinExamples[[i, 2]],
          β = robinExamples[[i, 3]],
          has0},
    Print["\n--- ", label, " ---"];
    If[α === Infinity || β === Infinity,
      Print["(Infinity treated as 'Dirichlet-like': expect NO zero mode.)"];
      Print["Zero mode present?  NO"];
      ,
      Print["Condition evaluates to: ", robinZeroModeCondition[LwVal, α, β] // Simplify];
      has0 = robinHasZeroMode[LwVal, α, β];
      Print["Zero mode present?  ", If[has0, "YES (k_w=0)", "NO"]];
    ];
  ],
  {i, Length[robinExamples]}
];

Print["\nInterpretation (Robin):"];
Print["  • A true zero mode is NOT generic under Robin BCs."];
Print["  • Neumann-Neumann (α=β=0) is the clean, natural way to guarantee k_w=0."];
Print["  • If your bulk confinement is 'reflective' (vanishing normal derivative), you get a gapless branch."];
Print["  • If confinement is 'hard wall' (Dirichlet-like), you generally do not."];

(*"
Output:

--- Setup ---
Eigenproblem:  -χ''(w) = λ χ(w),  λ = k_w^2,  w ∈ [0, Lw]
Using Lw = 1, printing up to mMax = 6 modes.

--- Dirichlet-Dirichlet  (χ(0)=0, χ(L)=0) ---
k_w values (first 6): {3.141592653589793, 6.283185307179586, 9.42477796076938, 12.566370614359172, 15.707963267948966, 18.84955592153876}
λ = k_w^2 (first 6): {9.869604401089358, 39.47841760435743, 88.82643960980423, 157.91367041742973, 246.74011002723395, 355.3057584392169}
Zero mode present?  NO

--- Neumann-Neumann      (χ'(0)=0, χ'(L)=0) ---
k_w values (first 6): {0., 3.141592653589793, 6.283185307179586, 9.42477796076938, 12.566370614359172, 15.707963267948966}
λ = k_w^2 (first 6): {0., 9.869604401089358, 39.47841760435743, 88.82643960980423, 157.91367041742973, 246.74011002723395}
Zero mode present?  YES (k_w=0)

--- Periodic             (χ(0)=χ(L), χ'(0)=χ'(L)) ---
k_w values (first 6): {0., 6.283185307179586, 12.566370614359172, 18.84955592153876, 25.132741228718345, 31.41592653589793}
λ = k_w^2 (first 6): {0., 39.47841760435743, 157.91367041742973, 355.3057584392169, 631.6546816697189, 986.9604401089358}
Zero mode present?  YES (k_w=0)

--- Mixed D-N            (χ(0)=0, χ'(L)=0) ---
k_w values (first 6): {1.5707963267948966, 4.71238898038469, 7.853981633974483, 10.995574287564276, 14.137166941154069, 17.27875959474386}
λ = k_w^2 (first 6): {2.4674011002723395, 22.206609902451056, 61.68502750680849, 120.90265391334464, 199.8594891220595, 298.5555331329531}
Zero mode present?  NO

--- Mixed N-D            (χ'(0)=0, χ(L)=0) ---
k_w values (first 6): {1.5707963267948966, 4.71238898038469, 7.853981633974483, 10.995574287564276, 14.137166941154069, 17.27875959474386}
λ = k_w^2 (first 6): {2.4674011002723395, 22.206609902451056, 61.68502750680849, 120.90265391334464, 199.8594891220595, 298.5555331329531}
Zero mode present?  NO

Key takeaway so far:
  • Dirichlet-Dirichlet and mixed BCs forbid k_w=0 (no true zero mode).
  • Neumann-Neumann and Periodic BCs allow k_w=0 (true zero mode exists).

=======================================================
 Robin BC zero-mode test
=======================================================
Robin BC: χ'(0)=α χ(0), χ'(L)=β χ(L)
Zero-mode condition (k_w=0):  α = β (1 + α L)

--- Neumann (α=0, β=0) ---
Condition evaluates to: True
Zero mode present?  YES (k_w=0)

--- Symmetric Robin (α=1, β=1) ---
Condition evaluates to: False
Zero mode present?  NO

--- One-sided Robin (α=0, β=1) ---
Condition evaluates to: False
Zero mode present?  NO

--- Dirichlet-like limit (α→∞, β→∞) ---
(Infinity treated as 'Dirichlet-like': expect NO zero mode.)
Zero mode present?  NO
"*)
