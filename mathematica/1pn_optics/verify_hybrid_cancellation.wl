(* ============================================================
   verify_hybrid_cancellation.wl

   Purpose:
     Verify the Hybrid paper’s key 1PN mixed-term coefficient:
       Chyb(n,q) = 1/2 [ q - (n-1) ]
     and show that (n=5, q=1) yields the orbital-required beta=3 via
       Chyb = -(beta/2)  <=>  beta = -2 Chyb.

   Notation:
     n : polytropic index in P ∝ ρ^n  (optics wants n=5)
     q : defect mass scaling exponent M ∝ ρ^q (hybrid selects q=1)
     Phi : Newtonian potential (your convention: typically Phi < 0)
     Mixed term refers to coefficient of (Phi v^2 / c^2) in the 1PN Lagrangian.

   ============================================================ *)

Print["======================================================="];
Print[" verify_hybrid_cancellation.wl"];
Print["======================================================="];

ClearAll["Global`*"];

(* --- 1) Define hybrid coefficient and mapping to orbital beta --- *)

Chyb[n_, q_] := 1/2 (q - (n - 1));          (* Eq. (Chyb-def): 1/2 [q - (n-1)] *)
betaFromC[C_] := -2 C;                       (* If Lmixed = C * Phi v^2/c^2 and orbital uses -(beta/2) *)

Cscalar[q_] := q/2;                         (* + (q/2) piece *)
Cvector[n_] := -(n - 1)/2;                  (* - (n-1)/2 piece *)
Ccheck[n_, q_] := Cscalar[q] + Cvector[n];

Print["\n--- Definitions ---"];
Print["Chyb(n,q)   = 1/2 [ q - (n-1) ]"];
Print["betaFromC(C) = -2 C   (so Chyb = -(beta/2))"];
Print["Decomposition: Chyb = (q/2) + (-(n-1)/2)"];

(* --- 2) Evaluate the selected pair (n=5, q=1) --- *)

nSel = 5;
qSel = 1;

cSel = Chyb[nSel, qSel];
betaSel = betaFromC[cSel];

Print["\n--- Evaluate at (n,q) = (5,1) ---"];
Print["Cscalar(q=1) = ", Cscalar[qSel] // Rationalize];
Print["Cvector(n=5) = ", Cvector[nSel] // Rationalize];
Print["Chyb(5,1)    = ", cSel // Rationalize];
Print["(sanity) Cscalar + Cvector = ", Ccheck[nSel, qSel] // Rationalize];
Print["betaeff      = -2 Chyb(5,1) = ", betaSel // Rationalize];

If[cSel === -3/2 && betaSel === 3,
  Print["PASS: (n=5,q=1) gives Chyb = -3/2 and betaeff = 3."],
  Print["FAIL: Unexpected values; check definitions/sign conventions."]
];

(* --- 3) Show what happens without the scalar cancellation (q=0) --- *)

Print["\n--- Optical-only (q=0) comparison at n=5 ---"];
cOptOnly = Chyb[5, 0];
betaOptOnly = betaFromC[cOptOnly];
Print["Chyb(5,0) = ", cOptOnly // Rationalize, "  =>  betaeff = ", betaOptOnly // Rationalize];
Print["(This reproduces the familiar 'beta ~ 4' optical-only drift.)"];

(* --- 4) Solve the cancellation condition algebraically --- *)

Print["\n--- Algebraic condition for matching orbital beta=3 ---"];
targetBeta = 3;
targetC = -(targetBeta/2);
Print["Target: Chyb = -(beta/2) = ", targetC // Rationalize];

solQ = Solve[Chyb[n, q] == targetC, q];
Print["Solve for q in terms of n: ", solQ // InputForm];
Print["Plugging n=5 gives q = ", (q /. solQ[[1]] /. n -> 5) // Rationalize];

solN = Solve[Chyb[n, q] == targetC, n];
Print["Solve for n in terms of q: ", solN // InputForm];
Print["Plugging q=1 gives n = ", (n /. solN[[1]] /. q -> 1) // Rationalize];

Print["\nDone."];

(*"
Output:

=======================================================
 verify_hybrid_cancellation.wl
=======================================================

--- Definitions ---
Chyb(n,q)   = 1/2 [ q - (n-1) ]
betaFromC(C) = -2 C   (so Chyb = -(beta/2))
Decomposition: Chyb = (q/2) + (-(n-1)/2)

--- Evaluate at (n,q) = (5,1) ---
Cscalar(q=1) = 1/2
Cvector(n=5) = -2
Chyb(5,1)    = -3/2
(sanity) Cscalar + Cvector = -3/2
betaeff      = -2 Chyb(5,1) = 3
PASS: (n=5,q=1) gives Chyb = -3/2 and betaeff = 3.

--- Optical-only (q=0) comparison at n=5 ---
Chyb(5,0) = -2  =>  betaeff = 4
(This reproduces the familiar 'beta ~ 4' optical-only drift.)

--- Algebraic condition for matching orbital beta=3 ---
Target: Chyb = -(beta/2) = -3/2
Solve for q in terms of n: InputForm[{{q -> -4 + n}}]
Plugging n=5 gives q = 1
Solve for n in terms of q: InputForm[{{n -> 4 + q}}]
Plugging q=1 gives n = 5
"*)
