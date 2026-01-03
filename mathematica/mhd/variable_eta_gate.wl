(* ::Package:: *)

(*
Step 2: Variable-η gate operator + dissipation identity (CLI-only)

Purpose (Paper 1):
  (i) Show the correct variable-diffusivity diffusion operator:
        Div(η Grad ψ) = η Laplacian ψ + (Grad η)·(Grad ψ)
      so if you implement only η ∇^2 ψ you are missing a drift term.
 (ii) Show why Div(η Grad ψ) is the "safe" choice: it is self-adjoint and
      guarantees monotone decay of ∫ |ψ|^2 (modulo boundary flux),
      whereas η ∇^2 ψ generally does not when η varies in space.

Notes:
 - This derivation is purely about the diffusion / gate term.
 - It assumes x,y,z,t are real variables and η(x,y,z,t) is real-valued.
 - ψ may be complex-valued.

Run:
  wolframscript -file step2_variable_eta_gate.wl
*)

ClearAll["Global`*"];

(* ---------- Controls ---------- *)
Verbosee = True;      (* print intermediate expansions *)
PrintTeX = False;    (* if True, prints TeXForm[...] to stdout; no files *)

(* ---------- Variables & fields ---------- *)
{x, y, z, t} = {x, y, z, t};

(* Use Greek symbols as function heads for readability *)
(* ψ(x,y,z,t), η(x,y,z,t) *)
psiSym = \[Psi][x, y, z, t];
etaSym = \[Eta][x, y, z, t];

(* ---------- Vector calculus helpers (3D) ---------- *)
Grad3[f_] := {D[f, x], D[f, y], D[f, z]};
Div3[v_] := D[v[[1]], x] + D[v[[2]], y] + D[v[[3]], z];
Lap3[f_] := Div3[Grad3[f]];

(* ---------- Assumptions ---------- *)
$Assumptions =
  Element[{x, y, z, t}, Reals] &&
  Element[etaSym, Reals] &&
  etaSym > 0;

(* ---------- Pretty print helpers ---------- *)
hdr[s_] := (Print["\n============================"]; Print[s]; Print["============================\n"]);
subhdr[s_] := (Print["\n--- ", s, " ---\n"]);

printMaybeTeX[label_, expr_] := Module[{},
  If[PrintTeX,
    Print[label, " (TeX):\n", TeXForm[expr], "\n"],
    Null
  ];
];

(* ====================================================================== *)
hdr["STEP 2A: Operator identity for variable diffusivity"];

subhdr["Define diffusion operators"];
opDiv = Expand[Div3[etaSym*Grad3[psiSym]]];
opNaive = Expand[etaSym*Lap3[psiSym]];

Print["Operator (recommended):  Dψ = Div(η Grad ψ)"];
Print["Operator (naive):        Dψ = η Laplacian ψ"];

If[Verbosee,
  Print["\nExpanded Div(η Grad ψ) =\n", opDiv, "\n"];
  Print["Expanded η Laplacian ψ =\n", opNaive, "\n"];
];

subhdr["Compute the missing term"];
missing = FullSimplify[opDiv - opNaive];
Print["Div(η Grad ψ) - η Laplacian ψ =\n", missing, "\n"];
Print["(This should equal (Grad η)·(Grad ψ).)"];

checkMissing = FullSimplify[missing - (Grad3[etaSym].Grad3[psiSym])];
Print["Check (should be 0): ", checkMissing];

printMaybeTeX["Missing term", missing];

(* ====================================================================== *)
hdr["STEP 2B: Local dissipation identity for Div(η Grad ψ)"];

subhdr["Key product-rule identity (pointwise)"];
(* Identity: Conj[ψ] Div(η Grad ψ) = Div(η Conj[ψ] Grad ψ) - η (Grad Conj[ψ])·(Grad ψ) *)
lhs1 = Conjugate[psiSym]*Div3[etaSym*Grad3[psiSym]];
rhs1 = Div3[etaSym*Conjugate[psiSym]*Grad3[psiSym]] - etaSym*(Grad3[Conjugate[psiSym]].Grad3[psiSym]);

id1 = FullSimplify[Expand[lhs1 - rhs1]];
Print["Conj[ψ] Div(η Grad ψ) - (Div(η Conj[ψ] Grad ψ) - η Grad Conj[ψ]·Grad ψ) =\n", id1];
Print["Check (should be 0): ", FullSimplify[id1 == 0]];

subhdr["From the pointwise identity to norm decay"];
(* For diffusion-only evolution: ∂t ψ = Div(η Grad ψ) *)
dNormHalf = D[1/2*Conjugate[psiSym]*psiSym, t] /. {
    D[psiSym, t] -> Div3[etaSym*Grad3[psiSym]],
    D[Conjugate[psiSym], t] -> Conjugate[Div3[etaSym*Grad3[psiSym]]]
};

(* With η real, Conjugate[Div(η Grad ψ)] = Div(η Grad Conjugate[ψ]) *)
dNormHalf2 = FullSimplify[
  dNormHalf /. {Conjugate[etaSym] -> etaSym},
  $Assumptions
];

flux = etaSym*Re[Conjugate[psiSym]*Grad3[psiSym]];
gradSq = (Grad3[Conjugate[psiSym]].Grad3[psiSym]);  (* = |Grad ψ|^2, real and >=0 *)

target = Div3[flux] - etaSym*gradSq;

delta = FullSimplify[Expand[dNormHalf2 - target], $Assumptions];

Print["d/dt (|ψ|^2/2) - (Div(η Re(Conj[ψ] Grad ψ)) - η |Grad ψ|^2) =\n", delta];
Print["Check (should be 0): ", FullSimplify[delta == 0, $Assumptions]];

Print[
  "\nResult (local form):\n",
  "  ∂t(|ψ|^2/2) = ∇·( η Re(ψ* ∇ψ) ) - η |∇ψ|^2\n"
];
Print[
  "Integrated (periodic or decaying boundary so the flux integral vanishes):\n",
  "  d/dt ∫(|ψ|^2/2) dV = - ∫ η |∇ψ|^2 dV ≤ 0\n",
  "Equivalently: d/dt ∫|ψ|^2 dV = -2 ∫ η |∇ψ|^2 dV.\n"
];

printMaybeTeX["Local dissipation identity", HoldForm[D[Abs[\[Psi]]^2/2, t] == Div3[flux] - etaSym*gradSq]];

(* ====================================================================== *)
hdr["STEP 2C: Why the naive η Laplacian ψ form is dangerous when η varies"];

subhdr["Pointwise identity for the naive form"];
lhs2 = Conjugate[psiSym]*(etaSym*Lap3[psiSym]);

(* Using the same flux term, the naive operator leaves an extra drift term involving Grad η *)
rhs2 = Div3[etaSym*Conjugate[psiSym]*Grad3[psiSym]] - etaSym*(Grad3[Conjugate[psiSym]].Grad3[psiSym]) - Conjugate[psiSym]*(Grad3[etaSym].Grad3[psiSym]);

id2 = FullSimplify[Expand[lhs2 - rhs2]];
Print["Conj[ψ] η Laplacian ψ - (Div(η Conj[ψ] Grad ψ) - η |Grad ψ|^2 - Conj[ψ] (Grad η·Grad ψ)) =\n", id2];
Print["Check (should be 0): ", FullSimplify[id2 == 0]];

subhdr["Corresponding local balance shows the extra term explicitly"];
dNormHalfNaive = D[1/2*Conjugate[psiSym]*psiSym, t] /. {
    D[psiSym, t] -> etaSym*Lap3[psiSym],
    D[Conjugate[psiSym], t] -> Conjugate[etaSym*Lap3[psiSym]]
} // Expand;

dNormHalfNaive2 = FullSimplify[dNormHalfNaive /. {Conjugate[etaSym] -> etaSym}, $Assumptions];

extraTerm = -Re[Conjugate[psiSym]*(Grad3[etaSym].Grad3[psiSym])];

targetNaive = Div3[flux] - etaSym*gradSq + extraTerm;

deltaNaive = FullSimplify[Expand[dNormHalfNaive2 - targetNaive], $Assumptions];

Print["d/dt (|ψ|^2/2) - (Div(η Re(Conj[ψ] Grad ψ)) - η |∇ψ|^2 - Re(Conj[ψ] (∇η·∇ψ))) =\n", deltaNaive];
Print["Check (should be 0): ", FullSimplify[deltaNaive == 0, $Assumptions]];

Print[
  "\nTakeaway:\n",
  "  If η varies in space and you use ∂tψ = η ∇^2ψ, you introduce an extra term\n",
  "    -Re(ψ* (∇η·∇ψ))\n",
  "  which is not sign-definite. So ∫|ψ|^2 need not decay monotonically.\n",
  "  Using ∂tψ = ∇·(η∇ψ) removes that term and restores guaranteed dissipation.\n"
];

(* ====================================================================== *)
hdr["STEP 2D: Implementation note (what to put in code)"];

Print[
  "To implement the recommended gate term in finite differences / spectral code:\n",
  "  ∇·(η∇ψ) = η ∇^2ψ + (∇η)·(∇ψ)\n",
  "So you need BOTH pieces when η=η(x,t) is localized.\n",
  "If η is time-only (η(t)), then ∇η=0 and η∇^2ψ is sufficient.\n"
];

Print["Done."];

(*"
Output:

============================
STEP 2A: Operator identity for variable diffusivity
============================

--- Define diffusion operators ---

Operator (recommended):  Dψ = Div(η Grad ψ)
Operator (naive):        Dψ = η Laplacian ψ

Expanded Div(η Grad ψ) =
Derivative[0, 0, 1, 0][η][x, y, z, t]*Derivative[0, 0, 1, 0][ψ][x, y, z, t] + η[x, y, z, t]*Derivative[0, 0, 2, 0][ψ][x, y, z, t] + Derivative[0, 1, 0, 0][η][x, y, z, t]*Derivative[0, 1, 0, 0][ψ][x, y, z, t] + η[x, y, z, t]*Derivative[0, 2, 0, 0][ψ][x, y, z, t] + Derivative[1, 0, 0, 0][η][x, y, z, t]*Derivative[1, 0, 0, 0][ψ][x, y, z, t] + η[x, y, z, t]*Derivative[2, 0, 0, 0][ψ][x, y, z, t]

Expanded η Laplacian ψ =
η[x, y, z, t]*Derivative[0, 0, 2, 0][ψ][x, y, z, t] + η[x, y, z, t]*Derivative[0, 2, 0, 0][ψ][x, y, z, t] + η[x, y, z, t]*Derivative[2, 0, 0, 0][ψ][x, y, z, t]


--- Compute the missing term ---

Div(η Grad ψ) - η Laplacian ψ =
Derivative[0, 0, 1, 0][η][x, y, z, t]*Derivative[0, 0, 1, 0][ψ][x, y, z, t] + Derivative[0, 1, 0, 0][η][x, y, z, t]*Derivative[0, 1, 0, 0][ψ][x, y, z, t] + Derivative[1, 0, 0, 0][η][x, y, z, t]*Derivative[1, 0, 0, 0][ψ][x, y, z, t]

(This should equal (Grad η)·(Grad ψ).)
Check (should be 0): 0

============================
STEP 2B: Local dissipation identity for Div(η Grad ψ)
============================


--- Key product-rule identity (pointwise) ---

Conj[ψ] Div(η Grad ψ) - (Div(η Conj[ψ] Grad ψ) - η Grad Conj[ψ]·Grad ψ) =
0
Check (should be 0): True

--- From the pointwise identity to norm decay ---

d/dt (|ψ|^2/2) - (Div(η Re(Conj[ψ] Grad ψ)) - η |Grad ψ|^2) =
-I*Conjugate[Derivative[0, 0, 1, 0][ψ][x, y, z, t]]*Im[Derivative[0, 0, 1, 0][η][x, y, z, t]]*ψ[x, y, z, t] - I*Conjugate[Derivative[0, 1, 0, 0][ψ][x, y, z, t]]*Im[Derivative[0, 1, 0, 0][η][x, y, z, t]]*ψ[x, y, z, t] - I*Conjugate[Derivative[1, 0, 0, 0][ψ][x, y, z, t]]*Im[Derivative[1, 0, 0, 0][η][x, y, z, t]]*ψ[x, y, z, t] + (η[x, y, z, t]*((Conjugate[Derivative[0, 0, 2, 0][ψ][x, y, z, t]] + Conjugate[Derivative[0, 2, 0, 0][ψ][x, y, z, t]] + Conjugate[Derivative[2, 0, 0, 0][ψ][x, y, z, t]])*ψ[x, y, z, t] - 2*Derivative[1][Conjugate][ψ[x, y, z, t]]*((-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 0, 1, 0][ψ][x, y, z, t]])*Derivative[0, 0, 1, 0][ψ][x, y, z, t]^2 + (-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 1, 0, 0][ψ][x, y, z, t]])*Derivative[0, 1, 0, 0][ψ][x, y, z, t]^2 + (-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[1, 0, 0, 0][ψ][x, y, z, t]])*Derivative[1, 0, 0, 0][ψ][x, y, z, t]^2) + Conjugate[ψ[x, y, z, t]]*((1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 0, 1, 0][ψ][x, y, z, t]])*Derivative[0, 0, 2, 0][ψ][x, y, z, t] + (1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 1, 0, 0][ψ][x, y, z, t]])*Derivative[0, 2, 0, 0][ψ][x, y, z, t] + (1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[1, 0, 0, 0][ψ][x, y, z, t]])*Derivative[2, 0, 0, 0][ψ][x, y, z, t])))/2
Check (should be 0): (-2*I)*Conjugate[Derivative[0, 0, 1, 0][ψ][x, y, z, t]]*Im[Derivative[0, 0, 1, 0][η][x, y, z, t]]*ψ[x, y, z, t] - (2*I)*Conjugate[Derivative[0, 1, 0, 0][ψ][x, y, z, t]]*Im[Derivative[0, 1, 0, 0][η][x, y, z, t]]*ψ[x, y, z, t] - (2*I)*Conjugate[Derivative[1, 0, 0, 0][ψ][x, y, z, t]]*Im[Derivative[1, 0, 0, 0][η][x, y, z, t]]*ψ[x, y, z, t] + η[x, y, z, t]*((Conjugate[Derivative[0, 0, 2, 0][ψ][x, y, z, t]] + Conjugate[Derivative[0, 2, 0, 0][ψ][x, y, z, t]] + Conjugate[Derivative[2, 0, 0, 0][ψ][x, y, z, t]])*ψ[x, y, z, t] - 2*Derivative[1][Conjugate][ψ[x, y, z, t]]*((-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 0, 1, 0][ψ][x, y, z, t]])*Derivative[0, 0, 1, 0][ψ][x, y, z, t]^2 + (-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 1, 0, 0][ψ][x, y, z, t]])*Derivative[0, 1, 0, 0][ψ][x, y, z, t]^2 + (-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[1, 0, 0, 0][ψ][x, y, z, t]])*Derivative[1, 0, 0, 0][ψ][x, y, z, t]^2) + Conjugate[ψ[x, y, z, t]]*((1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 0, 1, 0][ψ][x, y, z, t]])*Derivative[0, 0, 2, 0][ψ][x, y, z, t] + (1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 1, 0, 0][ψ][x, y, z, t]])*Derivative[0, 2, 0, 0][ψ][x, y, z, t] + (1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[1, 0, 0, 0][ψ][x, y, z, t]])*Derivative[2, 0, 0, 0][ψ][x, y, z, t])) == 0

Result (local form):
  ∂t(|ψ|^2/2) = ∇·( η Re(ψ* ∇ψ) ) - η |∇ψ|^2

Integrated (periodic or decaying boundary so the flux integral vanishes):
  d/dt ∫(|ψ|^2/2) dV = - ∫ η |∇ψ|^2 dV ≤ 0
Equivalently: d/dt ∫|ψ|^2 dV = -2 ∫ η |∇ψ|^2 dV.


============================
STEP 2C: Why the naive η Laplacian ψ form is dangerous when η varies
============================


--- Pointwise identity for the naive form ---

Conj[ψ] η Laplacian ψ - (Div(η Conj[ψ] Grad ψ) - η |Grad ψ|^2 - Conj[ψ] (Grad η·Grad ψ)) =
0
Check (should be 0): True

--- Corresponding local balance shows the extra term explicitly ---

d/dt (|ψ|^2/2) - (Div(η Re(Conj[ψ] Grad ψ)) - η |∇ψ|^2 - Re(Conj[ψ] (∇η·∇ψ))) =
-I*Conjugate[Derivative[0, 0, 1, 0][ψ][x, y, z, t]]*Im[Derivative[0, 0, 1, 0][η][x, y, z, t]]*ψ[x, y, z, t] - I*Conjugate[Derivative[0, 1, 0, 0][ψ][x, y, z, t]]*Im[Derivative[0, 1, 0, 0][η][x, y, z, t]]*ψ[x, y, z, t] - I*Conjugate[Derivative[1, 0, 0, 0][ψ][x, y, z, t]]*Im[Derivative[1, 0, 0, 0][η][x, y, z, t]]*ψ[x, y, z, t] + (η[x, y, z, t]*((Conjugate[Derivative[0, 0, 2, 0][ψ][x, y, z, t]] + Conjugate[Derivative[0, 2, 0, 0][ψ][x, y, z, t]] + Conjugate[Derivative[2, 0, 0, 0][ψ][x, y, z, t]])*ψ[x, y, z, t] - 2*Derivative[1][Conjugate][ψ[x, y, z, t]]*((-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 0, 1, 0][ψ][x, y, z, t]])*Derivative[0, 0, 1, 0][ψ][x, y, z, t]^2 + (-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 1, 0, 0][ψ][x, y, z, t]])*Derivative[0, 1, 0, 0][ψ][x, y, z, t]^2 + (-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[1, 0, 0, 0][ψ][x, y, z, t]])*Derivative[1, 0, 0, 0][ψ][x, y, z, t]^2) + Conjugate[ψ[x, y, z, t]]*((1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 0, 1, 0][ψ][x, y, z, t]])*Derivative[0, 0, 2, 0][ψ][x, y, z, t] + (1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 1, 0, 0][ψ][x, y, z, t]])*Derivative[0, 2, 0, 0][ψ][x, y, z, t] + (1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[1, 0, 0, 0][ψ][x, y, z, t]])*Derivative[2, 0, 0, 0][ψ][x, y, z, t])))/2
Check (should be 0): (-2*I)*Conjugate[Derivative[0, 0, 1, 0][ψ][x, y, z, t]]*Im[Derivative[0, 0, 1, 0][η][x, y, z, t]]*ψ[x, y, z, t] - (2*I)*Conjugate[Derivative[0, 1, 0, 0][ψ][x, y, z, t]]*Im[Derivative[0, 1, 0, 0][η][x, y, z, t]]*ψ[x, y, z, t] - (2*I)*Conjugate[Derivative[1, 0, 0, 0][ψ][x, y, z, t]]*Im[Derivative[1, 0, 0, 0][η][x, y, z, t]]*ψ[x, y, z, t] + η[x, y, z, t]*((Conjugate[Derivative[0, 0, 2, 0][ψ][x, y, z, t]] + Conjugate[Derivative[0, 2, 0, 0][ψ][x, y, z, t]] + Conjugate[Derivative[2, 0, 0, 0][ψ][x, y, z, t]])*ψ[x, y, z, t] - 2*Derivative[1][Conjugate][ψ[x, y, z, t]]*((-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 0, 1, 0][ψ][x, y, z, t]])*Derivative[0, 0, 1, 0][ψ][x, y, z, t]^2 + (-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 1, 0, 0][ψ][x, y, z, t]])*Derivative[0, 1, 0, 0][ψ][x, y, z, t]^2 + (-1 + Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[1, 0, 0, 0][ψ][x, y, z, t]])*Derivative[1, 0, 0, 0][ψ][x, y, z, t]^2) + Conjugate[ψ[x, y, z, t]]*((1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 0, 1, 0][ψ][x, y, z, t]])*Derivative[0, 0, 2, 0][ψ][x, y, z, t] + (1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[0, 1, 0, 0][ψ][x, y, z, t]])*Derivative[0, 2, 0, 0][ψ][x, y, z, t] + (1 - 2*Derivative[1][Re][Conjugate[ψ[x, y, z, t]]*Derivative[1, 0, 0, 0][ψ][x, y, z, t]])*Derivative[2, 0, 0, 0][ψ][x, y, z, t])) == 0

Takeaway:
  If η varies in space and you use ∂tψ = η ∇^2ψ, you introduce an extra term
    -Re(ψ* (∇η·∇ψ))
  which is not sign-definite. So ∫|ψ|^2 need not decay monotonically.
  Using ∂tψ = ∇·(η∇ψ) removes that term and restores guaranteed dissipation.


============================
STEP 2D: Implementation note (what to put in code)
============================

To implement the recommended gate term in finite differences / spectral code:
  ∇·(η∇ψ) = η ∇^2ψ + (∇η)·(∇ψ)
So you need BOTH pieces when η=η(x,t) is localized.
If η is time-only (η(t)), then ∇η=0 and η∇^2ψ is sufficient.
"*)
