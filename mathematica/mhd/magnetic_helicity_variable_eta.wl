(* magnetic_helicity_variable_eta_safe.wl
   STEP 4: Magnetic helicity balance with variable resistivity η(x,t)
   Uses custom operators to avoid protected symbols (Grad/Div/Curl/E).
*)

ClearAll["Global`*"];

(* ---------------------------
   Safe vector calculus (component-wise)
   --------------------------- *)
vGrad[f_] := {D[f, x], D[f, y], D[f, z]};

vDiv[v_List] /; Length[v] == 3 :=
  D[v[[1]], x] + D[v[[2]], y] + D[v[[3]], z];

vCurl[v_List] /; Length[v] == 3 := {
  D[v[[3]], y] - D[v[[2]], z],
  D[v[[1]], z] - D[v[[3]], x],
  D[v[[2]], x] - D[v[[1]], y]
};

vCross[u_List, v_List] /; Length[u] == 3 && Length[v] == 3 := {
  u[[2]] v[[3]] - u[[3]] v[[2]],
  u[[3]] v[[1]] - u[[1]] v[[3]],
  u[[1]] v[[2]] - u[[2]] v[[1]]
};

banner[s_] := (
  Print["\n============================"];
  Print[s];
  Print["============================\n"];
);

(* Assumptions: real scalars where appropriate *)
$Assumptions =
  Element[{x, y, z, t}, Reals] &&
  Element[η[x, y, z, t], Reals] &&
  Element[phi[x, y, z, t], Reals];

banner["STEP 4: Magnetic helicity balance with variable resistivity η(x,t)"];

(* ---------------------------
   Fields / definitions
   --------------------------- *)
A = {Ax[x, y, z, t], Ay[x, y, z, t], Az[x, y, z, t]};     (* vector potential *)
v = {vx[x, y, z, t], vy[x, y, z, t], vz[x, y, z, t]};     (* velocity field *)
ηf = η[x, y, z, t];
phif = phi[x, y, z, t];

B = vCurl[A];      (* magnetic field *)
J = vCurl[B];      (* current density (μ0=1 units) *)

Print["B = vCurl(A): ", B];
Print["Check vDiv(B) (should be 0): ",
  Simplify[vDiv[B], Assumptions -> $Assumptions]
];

banner["STEP 4A: Define eField, A_t and prove the local helicity identity"];

(* Ohm's law (resistive MHD):  E + v×B = η J  =>  E = -v×B + η J *)
eField = -vCross[v, B] + ηf*J;

(* Electromagnetic definition: E = -∂t A - ∇phi  =>  ∂tA = -E - ∇phi *)
At = -eField - vGrad[phif];

(* Since B = Curl(A), we have ∂tB = Curl(∂tA) *)
Bt = vCurl[At];

(* Faraday consistency: ∂tB + Curl(E) = 0 *)
Print["Check: Bt + vCurl(eField) (should be {0,0,0}): ",
  Simplify[Bt + vCurl[eField], Assumptions -> $Assumptions]
];

(* Helicity density h = A·B *)
h = A.B;

(* Time derivative using At and Bt: ∂t(A·B) = (∂tA)·B + A·(∂tB) *)
ht = At.B + A.Bt;

(* Local identity (exact for B = Curl(A)):
   ∂t(A·B) + ∇·( E×A + phi B ) = -2 E·B
*)
flux = vCross[eField, A] + phif*B;

resHel = Simplify[ht + vDiv[flux] + 2*(eField.B), Assumptions -> $Assumptions];
Print["Check: ht + vDiv(flux) + 2 (eField·B) (should be 0): ", resHel];

Print["\nLOCAL HELICITY BALANCE (verified):"];
Print["  ∂t(A·B) + ∇·( E×A + phi B ) = -2 E·B"];
Print["  where B=vCurl(A), E=-v×B + η vCurl(B)."];

banner["STEP 4B: Reduce RHS using E = -v×B + ηJ"];

(* Since (v×B)·B = 0, we expect E·B = η (J·B) *)
resEB = Simplify[eField.B - ηf*(J.B), Assumptions -> $Assumptions];
Print["Check: (eField·B) - η (J·B) (should be 0): ", resEB];

Print["\nRHS form:"];
Print["  E·B = η (J·B)  (since (v×B)·B = 0)"];
Print["So:  ∂t(A·B) + ∇·( E×A + phi B ) = -2 η (J·B)"];

banner["STEP 4C: Purely resistive limit (set v = 0)"];

v0 = {0, 0, 0};

eField0 = -vCross[v0, B] + ηf*J;  (* = η J *)
At0 = -eField0 - vGrad[phif];
Bt0 = vCurl[At0];
ht0 = At0.B + A.Bt0;
flux0 = vCross[eField0, A] + phif*B;

res0 = Simplify[ht0 + vDiv[flux0] + 2*(eField0.B), Assumptions -> $Assumptions];
Print["Check (v=0): ht + vDiv(flux) + 2 (E·B) (should be 0): ", res0];

Print["In this limit: E = η J, so helicity changes only where η and (J·B) overlap."];

banner["STEP 4D: Integrated statement (periodic/decaying boundaries)"];

Print["If the surface flux integral vanishes (periodic domain or decaying fields):"];
Print["  d/dt ∫ (A·B) dV = -2 ∫ η (J·B) dV."];
Print["NOTE: J·B is sign-indefinite, so helicity need not decay monotonically (unlike magnetic energy)."];

Print["\nDone."];

(*"
Output:


============================
STEP 4: Magnetic helicity balance with variable resistivity η(x,t)
============================

B = vCurl(A): {-Derivative[0, 0, 1, 0][Ay][x, y, z, t] + Derivative[0, 1, 0, 0][Az][x, y, z, t], Derivative[0, 0, 1, 0][Ax][x, y, z, t] - Derivative[1, 0, 0, 0][Az][x, y, z, t], -Derivative[0, 1, 0, 0][Ax][x, y, z, t] + Derivative[1, 0, 0, 0][Ay][x, y, z, t]}
Check vDiv(B) (should be 0): 0

============================
STEP 4A: Define eField, A_t and prove the local helicity identity
============================

Check: Bt + vCurl(eField) (should be {0,0,0}): {0, 0, 0}
Check: ht + vDiv(flux) + 2 (eField·B) (should be 0): 0

LOCAL HELICITY BALANCE (verified):
  ∂t(A·B) + ∇·( E×A + phi B ) = -2 E·B
  where B=vCurl(A), E=-v×B + η vCurl(B).

============================
STEP 4B: Reduce RHS using E = -v×B + ηJ
============================

Check: (eField·B) - η (J·B) (should be 0): 0

RHS form:
  E·B = η (J·B)  (since (v×B)·B = 0)
So:  ∂t(A·B) + ∇·( E×A + phi B ) = -2 η (J·B)

============================
STEP 4C: Purely resistive limit (set v = 0)
============================

Check (v=0): ht + vDiv(flux) + 2 (E·B) (should be 0): 0
In this limit: E = η J, so helicity changes only where η and (J·B) overlap.

============================
STEP 4D: Integrated statement (periodic/decaying boundaries)
============================

If the surface flux integral vanishes (periodic domain or decaying fields):
  d/dt ∫ (A·B) dV = -2 ∫ η (J·B) dV.
NOTE: J·B is sign-indefinite, so helicity need not decay monotonically (unlike magnetic energy).
"*)
