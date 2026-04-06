(* ============================================================
   verify_length_scaling_from_energy_balance.wl

   Fast version: no Solve/FullSimplify bottlenecks.
   Derives L*(ρ) scaling analytically by hand from dE/dL=0.

   Model:
     E(L;ρ) = A ρ^a L^p  +  B ρ^csExp L^(-r),
     csExp = (n-1)/2 from c_s ∝ ρ^((n-1)/2).

   Result:
     L*(ρ) ∝ ρ^s,  with  s = (csExp - a)/(p + r).

   LPI / universality condition (from your earlier script):
     s_required = csExp - q.

   ============================================================ *)

ClearAll["Global`*"];

(* calibrated defaults *)
nVal = 5;
qVal = 1;

csExp[n_] := (n - 1)/2;
sReq[n_, q_] := csExp[n] - q;

sFromExponents[n_, a_, p_, r_] := (csExp[n] - a)/(p + r);

Print["\n--- Calibrated inputs ---"];
Print["n = ", nVal, "   q = ", qVal];
Print["csExp = (n-1)/2 = ", csExp[nVal] // Rationalize];
Print["s_required = csExp - q = ", sReq[nVal, qVal] // Rationalize];

Print["\n--- General derivation (no CAS solve needed) ---"];
Print["dE/dL=0 gives:  A ρ^a p L^(p-1) = B ρ^csExp r L^(-r-1)"];
Print["=> L^(p+r) ∝ ρ^(csExp-a)  =>  s = (csExp-a)/(p+r)."];

Print["\n--- Plug in n=5 (leave a,p,r symbolic) ---"];
Print["s(n=5) = ", sFromExponents[5, a, p, r] // TraditionalForm];

(* Solve algebraically for a_required to enforce universality *)
targetS = sReq[nVal, qVal];
aReq[p_, r_] := Simplify[csExp[nVal] - (p + r) targetS];

Print["\n--- Universality condition expressed as a_required(p,r) ---"];
Print["Require s = s_required => a_required(p,r) = ", aReq[p, r] // TraditionalForm];

(* For (n,q)=(5,1), csExp=2 and s_required=1 -> a_required = 2 - (p+r) *)
Print["\n--- For (n,q)=(5,1): explicit simplification ---"];
Print["csExp=2, s_required=1 => a_required(p,r) = 2 - (p+r)."];

(* Table of simple integer cases *)
pList = {1, 2, 3};
rList = {1, 2, 3};

tab = Table[
  {p0, r0, aReq[p0, r0] // Rationalize},
  {p0, pList}, {r0, rList}
] // Flatten[#, 1] &;

Print["\nColumns: {p, r, a_required}"];
Print[tab // TableForm];

Print["\n--- Notable cases ---"];
Print["p=1,r=1 -> a_required = ", aReq[1, 1] // Rationalize,
      "  (binding prefactor independent of ρ)"];
Print["p=1,r=2 -> a_required = ", aReq[1, 2] // Rationalize];
Print["p=2,r=1 -> a_required = ", aReq[2, 1] // Rationalize];

Print["\n--- Stability check (2nd derivative) ---"];

Clear[A, B, a, p, r, ρ, L, EE, dE, d2E, LstarNice, d2Nice1, d2Nice2];

EE  = A ρ^a L^p + B ρ^2 L^(-r);
dE  = D[EE, L];
d2E = D[EE, {L, 2}];

(* Closed-form stationary point from dE/dL = 0:
   A ρ^a p L^(p-1) = B ρ^2 r L^(-r-1)
   => L^(p+r) = (B r)/(A p) ρ^(2-a)
*)
LstarNice = ((B r)/(A p))^(1/(p + r)) ρ^((2 - a)/(p + r));

(* Second derivative at L*:
   d2E = A ρ^a p(p-1)L^(p-2) + B ρ^2 r(r+1)L^(-r-2)
   Using the stationarity relation, it simplifies to either form below:
*)
d2Nice1 = A ρ^a p (p + r) LstarNice^(p - 2);
d2Nice2 = B ρ^2 r (p + r) LstarNice^(-r - 2);

(* --- Referee-friendly printing --- *)
Clear[K, s, Lexpr, Lsym, d2sym1, d2sym2, pp];

K = (B r)/(A p);
s = (2 - a)/(p + r);

(* Define L*(ρ) with explicit 1/(p+r) so it prints unambiguously *)
Lexpr = K^(1/(p + r)) ρ^s;

(* Use a symbolic placeholder so d2 lines don’t expand Lexpr inside powers *)
Lsym = Lstar[ρ];

d2sym1 = A p (p + r) ρ^a Lsym^(p - 2);
d2sym2 = B r (p + r) ρ^2 Lsym^(-r - 2);

pp[label_, expr_] := Print[label, ToString[expr, InputForm]];

(* Define Lexpr explicitly *)
Lexpr = ((B r)/(A p))^(1/(p + r)) ρ^((2 - a)/(p + r));

Print["L*(ρ) (TeX) = ", ToString[Lstar[ρ] == Lexpr, TeXForm]];
Print["d^2E/dL^2|_{L*} = ", ToString[A*p*(p + r)*ρ^a*Lstar[ρ]^(p - 2), InputForm]];
Print["(Equivalent form) = ", ToString[B*r*(p + r)*ρ^2*Lstar[ρ]^(-r - 2), InputForm]];

Print["Positivity: for A,B,ρ,p,r>0 and p+r>0, L*(ρ)>0 and d^2E/dL^2|_{L*}>0 (stable minimum)."];


(*"
Output:

--- Calibrated inputs ---
n = 5   q = 1
csExp = (n-1)/2 = 2
s_required = csExp - q = 1

--- General derivation (no CAS solve needed) ---
dE/dL=0 gives:  A ρ^a p L^(p-1) = B ρ^csExp r L^(-r-1)
=> L^(p+r) ∝ ρ^(csExp-a)  =>  s = (csExp-a)/(p+r).

--- Plug in n=5 (leave a,p,r symbolic) ---
s(n=5) = TraditionalForm[(2 - a)/(p + r)]

--- Universality condition expressed as a_required(p,r) ---
Require s = s_required => a_required(p,r) = TraditionalForm[2 - p - r]

--- For (n,q)=(5,1): explicit simplification ---
csExp=2, s_required=1 => a_required(p,r) = 2 - (p+r).

Columns: {p, r, a_required}
TableForm[{{1, 1, 0}, {1, 2, -1}, {1, 3, -2}, {2, 1, -1}, {2, 2, -2}, {2, 3, -3}, {3, 1, -2}, {3, 2, -3}, {3, 3, -4}}]

--- Notable cases ---
p=1,r=1 -> a_required = 0  (binding prefactor independent of ρ)
p=1,r=2 -> a_required = -1
p=2,r=1 -> a_required = -1

--- Stability check (2nd derivative) ---
L*(ρ) (TeX) = \text{Lstar}(\rho )=\rho ^{\frac{2-a}{p+r}} \left(\frac{B r}{A p}\right)^{\frac{1}{p+r}}
d^2E/dL^2|_{L*} = A*p*(p + r)*ρ^a*Lstar[ρ]^(-2 + p)
(Equivalent form) = B*r*(p + r)*ρ^2*Lstar[ρ]^(-2 - r)
Positivity: for A,B,ρ,p,r>0 and p+r>0, L*(ρ)>0 and d^2E/dL^2|_{L*}>0 (stable minimum).
"*)
