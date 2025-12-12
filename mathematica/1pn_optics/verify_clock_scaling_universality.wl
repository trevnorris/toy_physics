(* ============================================================
   verify_clock_scaling_universality.wl

   Purpose:
     Compare 1PN redshift scaling of:
       (A) atomic clock: ω_atom ∝ m(ρ) ∝ ρ^q
       (B) photonic cavity clock: ω_cav ∝ c_s(ρ)/L(ρ)

     Let:
       ρ = ρ0 (1 + δ), δ = Φ/c^2  (negative near mass)
       c_s ∝ ρ^((n-1)/2)  (polytrope)
       m ∝ ρ^q
       L ∝ ρ^s  (UNKNOWN; replace "Bohr scaling" by a general exponent)

     Then:
       ω_atom/ω0 = 1 + q δ
       ω_cav/ω0  = 1 + [(n-1)/2 - s] δ

     Universality requires:
       [(n-1)/2 - s] = q  ⇒  s = (n-1)/2 - q
   ============================================================ *)

ClearAll["Global`*"];

delta[Phi_, c_] := Phi/c^2;

(* Exponents *)
atomExp[q_] := q;
cavExp[n_, s_] := (n - 1)/2 - s;

Print["\n--- Exponent summary ---"];
Print["Atomic exponent:  q"];
Print["Cavity exponent: (n-1)/2 - s   where L ∝ ρ^s"];

(* Universality condition *)
solS = Solve[cavExp[n, s] == atomExp[q], s][[1]];

Print["\n--- Universality condition ---"];
Print["Solve (n-1)/2 - s = q  ⇒  s = ", (s /. solS) // InputForm];

(* Evaluate the model’s preferred pair *)
nSel = 5;
qSel = 1;

sNeeded = (s /. solS) /. {n -> nSel, q -> qSel} // Simplify;

Print["\n--- Evaluate at (n,q)=(5,1) ---"];
Print["Required s for universality: s = ", sNeeded, "   (i.e., L ∝ ρ^1)"];

(* Compare common assumptions *)
Print["\n--- Compare three length-scaling choices ---"];

sBohr = -qSel;  (* if L ∝ a0 ∝ 1/m ∝ ρ^{-q} *)
Print["Bohr-like rods: s = -q = ", sBohr,
      "  ⇒ cavity exponent = ", cavExp[nSel, sBohr] // Simplify];

sFixed = 0;
Print["Fixed rods:     s = 0",
      "  ⇒ cavity exponent = ", cavExp[nSel, sFixed] // Simplify];

sUni = sNeeded;
Print["Universal rods: s = ", sUni,
      "  ⇒ cavity exponent = ", cavExp[nSel, sUni] // Simplify];

Print["\nInterpretation:"];
Print["  • With (n,q)=(5,1), atomic exponent = 1."];
Print["  • If you assume Bohr scaling (s=-1), cavity exponent = 3 (the '3×' issue)."];
Print["  • If your cavity length instead scales as L ∝ ρ^(+1), cavity exponent = 1 and universality is restored."];

(*"
Output:

--- Exponent summary ---
Atomic exponent:  q
Cavity exponent: (n-1)/2 - s   where L ∝ ρ^s

--- Universality condition ---
Solve (n-1)/2 - s = q  ⇒  s = InputForm[(-1 + n - 2*q)/2]

--- Evaluate at (n,q)=(5,1) ---
Required s for universality: s = 1   (i.e., L ∝ ρ^1)

--- Compare three length-scaling choices ---
Bohr-like rods: s = -q = -1  ⇒ cavity exponent = 3
Fixed rods:     s = 0  ⇒ cavity exponent = 2
Universal rods: s = 1  ⇒ cavity exponent = 1
"*)
