(* ============================================================
   verify_bulk_guided_photon_index.wl

   Purpose:
     Model a photon as a bulk-guided wave excitation (NOT a throat),
     using a simple 4D wave equation with variable wave speed c_s(r).

     Derive:
       - dispersion: ω^2 = c_s(r)^2 (k^2 + k_w^2)
       - group speed in brane directions
       - effective refractive index N_eff = c∞ / v_g
       - show that k_w=0 (zero mode) gives achromatic lensing index
         N(r)=c∞/c_s(r), while k_w≠0 introduces chromatic corrections.

   Conventions:
     Φ < 0 near mass, and (to 1PN) δ ≡ Δρ/ρ0 = Φ/c∞^2 (κρ=1).
     EOS: c_s^2 ∝ ρ^(n-1) so c_s ≈ c∞ [1 + (n-1)/2 δ].
   ============================================================ *)

ClearAll["Global`*"];

(* ---- 1) 1PN density and sound-speed structure ---- *)
delta[Phi_, c_] := Phi/c^2;  (* δ = Δρ/ρ0 = Φ/c^2, so δ is negative near mass *)

cs[Phi_, n_, c_] := c * (1 + (n - 1)/2 * delta[Phi, c]);  (* linearized *)

N0[Phi_, n_, c_] := c/cs[Phi, n, c];  (* optical index for k_w=0 *)

Print["\n--- 1PN structure ---"];
Print["δ = Φ/c^2"];
Print["c_s ≈ c [1 + (n-1)/2 δ]"];
Print["N0 = c/c_s ≈ 1 - (n-1)/2 δ"];

(* ---- 2) Bulk-guided dispersion and group velocity ---- *)
omega[k_, kw_, Phi_, n_, c_] := cs[Phi, n, c] * Sqrt[k^2 + kw^2];

vg[k_, kw_, Phi_, n_, c_] := D[omega[k, kw, Phi, n, c], k];

Neff[k_, kw_, Phi_, n_, c_] := c/vg[k, kw, Phi, n, c];

Print["\n--- Dispersion / group speed ---"];
Print["ω = c_s sqrt(k^2 + k_w^2)"];
Print["v_g = ∂ω/∂k = c_s k/sqrt(k^2 + k_w^2)"];
Print["N_eff = c/v_g = (c/c_s) sqrt(1 + (k_w/k)^2)"];

(* ---- 3) Series expansions: small Φ/c^2 and small (kw/k)^2 ---- *)
exprSeries = Series[
  Neff[k, kw, Phi, n, c] /. {kw -> eps*kw},
  {Phi, 0, 1}, {eps, 0, 2}
] // Normal // Simplify;

Print["\n--- N_eff expanded to O(Φ) and O((k_w/k)^2) ---"];
Print[exprSeries // InputForm];

(* ---- 4) Specialize to the lensing choice n=5 ---- *)
exprN5 = exprSeries /. n -> 5 // Simplify;

Print["\n--- Specialize to n=5 ---"];
Print[exprN5 // InputForm];

Print["\nInterpretation:"];
Print["  • If k_w = 0 (true zero mode), N_eff = N0 = 1 - 2 Φ/c^2 = 1 + 2 GM/(rc^2)."];
Print["  • If k_w ≠ 0, N_eff has a + (1/2)(k_w/k)^2 correction ⇒ chromatic / dispersive lensing."];

(*"
Output:

--- 1PN structure ---
δ = Φ/c^2
c_s ≈ c [1 + (n-1)/2 δ]
N0 = c/c_s ≈ 1 - (n-1)/2 δ

--- Dispersion / group speed ---
ω = c_s sqrt(k^2 + k_w^2)
v_g = ∂ω/∂k = c_s k/sqrt(k^2 + k_w^2)
N_eff = c/v_g = (c/c_s) sqrt(1 + (k_w/k)^2)

--- N_eff expanded to O(Φ) and O((k_w/k)^2) ---
InputForm[((2*k^2 + eps^2*kw^2)*(2*c^2 + Phi - n*Phi))/(4*c^2*k*Sqrt[k^2])]

--- Specialize to n=5 ---
InputForm[((2*k^2 + eps^2*kw^2)*(c^2 - 2*Phi))/(2*c^2*k*Sqrt[k^2])]
"*)
