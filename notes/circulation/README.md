# Circulation Pair-Force Derivation Package

This package contains six derivation/audit steps for the 3D finite-mouth circulation / pair-force problem.

Each step has:

- a Markdown derivation note;
- a SymPy audit script;
- a captured output file from running the script.

## Step list

1. `step_01_fluxoid_firewall.md` / `step_01_fluxoid_firewall_sympy.py`  
   Audits the gauge-invariant loop integrand, derives integer phase winding from single-valued `psi`, and verifies that a single-valued gauge function has zero loop winding. The electric charge sign `eta_Q` and circulation integer `n` are kept as separate bookkeeping labels.

2. `step_02_no_universal_force.md` / `step_02_no_universal_force_sympy.py`  
   Shows that the fluxoid constraint alone is not a 3D radial-force law. It derives the current-closure force from `M(d)=B/d^3` and derives the fixed-current potential by a Legendre transform.

3. `step_03_current_loop_closure_mutual_inductance.md` / `step_03_current_loop_closure_mutual_inductance_sympy.py`  
   Adds an explicit fixed-current Maxwell/current-loop closure, derives the reduced Neumann integral from the double line integral, and expands the finite-mouth far-field force.

4. `step_04_facing_mouth_swirl_sign.md` / `step_04_facing_mouth_swirl_sign_sympy.py`  
   Translates global current sign into local swirl labels by computing a tangential swirl circulation and oriented loop area vectors. For two facing mouths, opposite local swirl labels attract under the fixed-current closure.

5. `step_05_3d_orientation_and_finite_size.md` / `step_05_3d_orientation_and_finite_size_sympy.py`  
   Builds the full 3D dipole-orientation law from vector dipoles and recovers the Step-3 coaxial limit. The finite-size expansion is treated as asymptotic, not a convergence proof.

6. `step_06_mixed_sector_plumbing_closure.md` / `step_06_mixed_sector_plumbing_closure_sympy.py`  
   Uses the mixed-sector transfer factor as a nonnegative magnitude in the current-like closure, leaves `Lambda_A` real, and identifies the remaining PDE-side sign gap.

## Main result of this pass

The derivation supports two distinct statements.

First, without a mixed-sector/current-like plumbing law,

\[
\boxed{\text{no universal radial attraction/repulsion sign follows from }n_1n_2\text{ alone.}}
\]

Second, under an identical passive fixed-current / Maxwell-like closure, the leading coaxial far-field force for two **facing** finite mouths is

\[
F_d
=+\frac{3\mu_0\pi R_1^2R_2^2}{2d^4}
\Lambda_1\Lambda_2 N_0 n_1n_2
\left[1-\frac{5}{2}\frac{R_1^2+R_2^2}{d^2}+\cdots\right],
\]
where `F_d < 0` means attraction. Therefore, if `Lambda1*Lambda2*N0 > 0`,

\[
\boxed{n_1n_2<0\quad\Longrightarrow\quad\text{attraction for facing mouths}.}
\]

In words: in the facing-mouth local convention, opposite local swirl labels attract under the magnetism-like current closure only after the plumbing branch fixes `Lambda1*Lambda2*N0 > 0`.

## Run all scripts

From this directory:

```bash
for f in *_sympy.py; do
  echo "=== $f ==="
  python "$f"
done
```

## Status

The remaining mathematical task is to derive the actual plumbing coefficient `Lambda_A` and its sign from the full moving-throat PDE/mixed-sector branch, rather than assuming it as a closure.
