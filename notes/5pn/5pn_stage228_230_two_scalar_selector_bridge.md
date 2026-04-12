# 5PN continuation notes — stages 228 through 230

These stages continue directly from the Stage-224/226 overlap-image branch-selection result.

The previous verdict was:

- the actual coherent moving-throat branch laws do **not** kill the compensated overlap image,
- they do **not** force it either,
- and at the raw overlap-state level they leave the image underdetermined by a four-equation selector quartet.

The point of this continuation was to sharpen that further by working **inside the coherent branch-law manifold itself** and solving out the latent support pair.

The outcome is a real reduction of the theorem gate:

> once the actual branch is already coherent-law clean, overlap-image realizability is no longer a four-selector problem. It is an exact **two-scalar** problem.

---

## Stage 228 — exact branch-law reduction of the overlap-image test

Files:
- `5pn_stage228_branch_law_two_scalar_reduction.py`

Start from the full Stage-226 overlap-image residual vector
\[
(R_K,
 R_M,
 R_{\Omega_U},
 R_{\Omega_W},
 R_{\epsilon_\eta},
 R_{\epsilon_W},
 R_{\delta_U}).
\]

Restrict to the actual coherent moving-throat branch-law manifold
\[
B_\eta=0,
\qquad
B_{\rm tr}=0,
\qquad
B_{\rm nt}=0.
\]

On the constructive branch with
\[
\chi_0=\frac29,
\qquad
\delta_U=1,
\qquad
E_*=\frac14,
\qquad
F_*=\frac56,
\]
these branch laws imply
\[
\delta\ln\epsilon_\eta=0,
\qquad
\delta\ln\delta_U=-\frac{18}{11}\,\delta\ln\chi_0,
\]

and
\[
\delta\ln\Omega_W
=
\frac{1}{2}
\Bigl(
\delta\ln Z_W+
E_*\,\delta\ln\epsilon_W-
F_*\,\delta\ln\delta_U
\Bigr).
\]

After substitution, the seven overlap-image residuals collapse exactly to:

- one support residual triple
  \[
  (R_K,
   R_M,
   R_{\Omega_U}),
  \]
- one orbit/shape selector
  \[
  S_{\rm shape}
  :=
  \delta\ln\epsilon_W-
  \bigl(2\,\delta\ln\chi_0+\delta\ln Z_W\bigr),
  \]
- and the exact identities
  \[
  R_{\epsilon_\eta}=0,
  \qquad
  R_{\delta_U}=0,
  \qquad
  R_{\Omega_W}=\frac18 S_{\rm shape},
  \qquad
  R_{\epsilon_W}=S_{\rm shape}.
  \]

So on the coherent branch-law manifold, overlap-image membership is already equivalent to
\[
R_K=R_M=R_{\Omega_U}=0,
\qquad
S_{\rm shape}=0.
\]

That is already sharper than the Stage-226 selector quartet.

---

## Stage 229 — exact support-plane compiler

Files:
- `5pn_stage229_support_plane_compiler.py`
- `5pn_stage229_support_plane_compiler_output.txt`

The remaining latent support coordinates are
\[
(\delta\ln C,
\delta\ln\varpi).
\]

Inside the coherent branch-law manifold, the support residual pair
\[
R_K=0,
\qquad
R_M=0
\]
solves **uniquely** for that latent support pair in terms of the observed support-side pair
\[
(\delta\ln K,
\delta\ln M)
\]
and the orbit coordinates
\[
(\delta\ln\chi_0,
\delta\ln Z_W).
\]

Numerically, on the constructive branch,
\[
\delta\ln C
\approx
-22.0706526877\,\delta\ln K
+287.773569763\,\delta\ln M
+49.3529290037\,\delta\ln Z_W
+364.327814331\,\delta\ln\chi_0,
\]
\[
\delta\ln\varpi
\approx
-22.5706526877\,\delta\ln K
+191.654233426\,\delta\ln M
+32.3559287390\,\delta\ln Z_W
+234.769642713\,\delta\ln\chi_0.
\]

After that solve, the whole support-side realizability question collapses to one exact codimension-1 plane:
\[
S_{\rm support}
:=
\delta\ln\Omega_U-
\delta\ln\Omega_U^{({\rm pred})}=0,
\]
with
\[
\delta\ln\Omega_U^{({\rm pred})}
\approx
1.95241347905\,\delta\ln K
-12.3996352800\,\delta\ln M
-1.99787596774\,\delta\ln Z_W
-14.3144730533\,\delta\ln\chi_0.
\]

So the support residual triple is no longer a three-condition problem once the latent support pair is allowed to adjust. It is one exact scalar plane.

---

## Stage 230 — reduced two-scalar realizability tester and restoration map

Files:
- `5pn_stage230_reduced_branch_realizability_tester.py`
- `5pn_stage230_reduced_branch_realizability_tester_output.txt`

After imposing the coherent branch laws and solving out the latent support pair from
\[
(\delta\ln K,\delta\ln M),
\]
the full seven-component residual vector collapses exactly to
\[
(0,
 0,
 S_{\rm support},
 S_{\rm shape}/8,
 0,
 S_{\rm shape},
 0).
\]

So the entire overlap-image membership problem on the coherent branch-law manifold is now:
\[
\boxed{
S_{\rm shape}=0,
\qquad
S_{\rm support}=0.
}
\]

Equivalently, the actual moving-throat tangent no longer needs to be tested in the full 11-component overlap state.
It only needs the reduced observable sextuple
\[
(\delta\ln K,
\delta\ln M,
\delta\ln\Omega_U,
\delta\ln\chi_0,
\delta\ln Z_W,
\delta\ln\epsilon_W),
\]
and then the two exact selectors above.

### Exact restoration map

The minimal restoration back to the compensated overlap image is now explicit:

1. fix the latent support pair to the Stage-229 target values
   \[
   (\delta\ln C,\delta\ln\varpi),
   \]
2. shift
   \[
   \delta\ln\Omega_U \to \delta\ln\Omega_U-S_{\rm support},
   \]
3. shift the coherent shape pair by
   \[
   \delta\ln\epsilon_W \to \delta\ln\epsilon_W-S_{\rm shape},
   \qquad
   \delta\ln\Omega_W \to \delta\ln\Omega_W-\frac18 S_{\rm shape},
   \]
which preserves the coherent branch-law manifold and kills the remaining overlap residuals.

So there is now an exact reduced restoration operator on the coherent branch-law manifold.

---

## What this changes in the 5PN roadmap

The next theorem gate is now much smaller than before.

Before these stages, the remaining question after Stage 226 was roughly:

> does the actual moving-throat branch land in the full compensated overlap image?

After Stages 228–230, that has sharpened to:

> on the coherent branch-law manifold, does the actual moving-throat branch satisfy the two exact selectors
> \[
> S_{\rm shape}=0,
> \qquad
> S_{\rm support}=0?
> \]

So the live reduced task is now:

1. extract the reduced observable sextuple
   \[
   (\delta\ln K,
   \delta\ln M,
   \delta\ln\Omega_U,
   \delta\ln\chi_0,
   \delta\ln Z_W,
   \delta\ln\epsilon_W)
   \]
   from the actual moving-throat branch,
2. evaluate the two selectors above,
3. and check whether both vanish.

That is a real tightening of the branch-selection bottleneck.
