# Moving-throat PDE PathA-32 grouped-P2 isotropy feed note

Computed verdict: `ISOTROPY_CALIBRATED` (`ISOTROPY_CALIBRATED`).

The SymPy and Mathematica engines independently computed the real l=2 angular basis, the isotropic Gram matrix `I5`, the per-harmonic `-Delta_S2` eigenvalue `6`, and the grouped conservative response assembled per lane from each harmonic self-overlap.

On the isotropic calibrated reference, all raw-D lane defects vanish exactly for orders `0,2,4`; the normalized `a2,b2,a4,b4` cross-check also vanishes. The verdict is `ISOTROPY_CALIBRATED`, not `ISOTROPY_PASS`, because `beta2(w)`, `T_Omega`, wall stiffnesses, and the `Btilde/Ztilde/Ktilde/Mtilde` radial data remain calibration inputs.

The able-to-fail probes now use computed gate flags. The remediated probes each fail with their mutation and stop failing under their self-ablation; neutering any probe's computed flag makes `able_to_fail_ok` false.

Dimensional retrofit: the engines walk the explicit `a^2 dw dOmega` wall measure, the `M2`/`K2` integrand terms, and the actual grouped `Mtilde` / `Ktilde + 6*TomegaTilde` expressions. Corrupting sourced `T_Omega` flips the recomputed verdict to `FAIL_DIMENSIONAL` and the self-ablation restores `ISOTROPY_CALIBRATED`.

Engine agreement status: `pass`, max numeric delta `2.220446049250313e-16`, grouped-lane D max delta `0.0`.
