# S11c-d finite scattering pilot

Implementation assembles all five fields and80 nonlocal rows from accepted
source factors, including the actual source derivatives in all35 generic source
amplitudes. It forms all four incident columns with the original physical current
matrices and full outward bases. Exact local operators and full symbolic source
operands are retained; the numerical adapter changes no upstream physics.

Focused boundary checks find five independent outward directions at each end:
two propagating transverse directions and three evanescent thickness directions.
Both trace-basis condition numbers are below5.94. Modal value/derivative map
residuals are below4.47e-16. Independent polynomial checks through derivative
order3 differ by at most3.56e-15. These checks validate the finite boundary
instrument; they do not prove a transparent nonlocal boundary condition.

The initial pilot (`efc61982`, launch `e7e81e01`) stopped after89 seconds
before momentum integration: native SymPy integer derivative orders produced
NumPy object arrays. Their exact conversion to Python integers restores float64
basis matrices. All native orders0–3 now pass independent differentiation checks;
the complete local low-degree action differs by at most1.34e-14. Negative and
nonintegral orders are rejected. The binding constructor and consumed helpers
are unchanged; all35 source jets,70 Gaussian checks and both boundary maps are
reused byte-for-byte with97 source and6 operand joins. Original logs and packets
remain under `pilot/`; the repair evidence is in
`S11c_d_finite_scattering_dtype_repair.json`.

The repaired retry completed in13.24 seconds (constructor11.91 seconds), with
empty stderr,17944 momentum nodes and216.3 MiB peak RSS (221500 KiB). Its325x325
system has rank325 and balanced condition number1600; all four incident columns
were solved. Maximum raw/scaled equation residuals are1.97e-12/1.40e-15,
boundary residuals4.48e-14 and matrix/direct-action residual3.54e-17. An independent
unscaled direct solve agrees in coefficients to9.19e-14. Every80 row,35 source,
70 Gaussian record,97 source pin,6 input operand and saved artifact hash passes.
Full matrices, fields, modal amplitudes and physical currents are preserved.

Computed outgoing/incident current ratios range0.9999630–0.9999931, with small
reflected amplitudes in this coarse pilot. Neither the apparent small loss nor
reflection is yet resolved physically. Finite regulator, cutoffs, modest
quadrature and approximate modal boundaries remain explicit. No continuum
expansion or converged physical scattering result is claimed. Next use a small
resolution comparison on the actual response, preserving domain and regulator.

The completed adaptive sequence is published/annex-verified at4ba95bcf/a6bfb7a6.
Future work follows exploratoryAcceptanceV1 and the focused completion plan;
no broad tail/regulator refinement campaign is queued.
