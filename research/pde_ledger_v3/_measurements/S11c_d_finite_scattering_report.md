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

The retry retains325 unknowns, four incoming columns, the15-minute budget, one
native thread and2 GiB ceiling. No matrix solve or physical scattering result is
accepted yet. Inspect saved matrices, action/measure controls, rank, residuals,
conditioning and amplitudes before increasing resolution. Finite regulator,
cutoffs, modest quadrature and approximate modal boundaries remain explicit.

The completed adaptive sequence is published/annex-verified at4ba95bcf/a6bfb7a6.
Future work follows exploratoryAcceptanceV1 and the focused completion plan;
no broad tail/regulator refinement campaign is queued.
