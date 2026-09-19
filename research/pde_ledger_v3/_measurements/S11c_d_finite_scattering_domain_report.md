# Finite scattering boundary and regulator comparison

The next selected set follows the accepted response resolution study at
`c1d466e8`. It compares a129-coefficient matching baseline, a larger finite
boundary/source interval, and a separate regulator change. Physical parameters,
all80 nonlocal rows, four incident directions and the native engine are retained.

Focused reverse AST joins preserve the complete finite constructor, cutoff
adapter and numerical inspector except the declared parameter sites. Both
80-row cutoff rebindings have zero coefficient reverse residuals. The retained
16-node single-momentum layout and common-grid fields reproduce exactly.
Direct modal plane-wave evaluation/projection checks the common-origin phase
maps to4.51e-16; wrong-sign phase controls respond. An actual regulator-dependent
coefficient changes by6.54e-5 under the selected regulator change.

The initial focused checker compared SymPy integer cutoffs with Python floats;
its corrected exact-rational comparison passes. Original logs are retained.
No constructor or physical equation changed for that checker correction.

The selected production completed in170.14 seconds with three clean sequential
children and empty stderr. Every645-unknown system is full rank; balanced
conditions range3773–4241. Scaled equation residuals are below6.60e-16,
boundary residuals below1.05e-13, matrix/direct-action residuals below1.39e-16,
and independent unscaled solves agree in coefficients to1.68e-13.

| Comparison | Largest amplitude change | Largest total-current change |
|---|---:|---:|
|97→129 coefficients at interval48, regulator0.2|1.2813e-7|1.8494e-7|
|Interval/source48→64 at regulator0.2|2.1422e-6|2.6116e-7|
|Regulator0.2→0.1 at interval64|5.100e-12|1.0950e-11|

Open amplitudes and common-grid fields use the same incident phase origin.
The largest field change is6.00e-6 in the declared field unit frames. Actual
rephased current forms, all80 rows,35 generic sources,70 Gaussian binding
checks, six profiles, three layouts per case, ordered limits, measures, exact
cutoff reverse joins and all saved comparisons/source/packet hashes pass.
The acceptance checker corrected its JSON tuple/list and native limits-key
lookups; original checker/logs are preserved. No integrations were repeated.

This supports the dominant finite response within the predeclared1% and1e-4
amplitude reporting targets. The final total-current ratios range
0.9999999308–1.0000001063. Their tiny departures from one and small reflected
signals remain unresolved at the1e-6 current reporting resolution. These are
empirical spreads, not rigorous error bounds, transparent-boundary results or
an Abel limit. Finite modal boundaries and retained operator truncation remain.

Next construct the independent continuum grade operands and required physical
controls from the saved symbolic sources. The finite response has not supplied
that expansion, a pole calculation or final exports. No broader quadrature
campaign is queued.
