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

Production (`df37ef94`) is running with a900-second total child budget, one native thread
and2 GiB per child. Preserve every completed case and partial. The unchanged
reporting targets are1% for resolved outputs, amplitude resolution1e-4 and
normalized-current resolution1e-6. Actual domain/regulator results are pending;
no transparent boundary, Abel limit, tiny reflection/loss, continuum expansion
or pole result is established by focused instrument checks.
