# Analog Four-Force Phenomenology Visualizer

This package visualizes the reduced, calibrated 2D phenomenology specified in
`notes/build_spec.md`. It does **not** solve the full 4D nonlinear PDE and does
not claim to derive or dynamically emerge the four forces.

The numerical core is under `physics/` and has no rendering dependency. The
Matplotlib layer is under `scenes/`; the GIF command selects the Agg backend
and writes GIFs with Pillow. Every scene exposes its sector's characterized
departure by default. `--hide-departures` is the explicit counterfactual toggle.

Run the live interactive app locally from the repository's `software/`
directory (a graphical desktop/display is required):

```bash
python -m force_visualizer.app
```

The sector selector switches among gravity, light, charge, and magnetism. Each
sector has live sliders and its own departures checkbox. In a headless session,
the command exits with a clear display-unavailable message; use the GIF command
below for headless rendering.

Each front end also samples the verified core on a spatial grid: gravity shows
one-way inward medium tracers, charge shows a non-flowing Coulomb field-vector
overlay, magnetism shows the circulating brane field of each current, and the
light wave is labeled as the brane shear field. Charge colors update red for
positive and blue for negative (the default pair is opposite-sign); magnetic
throats are colored by current direction and drawn end-on (⊙ out of page, ⊗
into page). Their core-force-driven transverse motion attracts for parallel
currents and repels for antiparallel currents, while the field grid follows
their positions. Particle/throat markers breathe subtly as a qualitative
dynamical cue.

Run from the repository root:

```bash
python -m software.force_visualizer.render_all --output-dir /tmp/force-scenes
python -m software.force_visualizer.report
python -m pytest software/force_visualizer/tests -q
```

The numeric report includes core-derived direction guards for inward gravity
flow, positive/negative electric fields, and current-signed magnetic
circulation, in addition to the force-law checks.

`--include-25pn-benchmark` adds the standard Burke–Thorne benchmark to the
gravity orbit. This is labeled calibrated because the cited 2.5PN paper leaves
the native toy-model response normalization genuinely blocked.

Source map:

- Gravity localization/departure: `pathA_29_brane_bulk_return.md` and results;
  conservative EIH 1PN: `research/4d_1pn_full/paper/4d_1pn_full.tex`;
  conditional 2.5PN benchmark: `research/4d_2_5pn/paper/4d_2_5pn.tex`.
- Light: `pathA_36_c5_phase_potential.md` and results, ledger stages 003/005;
  lensing: `research/1pn_optics/paper/1pn_optics.tex`.
- Charge: `pathA_38_throat_body_electric_localization.md` and results.
- Magnetism: `pathA_39_magnetic_force.md` and results, stage-4 field
  classification, and scalar-admixture screen.
- Cross-sector status: pathA_40 cone-lock and pathA_41 second-medium drift.
