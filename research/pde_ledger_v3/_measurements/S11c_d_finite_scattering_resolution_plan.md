# First finite-response resolution comparison

The validated pilot (`f01636e9`) computes the complete finite four-incident
response in13.24 seconds. It is now affordable to examine the response itself.
This implements exploratoryAcceptanceV1, with no new broad quadrature campaign.

Keep the physical inputs, modal boundary maps, position/source bound48,
momentum bound4 and regulator0.2 fixed. Reuse the verified35 source-jet bindings
and both boundary maps. Use the unchanged finite constructor sequentially:

1. At65 coefficients per field and momentum8/2/2, change source/profile
   orders128/128 to256/512. The previous actual transfer-range tests motivate
   these orders. This measures the combined source/profile-rule effect; it
   does not isolate the two transforms from each other.
2. Keep those rules and change65 to97 coefficients per field.
3. At97 coefficients and source/profile256/512, change momentum8/2/2 to16/4/4.
   This is one combined momentum-resolution comparison, not separate evidence
   for convergence of each nested order.

Compare all complex open-channel amplitudes, each end's full-current quadratic
form, total outgoing/incident current, evanescent amplitudes and fields evaluated
on one common position grid. Retain the actual current cross terms. Record
equation/boundary residuals, rank, conditioning and an independent direct solve.
The reference includes all80 native nonlocal rows and all four incident columns.

Before these comparisons, set a1% working target for resolved nonzero outputs.
For this initial diagnostic use amplitude resolution1e-4 and normalized-current
resolution1e-6 in their declared dimensionless frames. Below those floors,
record absolute changes and empirical resolution; do not claim a small reflected
signal, a loss, a sign or an absence. A stable total near one cannot resolve its
small deficit. These are reporting goals, not the withheld physical criterion.
They do not certify a boundary approximation or replace continuum bookkeeping.

Measure each case independently, retaining all original packets. Based on the
13-second coarse pilot this selected set is expected to take minutes; that is
a rough estimate, not a benchmark for the higher orders. Impose a900-second
total child-computation budget, one native numerical thread and2 GiB per child.
No expensive preflight is needed for unchanged construction. Save each completed
case before continuing. On budget exhaustion retain it and reassess the remaining
comparison rather than restarting it or filling missing entries.

After this one set, inspect which change controls the reported response. If
resolution is sufficient for the proposed precision, prioritize a targeted
boundary/domain check and a regulator comparison. If it is not, target the
observed cause. No automatic Cartesian sweep, rigorous tail/Abel certificate,
physical scattering acceptance or full S11c-d completion follows here.
