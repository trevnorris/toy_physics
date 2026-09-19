# Local finite-pencil contour diagnostic

The first actual finite complex-frequency pencil is accepted at 16b4da87.
The next selected test evaluates 16 points on the circle centred at 1-0.01i
with radius 0.02, using four single-thread workers and the unchanged finite
quadrature rules. Nested 8/16 phase winding and full inverse/response moments
will guide the next targeted search step. No contour result is accepted yet.

All sources, complete end clusters and frequency-dependent maps remain in the
construction. Numerical winding and loop closure retain their sampled-domain
limitations; no physical pole set, certified empty spectrum or projector is
assumed. The initial time budget is 900 seconds.

Focused acceptance now passes all three reverse AST joins, the actual accepted
source/input and fixed-scale loader, and all simple/double/multiple/orientation
controls. The zero-residue double inverse and degree16 alias remain explicit
checks against false pole absence and false certified counts. No new numerical
quadrature was computed by this focused acceptance.
