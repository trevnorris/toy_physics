# Finite continuum coefficient response

The constructor combines the accepted four interior matrices with actual
boundary/current coefficients. It solves both incident ends, all open incoming
directions and all outgoing/evanescent matching directions, including the full
mixed inverse terms. No new quadrature or mode computation is required.

Focused noncommuting matrix controls pass: coefficient recovery1.14e-15,
independent solve3.16e-15 and positive current square-root recovery8.89e-15.
Mixed-forcing and matrix-order mutations respond; the field-coordinate phase
control rejects applying a phase in the wrong basis. These are instrument
checks, not an accepted physical response.

Production is running under the local completion/error hook with one native
thread,2GiB and900seconds. It saves
all systems, solutions, currents, phase/normalization maps and formal truncation
diagnostics before output replay. Approximate modal boundaries and finite
regulator remain explicit. Open-current bookkeeping is distinct from thickness
conversion and bulk escape, which remain work alongside physical controls,
frequency poles, remaining cases and exports. No broad quadrature is queued.
