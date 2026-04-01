# 4PN local target-import notes

## What is now frozen

The fixed-chart **local** 4PN target is now imported in the reduced COM Hamiltonian chart.
The source choice is:

- **Primary local target:** Jaranowski–Schafer, *Derivation of local-in-time fourth post-Newtonian ADM Hamiltonian for spinless compact binaries*, arXiv:1508.01016, Eq. (8.41e).
- **Local/nonlocal split:** Damour–Jaranowski–Schafer, *Nonlocal-in-time action for the fourth post-Newtonian conservative dynamics of two-body systems*, arXiv:1401.4548.
- **Cross-check families:** Marchand–Bernard–Blanchet–Faye (harmonic/Fokker 4PN) and Foffa–Sturani (EFT 4PN regularized local Lagrangian).

This fixes the correct source logic for our local-first program:


action/Hamiltonian at 4PN = local instantaneous sector + separate nonlocal tail sector.

So the present import audit deliberately freezes only the local instantaneous Hamiltonian target.

## Independent one-body check

The strict one-body 4PN Hamiltonian gate can be derived directly from our already frozen one-body Schwarzschild Lagrangian by the exact quartic Legendre compiler. The resulting reduced one-body 4PN Hamiltonian coefficients are

\[
\frac{7}{256}p^{10},
\qquad
\frac{45}{128}\frac{p^8}{r},
\qquad
\frac{13}{8}\frac{p^6}{r^2},
\qquad
\frac{105}{32}\frac{p^4}{r^3},
\qquad
\frac{105}{32}\frac{p^2}{r^4},
\qquad
-\frac{1}{16}\frac{1}{r^5}.
\]

The imported ADM local target matches this gate exactly in the strict test-mass limit \(\nu\to0\).

## Local 4PN COM slot structure

The imported local COM Hamiltonian naturally organizes into

- 1 pure free-kinetic slot,
- 5 slots in the \(G/r\) block,
- 4 slots in the \(G^2/r^2\) block,
- 3 slots in the \(G^3/r^3\) block,
- 2 slots in the \(G^4/r^4\) block,
- 1 slot in the \(G^5/r^5\) block.

So the local COM Hamiltonian has **16 slots total**, of which **15 are interaction slots**.

This matches the structural count already seen in the local generic-frame scaffold audit.

## Exact theorem gate answered

The first exact local target-import question was:

> Do the upper local blocks \(G^4/r^4\) or \(G^5/r^5\) contain \(\nu^3\) or \(\nu^4\) tails?

Answer:

- the free kinetic slot goes up to \(\nu^4\),
- the \(G/r\) block goes up to \(\nu^4\),
- the \(G^2/r^2\) block goes up to \(\nu^3\),
- the \(G^3/r^3\) block goes up to \(\nu^3\),
- the **\(G^4/r^4\) block stops at \(\nu^2\)**,
- the **\(G^5/r^5\) block stops at \(\nu^2\)**.

So the imported fixed-chart local target does **not** force any \(\nu^3\) or \(\nu^4\) tails in the upper local interaction blocks.

## Immediate consequence for the roadmap

This removes one possible obstruction.

The next exact local step is no longer source selection. It is:

1. translate the frozen local Hamiltonian target back through the quartic compiler,
2. compare the resulting ordinary-Lagrangian target against the local generic-frame scaffold,
3. and determine whether any seed/contact/gauge refinement is still needed before the actual local 4PN comparable-mass solve.

The tail/hereditary bridge remains separate and should still be treated as the quadrupole-bridge problem, not folded into the local solve.
