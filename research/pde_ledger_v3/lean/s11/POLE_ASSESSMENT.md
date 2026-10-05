# Question 3: nonlinear-pencil residues, counts and projections

VC1–VC4 is complete at `9865050f`. The user authorized committing it after
both reviews cleared and proceeding to question 3. This assessment begins that
work under [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).

The resulting bounded NP1–NP4 increment is now complete: local verification
passes and both independent fidelity reviews are CLEAR. See
[POLE_FIDELITY_REVIEW.md](POLE_FIDELITY_REVIEW.md). The general analytic and
physical application hypotheses below remain separate from that finite proof.

The calculation session has already corrected the governing contract in
`directives/S11c_d_NONLINEAR_POLE_CONTRACT.md`, schema `nonlinearPoleV2`.
Its acceptance record is `S11c_d_nonlinear_pole_repair_checkpoint.json`;
the exact synthetic check records 89 satisfied checks. The original v10 shared
specification remains pinned during production. Its unrestricted §3b projector
wording is superseded by the addendum for new work. This Lean increment must
support that correction, not edit the pinned specification or repeat production.

## Which hypotheses support which conclusion?

| Object/conclusion | Required hypotheses and meaning |
|---|---|
| Meromorphic physical inverse near a finite-matrix zero | A holomorphic regular pencil on an open frequency domain, with determinant not identically zero. Work on a fixed sheet and away from the contour's singularities. |
| Meromorphic inverse for the actual nonlocal pencil | Fixed complex domain/codomain spaces and topology, an analytic Fredholm family of index zero and an invertible point on the connected domain, with finite-type singularities. An outgoing prescription alone supplies none of these domain assertions. |
| Semisimple residue | At one isolated zero, use full right and left nullspace bases V and W and an invertible full derivative pairing D=W L'(z*) V. Then R=V D^-1 W. A scalar pairing suffices only for nullity one. |
| Conditional local projections | P_X=R L'(z*) acts on fields X; P_Y=L'(z*) R acts on sources Y. These are idempotent under the full-pairing hypotheses. R:Y->X itself is not generally a projection. |
| Logarithmic-derivative contour | J=(2 pi i)^-1 integral L^-1 L' is an X->X map. For a positively oriented contour of winding one enclosing only that semisimple zero, J=P_X. For general nonlinear poles or multiple enclosed zeros it need not be idempotent. |
| Algebraic count | For regular analytic finite matrices, trace J counts enclosed determinant zeros with multiplicity. In infinite dimensions a finite-type reduction or an applicable characteristic-multiplicity theorem is needed; a pointwise operator trace need not exist. |
| Riesz projector | Supply a genuine fixed-operator spectral realization and an admissible resolvent contour. Its state space, domains and physical reconstruction maps are part of the claim. A compressed physical map need not remain idempotent. |
| Perturbative promotion | Analytic/Fredholm homotopy and a uniform contour inverse-error norm below one preserve the applicable total algebraic count. They do not alone preserve distinct-root count, nullity, semisimplicity or a linear pole-displacement rate. |

The finite-matrix full-pairing condition and residue formula are established in
[Schumacher, §7, Corollary 7.5 and Proposition 7.6](https://arxiv.org/html/2412.15985v1#S7).
The analytic existence theorem is a literature input to this assessment; it
will not be described as a proved Lean theorem merely because its algebraic
consequences compile. For operator pencils, see the fixed-space Fredholm
framework of [Beyn, Latushkin and Rottmann-Matthes](https://arxiv.org/abs/1210.3952).
The actual S11c realization still has to satisfy its hypotheses.

## Distinctions that need machine-checked controls

Write a local inverse principal part as

```
L(z)^-1 = sum_{j=1}^p C_-j (z-z*)^-j + H(z),  H holomorphic.
```

The residue is C_-1. It is only one coefficient. With the indicated analytic
expansions, multiplication gives the ordered local logarithmic coefficient

```
J = sum_{j=1}^p C_-j L^(j)(z*)/(j-1)!.
```

The order of the maps cannot be reversed. For observation O:X->V_o and forcing
B:U->Y, the physical response is O L^-1 B, not P_X applied to a Y-valued force.
At a double pole, even affine O and B give the residue

```
O0 C_-1 B0 + O1 C_-2 B0 + O0 C_-2 B1.
```

Thus freezing the forcing and observation at the pole loses genuine response.
These coefficient identities are algebraic consequences of expansions; their
use for an actual contour additionally needs those expansions and a suitable
holomorphic remainder. Nonlinear contour methods keep inverse principal parts
distinct from linear spectral projections; see [Beyn](https://arxiv.org/abs/1003.1580).

The existing calculation-side examples expose separate failure modes:

- L=z^2: inverse pole order two, residue zero, geometric multiplicity one and
  algebraic multiplicity two. The logarithmic contour is 2, whose square is 4.
- L=z^2-1: each single-root contour has scalar logarithmic count 1; the contour
  containing both simple roots gives 2, not a scalar projection. Semisimplicity
  at each root does not justify summing nonlinear modal projections globally.
- L=zI-N with N=[[0,1],[0,0]]: inverse z^-1 I+z^-2 N. Its residue and state
  Riesz projection are I even though the eigenvalue is defective. The nilpotent
  coefficient remains necessary; defectiveness does not make every contour
  object non-idempotent.
- The physical scalar transfer e1^T(zI-N)^-1 e2 is z^-2. Its residue vanishes
  despite a nonzero singular response. With O=2+5z and B=1+3z, the response
  residue is 11, while the frozen-map expression gives zero.
- Perturbing z^2 to z^2-epsilon splits the double root into +/-sqrt(epsilon).
  The total count stays two on a fixed admissible contour; the number of
  distinct roots changes and displacement need not be linear in epsilon.

The last item already has exact native evidence. This first Lean increment
will not add a general analytic homotopy or pole-displacement theorem.

## First bounded Lean increment

[POLE_COVERAGE.md](POLE_COVERAGE.md) specifies NP1–NP4: typed full-pairing
projection algebra and a conditional residue identification from leading
inverse-coefficient equations; finite Laurent/circle-moment identities; actual
small-pencil contour and response counterexamples; compact source matching and
meaningful controls followed by two independent fidelity reviews.

This is intentionally a finite mathematical core. It does not establish the
analytic Fredholm/Keldysh theorem, all root chains or partial multiplicities,
the argument principle for arbitrary operator pencils, a physical outgoing
realization, a pole search or the S11c bound-state tests. Those hypotheses and
their application remain explicit. The supplied full derivative-pairing inverse
cannot be presented as a proof that such an inverse exists for the physical
pencil. Likewise, coefficient data cannot be silently called an actual residue
without an inverse expansion or an explicit contour calculation.

All S11c sources, addendum and diagnostic records are read-only inputs. The
compact fidelity link will identify selected synthetic pencils and typed maps
already used in the correction, not reverify every native output. No previous
review-packet permission applies to the new NP1–NP4 packet.
