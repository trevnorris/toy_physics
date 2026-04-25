# Stage V2-30 - Electromagnetic Ontology and Status

## Purpose

This note records what the PDE/audit stack currently says about
electromagnetism, charge, circulation, and puncture orientation.

The short version:

```text
The program does contain an electromagnetic sector: a localized Maxwell field,
gauge-localization audit, mixed brane-bulk gauge invariants, charge variables,
and a circulation/fluxoid sector.  But electric charge and circulation are not
the same variable in the corrected ontology.
```

This matters for paper writing because older or informal language can make it
sound as if "charge = circulation" or "the puncture itself = electric charge."
The current stack is more precise.

---

## 1. What EM structure is present

The parent/reduced PDE materials include a Maxwell field

```text
A_M = (A_0, A_i, A_w),
F_MN = partial_M A_N - partial_N A_M,
```

with a localized kinetic weight `Z(w)` and source coupling to `J^M`.  The V2-02
audit checks the gauge-localization issue:

```text
S_EM ~ integral [-Z(w) F_MN F^MN/(4 mu0) - A_M J^M + gauge fixing].
```

The V2-09 audit then checks a reduced Maxwell/mixed block with variables

```text
Q, U, W.
```

Here:

- `Q` is the wall/worldtube amplitude;
- `U` is the localized brane-like Maxwell coordinate;
- `W` is the mixed `A_w/F_{mu w}/J^w` active coordinate.

The mixed gauge-invariant observables are

```text
E_w = -partial_t A_w - partial_w A_0,
C_a = partial_a A_w - partial_w A_a.
```

These are not gauge artifacts.  They are the tensor slots where brane-bulk
electromagnetic exchange, mixed-sector plumbing, and hidden-port response can
enter.

---

## 2. Electric charge variables

The corrected charge dictionary uses

```text
eta_Q = +/- 1,
q_* = eta_Q e_*,
q_eff = q_*/sqrt(Z_int).
```

Interpretation:

- `eta_Q` is the electric branch sign, tied to puncture orientation in the
  carried ontology;
- `q_*` is the microscopic branch charge;
- `q_eff` is the brane-observed/localization-dressed charge after zero-mode
  normalization;
- `Z_int = integral Z(w) dw` is the localization normalization entering the
  effective brane charge.

So the paper-safe phrasing is:

```text
The electric charge sign is carried by the oriented puncture branch eta_Q, and
the observed charge is the localization-dressed q_eff.
```

Do not write the stronger shorthand:

```text
The puncture is electric charge.
```

The puncture/throat is the physical defect.  Electric charge is one branch
label and coupling carried by that defect.

---

## 3. Circulation and magnetism

The corrected ontology keeps circulation separate from electric charge.

The circulation/fluxoid sector is the topological law for loops surrounding the
mouth.  In the charged superfluid notation it has the schematic form

```text
integral_C (partial_i theta - q_* A_i/hbar) dl^i = 2 pi n.
```

This law quantizes a tangential vortical/holonomy class.  It belongs to the
magnetic/vortical sector, not to the electric-charge dictionary.

Paper-safe phrasing:

```text
Circulation around or through the throat belongs to the magnetic/vortical
sector.  It is the natural place to encode magnetic/fluxoid topology, while
electric charge remains the separate puncture-orientation branch eta_Q with
microscopic charge q_*.
```

Avoid the stronger claim:

```text
Circulation is electric charge.
```

Also avoid claiming that the audit already derives all of magnetism from throat
intake circulation.  The files contain the Maxwell field, fluxoid/circulation
law, and mixed brane-bulk slots, but the full recirculation/plumbing law is not
closed yet.

---

## 4. What is reduced versus still open

The following pieces are present and audit-backed:

- localized Maxwell action and zero-mode normalization issues;
- finite gauge-localization warnings and patches;
- mixed gauge invariants `E_w` and `C_a`;
- reduced Maxwell/mixed kernel stability and transfer gates;
- conservative Maxwell/mixed moments `Z0,Z2,Z4`;
- charge dictionary `eta_Q`, `q_*`, `q_eff`;
- circulation/fluxoid sector as magnetic/vortical, not electric-charge
  defining.

The following pieces remain open:

- a complete moving-throat derivation of the source current `J^M`;
- a closed recirculation/plumbing law connecting exterior circulation, throat
  intake, mixed `A_w/F_{mu w}/J^w` transport, and brane magnetic fields;
- a derivation of the charge magnitude `e_*` rather than treating it as a
  branch parameter;
- a proof that the localized Maxwell sector itself emerges from the superfluid
  variables rather than being included as a parent sector;
- a completed same-charge internal twist/spin-like discretizer, if the lepton
  branch needs one.

---

## 5. Recommended paper paragraph

```text
The electromagnetic sector in the current PDE stack is represented by a
localized Maxwell field A_M with transverse localization profile Z(w), together
with mixed brane-bulk gauge-invariant components such as E_w and C_a.  Electric
charge is not identified with circulation.  The corrected charge dictionary
assigns the electric branch sign to the oriented puncture label eta_Q, with
microscopic charge q_* = eta_Q e_* and brane-observed charge
q_eff = q_*/sqrt(Z_int).  Circulation and fluxoid quantization around the mouth
belong to the magnetic/vortical sector and must be kept distinct from radial
throughput and mouth-output flux.  The present audit verifies the reduced
Maxwell/mixed algebra and gauge-localization consistency, but a full
moving-throat recirculation law deriving all electromagnetic response from the
superfluid branch remains future work.
```

---

## 6. Minimal claim boundary

Supported:

```text
The PDE program includes a localized Maxwell/mixed sector and a corrected
charge/circulation ontology.
```

Conditionally supported:

```text
Puncture orientation carries the electric charge sign eta_Q.
Circulation/fluxoid data belong to the magnetic/vortical sector.
```

Not yet supported as a completed theorem:

```text
The current PDE derives the whole of electromagnetism from superfluid
circulation entering the throat.
```
