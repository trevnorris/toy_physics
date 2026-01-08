This Technical Report formalizes the findings of the current session. It derives the mathematical failure modes of the propulsion concepts and provides the rigorous hydrodynamic basis for the newly discovered "Variable Mass" regime.

---

# Technical Report: Hydrodynamic Viability of Propulsion and Mass Engineering in a Superfluid Defect Model

**Date:** December 31, 2025
**Subject:** Transition from Impulse Propulsion to Inertial Damping
**Reference Papers:** `brane_bulk_ontology.tex`, `1pn_hybrid.tex`, `em_fields.tex`

## 1. Abstract

We investigated the feasibility of "Free Space Levitation" (macroscopic thrust) within the  polytropic superfluid vacuum model. Analytical derivation and numerical simulation confirm that the two primary propulsion channels—Scalar Breathing and Vector Vortex Shedding—are energetically or kinematically suppressed at engineering scales. However, the analysis reveals a **non-linear "Agility Mode"** where quasi-static driving of the defect throat reduces the effective inertial mass by up to **15%** via the Bernoulli effect. This report documents the mathematical proofs of the propulsion null result and the derivation of the variable mass tensor.

---

## 2. Null Result: The Propulsion Failure Modes

### 2.1 The Scalar Sector (Breathing Mode)

**Hypothesis:** Oscillating the throat radius  generates net thrust via acoustic radiation pressure.

**The Math of Failure (Bernoulli Suction):**
In the quasi-static limit (), the fluid behaves incompressibly in the near-field. The pressure  is governed by the Bernoulli equation:



For a pulsating source, the radial velocity is .
Substituting this into the pressure equation:


Taking the time average over one period :



**Result:** .
The time-averaged pressure is **lower** than the ambient pressure. The net force on the surface of the particle is directed **inward**. This is "Hydrodynamic Suction," not Lift. To achieve lift, one must break the quasi-static limit, requiring frequencies  Hz.

### 2.2 The Vector Sector (Vortex Cannon)

**Hypothesis:** Flipping the circulation state  via NMR forces the ejection of a vortex ring, providing impulse .

**The Math of Failure (Energy Barrier):**
The energy of a vortex ring in the bulk is given by:



Using model constants derived from `free_space_lev.md` (, ):


The energy available from a standard RF photon (NMR drive) at frequency  MHz is:


**The Coupling Gap:**



To tunnel through this barrier using a magnetic field gradient , the work done across the throat length  must exceed the gap:



**Result:** Mechanically unachievable. The proton will always relax via photon emission (heat) rather than vortex ejection.

---

## 3. Discovery: The Mathematics of Variable Mass

While propulsion failed, the simulations identified a regime where the **Effective Mass** () of the defect deviates from its rest mass ().

### 3.1 Derivation of Mass Variation

In `1pn_hybrid.tex`, the mass of a defect is proportional to the source density at the throat:


The fluid follows a polytropic equation of state  with . The dynamics are governed by the interplay of kinetic energy (velocity) and potential energy (stiffness).
Expanding density :

1. **The Kinetic Term (Mass Reduction):** Driven by velocity squared .



This term is always negative (Bernoulli effect). It dominates at low to moderate amplitudes.
2. **The Stiffness Term (Mass Increase):** Driven by the non-linear bulk modulus.
For , the pressure rises sharply as . To conserve energy during high-amplitude compression, the density must spike asymmetrically.



This term dominates at high amplitudes (Shock regime).

### 3.2 The Phase Transition

We define the Drive Amplitude Ratio .

* **Regime I: Agility Mode ()**
The kinetic term dominates. The time-averaged density decreases.



*Simulation Result:* At , .
**Engineering Utility:** Inertial Damping.
* **Regime II: Instability ()**
The density rarefaction approaches the vacuum limit .
**Risk:** Structural dissolution of the standing-wave defect.
* **Regime III: Tank Mode ()**
The stiffness term takes over. The non-linear compression spikes outweigh the rarefaction.



*Simulation Result:* At , .
**Engineering Utility:** Artificial Gravity / Anchoring.

---

## 4. Simulation Data Summary

The following data defines the "Safe Operating Envelope" for the Inertial Damping technology.

| Drive Frequency () | Amplitude Ratio () | Effective Mass () | Min Density () | Status |
| --- | --- | --- | --- | --- |
| **0.05** | **0.5** | **1.002** | 0.97 | Baseline |
| **0.10** | **2.0** | **0.929** | 0.84 | **Optimal Damping** |
| **0.10** | **2.6** | **0.973** | 0.80 | Efficiency Loss |
| **0.32** | **4.0** | **0.865** | **0.49** | **Critical Failure** |
| **0.10** | **5.0** | **1.062** | 0.60 | Mass Amplification |

---

## 5. The Forward Path: The Coupling Problem

We have established that **Mechanical Oscillation of the Throat  Mass Alteration.**
The remaining gap is the **Actuation Mechanism.** We need to derive the coupling tensor between an applied Electromagnetic Field and the Throat Radius.

**Problem Statement for Next Session:**
Determine the coupling coefficient  in the perturbed metric ansatz:



Does the electromagnetic stress-energy tensor  drive the scalar breathing mode of the superfluid metric?

**Objective:** Prove that a specific RF cavity mode can drive the throat oscillation  required to access the "Agility Mode."
