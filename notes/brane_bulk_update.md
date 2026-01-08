### Instructions for Updating Paper V: `brane_bulk_ontology.tex`

**Goal:** Transform Paper V from a speculative ontology paper into a rigorous derivation of **Magnetostatics (The Biot-Savart Law)**. This paper must solve the "Vanishing B-Field" limitation identified in Paper IV by using the Brane-World topology.

#### **1. The Narrative Shift: Solving the "Missing Link"**
* **Current State:** Paper IV (`em_fields.tex`) successfully derived Electrostatics (Coulomb’s Law) and the Hierarchy Problem ($1/a^2$) but admitted that the irrotational bulk vacuum cannot support magnetic fields ($\nabla \times \mathbf{v} = 0$).
* **New Purpose:** Paper V must present the "Brane-World" scenario not just as an interesting idea, but as the **necessary topological patch** that allows magnetic fields to exist.
* **The Hook:** "While Gravity lives in the irrotational bulk, Magnetism lives on the vortical brane."

#### **2. Core Physics Updates**
You need to implement the following physical model to derive the Biot-Savart Law:

* **The "Vortex Sheet" Model:**
    * Define the vacuum not just as a 3D bulk, but as a bulk containing a 3-brane (hypersurface).
    * **The Bulk:** Remains irrotational ($\nabla \times \mathbf{v} = 0$) to preserve the scalar gravity results of Papers I-III.
    * **The Brane:** Acts as a "Vortex Sheet" or phase boundary that *can* support non-zero vorticity (surface currents).

* **The Mechanism (Drag = B-field):**
    * Model a defect (particle) as an intersection or topological defect constrained to this brane.
    * When the defect moves with velocity $\mathbf{u}$ along the brane, it drags the local texture/phase of the brane surface.
    * **Derivation:** Show that this drag creates a "wake" of vorticity on the surface.
    * **The Result:** Identify this induced surface vorticity with the magnetic field $\mathbf{B}$. Show that for a steady current (moving line of defects), this vorticity creates a field that falls off as $1/r$, reproducing the **Biot-Savart Law**.

#### **3. Parameter Consistency Checks**
You must update the mathematical parameters in this paper to ensure it matches the successful calibration of Papers I–III.
* **Fix $\beta$:** Set $\beta = 3$ (matches 1PN perihelion precession).
* **Fix $\alpha^2$:** Set $\alpha^2 = 3/4$ (matches the EIH Lagrangian and vector gravity wake).
* *Note:* Do not use the generic or "wrong" values currently in the draft. The magnetic derivation must be compatible with these specific gravitational constants.

#### **4. Structural Changes to the Manuscript**
Please reorganize the sections to follow this logical flow:

* **Introduction:**
    * Recap the success of Paper IV (Charge is geometric).
    * State the problem: Bulk fluid cannot support B-fields.
    * Propose the solution: The Vacuum is a Brane.
* **Section 2: Brane Topology:**
    * Define the Bulk vs. Brane fluid properties.
* **Section 3: Induced Vorticity (The Meat):**
    * Perform the calculation of a defect moving on a superfluid membrane.
    * Derive the vorticity wake.
    * Map $\boldsymbol{\omega}_{surface} \to \mathbf{B}$.
* **Section 4: Recovering Biot-Savart:**
    * Show that the force between two parallel moving defects (currents) on the brane reproduces the Ampère force law via this induced vorticity.
* **Conclusion:**
    * Summarize that we now have Electrostatics (Bulk) and Magnetostatics (Brane).
    * Tease Paper VI (`1pn_hybrid.tex`) as the final step for Electrodynamics (Light Waves).

#### **5. Key "Takeaway" Sentence**
Ensure the abstract includes a sentence similar to this:
> "By restricting vorticity to a codimension-1 brane embedded in an irrotational bulk, we recover the Biot-Savart law and magnetic forces on the brane without disrupting the scalar PPN parameters of the gravitational sector."
