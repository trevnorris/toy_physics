# One contact followed through breakup at 0PN

**One modeled junction delivered 4.58643% of its entry prestrain energy to the ordered material before its last attachment broke. By that instant 3.93242% had crossed a section 0.25 model lengths into the receiver.** The held radius-one mouth supplies the four-dimensional carrier and order fields at intake 0.25 and half the reference conversion drive. For the contact and ordered receiver we used the following local material hypothesis. [Earlier contact test](mouth_contact_tracking.md)

**What counts as a break.** A tagged wall site connects to the same ordered site at two depths 0.04 apart. The relaxed sideways offsets of the attachments are fixed at opposite values of magnitude \(\epsilon=0.001\). In a perpendicular direction each tagged feature lies one eighth of a pitch from the center of its partner’s capture window. The window width is the local order \(\chi\) times the pitch. As the window shrinks past the feature, that attachment disconnects. Its remaining elastic energy is recorded separately. The model leaves its subsequent destination open. [Capture rule](mouth_contact_crossing_force.md)

| Stage, finer meshes | Time since entry | Order at junction center | Zero-force sideways offset, in units of ε |
|---|---:|---:|---:|
| Both attachments engaged at entry | 0 | 0.9500 | −0.0261 |
| Half order | 2.20594 | 0.5000 | −0.2477 |
| First attachment lets go | 2.62026 | 0.3563 | −0.3138 just before, −1 just after |
| Last attachment lets go | 3.31716 | 0.1671 | −1 just before, force zero at every offset afterward |

At each reading, two displacements of a frozen copy bracketed zero holding force. On the undisturbed trajectory a uniform elastic strip represented the ordered material. Its endpoint velocity evolved under the contact force and neighboring material stress. Integrating their signed product \(F_{\rm contact}v_{\rm end}\) gave \(2.16448\times10^{-8}\) model energy per unit area at breakup. The receiver’s kinetic and elastic energy independently matched the integral. Declining stiffness removed 70.61555% of the entry elastic energy and severance removed 24.79802%. The propagating wave continued across the interior section after the contact stopped doing work. [Travel and work plot](figures/mouth_tagged_contact.png) · [Full signed history](../../../research/pde_ledger_v3/_measurements/mouth_tagged_contact_history.jsonl.gz)

**What determines the answer.** The main case uses live stiffness \(k_i=\chi_i/2\). Holding each live stiffness at 1/2 until breakup gave **25.08964%** of that preparation’s entry energy at the same two break times. With the original stiffness rule, moving the tag to a quarter pitch gave 10.71722%. A centered tag stayed connected through the last observation. Freezing conversion at entry gave zero work for the initially balanced junction. Measure stiffness by the force-versus-slip slope on surviving individual contacts and measure their capture width independently. The full mouth also requires a measured connection to the feeding brane and a source of entry strain for sustained work.

The junction solver accepts `--tagged-readouts DIRECTORY` and an optional `--registry-offset PITCHES` to select one tag. Two independent reviewers checked the main results against outgoing-wave solutions and direct integration of the raw histories. Working base `49dca23e`. Formal review and commit are pending. [Solver entry point](../tools/mouth_contact_junction.py) · [Measurements and review scopes](../../../research/pde_ledger_v3/_measurements/mouth_tagged_contact.json)
