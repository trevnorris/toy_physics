## Verdicts

### Corrected physics and measurements
**Pass.** From the parent laws, the Noether balance gives body force \(f_i^{V_{\rm conf}}=-n\partial_i V_{\rm conf}\) and sideways momentum flux \(\Pi_{xz}\) with carriage \(m n v_x v_z\) plus mixed density stress. On the periodic neutral slab, steady \(\partial_z\int\Pi_{xz}=0\) plus a flat feed with \(v_x=U\) and uniform \(n\) fix the feeding-side aperture arrival rate as incoming momentum \(L_x=mJU\), not a support reaction. The executed signed budgets match that identity to the stated tolerances (\(L_x/U\approx0.25\), aperture vs feed \(\sim10^{-8}\), layers \(\sim10^{-5}\), odd in \(U\)).

Common differentiable phase with continuous normal number current sends the carried sideways jump \(mj(v_x^--v_x^+)\) through the cut. Any extra localized force needs stress jump or a separate coupling. The conditional real surface energy recovers net force \(-\partial_X\int e\) after tangential divergence; the \(\lambda=1\) restoring piece is \(\sim3\times10^{-5}\) of \(mJ\). Parent confinement at fixed reference shape yields only the distributed dipole load \((4\pi/3)\int r^2 n_1\partial_r V_{\rm conf}\,dr\), collar work \(\to0\), and a separate wall-endpoint traction \((4\pi/3)T_w\partial_w q\); the archived intake solves do not supply that endpoint response. Restoring virtual work on the unequal pair is \(\ll\) the recorded incoming loads and is an imposed density-coupling force, already inside the feeding-side stress bookkeeping.

**Physical force vs incoming momentum:** \(L_x\) is arrival rate of sideways carrier momentum at the feeding face. Holding force \(H_x\) (or passed-on \(E_x\)) is a different object; none of the four scripts invents \(H_x\), stored interface momentum, or acceleration.

### Four scripts
**Pass, with prior mutation gaps closed.**

| Script | Role | Check |
|---|---|---|
| `held_crossflow_budget.py` | Signed \(L_x\) feed vs aperture | Oddness, rate, aperture\(\approx\)feed; does not claim hold |
| `mouth_confinement_projection.py` | Parent \(V_{\rm conf}\) + wall endpoint | Volume, reaction, \(n_1\)/\(n_3\) coeffs, collar, \(T_w\) traction |
| `mouth_transmission_work.py` | Phase-continuous cut + conditional \(e\) | \(j\) continuous, stress identity, \(\lambda\in\{0,1\}\), tiny restoring |
| `unequal_mouth_boundary_work.py` | Restoring \(-\partial_X E_R\) + mixed \(\delta j\) | Energy FD match; fixed-rate vs fixed-profile mixed laws; no single scalar \((u,\theta)\) energy |

M1–M5 now reject. M6 separates fixed-support stencil change from support change; full-plane adjoint equals direct; aperture uses the program mask; changing support alone keeps periodic summation by parts. The old M2 \(2U^2\) false pass and M6 support/stencil mix are corrected.

### Four prose sections
**Pass, scoped to feeding-side and conditional identities.**

1. **Mouth equations / restoring:** Restoring is the imposed \(g_X\) coupling force, continuum-IBP piece of brane stress, \(\lt0.01\%\) of paired incoming load; mixed \(\delta j\) rules out one scalar surface energy; receiving traction left open.
2. **Signed one-mouth test:** Derives and tables \(L_x=mJU\); states \(H_x=E_x-L_x\) and that \(E_x'\) (or hold slope) is the missing datum.
3. **Phase-continuous transmission:** Carried momentum passes; conditional \(e\) gives the small restoring numbers; physical functional remains open.
4. **Parent confinement:** Matches the projection formula, collar \(\to0\), endpoint \(T_w\) traction; archived \(U\) solves give no coefficient for that channel.

No section equates incoming momentum with mechanical hold or asserts a global receiver solution.

### Register requirement
**Pass as correctly open.** The single 0PN held-mouth response \(e_tD_t=-\partial_U H_x(a,0)\) remains undetermined. That matches `throat_surface_functional` / `OPEN_INPUT`: feeding-side \(L_x\), tiny restoring, and distributed confinement are fixed or bounded by the present evidence; receiving tangential momentum, rim/endpoint reaction, and signed hold readings are not.

### Material open findings
None against these corrected claims. The open item is the already-registered missing hold / surface-functional datum. Existing executed results and controls are enough for this scoped verdict; no further witness is required for what these four calculations assert.
