# Where does a converting contact’s work go?

**A contact made of two elastic attachments can shift its force-free position and do work even when each attachment keeps its original relaxed length. In this test all the delivered work came from elastic energy already stored at entry.** We used the saved flow and shear order of the held radius-one mouth at intake 0.25. The attachments and receiver are hypotheses. [Earlier conversion-energy test](mouth_conversion_contact_transfer.md)

**What we could track.** On the finer mesh at half the reference drive, two neighboring carrier markers traversed 95% to 1% shear order in 5.29667 model time units. Their spacing at fixed angular separation grew by 1.1153% as the wall widened. The archive tracks present spacing. It has no individual contact identity, unloaded position or breakup reading. [Carrier reconstruction](../tools/mouth_contact_tracking.py) · [Earlier contact-survival tests](mouth_contact_observation.md)

**A test without a prescribed change in relaxed length.** Suppose one junction contains two attachments carrying opposite initial shear. Their preferred sideways offsets stay fixed at \(p_i=\pm\epsilon\). Both obey the same stiffness law \(k_i=k_0\chi_i/2\). We placed them a trial depth 0.04 apart and read the combined zero-force offset from two trial holding forces. It is \(x_0=\sum k_i p_i/\sum k_i\). As one attachment softened faster, the balance shifted. On the fine half-drive run the shift was \(0.4740\epsilon\) in magnitude between 95% and \(10^{-6}\) order. Reversing the initial offsets reversed the shift.

**Work reaching an assumed ordered attachment.** At 95% order the junction started in force balance. With the wall held fixed, its other end could launch an outgoing wave into a uniform ordered medium of impedance \(Z=1\), with \(k_0=1\). The force on that medium was \(F=\sum k_i(p_i-x)=Z\dot x\). The table gives integrated work per unit area as a fraction of entry prestrain energy. [Force and work calculation](../tools/mouth_contact_junction.py)

| Preparation, finer mesh | Work delivered by 1% order | Work delivered by \(10^{-6}\) order |
|---|---:|---:|
| Half reference drive, attachments 0.04 apart | 1.86817% | 1.86881% |
| Zero external conversion drive, same attachments | 1.88218% | 1.88286% |
| Half reference drive, both soften together | 0% | 0% |

Removing the initial strain also gave zero work. Halving and doubling the assumed impedance gave 1.92427% and 1.41865% at the final reading. The energy check was \(\dot E=\tfrac12\sum_i\dot k_i(x-p_i)^2-F\dot x\). Both stiffnesses decreased. In the main case the residual energy could supply at most another \(1.022\times10^{-6}\) of entry energy under continued passive softening without reloading.

**What would settle it for the mouth.** Tag both sites until separation. On matched copies measure force versus slip through zero load before renewal, with a nonconverting control of equal age and prestrain. On the undisturbed path integrate contact force times ordered-end velocity until separation. A distant feeding-brane probe also needs a travel-time allowance and background subtraction. Sustained 0PN motion requires a source of fresh prestrain. [Measurements and independent checks](../../../research/pde_ledger_v3/_measurements/mouth_contact_tracking.json)

Exploratory result on base `49dca23e`. Two local reviewers independently reproduced the reported calculations within their recorded scopes. Formal review and commit remain pending.
