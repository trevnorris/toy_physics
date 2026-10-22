# Gravity attribution pass (2026-10-22): queries, sources and quoted support

This file grounds the section "ATTRIBUTION PASS — the gravity sector and S16 (the worldtube)" in
`docs/medium_requirements_and_prior_art.md`. Row IDs `G1`–`G13` are that section's table rows. A claim in the
section without an entry here should be discarded.

**Round 1 repair (2026-10-22).** Review round 0 (`gravity_prior_art_review_r0_disposition.md`, findings R1–R7)
changed the entries for [S2] (full text opened), [S5], [S9], [S13], [S15] and [S18], and one row of the leads
table. The G5 and G6 relabelling (R2) also uses quotes already filed under [S6] and [S8]. Each new quote below was printed from a file I downloaded myself; the reviewers' files were used only to
find the sources.

## How the material was gathered

- **Tools.** WebSearch (standard and extended) and WebFetch found the sources. Full texts were downloaded with
  `curl` and extracted with `pdftotext` (poppler); quotes were then printed with `grep -n` and `sed -n`. The
  Crossref REST API supplied bibliographic data and abstracts, and the Wikipedia API supplied wikitext at a fixed
  revision. Every command here only downloads, extracts or prints text; none computes anything.
- **Scholar Gateway** (`semanticSearch`): two queries (analogue-gravity sink forces; Bjerknes forces as analogue
  gravity). Both returned `ACCESS_DENIED: Access denied`, so it was not used.
- **Downloads** are kept outside the repository at `/var/projects/_review_scratch/gravity-prior-art-search/`, one
  subdirectory per source (`dl_*`). The raw search log is `search_log.md` there.
- **Quoting rules.** Quotes are verbatim from the retrieved text. OCR texts (archive.org, the Zenodo scan of
  Lagally, the Thorne–Hartle scan) keep their OCR errors; a reading in `[square brackets]` is mine. Line breaks
  inside a quote are joined with one space. The doc section shortens some quotes (marked "…"), rejoins words split
  by a line-break hyphen (`pon- derablen` → `ponderablen`, `Körper- welt` → `Körperwelt`), and restores a hyphen that extraction dropped
  (`oneparticle` → `one-particle`). It changes no other character.
- **WebFetch summaries** were used only to locate sources. Where only a summary was available, the entry says so,
  and the claim is marked "secondary only" or "not opened" in the section.

## Queries run

WebSearch (S = standard, X = extended):

| # | query | mode |
|---|---|---|
| Q1 | `Einstein Infeld Hoffmann 1938 "The gravitational equations and the problem of motion" Annals of Mathematics` | S |
| Q2 | `Lagally theorem force on source sink proportional density strength velocity ZAMM 1922` | S |
| Q3 | `G. I. Taylor 1928 "forces on a body placed in a curved or converging stream of fluid" Proceedings Royal Society` | S |
| Q4 | `Bjerknes pulsating spheres hydrodynamic analogy gravitation "Fields of Force" 1906` | S |
| Q5 | `Lagally "Berechnung der Kräfte und Momente, die strömende Flüssigkeiten auf ihre Begrenzung ausüben" Zeitschrift für Angewandte Mathematik und Mechanik 2 409` | S |
| Q6 | `Papapetrou 1951 "Spinning test-particles in general relativity" Proceedings Royal Society A 209 248` | S |
| Q7 | `Mathisson 1937 "Neue Mechanik materieller Systeme" Acta Physica Polonica republished General Relativity and Gravitation 2010` | S |
| Q8 | `Thorne Hartle 1985 "Laws of motion and precession for black holes and other bodies" Physical Review D 31 1815` | S |
| Q9 | `Dixon 1970 "Dynamics of extended bodies in general relativity" Proceedings Royal Society A 314 499 momentum angular momentum` | S |
| Q10 | `Damour "effacement" internal structure problem of motion Newtonian Einsteinian gravity "300 Years of Gravitation"` | S |
| Q11 | `Einstein-Infeld-Hoffmann method equations of motion "surface integrals" surrounding each singularity field equations alone` | S |
| Q12 | `Riemann 1853 "Gravitation und Licht" ether flows into ponderable atoms gravitation Naturphilosophie` | S |
| Q13 | `Riemann "Neue mathematische Principien der Naturphilosophie" "ponderable Atom" Stoffmenge verschwindet Gravitationskraft` | S |
| Q14 | `Riemann 1853 ether "enters" "ponderable atoms" "disappears" gravitation "proportional to" translation Collected Papers Baker Christenson Orde` | S |
| Q15 | `Lamb Hydrodynamics pulsating spheres Bjerknes attract "same phase" inverse square article` | S |
| Q16 | `Mechanics of variable-mass systems Part 1 balance of mass and linear momentum Irschik Holl Applied Mechanics Reviews 2004 control volume` | S |
| Q17 | `Kirkwood 1953 "physical basis of gravitation" Physical Review ether flow into matter` | S |
| Q18 | `Landweber Yih 1956 "Forces, moments, and added masses for Rankine bodies" Journal of Fluid Mechanics 1 319` | S |
| Q19 | `analogue gravity acoustic black holes "equations of motion" sinks attract each other Newtonian force superfluid drain` | X |
| Q20 | `force on a quantized vortex from momentum balance in a cylinder surrounding the core independent of core structure Sonin superfluid Magnus force derivation` | S |
| Q21 | `sink flow model of matter gravity hydrodynamics inverse square derivation two sinks attract fluid ether paper` | X |
| Q22 | `Newtonian equations of motion extended bodies integrate Euler equation over body center of mass multipole moments tidal moments self-force vanishes Poisson Will Gravity chapter 1` | S |
| Q23 | `Bjerknes force analogue gravity bubbles "analogue" Newtonian gravity hydrodynamic arXiv` | X |
| Q24 | `"problem of motion" OR "worldtube" analogue gravity acoustic sink moving vortex body equations of motion emergent gravity fluid` | X |
| Q25 | `Soldner 1801 deflection of light Berliner Astronomisches Jahrbuch 1804 0.84 arcseconds half Einstein value` | S |

Other lookups:
- **Crossref** `api.crossref.org/works?query.bibliographic=…` for Landweber–Yih, Irschik–Holl, Auton–Hunt–Prud'homme,
  and Taylor's 1928 added-mass paper. `api.crossref.org/works/<DOI>` for every DOI cited below.
- **archive.org** `advancedsearch.php` for `title:(hydrodynamics) AND creator:(lamb)` and
  `title:(fields of force) AND creator:(bjerknes)`.
- **Wikipedia API** (`action=query&prop=revisions&rvprop=content|ids|timestamp`) for the articles named below.

Reproducing a Crossref lookup:
```
curl -s 'https://api.crossref.org/works/<DOI>' | python3 -I -c 'import json,sys; m=json.load(sys.stdin)["message"]; print(m["title"], m.get("container-title"), m.get("volume"), m.get("page")); print(m.get("abstract",""))'
```

---

## Sources

### [S1] O. Reynolds, *Papers on Mechanical and Physical Subjects*, vol. 3, *The Sub-Mechanics of the Universe* (Cambridge, 1903)
- **Rows:** G1, G2.
- **Found by:** the Wikipedia article *Reynolds transport theorem* (revision 1376375962, 2026-09-23), whose
  reference reads `|year=1903 |title=Papers on Mechanical and Physical Subjects |volume=3, The Sub-Mechanics of the Universe |publisher=Cambridge University Press |location=Cambridge |pages=12–13`.
- **Opened:** `https://archive.org/download/papersonmechanic03reynrich/papersonmechanic03reynrich_djvu.txt` (OCR).
- **Retrieval:** `curl -sL <url> -o r.txt; sed -n 1594,1627p r.txt; sed -n 1146,1156p r.txt; sed -n 1696,1700p r.txt`
- **Quotes:**
  - l.1594: `14.     Moving  Surface.`
  - l.1601–1603: "If, however, the motion of every point of the surface be taken into account it is possible to sum
    the results of equations (7), (8), (9) over the space enclosed by a surface in any manner of continuous motion."
  - l.1605–1607: "Putting u, v, w for the component velocities of the surface at the point x, y, z, then the component
    motions of the entity represented by Q relative to the surface at this point are respectively u — u, v — v, w — w"
    [OCR lost the bar that marks the surface's velocity].
  - l.1625–1626: "The instantaneous rate of production within the surface is not altered by the continuous motion of
    the surface."
  - Running head, l.1505/1623: "THE GENERAL EQUATIONS OF MOTION OF ANY ENTITY." (pp. 11, 13).
  - Lineage (gravity in a granular ether), l.1146–1147 and l.1155–1156: "Thus the diameter of a grain is the parameter
    of the effort"; "The efforts to cause the negative inequalities to approach correspond, exactly, to gravitation,
    if matter represents negative mass."; l.1699–1700: "recognised in mechanical philsophy as the medium of space or
    the ether."

### [S2] H. Irschik & H. J. Holl, "Mechanics of variable-mass systems—Part 1: Balance of mass and linear momentum", *Appl. Mech. Rev.* 57(2), 145–160 (2004), DOI 10.1115/1.1687409
- **Rows:** G1, G2, G3.
- **Query:** Q16 (it did not surface the paper); found by a Crossref bibliographic query. The full text was located
  in round 1 from the Codex review leg's evidence file `ih_centres_excerpt.txt` (finding R1).
- **Opened:** `https://api.crossref.org/works/10.1115/1.1687409` (metadata and abstract); round 1, the full text at
  `https://www.mobt3ath.com/uplode/book/book-98207.pdf` (a third-party copy of the ASME PDF; its pages carry
  "Downloaded from http://asmedigitalcollection.asme.org/appliedmechanicsreviews/article-pdf/57/2/B12/5440527/145_1.pdf
  by guest on 26 August 2020").
- **Retrieval (round 1):** `curl -sL -A 'Mozilla/5.0' <url> -o ih.pdf; pdftotext ih.pdf ih.txt; sed -n 1353,1372p ih.txt;
  sed -n 1040,1060p ih.txt; sed -n 1080,1120p ih.txt; sed -n 1755,1880p ih.txt` (in `dl_irschik_holl/`).
- **Quotes (abstract):** "mass is not conserved when sources of mass are present or when the equations of balance
  are applied to a non-material volume."; "We first summarize the general balance law and review the Reynolds
  transport theorem for a non-material volume."; "the equations of balance of mass and linear momentum in the
  presence of sources of mass in the interior of a material volume"; "when mass is flowing through the surface of a
  non-material volume in the absence of sources of mass in the interior"; "Our review remains within the framework
  of the classical, non-relativistic continuum mechanics of solids and fluids."
- **Not in the abstract:** any centre-of-mass relation. Round 0's G3 therefore said the flat form was not located.
  ⇒ **Corrected in round 1 (R1):** the full text has it.
- **Quotes (full text, round 1; equation layout lost in extraction):**
  - l.1353, the heading: "4 VARIABLE MASS CONTAINED IN A NONMATERIAL VOLUME".
  - l.1366–1370: "Recall that the control surface s(t) enclosing the volume v (t) moves at a velocity u, which
    generally differs from the velocity of the material particles instantaneously located on this surface. In the
    present section, we do not take into account sources of mass within v (t), since such sources have been
    discussed in Section 3 of the present review."
  - l.623–624, the heading of §3: "3 SOURCES OF MASS IN THE INTERIOR OF A MATERIAL VOLUME"; §3 defines the centre at
    eq. (3.13) (l.1040–1048) and reaches (3.18) at l.1103–1116.
  - l.1760–1763 (p. 157): "In an analogy to Eq. 共3.13兲, we now introduce the center of mass of the non-material volume
    under consideration. The instantaneous position of the center of mass is given by a vector c u defined through the
    relation" [eq. (4.13), `c_u M_u = ∫_{v(t)} ρ p dv`, l.1764–1773].
  - l.1775: "Differentiating with respect to time gives"; eq. (4.15) carries the surface term
    `∫_S da·(u − ṗ) ρ p` (l.1777–1793, layout lost).
  - l.1799–1800: "But, for a material volume in the absence of sources of mass in the interior, we may write" [eq.
    (4.16): the material-volume integral of `ρ ṗ` is the momentum `P_u`].
  - l.1832–1837: "Putting Eqs. 共4.15兲–共4.17兲 into Eq. 共4.14兲, we arrive at the following relation for the linear
    momentum contained in the non-material volume: P u ⫽M u ċ u ⫹Ṁ u 共 c u ⫺ĉ u 兲 (4.18)".
  - l.1851: "The equation of balance of linear momentum for the nonmaterial volume, Eq. 共4.7兲, then can be
    rewritten as" [eq. (4.20), `M_u c̈_u = F + Ṁ_u (w_u − 2ċ_u + dĉ_u/dt) − M̈_u (c_u − ĉ_u)`, l.1855–1863].
  - l.1865–1869: "In the present context of a non-material volume, we thus have arrived at a relation for the motion
    of the center of mass in complete analogy to the result of Federhofer 关65兴 for a material volume with sources of
    mass in the interior, Eq. 共3.18兲."
  - l.1870–1872: "For the case of a rocket motion with internal flow, Eq. 共4.20兲 was derived in the expositions of
    Rankin 关74兴, and independently by Leitmann 关83兴." (Rankin, Leitmann and Federhofer were not opened.)
  - The doc quotes these with the extraction's bracket glyphs dropped: "we now introduce the center of mass of the
    non-material volume under consideration"; "in complete analogy to the result of Federhofer"; "for a material
    volume with sources of mass in the interior".

### [S3] H. Lamb, *Hydrodynamics*, 3rd ed. (Cambridge, 1906)
- **Rows:** G2 (Art. 346), G9 (Art. 144).
- **Query:** Q15, then archive.org advancedsearch → identifier `hydrodynamics00lambuoft` ("Hydrodynamics | 1906").
- **Opened:** `https://archive.org/download/hydrodynamics00lambuoft/hydrodynamics00lambuoft_djvu.txt` (OCR).
- **Retrieval:** `curl -sL <url> -o l.txt; sed -n 30,45p l.txt; sed -n 12657,12750p l.txt; sed -n 39499,39565p l.txt; sed -n 620,624p l.txt`
- **Edition:** title page: "THIRD EDITION … CAMBRIDGE : AT THE UNIVERSITY PRESS 1906".
- **G9, Art. 144** (begins l.12657, "144. We may also notice one or two problems of Kineto-Statics"; the page head at
  l.12741 is "192 Motion of Solids through a Liquid"):
  - l.12680–12684: "Another interesting case is that of a number of open tubes, so narrow as not sensibly to impede
    the motion of the fluid outside them. If streams be established through the[se tubes] … then as regards the
    external space the extremities will act as sources and sinks."
  - l.12749–12750: "Hence, two sources of like sign attract, and two of unlike sign repel, with forces varying
    inversely as the square of the distance*."
  - l.12743–12746: "The force apparently exerted by mx [m₁] on m2, tending to increase r12, is therefore," followed
    by a formula the OCR garbles (`p'    d     mlm2_      p     mxm2` / `4tt  drn'    rn   ~      An'   rn2`). My
    reading is (ρ/4π) d/dr₁₂ (m₁m₂/r₁₂) = −(ρ/4π) m₁m₂/r₁₂². It comes from the kinetic energy of the source system
    (eq. (10) above it).
- **G2, Art. 346, p. 598** (Art. 346 begins l.39499, "346. Reynolds, in a remarkable paper§"; page head l.39546
  "598 Viscosity"): l.39559–39563: "They express the rate of variation of the momentum contained in a fixed
  rectangular space δxδyδz, as a consequence partly of the forces acting on the substance which at the moment
  occupies this space, and partly of the flux of matter across the boundary, carrying its momentum with it."
- **Bjerknes omission,** l.622–624 (Additions and Corrections): "Page 149, Art. 116. A reference to the researches
  of C. A. Bjerknes on pulsating spheres has unfortunately been omitted. A full account of these is given by V.
  Bjerknes, Vorlesungen über hydrodynamische Fernkräfte, Leipzig, 1900-1902."
- **Search note:** a grep of this edition for a statement "force on a body = momentum flux through any enclosing
  surface" (`grep -n -i -E 'flux of momentum|momentum.*(surface|boundary)' l.txt`) returned only l.39562 above.

### [S4] É. Racine & É. É. Flanagan, "Post-1-Newtonian equations of motion for systems of arbitrarily structured bodies", *Phys. Rev. D* 71, 044010 (2005), DOI 10.1103/PhysRevD.71.044010; arXiv:gr-qc/0404101
- **Rows:** G3, G5.
- **Query:** Q22.
- **Opened:** `https://arxiv.org/pdf/gr-qc/0404101v1`. Crossref confirms the journal data.
- **Retrieval:** `curl -sL https://arxiv.org/pdf/gr-qc/0404101v1 -o rf.pdf; pdftotext rf.pdf rf.txt; sed -n 8,30p rf.txt`
- **Quotes (abstract):** "We give a surface integral derivation of post-1-Newtonian translational equations of
  motion for a system of arbitrarily structured bodies, including the coupling to all the bodies' mass and current
  multipole moments."; "The derivation requires only that the post-1-Newtonian vacuum field equations are satisfied
  in weak-field regions between the bodies; the bodies' internal gravity can be arbitrarily strong."; "The
  definition of each body's current and mass multipole moments and center-of-mass worldline in terms of the
  behavior of the metric in a weak-field region surrounding the body."; "The derivation of a formula, for any given
  body, of the second time derivative of its mass dipole moment in terms of its other multipole and tidal moments
  and their time derivatives."

### [S5] Mathisson–Papapetrou–Dixon: M. Mathisson (1937); A. Papapetrou (1951); W. G. Dixon (1970)
- **Rows:** G4.
- **Queries:** Q6, Q7, Q9.
- **Opened:**
  - Wikipedia *Mathisson–Papapetrou–Dixon equations*, revision 1313946975 (2025-09-28), via the API.
  - Crossref for `10.1098/rspa.1951.0200` (Papapetrou), `10.1098/rspa.1970.0020` (Dixon) and
    `10.1007/s10714-010-0939-y` (Mathisson, English republication).
- **Citations as recorded:**
  - Mathisson, "Neue Mechanik materieller Systeme", *Acta Phys. Polon.* 6, 163–209 (1937), per Wikipedia. A
    WebSearch summary (Q7) gave 163–200, so the page range is uncertain. English republication: "Republication of:
    New mechanics of material systems", *Gen. Relativ. Gravit.* 42(4), 1011–1048 (2010), per Crossref.
  - Papapetrou, "Spinning test-particles in general relativity. I", *Proc. R. Soc. A* 209(1097), 248–258 (1951).
  - Dixon, "Dynamics of extended bodies in general relativity. I. Momentum and angular momentum", *Proc. R. Soc. A*
    314(1519), 499–527 (1970).
- **Quotes:**
  - Wikipedia: "the vector V^μ is the four-velocity of some reference point X^μ in the body, and the skew-symmetric
    tensor S^{μν} is the angular momentum … of the body about this point."; "The equations must therefore be
    supplemented by three additional constraints which serve to determine which point in the body has velocity V^μ."
    [math markup simplified].
  - Papapetrou, abstract (Crossref): "A method for the derivation of the equations of motion of test particles in a
    given gravitational field is developed."
  - Dixon, abstract (Crossref): "Definitions are proposed for the total momentum vector pα and spin tensor Sαβ of an
    extended body in arbitrary gravitational and electromagnetic fields."; "used to support the definition pβSαβ = 0
    for the centre of mass."
- **Not opened:** Dixon II (1970) and III (1974); Mathisson's original German text; Papapetrou's and Dixon's full
  texts (abstracts only).
- **Round 1 (R4):** the two Wikipedia sentences above are Wikipedia's summary, not the papers' text. The section now
  marks them secondary only. Its only primary quote in G4 is Dixon's abstract ("the definition pβSαβ = 0 for the
  centre of mass").

### [S6] A. Einstein, L. Infeld & B. Hoffmann, "The Gravitational Equations and the Problem of Motion", *Ann. Math.* 39(1), 65–100 (1938), DOI 10.2307/1968714
- **Rows:** G5, G6 (the Helmholtz remark).
- **Queries:** Q1, Q11.
- **Opened:**
  - `https://mathshistory.st-andrews.ac.uk/Extras/Hoffmann_papers/`, which reproduces the paper's introduction.
  - `https://edition-open-sources.org/sources/10/17/index.html` for the citation.
  - Crossref `10.2307/1968714`.
- **Retrieval:** `curl -sL -A 'Mozilla/5.0' 'https://mathshistory.st-andrews.ac.uk/Extras/Hoffmann_papers/' -o h.html`,
  then strip the tags and print around "surface integral".
- **Quotes:**
  - "the well-known theory of Helmholtz on the motion of vortices in a non-viscous fluid gives an instance where the
    motion of line singularities is actually determined by partial differential equations alone"
  - "We shall show in this paper that the gravitational equations for empty space are in fact sufficient to determine
    the motion of matter represented as point singularities of the field."
  - "(2) It is shown that for two-dimensional spatial surfaces containing singularities certain surface integral
    conditions are valid which determine the motion."
  - Edition Open Sources: "Albert Einstein, Leopold Infeld and Banesh Hoffmann (1938). The Gravitational Equations
    and the Problem of Motion. The Annals of Mathematics, 39: 65–100."

### [S7] L. Blanchet, "Post-Newtonian theory for gravitational waves", *Living Rev. Relativ.*; arXiv:1310.1528v6 (2024)
- **Rows:** G5 (the EIH method and Damour's effacing principle).
- **Queries:** Q10, Q11.
- **Opened:** `https://arxiv.org/pdf/1310.1528` (v6).
- **Retrieval:** `curl -sL https://arxiv.org/pdf/1310.1528v6 -o b.pdf; pdftotext b.pdf b.txt; grep -n -E 'surface-integral|effacing principle|^Damour T \((1983a|1987b)\)' b.txt`
- **Quotes:**
  - "Einstein et al. (1938) obtained the 1PN corrections by means of their famous “surface-integral” method, in which
    the equations of motion are deduced from the vacuum field equations, and are therefore applicable to any compact
    objects (be they neutron stars, black holes, or, perhaps, naked singularities)."
  - "This is via the effacing principle for the internal structure of the compact bodies."
  - "Therefore, in that sense the internal structure of the two compact bodies is “effaced” and their dynamics and
    radiation depend only, in first approximation, on the masses and possibly the spins (Damour 1983a)."
  - References: "Damour T (1983a) Gravitational radiation and the motion of compact bodies. In: Deruelle N, Piran T
    (eds) Rayonnement Gravitationnel / Gravitational Radiation. North-Holland, Amsterdam, pp 59–144"; "Damour T
    (1987b) The problem of motion in Newtonian and Einsteinian gravity. In: Hawking SW, Israel W (eds) Three Hundred
    Years of Gravitation. Cambridge University Press, Cambridge; New York, pp 128–198".
- **Not opened:** Damour's chapters themselves.

### [S8] K. S. Thorne & J. B. Hartle, "Laws of motion and precession for black holes and other bodies", *Phys. Rev. D* 31(8), 1815–1837 (1985), DOI 10.1103/PhysRevD.31.1815
- **Rows:** G5.
- **Query:** Q8.
- **Opened:** a PDF scan at
  `https://blackholes.tecnico.ulisboa.pt/gritting/pdf/gravitational_waves/Thorne-Hartle_Laws-of-motion-and-precession-for-black-holes-and-other-bodies.pdf`
  (OCR). Crossref confirms the journal data.
- **Retrieval:** `curl -sL <url> -o th.pdf; pdftotext th.pdf th.txt; sed -n 14,24p th.txt; sed -n 455,470p th.txt; sed -n 1060,1080p th.txt`
- **Quotes:**
  - Abstract: "Laws of motion and precession are derived for a Kerr black hole or any other body which is far from
    all other sources of gravity {"isolated body" ) and has multipole moments that change slowly with time."
  - "The rates of change of momentum and angular momenturn of a body could be evaluated by integrating expressions
    for the momentum flux and angular momentum flux over a closed two surface in the buffer region."
  - "Now, only the values of the body's multipole moments are important for the gravitational field, not their
    source. In particular these values carry no information as to whether the source is compact or diffuse."
  - On the rates of M, P and angular momentum: "These can be expressed as surface integrals of a symmetric
    energy-momentum pseudotensor in the standard way"; "the surface integrals again are over the closed two-surface
    of constant r in the buffer region."

### [S9] D. J. Thouless, P. Ao & Q. Niu, "Transverse force on a quantized vortex in a superfluid", *Phys. Rev. Lett.* 76(20), 3758–3761 (1996), DOI 10.1103/PhysRevLett.76.3758; arXiv:cond-mat/9603196
- **Rows:** G6.
- **Query:** Q20.
- **Opened:** `https://arxiv.org/pdf/cond-mat/9603196`. Crossref confirms the journal data.
- **Retrieval:** `curl -sL https://arxiv.org/pdf/cond-mat/9603196v1 -o tan.pdf; pdftotext tan.pdf tan.txt; sed -n 11,18p tan.txt; sed -n 355,390p tan.txt`
- **Quotes:**
  - Abstract, l.13–14 (round 1, R2): "We have derived an exact expression for the total nondissipative transverse
    force acting on a quantized vortex moving in a uniform background."
  - "The derivation is valid for neutral boson or fermion superfluids, provided the order parameter is a complex
    scalar quantity. This force is determined by the oneparticle density matrix far away from the vortex core, and is
    found to be the Magnus force proportional to the superfluid density. We conclude that contributions of the
    localized core states do not change this force."
  - "there are no contributions to the integrals from the neighborhood of the vortex core, since we have chosen the
    boundaries of the regions of integration to be well away from the core"
  - "If fluid is streaming past a fixed vortex with velocity v, one can calculate the momentum balance on a very
    large cylinder centered on the vortex. The net pressure force per unit length on the interior of this cylinder is
    ρn K × v/2, the momentum flow out from the cylinder is −ρn K × v/2, and these must be balanced by a force −ρn K ×
    v acting on the vortex core."
  - "for any infinite, homogeneous, neutral superfluid with a scalar order parameter the transverse force has the
    form ρs K × vV"
  - l.368–369 (momentum density, for the NOT FOUND list): "For a neutral superfluid the integrand in eq. (10) is just
    the momentum density, equal, by definition, to ρs us + ρn un , where us and un are the local values of the
    superfluid and normal velocities."

### [S10] M. Lagally, "Berechnung der Kräfte und Momente, die strömende Flüssigkeiten auf ihre Begrenzung ausüben", *ZAMM* 2(6), 409–422 (1922), DOI 10.1002/zamm.19220020601
- **Rows:** G7.
- **Queries:** Q2, Q5.
- **Opened:**
  - Zenodo record `https://zenodo.org/records/1548357`; PDF from
    `https://zenodo.org/api/records/1548357/files/article.pdf/content` (a scan with an OCR layer).
  - Crossref `10.1002/zamm.19220020601`: "ZAMM … | 2 6 409-422 | 1922".
- **Retrieval:** `curl -sL <pdf url> -o lagally.pdf; pdftotext lagally.pdf l.txt; sed -n 112,160p l.txt`
- **Quotes (OCR):**
  - l.118–123: "Wenn eine stationLre StrSmnng das Innere einer geschlossenen b’l‘llche erfullt, so h a n g t d e r v
    o n d e r F l u s s i g k e i t auf d i e F l l c h e a u s g e u b t e D r n c k u n d s e i n M o m e n t
    ansschlie5lich von d e r Verteilung der Quellen und Wirbel i n d e r Flussigk e i t n n d v o n d e r a n i h
    r e m O r t h e r r s c h e n d e n G e s c h w i n d i g k e i t a b ; von der Gestalt der Begrenzung nicht
    unmittelbar, sondern nur inaofern, als die Verteilung der Geschwindigkeit durch die Begrenzung bedingt ist."
  - My reading of that sentence: if a steady flow fills the interior of a closed surface, the pressure the fluid
    exerts on the surface, and its moment, depend only on the distribution of sources and vortices in the fluid and
    on the velocity at their location. They depend on the boundary's shape only through the velocity distribution
    the boundary sets.
  - l.123–125: "Aehnliches gilt, wenn die Flussigkeit den AuEenranm erfiillt , jedoch mud in diesem Fall auch die
    analytische Fortsetzung der Stromung in den Innenraum bekannt sein." [the same holds when the fluid fills the
    exterior, but then the analytic continuation of the flow into the interior must be known].
  - l.130: "Die Stromung soll stationiir sein" [the flow shall be steady].

### [S11] L. Landweber & C.-S. Yih, "Forces, moments, and added masses for Rankine bodies", *J. Fluid Mech.* 1(3), 319–336 (1956), DOI 10.1017/S0022112056000184
- **Rows:** G7 (and G8, for the Kirchhoff–Taylor–Lagally triad).
- **Query:** Q18, then a Crossref bibliographic query. The Cambridge Core pages returned an error page to `curl`.
- **Opened:** Crossref (abstract).
- **Quotes (abstract):** "The dynamical theory of the motion of a body through an inviscid and incompressible fluid
  has yielded three relations: a first, due to Kirchhoff, which expresses the force and moment acting on the body in
  terms of added masses; a second, initiated by Taylor, which expresses added masses in terms of singularities within
  the bòdy; and a third, initiated by Lagally, which expresses the forces and moments in terms of these
  singularities."; "generalizations of the Taylor and Lagally theorems to include unsteady flow and arbitrary
  translational and rotational motion of the body"

### [S12] L. Landweber & T. Miloh, "Unsteady Lagally theorem for multipoles and deformable bodies", *J. Fluid Mech.* 96(1), 33–46 (1980), DOI 10.1017/S0022112080002005
- **Rows:** G7.
- **Query:** Q2.
- **Opened:** Crossref (abstract). The Cambridge Core page was reached through WebFetch, which returned a summary only.
- **Quotes (abstract):** "The Lagally theorem for unsteady flow expresses the forces and moments acting on a rigid
  body moving in an inviscid and incompressible fluid in terms of the singularities of the analytically continued
  flow within the body."; "Previous generalizations of the Lagally theorem, originally given by Lagally (1922) for
  steady flows, are due to Cummins (1957) and Landweber & Yih (1956)"; "to generalize the Lagally theorem so as to
  include the effects of multipoles in unsteady flow, and deformability of the body"

### [S13] X.-S. Wang, "Derivation of Newton's law of gravitation based on a fluidic continuum model of vacuum and a sink flow model of particles", arXiv:physics/0506062v5 (2005, revised 2017); *Progress in Physics* 4, 25–30 (2008)
- **Rows:** G7, G9, G11, G12, G13.
- **Queries:** Q2, Q21.
- **Opened:**
  - `https://arxiv.org/abs/physics/0506062`, whose listing gives "Journal reference: Progress in Physics, 2008,
    Vol. 4, 25-30".
  - `https://arxiv.org/pdf/physics/0506062` (served v5).
- **Retrieval:** `curl -sL https://arxiv.org/pdf/physics/0506062v5 -o w.pdf; pdftotext w.pdf w.txt`, then
  `grep -n -E 'Theorem 2 Suppose|FQ = −ρ Q|F21 = |universal damping|Assumption 6|increasing linearly|Lagally \[52\]|fixed in space|arbitrarily small|does not include the singularity|\[52\] M. Lagally' w.txt`
  and `sed -n 147,232p w.txt; sed -n 505,700p w.txt`.
- **Quotes:**
  - Abstract: "Thus, the inverse-square law of gravitation is derived by methods of hydrodynamics based on a sink flow
    model of particles."; "Newton’s law of gravitation is derived if we introduce an assumption that G and the masses
    of particles are changing so slowly that they can be treated as constants."
  - Theorem 2 conditions: "Suppose (1) there exists an ideal fluid (2) the ideal fluid is irrotational and
    barotropic, (3) the density ρ is homogeneous, that is ∂ρ/∂x=∂ρ/∂y=∂ρ/∂z=∂ρ/∂t = 0 , (4) there are no external
    body forces exerted on the fluid, (5)the fluid is unbounded and the velocity of the fluid at the infinity is
    approaching to zero. Suppose a source or sink is stationary and is immersed in the ideal fluid."
  - Theorem 2 result: "FQ = − ρQu0"; "u0 is the velocity of the fluid at the location of the source induced by all
    means other than the source itself."
  - The control volume: "We surround the source by two arbitrary spherical surfaces Sε and S centered at the origin
    of the coordinates."; "Note that the mass system τ does not include the singularity at the origin. According to
    Reynolds’ transport theorem"; "where V is the volume fixed in space which coincide with the mass system τ (t) at
    time t"; "Since the radius r can be arbitrarily small, the velocity u0 can be treated as a constant".
  - "Remark. Lagally [52], Landweber and Yih [53, 54], Faber [55] and Currie [56] obtained the same result of
    Theorem 2 for the special case where the velocity field is steady."; "[52] M. Lagally, Z. Angew. Math. Mech. 2,
    409 (1922)."
  - Theorem 3 (`wang_raw.txt` l.535–539): "Suppose the presuppositions (1), (2), (3), (4) and (5) in Theorem 2 are
    valid and a source or a sink is moving in the fluid with a velocity vs , then there is a force FQ = −ρ Q (uf −
    vs )". Its proof is "The velocity of the fluid relative to the source at the location of the source is uf − vs .
    Let us select the coordinates that is attached to the source and set the origin of the coordinates at the
    location of the source. Then (27) can be arrived following the same procedures in the proof of Theorem 2."
  - Corollary 4: "F21 = −ρQ2 Q1/(4πr²) r̂21 + ρQ2 v2" [fraction layout flattened].
  - Round 1 (R6), Corollary 4's conditions and symbols, from `pdftotext` without `-layout` (`wang_raw.txt`
    l.580–583 and l.615–618): "Suppose the presuppositions (1), (2), (3), (4) and (5) in Theorem 2 are valid and a
    source or a sink with strength Q2 is exposed to the velocity field of another source or sink with strength Q1";
    "r is the distance between the two singularities, v2 is the velocity of the source with strength Q2." So the
    inverse-square term alone is the force on a sink at rest (`v2 = 0`); the doc writes the term as
    `−ρQ₂Q₁/(4πr²) r̂₂₁` and the velocity term as `ρQ₂v₂`.
  - Mass from inflow: "Assumption 6 All the microscopic particles were made up of a kind of elementary sinks of Ω(0)
    substratum."; eq. (30) "dm/dt = −ρQ" [layout flattened]; "From (31), we also see that the mass of a monad or a
    particle is increasing linearly."
  - Drag: "Formula (32) shows that there exists a universal damping force".
  - Momentum density: "K = ∫∫∫ ρ udV" [eq. (10), layout flattened].
- **Assessment search:** Q21 returned the same author's related arXiv papers and a forum thread, but no independent
  peer assessment. G11 records this as "no independent assessment found".

### [S14] G. I. Taylor, "The forces on a body placed in a curved or converging stream of fluid", *Proc. R. Soc. A* 120(785), 260–283 (1928), DOI 10.1098/rspa.1928.0148
- **Rows:** G8.
- **Query:** Q3. Crossref for the DOI, and a bibliographic query for the companion paper.
- **Opened:** Crossref (abstract). Also the Royal Society referee-report record `https://makingscience.royalsociety.org/items/rr_38_147`,
  reached through the search summary only.
- **Quotes (abstract):** "The force which a body experiences when placed in a converging stream of fluid has a
  certain practical interest in aeronautics"; "if dp/dx is the pressure gradient, i. e., the gradient of static
  pressure in the channel, and V the volume of the body, the “horizontal buoyancy” is — V dp/dx."
- **Companion paper (Crossref):** G. I. Taylor, "The energy of a body moving in an infinite fluid, with an
  application to airships", *Proc. R. Soc. A* 120, 13–21 (1928), DOI 10.1098/rspa.1928.0131. This is the "Taylor"
  of S11's added-mass theorem. It was not opened beyond the metadata.

### [S15] T. R. Auton, J. C. R. Hunt & M. Prud'homme, "The force exerted on a body in inviscid unsteady non-uniform rotational flow", *J. Fluid Mech.* 197, 241–257 (1988), DOI 10.1017/S0022112088003246
- **Rows:** G8.
- **Query:** Crossref bibliographic query (no WebSearch).
- **Opened:** Crossref (abstract). Round 1 (R7): re-fetched with
  `curl -s 'https://api.crossref.org/works/10.1017/S0022112088003246' -o auton.json` (in `dl_crossref/`) and printed
  with the Crossref one-liner above.
- **Round 1 quote (R7), the sentence's second half:** "It is assumed that the radius is small compared with the scale
  over which the strain rate changes, though for the sphere it is also assumed that the changes in the ambient
  velocity field over the scale of the sphere are small compared with the velocity of the body relative to the
  flow."
- **Quotes (abstract):** "A general expression is derived for the fluid force on a body of simple shape moving with a
  velocity v through inviscid fluid in which there is an unsteady non-uniform rotational velocity field u0(x, t) in
  two or three dimensions." [JATS markup removed]; "It is assumed that the radius is small compared with the scale
  over which the strain rate changes"; "can be simply calculated from the added-mass coefficient".

### [S16] V. Bjerknes, *Fields of Force* (Columbia University Press / Macmillan, 1906)
- **Rows:** G9.
- **Query:** Q4, then archive.org advancedsearch → `fieldsofforce00bjeruoft` ("Fields of force | 1906").
- **Opened:** `https://archive.org/download/fieldsofforce00bjeruoft/fieldsofforce00bjeruoft_djvu.txt` (OCR).
- **Retrieval:** `curl -sL <url> -o f.txt; sed -n 150,156p f.txt; sed -n 1230,1345p f.txt; sed -n 2219,2226p f.txt`
- **Quotes:**
  - l.152–153: "The results which I shall bring before you were discovered originally by my late father. Professor
    C. A. Bjerknes, of Christiania."
  - l.1235–1237: "If its volume be constant, so that the displaced mass of water is constant, the force will be
    strictly periodic, with a mean value zero in the period."
  - l.1257–1260: "A pulsating body in a synchronously oscillating current is subject to the action of a resultant
    force, the direction of which is that of the acceleration in the current at the time when the pulsating body has
    its maximum volume."
  - l.1297–1301 (p. 32): "Between bodies puhating in the same phase there is nn apparent attraction ; between bodies
    pulsating in the opposite phase there is an apparent repulsion, the force being proportional to the product of
    the two intensities of pulsation, and proportional to the inverse square of t/ie distance." [puhating =
    pulsating; nn = an; t/ie = the]
  - l.1309–1314: "the law enunciated above has exactly the form of Coulomb's law for the action between two
    electrically charged particles, with one striking difference ; the direction of the force in the iiydrodynamic
    field is opposite to that of the corresponding force in the electric or magnetic field."
  - l.1341–1343 (p. 33): "Palliating bodies act upon each other as if they were electrically charged particles or
    magnetic poles, but vnth the difference that charges or poles of the same sign attract, and charges or poles of
    opposite sign repel each other." [Palliating = Pulsating; vnth = with]
  - l.2222–2225: "While searching for phenomena of hydrodynamics which should have the appearance of actions at a
    distance, he solved the problem of the simultaneous motion of any number of spherical bodies in a liquid."
- **Search note:** gravity is not the analogy drawn in these pages. `grep -n -i gravit f.txt` hits concern static
  buoyancy only (l.452, 455, 639, 1000).

### [S17] Wikipedia, *Mechanical explanations of gravitation*, revision 1360707870 (2026-06-23)
- **Rows:** G9 (pulsation as gravity), G10 (Newton 1675), G12 and G13 (Yarkovsky; criticism). This is a secondary
  source, and the section labels everything it supports as "secondary only".
- **Query:** follow-up to Q4/Q12.
- **Opened:** `https://en.wikipedia.org/w/api.php?action=query&prop=revisions&rvprop=content|ids|timestamp&rvslots=main&format=json&titles=Mechanical_explanations_of_gravitation`
  (wikitext, sections "Streams" and "Pulsation").
- **Quotes:**
  - "In a 1675 letter to [[Henry Oldenburg]], and later to [[Robert Boyle]], Newton wrote the following: [Gravity is
    the result of] “a condensation causing a flow of ether with a corresponding thinning of the ether density
    associated with the increased velocity of flow.”"
  - "[[Bernhard Riemann]] assumed in 1853 that the gravitational aether is an [[Incompressible flow|incompressible
    fluid]] and normal matter represents sinks in this aether."; "Riemann speculated that the absorbed aether is
    transferred into another world or dimension."
  - "Another attempt to solve the energy problem was made by [[Ivan Osipovich Yarkovsky]] in 1888. Based on his
    aether stream model, which was similar to that of Riemann, he argued that the absorbed aether might be converted
    into new matter, leading to a mass increase of the celestial bodies."
  - "''Criticism'': As in the case of Le Sage's theory, the disappearance of energy without explanation violates the
    [[energy conservation law]]. Also some drag must arise, and no process which leads to a creation of matter is
    known."
  - "[[William Thomson, 1st Baron Kelvin|Lord Kelvin]] (1871) and [[Carl Anton Bjerknes]] (1871) assumed that all
    bodies pulsate in the aether."; "To explain universal gravitation, one is forced to assume that all pulsations in
    the universe are in phase—which appears very implausible." (cited there to "Zenneck (1903)", not opened).

### [S18] B. Riemann, "Gravitation und Licht" and "Neue mathematische Principien der Naturphilosophie", in *Bernhard Riemann's Gesammelte mathematische Werke und wissenschaftlicher Nachlass*, ed. H. Weber with R. Dedekind (Teubner, Leipzig, 1876), pp. 496–506
- **Rows:** G7 (Riemann's question), G10, G12.
- **Queries:** Q12, Q13, Q14. Galindo & Cervantes-Cota [S19] give the archive.org link and the page numbers.
- **Opened:** `https://archive.org/stream/bernardrgesamm00riemrich/bernardrgesamm00riemrich_djvu.txt` (OCR).
- **Retrieval:** `curl -sL <url> -o rw.txt`, then `sed -n 2974,2985p rw.txt` (title page), `sed -n 37808,37965p rw.txt`
  (pp. 496–498), `sed -n 38312,38440p rw.txt` (pp. 502–504) and `sed -n 38500,38616p rw.txt` (pp. 505–506).
- **Edition (title page):** "HERAUSGEGEBEN UNTER MITWIRKUNG VON R. DEDEKIND [VON] H. WEBER. … 1876."
- **Round 1 (R3): the two fragments are separate.** Fragment 2, "Gravitation und Licht", runs from its heading at
  l.37838 (p. 496) to l.38352 (p. 502); fragment 3, "Neue mathematische Principien", runs from its heading at l.38353
  (p. 502) to the end of the text at l.38616 (p. 506), before "Bernhard Riemann's Lebenslauf". The page heads are
  l.37808 (496), 38312 (502), 38433 (504) and 38556 (506). Inside fragment 3, `grep -n '^[0-9]\+\. \|^[IVX]\+\. '`
  finds only its own heading and page heads, so pp. 502–506 are one fragment. Only fragment 3 carries the date
  footnote. Retrieval: `sed -n 37838,38352p rw.txt > frag2.txt; sed -n 38353,38620p rw.txt > frag3.txt`.
- **Search statement (R3):** `grep -c -i` over `frag2.txt` gives `verschwind` 0, `Atom` 0, `einström` 1 (the
  conservation passage below); over `frag3.txt`, `verschwind` 4, `Atom` 7. No inflow into atoms is stated in
  fragment 2.
- **"2. Gravitation und Licht"** (l.37838; p. 496; no date):
  - l.37885–37892: "Die nach Grösse und Richtung bestimmte Ursache (beschleunigende Schwerkraft), welche nach 3. in
    jedem Punkte des Raumes stattfindet, suche ich in der Bewegungsform eines durch den ganzen unendlichen Raum stetig
    verbreiteten Stoffes, und zwar nehme ich an, dass die Richtung der Bewegung der Richtung der aus ihr zu
    erklärenden Kraft ‚gleich, und ihre Geschwindigkeit der Grösse der Kraft proportional sei. Dieser Stoff kann also
    vorgestellt werden als ein physischer Raum, dessen Punkte sich in dem geometrischen bewegen."
  - l.37877–37881 (p. 497, Newton's item 4, the force law the motion must reproduce): "eine dem Quadrat der Entfernung
    von diesem ponderablen Punkte umgekehrt und seiner Schwerkraft direet propor- tionale beschleunigende Kraft".
  - l.37902–37905: "Ich nehme nun an, dass die wirkliche Bewegung des Stoffes im leeren Raum zusammengesetzt ist aus
    der Bewegung, welche zur Erklärung der Gravitation, und aus der, welche zur Erklärung des Lichtes angenommen werden
    muss."
  - l.37960–37963: "behält der Stoff, wenn er in Einem Zeitpunkte überall gleich dicht ist, stets allenthalben
    dieselbe’Dichtigkeit".
  - l.37956–37959 (round 1): "der Bewegung, welche die Gravitationserscheinungen verursacht, durch ,, %,, U, der
    Bewegung, welche die Lichterscheinungen verursacht, durch w,, 20,, %0,, der wirklichen Bewegung durch ®,, %,, %,,
    so dass v=u-+ w." [OCR; the reading is v = u + w].
  - l.38287–38296 (p. 501, round 1): "Die Gleichung (I) beweist unsere frühere Behauptung, dass bei der
    Stoffbewegung die Dichtigkeit ungeändert bleibe; denn … drückt die in das Raum- element dx, dx,dx, im
    Zeitelement dt einströmende Stoffmenge aus, und die in ihm enthaltene Stoffimenge bleibt daher constant."
- **"3. Neue mathematische Principien der Naturphilosophie"** (pp. 502–506; dated):
  - footnote l.38372: "*) Gefunden am 1. März 1853."
  - l.38400–38403 (p. 503): "Von dieser Thatsache geleitet, mache ich die Hypothese, dass der Weltraum mit einem Stoff
    erfüllt ist, welcher fortwährend in die pon- derablen Atome strömt und dort aus der Erscheinungswelt (Körper-
    welt) verschwindet."
  - l.38405–38406: "dass in allen ponderablen Atomen beständig Stoff aus der Körperwelt in die Geisteswelt eintritt."
  - l.38420–38422 with l.38436: "Nimmt man nun an, dass der raumerfüllende Stoff eine incom- pressible homogene
    Flüssigkeit ohne Trägheit sei, und dass in jedes ponderable Atom in gleichen Zeiten stets gleiche, seiner Masse
    pro- … portionale Mengen einströmen, so wird offenbar der Druck, den das ponderable Atom erfährt, (der
    Geschwindigkeit der Stoffbewegung an dem Orte des Atoms proportional sein?)"
    [The word "pro- portionale" is split by the page break and the footnote; the doc quotes it rejoined as "seiner
    Masse proportionale Mengen".]
  - footnote l.38425–38426: "In jedes ponderable Atom tritt in jedem Augenblick eine bestimmte, der
    Gravitationskraft proportionale Stoffmenge ein und verschwindet dort."
    [This footnote's inflow law (∝ the gravitational force) differs from the text's (∝ the mass). The doc gives both
    (R3).]
  - l.38440–38443 (p. 504): "Es kann also die Wirkung der allgemeinen Gravitation auf ein ponderables Atom durch den
    Druck des raumerfüllenden Stoffes in der unmittelbaren Umgebung desselben ausgedrückt und von demselben abhängig
    gedacht werden."
  - l.38445–38447 (p. 504, round 1): "Aus unserer Hypothese folgt nothwendig, dass der raumerfüllende Stoff die
    Schwingungen fortpflanzen muss, welche wir als Licht und Wärme wahrnehmen."
  - l.38512–38516 (p. 505): "dass der Einfluss einer früheren Form (caeteris paribus) desto geringer wird, je länger
    vor £ [t] sie stattfand"
  - l.38564–38565 (p. 506): "Es ist nun kein Grund vorhanden, apzunelimen [anzunehmen], das die Wirkungen beider
    Ursachen nach demselben Gesetz mit der Zeit sich änderten"
  - l.38582–38586: "Wie müssen nun die Funetionen % und @ beschaffen sein, damit Gra- vitation, Licht und strahlende
    Wärme durch den Raumstoff vermittelt werde ?" [the two function symbols are lost in OCR]
  - l.38607–38616: "1) in den Widerstand, mit welchem ein Stofftheilchen einer Yu [Volum-] änderung, und 2) in den
    Widerstand, mit welchem ein physisches Linienelement einer Längenänderung widerstrebt. Auf dem ersten Theil beruht
    die Gravitation und die eleetrostatische Anziehung und Abstossung, auf dem zweiten die Fortpflanzung des Lichts
    und der Wärme [und] die electrodynamische oder magnetische Anziehung und Abstossung."
- **Search note:** `grep -c -i 'Ablenk' rw.txt` returned `0`, and the same grep over l.37600–38620 also returned `0`.
  I found no light-deflection computation in the fragments.

### [S19] S. Galindo & J. L. Cervantes-Cota, "Clifford's attempt to test his gravitation hypothesis", *Rev. Mex. Fís. E* 64, 162 (2018); arXiv:1807.09230
- **Rows:** G10 (the English rendering, and the page pointers to S18).
- **Queries:** Q12, Q13.
- **Opened:** `https://arxiv.org/abs/1807.09230` and `https://arxiv.org/pdf/1807.09230`.
- **Retrieval:** `curl -sL https://arxiv.org/pdf/1807.09230 -o g.pdf; pdftotext g.pdf g.txt; sed -n 95,107p g.txt; sed -n 399,406p g.txt`
- **Quotes:**
  - "in an unfinished manuscript dated March 1st 1853 … he proposed the very audacious idea that the universe is filled
    with some kind of ether which he called “Stoff”, that ran into atoms to fade away, “I make the hypothesis that
    space is filled with a substance which continually flows into the ponderable atoms where it disappears from the
    world of phenomena “[7]."
  - "“I seek the cause [of gravitation] in the state of motion of the […] continuous substance (Stoff) spread
    throughout the entire infinite space. […] Thus this substance may be thought of as a physical space whose points
    move in the geometrical space.”[8]"
  - References [7] and [8] give Werke vol. II "p. 502" and "p. 496" respectively, with the archive.org link.

### [S20] U. Bottazzini & R. Tazzioli, "Naturphilosophie and its role in Riemann's mathematics", *Revue d'histoire des mathématiques* 1, 3–38 (1995)
- **Rows:** G10 (context only).
- **Queries:** Q12, Q13.
- **Opened:** the publisher's sample PDF `https://smf.emath.fr/sites/default/files/2024-11/smf_rhm_1_3-38__sample.pdf`
  (pp. 3–7 only).
- **Retrieval:** `curl -sL <url> -o bt.pdf; pdftotext -layout bt.pdf bt.txt; sed -n 1,40p bt.txt; sed -n 86,91p bt.txt`
- **Quotes:** "Revue d’histoire des mathématiques, 1 (1995), p. 3–38."; "Riemann’s writings on Naturphilosophie can be
  regarded as the result of his attempt to find a unified, mathematical explanation of various physical phenomena
  such as gravitation, electricity, magnetism and light."; "Instead of this, in his later lectures on gravitation,
  electricity and magnetism Riemann preferred to use the term ether."
- **Not verified:** what Bottazzini & Tazzioli say about the open memory-function question. Only pp. 3–7 were
  available. The "open question" in G10 is quoted from Riemann himself [S18].

### [S21] A. J. S. Hamilton & J. P. Lisle, "The river model of black holes", *Am. J. Phys.* 76(6), 519–532 (2008), DOI 10.1119/1.2830526
- **Rows:** G11 (already cited in the 4D-bulk table; re-checked here).
- **Opened:** Crossref (abstract).
- **Quote:** "In a spherical black hole, the river of space falls into the black hole at the Newtonian escape velocity,
  hitting the speed of light at the horizon."

### [S22] B. Kulangiev, "Emergent Newtonian Gravity from the Collective Acoustic Metric of a Logarithmic Superfluid Vacuum", Qeios preprint v3, DOI 10.32388/INA4A7.3 (2026-07-28)
- **Rows:** G11 (ADJACENT).
- **Query:** Q24.
- **Opened:** `https://www.qeios.com/read/INA4A7.3/pdf`.
- **Retrieval:** `curl -sL <url> -o k.pdf; pdftotext k.pdf k.txt; sed -n 1,30p k.txt; sed -n 2540,2560p k.txt`
- **Quotes:** "28 July 2026, Preprint v3"; "Boris Kulangiev1 1. Independent researcher"; "the Painlevé-Gullstrand
  (PG) flow"; "Each soliton core … acts as a sink for the background condensate."; "so the equality can hold at one
  radius only."; "The continuity constraint is discussed as an open problem."
  - l.2545–2554 (the context of "one radius only"; symbols lost in extraction): "Each soliton core ( at ) acts
    as a sink for the background condensate. The absorption rate per particle, summed over particles, could balance
    the PG inflow."; "A sustained balance across a range of radii would require the absorption rate itself to vary
    with position".

### [S23] Wikipedia, *Johann Georg von Soldner*, revision 1342726980 (2026-03-10)
- **Rows:** none in the table. It is used only for the correction of a lead in the final report.
- **Query:** Q25.
- **Quote:** "In a paper written in 1801 and published in 1804, he calculated the amount of deflection of a light ray by
  a star"; "a combined prediction that was twice as high as the earlier predictions."

---

## Leads from the two earlier notes: what this pass found

| lead (as given) | finding |
|---|---|
| EIH 1938, *Ann. Math.* 39: motion from surface integrals around each body | **Verified** [S6, S7]: pp. 65–100. The motion follows from the vacuum ("empty space") field equations, for point singularities. |
| Mathisson 1937; Papapetrou 1951, *Proc. R. Soc. A* 209: integrate over a worldtube, motion in moments (pole–dipole) | **Citations verified** [S5]. Mathisson's page range is 163–209 (Wikipedia) or 163–200 (search summary). Papapetrou's abstract says "test particles in a given gravitational field". "Worldtube" and "pole–dipole" are **not verified** from an opened text. |
| Dixon (1970–74) made it exact | **Partly verified** [S5]: Dixon I (1970) defines momentum, spin and the centre of mass (`pβSαβ = 0`). Dixon II/III were **not opened**, and "made it exact" is not a quoted claim. |
| Thorne–Hartle 1985, *PRD* 31: surface integrals in a buffer zone | **Verified** [S8]: pp. 1815–1837; the paper's term is "buffer region". |
| Damour's "effacement" | **Verified** through Blanchet [S7]: the "effacing principle" (Damour 1983a, Les Houches proceedings, pp. 59–144). The 1987 chapter title is verified from a reference list; the chapter itself was not opened. |
| Control-volume momentum theorem (Lamb): force on a body = momentum flux through any enclosing surface | **Corrected / partly verified** [S3]: Lamb (1906) Art. 346, p. 598 states the momentum balance of a **fixed element** (forces plus "flux of matter across the boundary, carrying its momentum with it"). The "any enclosing surface" body-force form was **not found in Lamb**. It is stated in TAN's classical cylinder argument [S9] and used in Wang's proof [S13]. |
| Lagally (1922, ZAMM): force on a source/sink = ρ × strength × outside flow | **Verified with conditions** [S10, S13]: *ZAMM* 2(6), 409–422. Lagally's paper is for **steady** flow and states the general dependence (sources, vortices, and the velocity at their location). The explicit `−ρQu₀` form is Wang's restatement, which Wang credits to Lagally "for the special case where the velocity field is steady". Unsteady extensions: [S11, S12]. |
| Two steady sinks attract by an inverse square ∝ the product of drain rates | **Verified** [S3, S13]: Lamb Art. 144, p. 192 (like sources attract), and Wang Corollary 4. Lamb's sources are the ends of narrow tubes, so fluid is conserved overall. |
| C. A. and V. Bjerknes (1870s–1906, *Fields of Force*): pulsating spheres as a **gravity** analogue | **Corrected** [S16, S17]: *Fields of Force* (1906) draws the analogy to **electric charges and magnetic poles with the sign reversed** (like phases attract), ⛔ not to gravity. The pulsation theory of **gravity** is credited to Kelvin (1871) and C. A. Bjerknes (1871) only by Wikipedia (citing Zenneck 1903), which is secondary only. |
| G. I. Taylor (1928): forces on a body in a converging stream | **Verified** [S14]: *Proc. R. Soc. A* 120, 260–283; the "horizontal buoyancy" `−V dp/dx`. Note that the "Taylor theorem" of the Lagally–Taylor–Kirchhoff triad is a different 1928 paper (*Proc. R. Soc. A* 120, 13–21) [S11, S14]. |
| Riemann 1853, two fragments, printed 1876 by Dedekind and Weber | **Verified** [S18]: pp. 496–506; edited by H. Weber "unter Mitwirkung von R. Dedekind". "Neue mathematische Principien" carries the footnote "Gefunden am 1. März 1853". |
| Riemann: "Gravitation und Licht" tied to his letter of 28 Dec 1853 | **Unverified**: no opened source gives it. |
| Riemann: the substance vanishes into "another world" | **Corrected** [S18, S17]: Riemann's word is **"Geisteswelt"** (the mental or spirit world, from Herbart's psychology). "Another world or dimension" is Wikipedia's retelling. |
| Riemann: inflow ∝ mass; incompressible, inertia-free; speed ∝ force (so v ∝ 1/r², not the river model's escape speed) | **Verified, from two fragments** [S18]: the drain, "incompressible homogene Flüssigkeit ohne Trägheit" and "seiner Masse proportionale Mengen" are in fragment 3 (dated), which also has a footnote with inflow ∝ the gravitational force; "Geschwindigkeit der Grösse der Kraft proportional" and the inverse-square force he must reproduce (his item 4) are in fragment 2 (undated), whose substance is conserved in each volume element. The river model's escape speed: [S21]. |
| Riemann: volume resistance → gravity and electrostatics; line-element resistance → light, heat, magnetism | **Verified** [S18]. The text says "electrodynamische oder magnetische" attraction and repulsion. |
| Riemann: stops on how to choose two memory kernels (Bottazzini & Tazzioli) | **Verified from Riemann himself** [S18]: the question near the end (l.38582–38586) is about the two functions weighting a substance element's earlier forms; the fragment then closes on the two-resistance split (l.38607–38616) without choosing them. The Bottazzini & Tazzioli attribution is **not verified** (sample pages only) [S20]. |
| Riemann: never computed it, no light-bending number | **Search statement only** [S18]: no deflection term found in pp. 494–506. |
| English translation: Baker–Christenson–Orde, *Collected Papers* (2004) | **Not opened.** It appears only in search summaries (Q14), and the MAA review page returned HTTP 521. |
| Newton 1675 letter to Oldenburg; 1679 letter to Boyle | 1675 to Oldenburg, "and later to Robert Boyle": **secondary only** [S17]. The year 1679 is **unverified**. |
| Soldner (1801): Newtonian deflection, half the GR value | **Corrected** [S23]: written 1801, **published 1804**; GR's prediction is "twice as high". |

## Checked, not cited

- **R. L. Kirkwood**, "Gravitational Field Equations", *Phys. Rev.* 95, 1051 (1954), DOI 10.1103/PhysRev.95.1051
  (APS abstract page, Q17). It gives ether equations of motion from Newton's law; the abstract shows no sink. The
  1953 paper named in Q17 was not found.
- **I. Simaciu, G. Dumitrescu, Z. Borsos**, "Acoustic force of the gravitational type", arXiv:1905.03622 (2019),
  Q23. Abstract: "under certain restrictive conditions, Bjerknes secondary forces are attractive"; "proportionate to
  the product of the virtual masses of the two bubbles." Not used (bubble acoustics, no drain).
- **Wikipedia, *Bjerknes force*** (Q4): superseded by the primary text [S16].
- **St Andrews, C. A. Bjerknes biography** (Q4): the WebFetch summary reports no mention of gravitation. Not used.
- **Papadopoulos**, "Physics in Riemann's mathematical papers", arXiv:1711.01755: cites "Gravitation und Licht" at
  "p. 532–538" of another edition. Not used; the pages here are from the 1876 edition [S18].
- **Analogue-gravity results** for Q19 and Q24 (acoustic black holes, the draining bathtub, analogue lensing): none
  treats the motion of the sink or vortex as a body. This is recorded in the section's NOT FOUND list.
