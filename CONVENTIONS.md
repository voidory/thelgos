# CONVENTIONS.md — LGO Method Explorer

**Version: 1.6 — not frozen (applied 2026-11-07 on the owner's acceptance of `drafts/CONVENTIONS-v1.6-draft.md`; owner sign-off B8 of this text pending)**  
**Owner / change authority: the project owner**  
**Date: 2026-11-07, Singapore**  
**Release 1 blockers: none. CrystOD comparison deferred until the programme is finished (B1). Diatomic fixtures (periods 2 and 3) and browser acceptance tests (C4) in place. v1.5 adds C3v, NH3, H2S, SiH4 and the hydrides NH, OH, HF.**

*v1.3 record:* the owner decided B1 and B4 and the ring orientation in B3, delegated B2, B3 (in-plane axis) and B5 to Claude, and gave standing permission (2026-11-06) for Claude to apply the drafted changes to this file. Every decision and its reasons are in `DECISIONS.md`. Checker changes are in `CHECKER-CHANGES.md` (K1–K19, plus the M7 data change).

## Review and evidence record

v1.0 was AI-prepared and owner-reviewed. Its “AI-checked” status preceded the numerical review that found the water-ordering problem, the SF6 construction-stage degeneracy and malformed Markdown. It did not certify numerical correctness. v1.1 corrected those problems, but its hardcoded-data audit and character-matched D2h PASS were insufficient. In particular, raw MolSym 1.3.1 associates x with B2u and y with B3u in the generated D2h matrices; the prior comparison concealed that swap. The old D2h PASS is retracted.

The attached v1.2 TODO records Claude's independent hand review of the physics, for the owner. This revision records the agent's executed engineering corrections; it does not turn that earlier review into owner approval of this new bundle. v1.1 was unavailable as a local file in this turn; this candidate reconstructs the accepted contract from the conversation and the original attachments, with explicit retained/removed content below.

New evidence: offline MolSym generation, geometry-first class mapping, explicit D2h label repair, data-file-driven Hamiltonians and oracle substitution, committed spectrum guards, a corrupt-pair negative test, schema validation, reproducibility via temporary regeneration and diffs, exact-value snapping, open-shell warnings, and a complete hash-locked checker environment. Twelve current fixtures execute, including neutral H3. Raw group-source snapshots and install evidence are packaged. No claim of new blind derivation: this is a regression audit of previously reviewed analytic results.

**Tech stack:** `CPython 3.12.7 + NumPy 2.0.2 engine; Pyodide 0.27.7 in a Web Worker; vanilla ES2022 + HTML/CSS/SVG + WebGL2 frontend; static hosting; no backend.`

Pin evidence: [Pyodide package list](https://pyodide.org/en/0.27.7/usage/packages-in-pyodide.html), [changelog](https://pyodide.org/en/0.27.7/project/changelog.html), [distribution lock](https://cdn.jsdelivr.net/pyodide/v0.27.7/full/pyodide-lock.json), checked 2026-11-07. The lock reports Python 3.12.7 and NumPy 2.0.2. This audit's host is **CPython 3.12.8** with NumPy 2.0.2: report that actual runtime, rather than claiming it is the production 3.12.7/browser test. Browser tests remain open.

## 1. Units, signs and angular factors

Use eV, Å, degrees at input and radians internally. On-sites are approximate isolated-orbital energies on an approximate vacuum reference; molecular eigenvalues are model results, not measured ionisation energies. The basis is orthonormal; overlap is the identity. For each stored pair a<b, the direction vector points from a to b. Compute that pair block once and set the reverse block to its transpose. Do not reconstruct a reverse block with an unreversed direction vector. Stabilising in-phase s–s hopping is negative; directional-orbital contact phases include AO signs, not just coefficients.

Reference: [Slater & Koster, Phys. Rev. 94, 1498 (1954)](https://doi.org/10.1103/PhysRev.94.1498). The +z-bond table is authoritative. Markdown bars are escaped inside code spans.

| Matrix element (+z bond) | Value |
|---|---|
| `<s_a\|H\|s_b>` | \(V_{ss\sigma}\) |
| `<s_a\|H\|p_z,b>` | \(+V_{sp\sigma}\) |
| `<p_z,a\|H\|s_b>` | \(-V_{sp\sigma}\) |
| `<p_z,a\|H\|p_z,b>` | \(V_{pp\sigma}\) |
| `<p_x,a\|H\|p_x,b>`, likewise y | \(V_{pp\pi}\) |
| `<s_a\|H\|d_z2,b>`, reverse types | \(V_{sd\sigma}\) |
| `<p_z,a\|H\|d_z2,b>` | \(+V_{pd\sigma}\) |
| `<d_z2,a\|H\|p_z,b>` | \(-V_{pd\sigma}\) |
| `<p_x,a\|H\|d_xz,b>`, likewise y/yz | \(+V_{pd\pi}\) |
| `<d_xz,a\|H\|p_x,b>`, likewise yz/y | \(-V_{pd\pi}\) |

Test-only reference: [pysktb 0.5.6](https://pypi.org/project/pysktb/0.5.6/), `pysktb._params.get_hop_int`. Map its `s,px,py,pz,dxy,dyz,dxz,dx2-y2,dz2,S` order by indices `[0,1,2,3,8,6,5,7,4]`. Independently expressed vector/tensor contractions agree for 124 deterministic directions; reversal-transpose parity also passes. s/p, s–d and p–d checked; d–d not checked. This reference is not an application dependency.

## 2. Axes, numbering and geometries

Generate exact symmetric positions from unit directions, not rounded coordinate tables.

| Case | Positions and atom order |
|---|---|
| Diatomic | z axis; atom 1 at −R/2, atom 2 at +R/2 |
| H2O | O=1 at origin, yz plane, z bisector; H2 on +y, H3 on −y, both +z |
| CH4 | C=1; H2–5 directions (1,1,1),(1,−1,−1),(−1,1,−1),(−1,−1,1), each divided by √3 |
| BH3 | B=1; H2–4 in xy at 0°,120°,240° from +x |
| NH3 | N=1 at origin, C3 along z, lone pair towards +z; H2–4 at azimuth 0°,120°,240° from +x, each at angle β from −z with cos(HNH)=cos²β−½sin²β (β=67.8549° for 106.67°); H2 in the xz mirror plane σv(0) |
| H2S | As H2O, with S=1 |
| SiH4 | As CH4, with Si=1 |
| NH, OH, HF | Heavy atom 1 at origin, H 2 on +z (C2v finite subgroup of C∞v, as for CO) |
| SF6 | S=1; F2–7 at +x,−x,+y,−y,+z,−z |
| Six H Oh | Same six directions, no centre; H1–6 |
| Cyclic H3 | xy plane, centroid origin; H1–3 at 0°,120°,240° |
| Allyl | yz plane, C2 at origin, C1 on +y and C3 on −y, both +z; π=px |
| Benzene, later release | xy plane; atoms k=1–6 at R(cos((k−1)π/3),sin((k−1)π/3),0). Atoms lie on ±x. C2-prime axes pass through opposite atoms (0°,60°,120°); C2-double-prime pass through opposite bond midpoints (30°,90°,150°). B3 decided: rings flat and parallel, above and below the metal, ring planes perpendicular to z (owner); one ring atom on +x (delegated). Same rule for the eight-membered rings in D8h. |

| Fixture | Pinned geometry | Source / qualification |
|---|---|---|
| H2O | O–H 0.958 Å; HOH 104.4776° | Retained [NIST input](https://cccbdb.nist.gov/expgeom2x.asp?casno=7732185) |
| CH4 | C–H 1.087 Å | Retained [NIST input](https://cccbdb.nist.gov/exp2x.asp?casno=74828) |
| BH3 | B–H 1.190 Å | [NIST r0 input](https://cccbdb.nist.gov/exp2x.asp?casno=13283313&charge=0), not re |
| SF6 | S–F 1.561 Å | Retained [NIST input](https://cccbdb.nist.gov/exp2x.asp?casno=2551624) |
| H2 | H–H 1.000 Å | Schematic test geometry, not equilibrium |
| H3+/neutral H3 | Edge 1.000 Å; radius 1/√3 Å | Schematic test geometry |
| Six H Oh | Radius 1.000 Å | Schematic test geometry |
| Four H Td | Radius 1.087 Å; edge R√(8/3) | CH4 ligand-only diagnostic |
| Allyl π | C–C 1.400 Å; angle 120° | Idealised fixture |
| B2 | B–B 1.5900 Å | NIST CCCBDB r_e (Gurvich); `evidence/cccbdb_B2.html` |
| C2 | C–C 1.2425 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_C2.html` |
| N2 | N–N 1.0977 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_N2.html`; owner decision B4 |
| O2 | O–O 1.2075 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_O2.html` |
| F2 | F–F 1.4119 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_F2.html` |
| CO | C–O 1.1282 Å; C atom 1 | NIST CCCBDB r_e (NIST diatomic); `evidence/cccbdb_CO.html` |
| Al2 | Al–Al 2.7010 Å | NIST CCCBDB r_e (Fu, Lemire et al. 1990); `evidence/cccbdb_Al2.html`; owner request M5 |
| Si2 | Si–Si 2.2460 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_Si2.html`; owner request M5 |
| P2 | P–P 1.8934 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_P2.html`; owner request M5 |
| S2 | S–S 1.8892 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_S2.html`; owner request M5 |
| Cl2 | Cl–Cl 1.9879 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_Cl2.html`; owner request M5 |
| NH3 | N–H 1.0124 Å; HNH 106.67° | NIST CCCBDB (reference 1966Herzberg); `evidence/cccbdb_NH3.html`; v1.5 |
| H2S | S–H 1.3356 Å; HSH 92.11° | NIST CCCBDB (reference 1975Coo/DeL); `evidence/cccbdb_H2S.html`; v1.5 |
| SiH4 | Si–H 1.4799 Å | NIST CCCBDB (reference 1955Boyd; listed 1.480); `evidence/cccbdb_SiH4.html`; v1.5 |
| NH | N–H 1.0362 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_NH.html`; v1.5 |
| OH | O–H 0.9697 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_OH.html`; v1.5 |
| HF | H–F 0.9168 Å | NIST CCCBDB r_e (NIST diatomic); `evidence/cccbdb_HF.html`; v1.5 |

Water positions are \((0,\pm R\sin(\theta/2),R\cos(\theta/2))\). The same formula defines C1/C3 for allyl. Geometry overrides live in presets, not duplicated checker constants.

The eleven diatomics B2, C2, N2, O2, F2, CO, Al2, Si2, P2, S2 and Cl2 are Release 1 fixtures with the geometries above (CCCBDB lists three decimals and the Cartesian distance to four; the four-decimal value is pinned). Do not interpret the H2 D2h subgroup as its exact point group: the exact linear homonuclear symmetry is D-infinity-h; D2h is the finite label-checking subgroup.

## 3. Orbital phases and local frames

Positive real normalisation constants multiply s=1; px=x, py=y, pz=z; dz2=2z²−x²−y²; dxz=xz; dyz=yz; dx2-y2=x²−y²; dxy=xy. These explicit signs override ambiguous phase-name shorthand.

\[
\mathbf f_d=\left((3n^2-1)/2,\sqrt3ln,\sqrt3mn,
\sqrt3(l^2-m^2)/2,\sqrt3lm\right),\qquad\sum_\mu f_{d,\mu}^2=1.
\]

Distinguish effective s-type donors from actual inward p sigma functions. For outward radial unit u, inward p sigma is −u·p; use actual p angular channels. To define local π axes, choose the Cartesian axis e least parallel to u (ties x,y,z), project/normalise it to v perpendicular to u, then w=u×v. Functions are v·p and w·p. Allyl uses common global px. Frame-chart changes are basis transformations, not physical jumps.

## 4. Colours and basis order

Red=positive amplitude, blue=negative. Lobe colour includes AO phase and MO coefficient. Fix-and-flip makes the lowest-numbered member positive. Contact lines use neither colour: solid stabilising, dashed destabilising. Basis order: atom index, then `s,px,py,pz,dz2,dxz,dyz,dx2-y2,dxy`; reduced/local donor functions are explicitly declared. Isolated eigenvector phase makes its first above-threshold coefficient positive. Display thresholds do not truncate computation. Degenerate spaces are compared by projectors, not individual solver vectors.

## 5. Group-data provenance and axis repair

Generator: [MolSym 1.3.1](https://github.com/NASymmetry/MolSym), `molsym.symtext.general_irrep_mats.pg_to_symels`. Independent reference: [Bilbao POINT](https://cryst.ehu.es/rep/point.html); fresh pages fetched 2026-11-07 and packaged as `evidence/bilbao_*.html`. Generation and all baseline checks are offline thereafter.

**Geometry first:** for each Bilbao column construct a Cartesian target operation from the rotation axis/angle or mirror normal. Find the generated conjugacy class containing that matrix. Only then compare characters in that class. Never choose a class by matching character values or its library display label. Every target matrix and matched operation index is saved in `checks/bilbao_*.json` and group metadata.

For D2h the target order is E,C2z,C2y,C2x,i,mz,my,mx; mirror subscripts are normals. Raw library data fail this mapping with maximum character difference 2. Repair the **irrep names**, not geometry: swap MolSym B2g↔B3g and B2u↔B3u, and normalise A1g/A1u to Ag/Au. `library_label_map` is explicit in the file. Corrected assertions: **x=B3u, y=B2u, z=B1u**. C2v assertions: x=B1,y=B2,z=A1. D3h explicitly contains C2 about +x. Thus the earlier raw D2h PASS is not repeated.

For the hexagonal table setting, record the chosen Cartesian alignment: Bilbao's 2_120 family maps to C2-prime through +x; 2_100 maps to the perpendicular C2-double-prime family in D6h. This is the fixed setting transformation used here, not a swap chosen to fit characters. The table's improper classes follow the corresponding geometric operations. The target matrices in JSON are authoritative. Benzene's atom/bond families are separately asserted against coordinates; no automatic renaming from character resemblance is allowed.

| Coverage | Supplied files | Status |
|---|---|---|
| Release 1 | Oh,Td,D3h,C2v,D2h | Geometry-mapped Bilbao comparison passes; D2h includes documented label repair |
| Later crystallographic | D4h,D2d,C2h,D6h | Geometry-mapped comparison passes; D6h uses the stated setting and B3 default |
| v1.5 | C3v | Geometry-mapped comparison passes (`evidence/bilbao_C3v.html`, point group 3m); Bilbao's m₁₋₁₀ is in hexagonal axes, so its target is the mirror with Cartesian normal (3/2,−√3/2,0), in the class of σv(0), the xz plane (checker patch K21) |
| Release 4 | D8h | Not crystallographic; closure and character orthogonality only, no Bilbao coverage claim |

JSON contains operation matrices, generating indices, multiplication tables, conjugacy classes, characters and irrep dimensions. Snap exact values to binary representations of 0,±1,±1/2,±√2/2,±√3/2, etc., including ±√2 character values. Assert maximum snapping residual <1e−12 and record it for matrices and characters. Do not merely round near-zero trigonometric artifacts into long decimal approximations.

Offline reproduction: `python scripts/generate_groups.py --out /tmp/regenerated`. Fresh retrieval: `python scripts/fetch_bilbao.py --out /tmp/bilbao`, then `python scripts/generate_groups.py --evidence /tmp/bilbao --out /tmp/regenerated`. Fresh results are candidates for review, not silent replacements of committed evidence.

## 6. Aliases and construction schedule

Keep `canonical_id` separate from source-backed `lecture_alias`. MO numbering is valence-only, ascending energy. Pairs are ordered by lowest atom index, bonding before antibonding. Exact lecture LGO numbering still needs the slides.

Canonical SALCs:

- Water H2,H3: L+=(1,1)/√2; L−=(1,−1)/√2.
- CH4 H2–5: a1=(1,1,1,1)/2; Lx=(1,1,−1,−1)/2; Ly=(1,−1,1,−1)/2; Lz=(1,−1,−1,1)/2.
- Oh +x,−x,+y,−y,+z,−z: a1g=(1,1,1,1,1,1)/√6; Lx=(1,−1,0,0,0,0)/√2; Ly=(0,0,1,−1,0,0)/√2; Lz=(0,0,0,0,1,−1)/√2; eg1=(−1,−1,−1,−1,2,2)/√12; eg2=(1,1,−1,−1,0,0)/2.

B5 decided under delegation (DECISIONS.md): skip a construction stage whose newly enabled coupling block is identically zero; show **“No ligand–ligand coupling in this model.”** Do not secretly turn on H–H edges for ligand stages while displaying the H–H-off molecular spectrum. SF6 retains its nonzero combine stage. Visual K is the number of effective stages after filtering, not an unconditionally retained combinatorial pair count. This is a convention awaiting UI implementation/testing.

## 7. Electron counting and warnings

Neutral-atom counting over the included basis for molecules; named oxidation-state counting for advanced complexes. Save the counting convention and occupations with every result. Counts: H2=2; water=8; CH4=8; BH3=6; SF6 reduced sigma basis=12; allyl cation/radical/anion=2/3/4; H3+/H3/H3−=2/3/4; NH3=8; H2S=8; SiH4=8; NH=6; OH=7; HF=8. Omitted SF6 F orbitals do not silently contribute electrons to the reduced basis.

Any level of degeneracy d>1 occupied by n with 0<n<2d triggers **occupation_not_unique**, **Jahn_Teller_warning**, and **spin_state_not_determined**. This conservative warning policy flags the limits of the one-electron model; it is not a proof of a particular spin ground state or distortion (especially for linear diatomics). Distinguish explicit Hund filling from equal-occupation ensembles. Tested warning fixtures: neutral H3 e-prime¹, six-H t1u⁴, four-H t2², and the real molecules O2 (πg²), C2 (πu²), NH (π²) and OH (π³). B2 does not fire it, because the model gives B2 a closed shell; that is reported as a Spin-state fidelity disagreement (A.6), not forced.

## 8. Bond orders and contact ledger

Counting bond order is (bonding−antibonding electron count)/2, shared over the declared equivalent-bond set; nonbonding electrons count zero. Default counting is not an observable.

\[
P_{\mu\nu}=\sum_i n_i C_{\mu i}C_{\nu i},\quad
p_{ab}=P_{ab},\quad W_{AB}=\sum_{\mu\in A,\nu\in B}P_{\mu\nu}^{2}.
\]

Use named Coulson p for one orbital per atom and Wiberg W for multiple orbitals. Report assumed sigma bonds separately for π-only calculations. Do not force these definitions to agree.

The contact ledger counts each unordered real-basis off-diagonal pair once with contribution 2 c_mu c_nu H_mu,nu, plus all diagonal c_mu² H_mu,mu terms. Its total is the eigenvalue for a normalised eigenvector. A raw contact weight c_mu c_nu t is half that unordered-pair energy contribution; label the distinction explicitly. See restored contractual Appendix A.7; browser ledger implementation remains untested.

## 9. Parameter and preset contract

### 9.1 On-sites

`params/onsite.csv`: `element,orbital,energy_eV,source`. Single table for every element: the orbital potential energies of the [teaching supplement](https://chem.libretexts.org/Bookshelves/Inorganic_Chemistry/Supplemental_Modules_and_Websites_(Inorganic_Chemistry)/Advanced_Inorganic_Chemistry_(Wikibook)/01%3A_Chapters/1.07%3A_Diatomic_Molecular_Orbitals), whose table is credited to Miessler, Fischer & Tarr 2014. Snapshot: `evidence/libretexts_diatomic_MOs.html`. Owner decision M7 (2026-11-07) replaced the earlier composite, which took Harrison's solid-state table for B, C, N, O, Si, P and S, because the owner's module uses Miessler. The teaching values are configuration-averaged ionisation energies, a few eV deeper than Harrison's computed term values. On-sites are not Harrison's; only the couplings are (§9.2).

| Element | s / eV | p / eV |
|---|---:|---:|
| H (n=1) | −13.61 | absent |
| B (n=2) | −14.05 | −8.30 |
| C (n=2) | −19.43 | −10.66 |
| N (n=2) | −25.56 | −13.18 |
| O (n=2) | −32.38 | −15.85 |
| F (n=2) | −40.17 | −18.65 |
| Al (n=3) | −11.32 | −5.98 |
| Si (n=3) | −15.89 | −7.78 |
| P (n=3) | −18.84 | −9.65 |
| S (n=3) | −22.71 | −11.62 |
| Cl (n=3) | −25.23 | −13.67 |
| He (n=1) | −24.59 | absent |
| Li (n=2) | −5.39 | not listed |
| Be (n=2) | −9.32 | not listed |
| Ne (n=2) | −48.47 | −21.59 |
| Na (n=3) | −5.14 | not listed |
| Mg (n=3) | −7.65 | not listed |
| Ar (n=3) | −29.24 | −15.82 |
| K (n=4) | −4.34 | not listed |
| Ca (n=4) | −6.11 | not listed |

All p directions share the p on-site. No main-group d parameters. Nitrogen was added by owner decision B4. *v1.5 erratum (V1):* the table signed in v1.4 still showed the earlier Harrison composite for B, C, N, O and S and omitted Al, Si, P and Cl; the values above are the ones in `params/onsite.csv`, which every result used. No computed value changed. *v1.6 (W1):* the rest of the same table (He to Ca) is added; the teaching table gives Li, Be, Na, Mg, K and Ca an s value only, so a molecule needing their p is "not built yet: element without parameters" (owner decision R3-1).

### 9.2 pairs.csv columns and lookup

`element_a,element_b,channel,eta,K_eV_A2,value_eV,units,distance_law,s_carrier,scale,provenance_kind,source,override_reason`

- `element_a,element_b`: unordered element pair; duplicates per channel are errors.
- `channel`: `ss_sigma`, `sp_sigma`, `pp_sigma`, `pp_pi` in Release 1. sp direction signs are handled geometrically, not by duplicating reversed rows.
- `eta,K_eV_A2`: dimensionless coefficient and 7.62 eV Å² conversion for `harrison_inverse_square`.
- `value_eV`: numeric fixed value for `fixed`; blank for inverse-square. `off` is an explicit zero, not missing data.
- `units`: `eV`.
- `distance_law`: enum `harrison_inverse_square`, `fixed`, `off`.
- `s_carrier`: `H` if hydrogen triggers the scale; otherwise `none`. This is a **scale carrier**, not a claim that the s orbital belongs to hydrogen.
- `scale`: 0.6 for H-containing pair blocks, otherwise 1. Apply once per block including H–H, never once per H atom.
- `provenance_kind`: enum `harrison`, `harrison_extrapolated_H`, `pedagogical_override`, `lecture_decoupling`.
- `source`: mandatory source URL or named decision provenance.
- `override_reason`: mandatory for an override; H rows explicitly record extrapolation and pedagogical scaling.

Precedence: **preset matrix-element override → element-pair/channel row → explicit channel default → error**. Defaults live in `params/channel_defaults.json`. *v1.6 (W2, W3):* the defaults are Harrison's law for the four channels with the coefficients below, scale 0.6 when the pair contains H (s_carrier H) and 1 otherwise, the same rule as every H row. A default applies only to an element pair with no row at all in `pairs.csv`, and only when both elements have on-sites; a pair that has rows keeps the rule that a missing channel is an error (O–H has no pp_pi), and an element without on-sites still errors. Every pair that had a row keeps it, so no earlier value changes. A no-match lookup must raise an error, never silently return zero. Missing orbital channels are not requested; absence is not a fit to zero. Reverse-pair Hermiticity follows §1.

For `harrison_inverse_square`, V=scale×eta×K/d². Coefficients remain ss=−1.40, sp=+1.84, pp sigma=+3.24, pp pi=−0.81: Harrison's 1980 book values as transcribed in `evidence/harrison.txt`. The same file lists a second set attributed to Shi & Papaconstantopoulos (2004), which is not used (DECISIONS.md M2). **All H-containing couplings extrapolate this scheme beyond the Harrison table's H coverage.** v1.5 adds Si–H ss and sp rows of the same form (V3). The 0.6 factor is a declared pedagogical correction, not a published coefficient or measured fit. The two strict ordering constraints are: BH3's empty pz (a2″) below the antibonding a1′ (lower bound); and in lecture water, the bonding 1b2 above the decoupled O 2s (upper bound). With the M7 on-sites they give 0.4234786093<sH<1.0312679542; with the earlier composite they gave 0.4914651002<sH<0.8938796155. Two further orderings hold at 0.6 and are kept on the owner's instruction to follow Miessler (DECISIONS.md B2): water's LUMO is 3a1*, below 2b2*, which holds for sH<0.666786; and CH4's LUMO is a1*, below t2*, which holds for sH<0.970429. Both are recorded as textbook orderings, from Miessler's figures as recalled, pending the owner's check against the book. They are not asserted by the checker. Together the four orderings allow 0.4235<sH<0.6668. B2 decided under delegation: 0.6 (DECISIONS.md). Because it was chosen to reproduce orderings, those orderings are not independent fidelity evidence. The audit reads the CSV rather than these prose numbers. The v1.5 molecules add no constraint: every ordering they show holds for any factor between 0.30 and 1.20, except NH3's two empty levels (a1* and e*), which swap near 0.80, outside the window.

Neighbour lists are explicit: H2O/CH4/BH3/NH3/H2S/SiH4 central-H on, H–H off; NH/OH/HF their one bond; SF6 central-F and all F–F on; ligand-only H diagnostics have their declared edges on. No distance cutoff silently re-enables absent edges.

### 9.3 presets.json schema

`params/presets.schema.json` is JSON Schema Draft 2020-12; validate before calculation. Top-level `schema_version` and `presets` are mandatory. Each preset requires `id`, `edge_policy` (`edges_on`, `edges_off`), `coupling_overrides`, `hydrogen_scale`, `geometry_overrides`, and `provenance_tag`; optional `note`. Unknown fields are rejected. Coupling overrides identify both element/orbital endpoints, `value_eV` and a provenance enum. Geometry overrides include bond length/radius/angle with range checks. Checks also reject a missing edge policy and a row/preset hydrogen-scale conflict. Array IDs must be unique.

| Lecture preset | Values / meaning |
|---|---|
| H2O | O 2s–H forced exactly zero; H scale 0.6; H–H off; O–H 0.958 Å, angle 104.4776° |
| SF6 | Equal to model; actual inward F p basis; cis −3.166242667913 eV, trans −2.532994134331 eV at S–F 1.561 Å |
| Six H Oh | Radius 1 Å; H scale 0.6; cis −3.200400 eV, trans −1.600200 eV |
| NH, OH, HF (v1.5, V5) | Central s–H 1s forced exactly zero, as in H2O; H scale 0.6; geometry as §2. NH3, H2S and SiH4 have no lecture preset (V5): their textbook diagrams are the model ones |

The SF6 table numbers are consequences of the row law and preset geometry, not duplicate mutable fixed overrides. With cis distance √2R and trans 2R, inward-p contraction gives tc=(Vpppi−Vppsigma)/2 and tt=−Vppsigma. The combine-two-axes symmetric-pair block has diagonal epsilon_F+tt and off-diagonal 2tc. Its split is 4|tc|=12.664970671653 eV. The v1.0 cis-zero/trans−5 override stays removed. No slide-exact calibration or aliases are asserted without slides.

`checks/unscaled_reference_couplings.csv` is generated from actual rows with H scaling cancelled. Includes O–H sp=15.277130068296 eV at 0.958 Å, O–H ss, C–H/B–H ss and sp, H2/H3, six-H cis/trans, four-H edge and allyl pp_pi. Uses `d_angstrom` and full channel names.

## 10. Tolerances, locks and reproducibility

Float64. Hamiltonian symmetry ≤1e−12 eV; eigen-residual norm ≤1e−10 eV; orthonormality maximum-entry error ≤1e−10; ledger error ≤1e−10 eV. Clean irrep weight >0.999. Degeneracy clustering uses a total-span tolerance of 1e−8 max(1 eV, ||H||2), never chained near-equalities. Current small fixtures use an absolute 1e−8 eV clustering threshold (stricter than that bound); do not mistake this baseline helper for the application's full clustering implementation.

`requirements-checker.lock` is a complete 50-distribution hash-pinned graph for CPython 3.12/Linux x86_64. Install with `python -m pip install --require-hashes -r requirements-checker.lock`. A clean environment was installed from hash-verified wheels and `pip check` passed. `evidence/hash_install.txt` records it. Bootstrap pip is tooling, not an application/checker dependency row. Tests run on host Python 3.12.8, not falsely labelled 3.12.7. Hash-locked means distribution integrity/reproducibility in this target environment, not support for every platform.

`make check-baseline` generates groups and checks into a temporary directory, runs `diff -r` against both committed directories, then `sha256sum -c SHA256SUMS`; any mismatch is a failure. It never overwrites committed groups/checks. Python bytecode writing is disabled; an additional tree-hash comparison checks nonmutation. The parameter mutation uses a temporary CSV copy; spectra are compared to committed numerical snapshots in `evidence/reference_spectra.json`, so changing the same data used by both assembly and analytic substitution cannot create a false pass.

The committed spectrum snapshot preserves the reviewed numbers; it was reconstructed from the analytic results, not an independent experimental dataset. Numeric source data are not “simulated measurements”. Every committed file except SHA256SUMS itself is hashed. Store the exact UTF-8 convention-file SHA-256 with every result; no self-referential embedded digest.

## 11. CrystOD and release gates

CrystOD must come from the owner's **original supplied zip**. Default extended-Hückel route only, no PySCF. No zip is accessible in the current thread, so version/archive name/archive hash are **unavailable**, installation and real comparisons **not run**. These are blocking unknowns, not placeholders that an agent may fill by guessing.

`compare_crystod.py` computes `checks/crystod_manifest.json`; it does not read a trusted PASS flag. It hashes original archive bytes, extracts an embedded version, and invokes an explicitly supplied adapter with geometry/basis/edges/electron count only. The adapter must run the actual archive's default EHT route; expected oracle eigenvalues/order are not sent. Each result must record full point group, finite-subgroup irrep mapping, orthogonal axis map, basis description/differences, and independently justified level identities. The script compares point group, shared-irrep degeneracy, mapped labels and qualitative order within each shared irrep. Numerical energies are used only to rank levels; no energy values, compositions or bond orders are compared. Distinguish exact H2 D-infinity-h from its D2h subgroup.

**Appendix A wins over CrystOD.** Log disagreements for the owner, never change oracles to match the external checker. No common irrep, ambiguous level identities, missing metadata or incomplete required cases block the comparison. Different basis sizes are recorded rather than silently called agreement.

The adapter protocol is implemented and unit-tested with explicitly mock inputs, but **no native CrystOD adapter can be claimed integrated or executed until the archive/API is inspected**. This is the narrowed outstanding part of A4/C1. The current script exits BLOCKED for the missing zip; `make check` consequently fails for the right reason. Filling a manifest manually cannot make it pass because the manifest is recomputed in a temporary directory. This is a provenance gate, not tamper-proof attestation of arbitrary user-written adapters.

**Comparison rules (from CHECKER-CHANGES K5–K11).** Each adapter level carries `contract_basis_weight` in [0,1] and a `descriptor` {character, dominant_fragment} computed from the checker's own eigenvectors by the K4 rule: character is the sign of 2Σc_i c_j H_ij over central–ligand pairs (all inter-atomic pairs when there is no central atom), averaged over the degenerate set, with |x| ≤ 1e−6 eV nonbonding; dominant_fragment is central when the central-atom weight exceeds 0.5. Levels with weight below 0.5 are excluded before counting. Within-irrep order is compared through descriptors, not adapter-chosen ids. Point-group, shared-count and degeneracy mismatches block (DISAGREEMENT_ORACLES_PREVAIL); descriptor-order differences are logged for the owner without failing (ORDER_REVIEW_FOR_OWNER). Linear molecules take D∞h/C∞v from geometry.

**Deferral (B1).** The CrystOD comparison is deferred until the programme is finished. `make check` runs the baseline and engine tests; the gate above runs as `make check-release` and stays BLOCKED until the original archive and a real adapter are supplied.

Release 1 gates (all closed except C1, which owner decision B1 defers):

- C1 (deferred by B1): original archive identity, tested real default-EHT adapter and executed qualitative comparisons.
- C2: closed. Source-backed B2, C2, N2, O2, F2, CO, Al2, Si2, P2, S2 and Cl2 geometries, explicit blocks (A.6) and regression rows (Appendix B) are in place.
- C4: closed 2026-11-07. `web/tests/web_test.js` runs in Pyodide 0.27.7 in Chromium and checks every case: both slider end points (s = 1 against `evidence/reference_spectra.json` to 1e−9 eV; s = 0 against the on-site table to 1e−12 eV, with multiplicities); the ledger identity at s = 1 and s = 0.5 (1e−10 eV); counting bond orders against A.5; and Watch mode playing to the finished diagram. It is not part of `make check`, because it needs the network.
- B7: closed 2026-11-07. The owner signed off v1.5 ("ok i sign off v1.5") after accepting `drafts/CONVENTIONS-v1.5-draft.md` and reviewing the applied text. The signed text is this file at commit 32d6376, sha256 `312ea75e7a92c46aadc24716ad4de8211b136223c15e1d402b42c4c3c2c7f62b`. The only edits made afterwards are this status record, the header and the changelog row.
- B6: closed 2026-11-07. The owner signed off v1.4 ("I sign off v1.4") after reading `drafts/signoff-summary-v1.4.md`, which lists the delegated decisions (B2, B3 in-plane axis, B5, M1–M4) and the points resting on weaker sources. The signed text is this file at commit 81cf760, sha256 `ca861a1d6f47f5fe9962eb229a3fbdf4313c22ed2e83083dcfa8fdbe71b98d36`. The only edits made afterwards are this status record, the header and the changelog row.

Slides/aliases are nonblocking. Later d-channel parameters remain outside Release 1. Spec-document edits in tier D are reserved for the reviewer after approval, not silently applied here.

## 12. Change log and retained content

| Version | Date | Change |
|---|---|---|
| 0.1 | 2026-11-07 | Original skeleton |
| 0.2 / 1.0 | 2026-11-07 | First pass and owner review; AI-check label preceded numerical audit |
| 1.1 | 2026-11-07 | H scale 0.6; SF6 lecture=model; Markdown repair; group/angular/stack evidence; incomplete CrystOD gate; reconstructed rather than bytewise-edited legacy file |
| 1.2 | 2026-11-07 | A1 geometry-first mappings and D2h g/u label repair with Cartesian assertions; A2 file-driven data and corrupt-pair negative test; A3 temp regeneration/diffs/hashes; A4 computed CrystOD manifest/protocol and explicit real-adapter blocker; A5 restore column contract and define replacement A.7 with legacy-text limitation below; A6 typed pair fields and error-on-missing lookup; A7 preset schema; A8 exact snapping/residuals; A9 neutral-H3 and diagnostic warnings; A10 complete hashed lock/clean install; A11 full unscaled-coupling table; A12 benzene default and geometry assertions; A13 header/review includes CrystOD, diatomic and browser blockers. Owner approval pending. |
| 1.3 | 2026-11-06 | Applied by Claude under the owner's standing permission. Owner decisions B1 (CrystOD deferred), B4 (N2 in), B3 orientation; delegated B2 (0.6), B3 in-plane axis, B5 (skip zero stages), M1–M4 (DECISIONS.md). Nitrogen on-sites and pairs; six diatomic geometries, A.6 blocks, Appendix B rows; CrystOD comparison rules K5–K11 written into §11; coefficient-set provenance in §9.2. Not frozen: owner sign-off B6 pending. |
| 1.3.1 | 2026-11-06 | Applied by Claude under the owner's standing permission. Period-3 diatomics Si2, P2, S2 at the owner's request (DECISIONS.md M5): Si and P on-sites and Si–Si, P–P pair rows from the same Harrison transcription; CCCBDB geometries and ground states; Appendix B rows; checker patch K18. Al2 and Cl2 wait for a sourced Al/Cl on-site. Not frozen: owner sign-off B6 pending. |
| 1.3.2 | 2026-11-07 | Applied by Claude under the owner's standing permission. Al2 and Cl2 added (DECISIONS.md M5, checker patch K19): Al and Cl on-sites from the teaching table (snapshot packaged), CCCBDB geometries and ground states, Appendix B rows. Al2's spin disagrees with experiment and is recorded as a fidelity warning. Not frozen: owner sign-off B6 pending. |
| 1.4 | 2026-11-07 | Applied by Claude on the owner's decision M7: all on-sites from the Miessler teaching table (§9.1); couplings unchanged (Harrison η, H scale 0.6). Every fixture involving B, C, N, O, Si, P or S has new Appendix B values. Unchanged: the H-only sets, F2, Al2, Cl2. No qualitative result changed: orderings, spin agreements and counting bond orders are as before. H-scale window updated. Not frozen: owner sign-off B6 pending. |
| 1.4 (frozen) | 2026-11-07 | Owner sign-off B6. Signed text: commit 81cf760, sha256 `ca861a1d6f47f5fe9962eb229a3fbdf4313c22ed2e83083dcfa8fdbe71b98d36`. No content change. |
| 1.5 | 2026-11-07 | Applied on the owner's acceptance of `drafts/CONVENTIONS-v1.5-draft.md` ("accept v1.5 as recommended"; V5 as revised by the owner). V1 erratum (§9.1 table now the Miessler values in use); C3v (K21); NH3, H2S, SiH4, NH, OH, HF fixtures with CCCBDB geometries; Si–H rows; lecture presets for NH, OH, HF; A.3, A.5, A.6 and Appendix B rows. No earlier fixture's values changed. Not frozen: owner sign-off B7 pending. |
| 1.6 | 2026-11-07 | Applied on the owner's acceptance of `drafts/CONVENTIONS-v1.6-draft.md` ("yes"). W1: on-site table completed from the same Miessler snapshot (He to Ca). W2/W3: universal-law channel defaults for element pairs with no row, both elements in the table (checker patch K22). No fixture value changes; `checks/` regenerate byte-identical. Not frozen: owner sign-off B8 pending. |
| 1.5 (frozen) | 2026-11-07 | Owner sign-off B7. Signed text: commit 32d6376, sha256 `312ea75e7a92c46aadc24716ad4de8211b136223c15e1d402b42c4c3c2c7f62b`. No content change. |

**A5 disposition:** the pairs.csv column specification is restored in §9.2. Exact legacy Appendix A.7 text was not present in the accessible attachments or earlier inline answer; it is not silently reconstructed as a quotation. The unknown legacy A.7 is **removed as an unrecoverable legacy section**, because its original content cannot be verified; it is replaced by the explicitly newly authored ledger/limit/invariance contract below. If the owner supplies the original, review whether it contains additional obligations before freezing. No claim that its exact old content has been recovered.

## Appendix A. Analytic oracles

### A.1 Small blocks

\[
D(a,b;g)=\begin{pmatrix}a&g\\g&b\end{pmatrix},\quad
E_\pm=(a+b)/2\pm\sqrt{((a-b)/2)^2+g^2}.
\]

\[
A(a,b,c;g,h)=\begin{pmatrix}a&0&g\\0&b&h\\g&h&c\end{pmatrix},
\]

with characteristic polynomial (E−a)(E−b)(E−c)−g²(E−b)−h²(E−a). Do not reduce a three-level block to bright/dark two-level form without the required central degeneracy or a vanishing coupling.

### A.2 Generic s-donor sigma blocks

Ligand–ligand edges off. v values are the effective channels. d rows are later-release analytic contracts, not executed advanced-complex fixtures.

| Block | Matrix/energy | Copies |
|---|---|---:|
| Oh a1g | D(eps_s,eps_L;√6 vss) | 1 |
| Oh t1u | D(eps_p,eps_L;−√2 vsp) | 3 |
| Oh eg | D(eps_d,eps_L;√3 vsd) | 2 |
| Oh t2g | eps_d | 3 |
| Td a1 | D(eps_s,eps_L;2 vss) | 1 |
| Td t2 | A(eps_p,eps_d,eps_L;−2vsp/√3,2vsd/√3) | 3 |
| Td e | eps_d | 2 |
| D4h a1g | A(eps_s,eps_d,eps_L;2vss,−vsd) | 1 |
| D4h b1g | D(eps_d,eps_L;√3vsd) | 1 |
| D4h eu | D(eps_p,eps_L;−√2vsp) | 2 |
| D4h b2g,eg | eps_d | 1+2 |
| D4h a2u | eps_p | 1 |

Td Lx couples px and dyz, cyclically. D4h ligand order +x,−x,+y,−y; a1g all-plus/2, b1g=(1,1,−1,−1)/2; eu axis differences. Signs change with SALC phase, eigenvalues do not.

### A.3 Release 1 molecular blocks

H–H off in H2O,CH4,BH3; H factors included in v.

| Case | Block |
|---|---|
| H2 | eps_H+t and eps_H−t, t<0 |
| CH4 a1 | D(eps_Cs,eps_H;2vss) |
| CH4 t2 | D(eps_Cp,eps_H;−2vsp/√3), three copies |
| BH3 a1-prime | D(eps_Bs,eps_H;√3vss) |
| BH3 e-prime | D(eps_Bp,eps_H;−√(3/2)vsp), two copies |
| BH3 a2-double-prime | eps_Bp, central pz |
| H2O a1 | A(eps_Os,eps_Op,eps_H;√2vss,−√2vsp cos(theta/2)) |
| H2O b2 | D(eps_Op,eps_H;−√2vsp sin(theta/2)) |
| H2O b1 | eps_Op, central px |

Lecture water sets s–H=0 explicitly, not by symmetry.

| Case (v1.5) | Block |
|---|---|
| NH3 a1 | A(eps_Ns,eps_Np,eps_H;√3vss,√3vsp cos β), β the angle of each N–H bond from −z |
| NH3 e | D(eps_Np,eps_H;−√(3/2)vsp sin β), two copies |
| H2S | The H2O blocks with S on-sites and theta=92.11° |
| SiH4 | The CH4 blocks with Si on-sites |
| NH, OH, HF σ | (s, pz, H 1s) with rows (eps_s,0,vss), (0,eps_p,−vsp), (vss,−vsp,eps_H) for H on +z; lecture vss=0 |
| NH, OH, HF π | eps_p, two copies, nonbonding | For actual inward F p donors in SF6, lambda_a=eps_F+4tc+tt; lambda_t=eps_F−tt; lambda_e=eps_F−2tc+tt. a1g is D(eps_Ss,lambda_a;−√6 Vsp_sigma); t1u is D(eps_Sp,lambda_t;−√2 Vpp_sigma), three copies; eg is lambda_e, two copies, fluorine-only. No sulfur d functions.

### A.4 Ligand-only and pi spectra

Four H Td: eps_H+3t (a1), eps_H−t (t2, three). Six H Oh: eps_H+4tc+tt (a1g), eps_H−tt (t1u, three), eps_H−2tc+tt (eg, two). Cyclic H3: eps_H+2t (a1-prime), eps_H−t (e-prime, two). Allyl: eps_C+√2t (b1), eps_C (a2), eps_C−√2t (b1), t<0; vectors (1,√2,1)/2, (1,0,−1)/√2, (1,−√2,1)/2.

Six-H order a1g<t1u<eg requires 2tc+tt<0 and tt>tc. Spectral trace equals the sum of on-sites. Neutral H3 has the same one-electron spectrum as H3+ but a different filling and warnings.

### A.5 Bond-order oracles

Allyl cation: counting 1/2 per C–C, Coulson 1/√2. H3+: counting 1/3 per H–H, Coulson 2/3. BH3,CH4,lecture water: counting 1 per bond under the stated occupations. Diatomics, counting: B2 1, C2 2, N2 3, O2 2, F2 1, CO 3, Al2 1, Si2 2, P2 3, S2 2, Cl2 1. SF6: counting 2/3 per S–F if four bonding and two ligand-only MOs are doubly occupied. v1.5, counting: NH3 4/3, H2S 3/2, SiH4 1 (model; for NH3 and H2S the very slightly bonding lone pair counts as bonding, and taking it as nonbonding, as textbooks do, gives 1, V5); NH, OH, HF 2 (model) and 1 (lecture). Wiberg values depend on mixing and on-sites; do not set them equal to counting values. pz being BH3 LUMO requires the actual parameter ordering, not symmetry alone.

### A.6 Diatomic blocks

Bond along z; atom 1 at −R/2, atom 2 at +R/2; basis s, px, py, pz on each atom; all channels from `pairs.csv` at R; no H scale. Homonuclear D2h finite subgroup (x=B3u, y=B2u, z=B1u):

| D∞h | D2h | Block | Copies |
|---|---|---|---:|
| σg | Ag | D(eps_s+V_sss, eps_p−V_pps; −V_sps) | 1 |
| σu | B1u | D(eps_s−V_sss, eps_p+V_pps; +V_sps) | 1 |
| πu | B2u+B3u | eps_p+V_ppp | 2 |
| πg | B2g+B3g | eps_p−V_ppp | 2 |

CO (C2v, C atom 1): σ (A1) is the 4×4 block in (s_C, pz_C, s_O, pz_O) with rows (eps_sC, 0, V_sss, V_sps), (0, eps_pC, −V_sps, V_pps), (V_sss, −V_sps, eps_sO, 0), (V_sps, V_pps, 0, eps_pO); π (B1+B2) is D(eps_pC, eps_pO; V_ppp), two copies. In the finite subgroup a π pair spans two irreps, so degenerate sets are labelled by their exact irrep decomposition and mapped to D∞h/C∞v names by a fixed lookup.

Reference orderings: B2, C2, N2 πu below σg(2p); O2, F2 σg(2p) below πu, with σg(2p) the higher of the two σg levels. Reference ground states from the CCCBDB pages: B2 ³Σg, C2 ¹Σg⁺, N2 ¹Σg, O2 ³Σg, F2 ¹Σg, CO ¹Σ, Al2 ³Σg, Si2 ³Σg, P2 ¹Σg, S2 ³Σg, Cl2 ¹Σg (CCCBDB lists Al2's lowest level with degeneracy 4; every candidate ground state is a triplet, so the spin comparison stands). No orbital-ordering reference is sourced for the period-3 diatomics, so only their spin state is compared. v1.5 hydrides, from the CCCBDB pages: NH ³Σ⁻, OH ²Π, HF ¹Σ; the model agrees for all three (π², π³ open; π⁴ closed). No ordering is compared. A disagreement is recorded as an Orbital-ordering or Spin-state fidelity warning, never tuned away. Asserted: O2 ends in a partly filled πg with the warnings; CO has counting bond order 3 and a σ HOMO with carbon weight above 0.5.

Result under the current scheme (DECISIONS.md M1): B2, C2 and N2 disagree on ordering; B2 and C2 also on spin state; O2, F2 and CO agree. Si2, P2, S2 and Cl2 agree on spin state; Al2 disagrees (model closed shell, experiment triplet). The model puts σg(3p) below πu in all five.

### A.7 Replacement ledger, limit and invariance obligations

This is newly authored replacement content, not a recovered legacy A.7. For normalised real c, diagonal sum plus unordered pair terms 2c_mu c_nu H_mu,nu equals c-transpose H c. Verify this along slider endpoints and every displayed ledger state. Zero coupling returns isolated on-sites; MO count equals basis size. Global eigenvector sign cannot change physical outputs; rotations within a degenerate space preserve energy and summed composition/projector. With water theta→180°, the a1 pz coupling vanishes and b2 reaches magnitude √2vsp. Smooth bending must preserve continuous subspaces even where individual eigenvectors are not unique. Later d-mode weak-coupling splittings require separate advanced fixtures. The current candidate does not claim executed browser versions of any of these tests.

## Appendix B. Regression energies and coverage

Energies in eV; repeated levels have multiplicities shown. Values are model calculations, not measurements. Machine precision, occupations and warnings: `checks/release1_energies.csv` and `release1_validation.json`.

| Fixture | Levels, lowest first |
|---|---|
| H2O lecture | a1(O2s) −32.380000; b2 −25.039236; a1 −22.746843; b1 −15.850000; a1 −6.713157; b2 −4.420764 |
| H2O model | a1 −37.115699; b2 −25.039236; a1 −20.181834; b1 −15.850000; a1 −4.542467; b2 −4.420764 |
| CH4 | a1 −27.738399; t2×3 −20.487453; a1 −5.301601; t2×3 −3.782547 |
| BH3 | a1-prime −21.661994; e-prime×2 −18.700003; a2-double-prime −8.300000; a1-prime −5.998006; e-prime×2 −3.209997 |
| SF6 model and lecture (two executions) | a1g −43.433587; t1u×3 −28.372629; eg×2 −14.850509; a1g −13.124378; t1u×3 +0.635623 |
| H2 | Ag −20.010800; B1u −7.209200 (D2h subgroup) |
| H3+ and neutral H3 (two executions) | a1-prime −26.411600; e-prime×2 −7.209200; neutral e-prime¹ warns |
| Allyl cation π | b1 −15.113474; a2 −10.660000; b1 −6.206526 |
| Six H Oh | a1g −28.011800; t1u×3 −12.009800 (4 electrons, warns); eg×2 −8.809400 |
| Four H Td | a1 −19.704354; t2×3 −11.578549 (2 electrons, warns) |
| B2 | σg −23.714683; σg −12.620840; σu −12.097903; πu×2 −10.741438; πg×2 −5.858562; σu +3.733425 |
| C2 | σg −35.579486; σg −17.412870; σu −16.327845; πu×2 −14.658041 (2 electrons, warns); πg×2 −6.661959; σu +9.140201 |
| N2 | σg −45.683529; σg −22.399546; σu −21.419414; πu×2 −18.302390; πg×2 −8.057610; σu +12.022489 |
| O2 | σg −46.458232; σu −28.219157; σg −26.021038; πu×2 −20.083170; πg×2 −11.616830 (2 electrons, warns); σu +4.238427 |
| F2 | σg −48.374432; σu −36.456980; σg −28.181952; πu×2 −21.746222; πg×2 −15.553778; σu −4.626636 |
| Al2 | σg −13.645089; σu −10.334969; σg −8.501361; πu×2 −6.826040; πg×2 −5.133960; σu −2.118580 |
| Si2 | σg −19.190314; σu −14.443620; σg −11.488644; πu×2 −9.003547 (2 electrons, warns); πg×2 −6.556453; σu −2.217422 |
| P2 | σg −23.894620; σu −16.942957; σg −14.457907; πu×2 −11.371691; πg×2 −7.928309; σu −1.684516 |
| S2 | σg −27.433710; σu −20.686488; σg −16.802717; πu×2 −13.349355; πg×2 −9.890645 (2 electrons, warns); σu −3.737084 |
| Cl2 | σg −29.274860; σu −23.322165; σg −18.572273; πu×2 −15.231892; πg×2 −12.108108; σu −6.630702 |
| CO | σ −45.350342; σ −26.597061; π×2 −18.754865; σ −17.057093 (HOMO, C weight 0.977); π×2 −7.755135; σ +10.684496 |
| NH3 | a1 −32.346972; e×2 −22.708244; a1 −15.354182 (HOMO); a1 −4.648846 (LUMO); e×2 −4.081756 |
| H2S | a1 −25.263728; b2 −17.519005; a1 −15.520417; b1 −11.620000 (HOMO, S 3px); b2 −7.710995 (LUMO); a1 −7.155856 |
| SiH4 | a1 −20.705332; t2×3 −16.002506 (HOMO); a1 −8.794668 (LUMO); t2×3 −5.387494 |
| NH model | σ −28.738041; σ −18.928837; π×2 −13.180000 (2 electrons, warns); σ −4.683121 |
| NH lecture | σ −25.560000 (N 2s, nonbonding); σ −21.232911; π×2 −13.180000 (2 electrons, warns); σ −5.557089 |
| OH model | σ −35.060791; σ −22.045407; π×2 −15.850000 (3 electrons, warns); σ −4.733802 |
| OH lecture | σ −32.380000 (O 2s, nonbonding); σ −23.746253; π×2 −15.850000 (3 electrons, warns); σ −5.713747 |
| HF model | σ −42.516934; σ −25.152810; π×2 −18.650000 (HOMO); σ −4.760256 (LUMO) |
| HF lecture | σ −40.170000 (F 2s, nonbonding); σ −26.451006; π×2 −18.650000 (HOMO); σ −5.808994 (LUMO) |

All 32 current fixtures satisfy file-driven analytic assembly, symmetry/projector labels and residual/orthonormality checks. Still open: real CrystOD comparison (deferred, B1). Browser gates (C4) closed. Later advanced d blocks, exact slide aliasing and spectral fits remain unclaimed.
