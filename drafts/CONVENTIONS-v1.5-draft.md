# CONVENTIONS v1.5 — draft for owner sign-off

**Status: draft, not applied.** `CONVENTIONS.md` v1.4 stays frozen until the owner signs this off. Prepared 2026-11-07 by Claude. Nothing in `params/`, `groups/`, `checks/` or the engine has changed. The numbers below come from a scratch prototype: the unchanged engine, run on a copy of `params/` and `groups/` with the additions listed here.

## One-page summary: what you are asked to decide

v1.5 adds the molecules Release 2 still needs: NH3 (for the BH3 → NH3 morph), H2S and SiH4 (for the H2O vs H2S and CH4 vs SiH4 comparisons). It also corrects one erratum in the signed text. Optionally it adds the diatomic hydrides NH, OH and HF, which your exam's Question 2(a) asks about.

| # | Decision | Recommendation |
|---|---|---|
| V1 | **Erratum in signed v1.4.** The §9.1 table still shows the old Harrison composite values for B, C, N, O and S, and omits Al, Si, P and Cl. The prose above it, `params/onsite.csv`, Appendix B and every result use Miessler. No computed number is affected (the checker reads the CSV), but the frozen document contradicts itself. | Replace the table with the Miessler values actually in use (below). |
| V2 | **Add NH3, H2S and SiH4** as fixtures, with CCCBDB geometries and the C3v group. | Yes: the spec's remaining Release 2 items need them. |
| V3 | **Si–H coupling.** No Si–H row exists. | Add Si–H ss and sp rows with the same Harrison η and hydrogen factor 0.6 as every other X–H row (flagged as extrapolated, like all H rows). |
| V4 | **Diatomic hydrides NH, OH, HF.** | Optional. Yes if you want your exam's NH question in the app. |
| V5 | **Lecture presets for the new hydrides.** In the model, the central atom's s mixes with H 1s, so the lowest σ (or a1) level counts as bonding. The counting bond order is then 2 for NH, OH and HF (textbook 1), 4/3 for NH3 and 3/2 for H2S (textbook 1). Model water already gives 3/2 for the same reason. Switching the central s–H coupling off, as lecture water does, gives 1 in all five. | Give NH3, H2S, NH, OH and HF a Model/Lecture pair, like water. Not SiH4: switching Si 3s off breaks it (bond order 3/4, wrong HOMO), so SiH4 stays model-only like CH4. |

The hydrogen factor stays 0.6. Every new ordering holds for any factor between 0.30 and 1.20, except NH3's two empty levels (a1* and e*), which swap at about 0.80, outside the signed window 0.4235–0.6668. So the new molecules neither move nor tighten the window.

The exam's hint for Q2(a) (H 1s "slightly higher" than N 2p) still disagrees with Miessler (H 1s −13.61, below N 2p −13.18). As decided in M8, the app follows Miessler. NH's answer is the same either way: two electrons in a nonbonding π pair, a triplet.

## Proposed changes, section by section

### Header and §12

Version 1.5, not frozen until owner sign-off B7. Change-log row: "V1 erratum (§9.1 table); C3v; NH3, H2S, SiH4 [and NH, OH, HF]; Si–H rows; lecture presets (V5); Appendix A.3 and B rows; checker patch K21."

### §2 Geometries (V2, V4)

| Case | Positions and atom order |
|---|---|
| NH3 | N=1 at origin, C3 along z, lone pair towards +z; H2–4 at azimuth 0°, 120°, 240° from +x, each at angle β from −z, with cos(HNH) = cos²β − ½ sin²β (β = 67.8549° for 106.67°). H2 lies in the xz mirror plane σv(0). |
| H2S | As H2O: S=1 at origin, yz plane, z bisector; H2 on +y, H3 on −y |
| SiH4 | As CH4 |
| NH, OH, HF | Heavy atom 1 at origin, H 2 on +z (C2v finite subgroup of C∞v, as for CO) |

| Fixture | Pinned geometry | Source |
|---|---|---|
| NH3 | N–H 1.0124 Å; HNH 106.67° | NIST CCCBDB (reference 1966Herzberg); `evidence/cccbdb_NH3.html` |
| H2S | S–H 1.3356 Å; HSH 92.11° | NIST CCCBDB (reference 1975Coo/DeL); `evidence/cccbdb_H2S.html` |
| SiH4 | Si–H 1.4799 Å | NIST CCCBDB (reference 1955Boyd; listed 1.480); `evidence/cccbdb_SiH4.html` |
| NH | N–H 1.0362 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_NH.html` |
| OH | O–H 0.9697 Å | NIST CCCBDB r_e (Huber & Herzberg 1979); `evidence/cccbdb_OH.html` |
| HF | H–F 0.9168 Å | NIST CCCBDB r_e (NIST diatomic); `evidence/cccbdb_HF.html` |

Four-decimal values are the Cartesian distances CCCBDB lists, as for the diatomics.

### §5 Group data (V2)

Add C3v to the generated groups and to the "Later crystallographic" coverage row. Generator: MolSym 1.3.1, as for the others. Bilbao reference: point group 3m (`sg=156&num=19`), snapshot `evidence/bilbao_C3v.html` (sha256 starts `70bf4bd94c5324ce`). Bilbao's mirror column m₁₋₁₀ is in hexagonal axes: direction [1, −1, 0] is Cartesian (3/2, −√3/2, 0), so the target is the mirror with that normal, which lies in the class of σv(0) (the xz plane). Prototype result: `PASS_GEOMETRY_MAPPED`, maximum character error 0, closure 1.1e−16. All ten existing group files regenerate byte-identical. Checker patch K21 adds the C3v job and target.

### §7 Electron counting

NH3 8, H2S 8, SiH4 8; NH 6, OH 7, HF 8 (neutral-atom counting). NH's π² and OH's π³ fire the partly-filled warnings, as O2's πg² does.

### §9.1 On-sites (V1 erratum)

Replace the table with the values in `params/onsite.csv` (unchanged by v1.5):

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

### §9.2 Pairs (V3)

Add `H,Si,ss_sigma,-1.4,…,H,0.6,harrison_extrapolated_H` and `H,Si,sp_sigma,1.84,…,H,0.6,harrison_extrapolated_H`, with the same source and override-reason fields as the other X–H rows. Neighbour lists: NH3, H2S, SiH4 central–H on, H–H off, as for H2O, CH4 and BH3 (B5).

### §9.3 Presets (V5)

`NH3_model`, `H2S_model`, `SiH4_model`, `NH_model`, `OH_model`, `HF_model`: geometry as §2, hydrogen scale 0.6. If V5 is accepted, also `NH3_lecture`, `H2S_lecture`, `NH_lecture`, `OH_lecture`, `HF_lecture`: the central s–H 1s coupling forced to zero, as in `H2O_lecture`.

### Appendix A.3 (new blocks)

NH3, with β the angle of each N–H bond from −z: a1 is A(ε_Ns, ε_Np, ε_H; √3 v_ss, √3 v_sp cos β); e is D(ε_Np, ε_H; −√(3/2) v_sp sin β), two copies. No a2 functions in this basis. Checked against the engine: largest difference 7e−15 eV. H2S has the H2O blocks with S on-sites and θ = 92.11°. SiH4 has the CH4 blocks with Si on-sites.

### Appendix B (new rows, model presets, hydrogen factor 0.6)

| Fixture | Levels, lowest first |
|---|---|
| NH3 | a1 −32.346972; e×2 −22.708244; a1 −15.354182 (HOMO); a1 −4.648846 (LUMO); e×2 −4.081756 |
| H2S | a1 −25.263728; b2 −17.519005; a1 −15.520417; b1 −11.620000 (HOMO, S 3px); b2 −7.710995 (LUMO); a1 −7.155856 |
| SiH4 | a1 −20.705332; t2×3 −16.002506 (HOMO); a1 −8.794668 (LUMO); t2×3 −5.387494 |
| NH | σ −28.738041; σ −18.928837; π×2 −13.180000 (2 electrons, warns); σ −4.683121 |
| OH | σ −35.060791; σ −22.045407; π×2 −15.850000 (3 electrons, warns); σ −4.733802 |
| HF | σ −42.516934; σ −25.152810; π×2 −18.650000 (HOMO); σ −4.760256 (LUMO) |

Lecture rows, if V5 is accepted, are generated when applied.

### A.5 bond-order oracles

Model counting: NH3 4/3, H2S 3/2, SiH4 1, NH 2, OH 2, HF 2. With V5, lecture counting: NH3, H2S, NH, OH, HF 1.

### Fidelity references (no new claims)

- H2S's occupied order (2b2 below 3a1 below 1b1) agrees with its photoelectron band order as recalled, not yet sourced: bands near 10.5, 13.4 and 15.5 eV. If a source is wanted, it would be added as a curated check, like CH4's.
- NH and OH fill a nonbonding π pair partly, matching their ³Σ⁻ and ²Π ground states, and HF is closed-shell.

None is used to tune anything. No orbital-ordering reference is sourced for the empty levels.

## What happens on sign-off

1. Apply the text above to `CONVENTIONS.md` as v1.5.
2. Generate `groups/C3v.json` (K21) and extend `params/` and the checks.
3. Add the cases to the engine.
4. Build the deferred Release 2 items: the BH3 → NH3 morph, H2O vs H2S and CH4 vs SiH4 comparisons, and the hydrides in Practice and Exam.
5. Run `make check` and the browser tests.
6. Freeze v1.5 only after you sign off the applied text, as with v1.4.
