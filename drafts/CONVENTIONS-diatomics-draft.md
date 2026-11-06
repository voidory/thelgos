> **Applied 2026-11-06** to `CONVENTIONS.md` v1.3 under the owner's standing permission. This draft is kept as the review record of what changed and why.

# Draft for the owner: diatomic additions to CONVENTIONS.md (gate C2)

**Status:** draft prepared by Claude for the owner to accept, edit or reject. Only the owner edits `CONVENTIONS.md`; nothing here is in force until copied there. The checker already runs these fixtures (`CHECKER-CHANGES.md` K14).

## 1. Proposed rows for the §2 fixture table

| Fixture | Pinned geometry | Source / qualification |
|---|---|---|
| B2 | B–B 1.5900 Å | NIST CCCBDB, r_e; primary reference Gurvich et al. (1989). Snapshot `evidence/cccbdb_B2.html` |
| C2 | C–C 1.2425 Å | NIST CCCBDB, r_e; Huber & Herzberg (1979). `evidence/cccbdb_C2.html` |
| N2 | N–N 1.0977 Å | NIST CCCBDB, r_e; Huber & Herzberg (1979). `evidence/cccbdb_N2.html`. Included by owner decision B4 |
| O2 | O–O 1.2075 Å | NIST CCCBDB, r_e; Huber & Herzberg (1979). `evidence/cccbdb_O2.html` |
| F2 | F–F 1.4119 Å | NIST CCCBDB, r_e; Huber & Herzberg (1979). `evidence/cccbdb_F2.html` |
| CO | C–O 1.1282 Å | NIST CCCBDB, r_e; NIST diatomic data. `evidence/cccbdb_CO.html`. C is atom 1 (−R/2), O atom 2 (+R/2) |

The CCCBDB pages list each bond length rounded to three decimals (for example `rOO 1.208`) and give the Cartesian distance to four (`1.2075`). The four-decimal value is pinned.

## 2. Proposed text replacing Appendix A.6

### A.6 Diatomic blocks

Bond along z; atom 1 at −R/2, atom 2 at +R/2; basis s, px, py, pz on each atom; all channels from `pairs.csv` at distance R; no H, so no hydrogen scale. With V written for the row-law values at R:

**Homonuclear (D2h finite subgroup; x = B3u, y = B2u, z = B1u):**

| D∞h | D2h | Block | Copies |
|---|---|---|---:|
| σg | Ag | D(ε_s + V_ssσ, ε_p − V_ppσ; −V_spσ) | 1 |
| σu | B1u | D(ε_s − V_ssσ, ε_p + V_ppσ; +V_spσ) | 1 |
| πu | B2u + B3u | ε_p + V_ppπ | 2 |
| πg | B2g + B3g | ε_p − V_ppπ | 2 |

SALCs: σg from (s1+s2)/√2 and (pz1−pz2)/√2; σu from (s1−s2)/√2 and (pz1+pz2)/√2; πu from (px1+px2)/√2, likewise y; πg from the differences. Signs follow §1; eigenvalues do not depend on them.

**CO (C2v finite subgroup, C∞v):** σ (A1) is the 4×4 block in basis (s_C, pz_C, s_O, pz_O):

```
[[ε_sC, 0,     V_ssσ,  V_spσ],
 [0,    ε_pC, −V_spσ,  V_ppσ],
 [V_ssσ, −V_spσ, ε_sO, 0    ],
 [V_spσ,  V_ppσ, 0,    ε_pO ]]
```

π (B1 + B2) is D(ε_pC, ε_pO; V_ppπ), two copies.

**Labelling.** In the finite subgroup a π pair spans two irreps, so degenerate sets are labelled by their exact irrep decomposition (for example `B2u+B3u`) and mapped to D∞h/C∞v names by a fixed lookup.

**Required checks (unchanged in substance from v1.2, plus N2):** reference orderings B2, C2, N2 πu below σg(2p); O2, F2 σg(2p) below πu, where σg(2p) is taken as the higher of the two σg levels. A disagreement is recorded as an Orbital-ordering fidelity warning, never tuned away. O2 must end in a partly filled πg with occupation and spin warnings. CO: counting bond order 3, σ HOMO with carbon weight above 0.5.

## 3. Proposed Appendix B rows

Energies in eV, lowest first; model calculations, not measurements. Electrons: B2 6, C2 8, N2 10, O2 12, F2 14, CO 10 (neutral-atom valence counting).

| Fixture | Levels, lowest first |
|---|---|
| B2 | σg −22.131569; σg −11.033954; σu −10.566590; πu×2 −9.081438; πg×2 −4.198562; σu +5.372113 |
| C2 | σg −33.782041; σg −15.610315; σu −14.450984; πu×2 −12.968041 (2 electrons, warns); πg×2 −4.971959; σu +10.863339 |
| N2 | σg −43.562631; σg −20.290444; σu −19.015983; πu×2 −16.592390; πg×2 −6.347610; σu +13.849058 |
| O2 | σg −43.746764; σu −25.133433; σg −23.772506; πu×2 −18.363170; πg×2 −9.896830 (2 electrons, warns); σu +6.112703 |
| F2 | σg −48.374432; σu −36.456980; σg −28.181952; πu×2 −21.746222; πg×2 −15.553778; σu −4.626636 |
| CO | σ −42.990247; σ −24.079019; π×2 −17.042804; σ −15.182634 (HOMO, C weight 0.977); π×2 −6.057196; σ +12.491900 |

I checked the F2 and O2 rows by hand from the blocks above. Counting bond orders: B2 1, C2 2, N2 3, O2 2, F2 1, CO 3. The N2 row matches a hand evaluation of the same blocks.

## 4. Results and the decision taken

Decided under delegation as M1 in `DECISIONS.md`: keep the parameters and show the red flags. The table and options below are kept for review.


| Molecule | Orbital ordering | Ground state | Model |
|---|---|---|---|
| B2 | **Disagrees**: σg(2p) below πu | **Disagrees**: reference ³Σg, model closed shell (σg² σg² σu²) | Predicts diamagnetic B2 |
| C2 | **Disagrees**: σg(2p) below πu | **Disagrees**: reference ¹Σg⁺, model partly fills πu | Predicts an open-shell C2 |
| N2 | **Disagrees**: σg(2p) below πu | Agrees: closed shell, ¹Σg⁺ | Bond order 3 |
| O2 | Agrees | Agrees: ³Σg, πg² open shell | — |
| F2 | Agrees | Agrees: closed shell | — |
| CO | n/a | Agrees: closed shell; HOMO carbon σ | Bond order 3 |

Ground states are as listed on the same CCCBDB pages.

The current law includes s–p mixing, but it is not strong enough to lift σg(2p) above πu for B2 and C2. Spec §7 says "B2 red without mixing, green with the R2 switch" and §6 plans the R2 switch as its demonstration. With these parameters that demonstration cannot work: B2 stays red with mixing on. Options for the owner:

1. **Accept and rewrite the spec (my recommendation for Release 1):** B2 and C2 show red Orbital-ordering and Spin-state chips with mixing on, which is an honest limit of the model; spec §6/§7 change to say so.
2. **Revisit the parameter scheme, as a stated scheme choice and not a fit to these two molecules:** for example the alternative η set in `evidence/harrison.txt` or a different on-site source for B and C. Any change moves every Appendix B row and the admissible H-scale window, so it needs the full regression rerun.

Not to do: tune B and C parameters until B2 and C2 pass. CONVENTIONS A.6 and the v1.2 TODO both rule this out.

## 5. Proposed §9.1 and §9.2 text for nitrogen and the coefficient set

- §9.1 on-site table: add `N (n=2) | −23.04 | −11.47`, from the same Harrison transcription (`evidence/harrison.txt`, line 87). Replace "nitrogen is not added without the B4 scope decision" with "nitrogen added by owner decision B4".
- §9.2: after the coefficient list, add: "These are the coefficients of Harrison's 1980 book (Dover reprint), as transcribed in `evidence/harrison.txt`. The same file lists a second set attributed to Shi & Papaconstantopoulos (2004); it is not used. See DECISIONS.md M2."
- §2 and §11: replace "N2 remains excluded pending B4" and "N2 only if the owner includes it" with N2 as a Release 1 fixture.
- §11: record that the CrystOD gate C1 is deferred by owner decision B1 and runs as `make check-release`.
