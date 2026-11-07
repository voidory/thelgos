# Molecule packs

A molecule pack adds a molecule to the app as data, with no engine code (spec section 14). Each pack is `packs/<id>.json`, validated against `pack.schema.json`, and listed in `index.json`. The browser loads the packs listed there.

A pack declares:
- **identity:** name, formula, charge, SMILES, aliases;
- **template and shape:** the template (T0–T8) and one of the shapes in the schema;
- **geometry:** bond lengths and angles from a cited source (owner decision R3-3: no estimated geometries); keep the source page in `evidence/`;
- **basis:** central s and p; terminal H 1s, or one inward σ p for other terminals;
- **electrons:** the count, by neutral-atom counting;
- **point group:** declared; a group missing from `groups/` falls back to C1 (F2.14);
- **parameters:** always "scheme", meaning Miessler on-sites and Harrison couplings with hydrogen factor 0.6 (R3-1). An element, orbital or pair with no parameter makes the gate fail, and the gate names it;
- **notes:** optional reading notes shown in the app.

The engine builds the molecule in its template's CONVENTIONS §2 orientation, so the stored group files apply.

To add one, write the file, then run:

```bash
make add-molecule PACK=packs/<id>.json PYTHON=.venv/bin/python
```

If every check passes, add the id to `index.json`. `make check` then runs the same gate on every listed pack (F2.15).
