# Article_004 — Reproducibility Package

This directory contains the Python reproducibility package for
*Golden Albert Algebras and the Integral Tits–Kantor–Koecher
Envelope of Type E₇ over Q(√5)* (Article_004).

The implementation uses Python 3.11 and the standard library only.
All arithmetic in O_K = Z[φ] is exact (no floating point).

## Modules

- `Z_phi.py` — Exact arithmetic in Z[φ] (and Q[φ]), with Galois
  conjugation, numerical evaluation at the two real embeddings
  σ_+ and σ_-, and integrality test.
- `quaternion.py` — Quaternion algebra H_K = (−1,−1)/K and the
  icosian ring I ⊂ H_K with ordered basis {1, i, h, g} (Tits
  1980 convention).
- `g0_j3.py` — The Cayley–Dickson double G_0 = I + I·ℓ, the
  hermitian module J_3(G_0), the ordinary and doubled Albert
  products, and the integrality test against the G_0-basis
  (icosian/ℓ-shifted icosian).

## Pipeline scripts

- `02_jordan_layer.py` — Jordan-layer certificates:
  G_0 Gram half-integer count (Lemma 4.4: 20 of 64), the
  Peirce-block partition 540 = 96 + 60 + 384 (Proposition 7.1),
  and the doubled-product closure (Theorem 8.2: 0 failures
  over the 729 ordered basis pairs).
- `99_master_manifest.py` — Consolidates all certificate JSONs
  into a single manifest with SHA-256 hashes.

## How to run

From the workspace root,

```
cd "2. Scripts/Article_004"
python 02_jordan_layer.py
python 99_master_manifest.py
```

The first script writes:
- Certificates to `3. Certificates/Article_004/*.json`.
- An immutable run snapshot to
  `7. Results/Article_004/run_jordan_layer_<timestamp>/` with a
  per-run `manifest.json` recording SHA-256 hashes and metadata.

The second script writes
`3. Certificates/Article_004/master_manifest.json`, listing every
JSON certificate with its SHA-256 hash.

## Smoke tests

Each module exposes a `smoke_test()` callable as
`python <module>.py`.

## Architecture

The pipeline reflects the manuscript section structure:

- Setup layer (in `Z_phi`, `quaternion`, `g0_j3`).
- Jordan layer (`02_jordan_layer.py` — covers Lemma 4.4,
  Proposition 7.1, Theorem 8.2; remaining Jordan-layer
  certificates from Article §10 are routine continuations of
  the same pattern and are listed in Article §20).
- Lie layer (forthcoming `03_lie_layer.py` — see Article §20
  for the planned certificate set).
- Identification layer (forthcoming `04_identifications.py`).

Three of the most load-bearing Jordan-layer certificates are
included in v1 of the reproducibility package; the remaining
certificates listed in the manuscript Table `tab:certs-summary`
follow the same coding pattern and will be added before the
Zenodo archive is published.

## Versioning

The certificate names carry a `_v1` suffix; an immutable per-run
folder under `7. Results/` records the exact timestamped output.
The master manifest collects SHA-256 hashes of every certificate.

## License

Released under the MIT licence; see the manuscript for the
mathematical context.
