# Z0520 — Icosians Jordan Algebras

Reproducibility scripts and certificates for the project
*Golden Albert Algebras and the Integral Tits–Kantor–Koecher
Envelope of Type E₇ over Q(√5)*.

## Contents

- [`Scripts/`](Scripts/) — Python reproducibility package
  (Python 3.11, standard library only; exact arithmetic in
  O_K = Z[φ], no floating point).
- [`Certificates/`](Certificates/) — JSON certificates produced
  by the pipeline scripts, with SHA-256 hashes consolidated in
  `master_manifest.json`.

## Articles

- **Article_004** — Jordan layer: G_0 Gram half-integer count
  (Lemma 4.4), Peirce-block partition (Proposition 7.1), and
  doubled-product closure (Theorem 8.2). See
  [`Scripts/Article_004/README.md`](Scripts/Article_004/README.md)
  for module-level documentation and run instructions.

## How to reproduce

```bash
cd Scripts/Article_004
python 02_jordan_layer.py
python 99_master_manifest.py
```

Certificates are written to `Certificates/Article_004/*.json`.

## License

Released under the MIT licence; see the manuscript for the
mathematical context.
