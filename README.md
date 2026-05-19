# Icosians TKK-Construction
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3110/)
[![stdlib only](https://img.shields.io/badge/dependencies-stdlib%20only-success)](https://docs.python.org/3.11/library/)
[![Exact arithmetic](https://img.shields.io/badge/arithmetic-exact%20%E2%84%A4%5B%CF%86%5D-orange)](https://en.wikipedia.org/wiki/Golden_ratio)
[![Last commit](https://img.shields.io/github/last-commit/DCorradetti/Z0520-Icosians-Jordan)](https://github.com/DCorradetti/Z0520-Icosians-Jordan/commits/main)
[![Repo size](https://img.shields.io/github/repo-size/DCorradetti/Z0520-Icosians-Jordan)](https://github.com/DCorradetti/Z0520-Icosians-Jordan)

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

## Related projects

- [**Non-crystallographic-Integers**](https://github.com/DCorradetti/Non-crystallographic-Integers)
  — Companion supplementary archive on non-crystallographic systems
  of integers over composition algebras (arXiv:2605.15075). Same
  exact-arithmetic, stdlib-only verification philosophy and the
  same icosian / Z[φ] integral backbone used here.

## License

Released under the MIT licence; see the manuscript for the
mathematical context.
