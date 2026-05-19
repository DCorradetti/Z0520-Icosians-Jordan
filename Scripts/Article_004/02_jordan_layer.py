"""
02_jordan_layer.py — Jordan-layer certificates for Article_004.

Generates:
  - jordan_peirce_failure_count.json  (the 540 = 96 + 60 + 384 count)
  - jordan_doubled_closure.json       (doubled product closes; 0 failures)
  - g0_gram_half_integer.json         (Lemma 4.4: 20 of 64)

Writes versioned JSON certificates under
  3. Certificates/Article_004/
and an immutable run snapshot under
  7. Results/Article_004/run_jordan_layer_<timestamp>/.

This script is intentionally compact — full pipeline scope (norm
integrality, LDL pivots, principal-lattice criterion) follows the
same pattern and is left to forthcoming work (see manuscript §20).
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import platform
import sys

import g0_j3
from Z_phi import Zphi


def _here() -> str:
    return os.path.dirname(os.path.abspath(__file__))


WORKSPACE_ROOT = os.path.abspath(os.path.join(_here(), "..", ".."))
CERT_DIR = os.path.join(WORKSPACE_ROOT, "3. Certificates", "Article_004")
RESULTS_DIR = os.path.join(WORKSPACE_ROOT, "7. Results", "Article_004")


def _sha256_of_file(p: str) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_certificate(name: str, payload: dict) -> str:
    os.makedirs(CERT_DIR, exist_ok=True)
    path = os.path.join(CERT_DIR, name)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
    return path


def _timestamp() -> str:
    return _dt.datetime.utcnow().strftime("%Y%m%d_%H%M%SZ")


def _run_folder(label: str) -> str:
    folder = os.path.join(RESULTS_DIR, f"run_{label}_{_timestamp()}")
    os.makedirs(folder, exist_ok=True)
    return folder


def cert_g0_gram(run_dir: str) -> str:
    half_int = g0_j3.g0_gram_half_integer_count()
    payload = {
        "certificate": "g0_gram_half_integer",
        "version": "v1",
        "conclusion": (
            "On the ordered basis (1,i,h,g,ell,i*ell,h*ell,g*ell) of "
            "G_0, the 8x8 Gram matrix of B(x,y) = Re(x conj(y)) has "
            "exactly 20 of its 64 entries in (1/2)O_K \\setminus O_K."
        ),
        "expected": 20,
        "observed": half_int,
        "passed": (half_int == 20),
        "reference": "Article_004 Lemma 4.4",
    }
    name = "g0_gram_half_integer_v1.json"
    path = _write_certificate(name, payload)
    return path


def cert_jordan_failure_count(run_dir: str) -> str:
    DD, DO_total, OO_same, OO_diff = g0_j3.jordan_failure_count()
    total = DD + DO_total + OO_same + OO_diff
    payload = {
        "certificate": "jordan_peirce_failure_count",
        "version": "v1",
        "conclusion": (
            "Of the 27x27 = 729 ordered basis pairs of J_3(G_0), the "
            "ordinary Albert product X o Y = (XY + YX)/2 has 540 "
            "integrality failures, partitioned as 96 + 60 + 384 "
            "(Peirce-incident D-O, same-slot O-O, different-slot O-O)."
        ),
        "expected": {
            "D_D": 0, "D_O_or_O_D": 96, "O_O_same": 60, "O_O_diff": 384,
            "total": 540,
        },
        "observed": {
            "D_D": DD, "D_O_or_O_D": DO_total,
            "O_O_same": OO_same, "O_O_diff": OO_diff,
            "total": total,
        },
        "passed": (DD == 0 and DO_total == 96 and OO_same == 60
                   and OO_diff == 384),
        "reference": "Article_004 Proposition 7.1",
    }
    name = "jordan_peirce_failure_count_v1.json"
    path = _write_certificate(name, payload)
    return path


def cert_doubled_closure(run_dir: str) -> str:
    failures = g0_j3.doubled_failure_count()
    payload = {
        "certificate": "jordan_doubled_closure",
        "version": "v1",
        "conclusion": (
            "Of the 729 ordered basis pairs, the doubled product "
            "{X, Y} = XY + YX has 0 integrality failures."
        ),
        "expected_failures": 0,
        "observed_failures": failures,
        "passed": (failures == 0),
        "reference": "Article_004 Theorem 8.2",
    }
    name = "jordan_doubled_closure_v1.json"
    path = _write_certificate(name, payload)
    return path


def main() -> int:
    run_dir = _run_folder("jordan_layer")
    print(f"Run folder: {run_dir}")
    cert_paths = []
    print("[1/3] G_0 Gram half-integer count ...")
    cert_paths.append(cert_g0_gram(run_dir))
    print("[2/3] Jordan Peirce failure count (729 pairs) ...")
    cert_paths.append(cert_jordan_failure_count(run_dir))
    print("[3/3] Doubled-product closure (729 pairs) ...")
    cert_paths.append(cert_doubled_closure(run_dir))

    # Collect SHA-256 hashes
    manifest = {
        "run": "jordan_layer_v1",
        "timestamp_utc": _timestamp(),
        "python": sys.version,
        "platform": platform.platform(),
        "certificates": [
            {
                "path": os.path.relpath(p, WORKSPACE_ROOT).replace("\\", "/"),
                "sha256": _sha256_of_file(p),
                "name": os.path.basename(p),
            }
            for p in cert_paths
        ],
    }
    out = os.path.join(run_dir, "manifest.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    print(f"Run manifest: {out}")

    # Print pass/fail summary
    for p in cert_paths:
        with open(p, encoding="utf-8") as f:
            doc = json.load(f)
        flag = "PASS" if doc.get("passed") else "FAIL"
        print(f"  {flag}: {doc['certificate']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
