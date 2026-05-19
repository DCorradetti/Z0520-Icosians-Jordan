"""
99_master_manifest.py — Consolidate the certificate package for
Article_004 into a single manifest with SHA-256 hashes.

Iterates over every JSON file in
  3. Certificates/Article_004/
records (name, sha256, size, modification time, principal
conclusion / observed result if present), and writes
  3. Certificates/Article_004/master_manifest.json.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import platform
import sys


def _here() -> str:
    return os.path.dirname(os.path.abspath(__file__))


WORKSPACE_ROOT = os.path.abspath(os.path.join(_here(), "..", ".."))
CERT_DIR = os.path.join(WORKSPACE_ROOT, "3. Certificates", "Article_004")


def _sha256_of_file(p: str) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if not os.path.isdir(CERT_DIR):
        print(f"No certificate directory at {CERT_DIR}", file=sys.stderr)
        return 1
    entries = []
    for name in sorted(os.listdir(CERT_DIR)):
        if not name.endswith(".json"):
            continue
        if name == "master_manifest.json":
            continue
        path = os.path.join(CERT_DIR, name)
        with open(path, encoding="utf-8") as f:
            doc = json.load(f)
        entries.append({
            "name": name,
            "sha256": _sha256_of_file(path),
            "size_bytes": os.path.getsize(path),
            "modified_utc": _dt.datetime.utcfromtimestamp(
                os.path.getmtime(path)).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "certificate": doc.get("certificate"),
            "passed": doc.get("passed"),
            "reference": doc.get("reference"),
        })
    manifest = {
        "manifest": "master_manifest",
        "version": "v1",
        "timestamp_utc": _dt.datetime.utcnow().strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "python": sys.version,
        "platform": platform.platform(),
        "certificate_dir": os.path.relpath(
            CERT_DIR, WORKSPACE_ROOT).replace("\\", "/"),
        "total_certificates": len(entries),
        "all_passed": all(e.get("passed") for e in entries),
        "certificates": entries,
    }
    out = os.path.join(CERT_DIR, "master_manifest.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    print(f"Master manifest: {out}")
    print(f"  Total certificates: {manifest['total_certificates']}")
    print(f"  All passed: {manifest['all_passed']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
