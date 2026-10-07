#!/usr/bin/env python3
"""Fetch public GEO series matrices and verify the frozen input byte hashes."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=root / "csv")
    parser.add_argument("--manifest", type=Path, default=root / "environment/geo_inputs.tsv")
    parser.add_argument("--reference-dir", type=Path,
                        help="Also compare with existing local matrices; never modifies them.")
    parser.add_argument("--verify-only", action="store_true", help="Check local inputs without network access.")
    parser.add_argument("--force-download", action="store_true", help="Download again even if local input passes.")
    parser.add_argument("--timeout", type=float, default=90)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--receipt", type=Path, help="JSON verification record; defaults to DATA_DIR/geo_fetch_receipt.json.")
    args = parser.parse_args()
    if args.verify_only and args.force_download:
        parser.error("--verify-only and --force-download are mutually exclusive")
    if args.retries < 1 or args.timeout <= 0:
        parser.error("--retries and --timeout must be positive")
    args.data_dir.mkdir(parents=True, exist_ok=True)
    entries = list(csv.DictReader(args.manifest.open(encoding="utf-8"), delimiter="\t"))
    records = []
    for entry in entries:
        name = entry["filename"]
        if Path(name).name != name:
            raise ValueError("Manifest filenames must be basenames")
        target = args.data_dir / name
        record = dict(entry)
        record["expected_bytes"] = int(entry["bytes"])
        try:
            current = sha256(target) if target.is_file() else None
            if args.verify_only:
                record["source"] = "local verification only"
            elif current == entry["sha256"] and not args.force_download:
                record["source"] = "previously verified local file"
            else:
                if not entry["url"].startswith("https://ftp.ncbi.nlm.nih.gov/"):
                    raise ValueError("Expected an HTTPS NCBI GEO download URL")
                with tempfile.TemporaryDirectory(prefix="geo_fetch_", dir=args.data_dir) as temp:
                    archive = Path(temp) / (name + ".gz")
                    unpacked = Path(temp) / name
                    for attempt in range(args.retries):
                        try:
                            request = urllib.request.Request(entry["url"], headers={"User-Agent": "CAVD-GSEA-reproducibility/1.0"})
                            with urllib.request.urlopen(request, timeout=args.timeout) as response, archive.open("wb") as output:
                                shutil.copyfileobj(response, output)
                            break
                        except (OSError, urllib.error.URLError):
                            if attempt + 1 == args.retries:
                                raise
                            time.sleep(min(2 ** attempt, 4))
                    record["downloaded_gzip_sha256"] = sha256(archive)
                    with gzip.open(archive, "rb") as compressed, unpacked.open("wb") as output:
                        shutil.copyfileobj(compressed, output)
                    observed = sha256(unpacked)
                    record["downloaded_uncompressed_sha256"] = observed
                    record["downloaded_uncompressed_bytes"] = unpacked.stat().st_size
                    if observed != entry["sha256"] or unpacked.stat().st_size != int(entry["bytes"]):
                        raise ValueError("Downloaded bytes differ from frozen input; existing input was not replaced")
                    unpacked.replace(target)
                    current = observed
                record["source"] = "fresh NCBI HTTPS download, gzip decompressed"
            record["observed_sha256"] = current
            record["observed_bytes"] = target.stat().st_size if target.exists() else None
            if current != entry["sha256"] or record["observed_bytes"] != int(entry["bytes"]):
                raise ValueError("Local input does not match frozen-input manifest")
            if args.reference_dir:
                reference = args.reference_dir / name
                record["reference_sha256"] = sha256(reference)
                record["reference_bytes"] = reference.stat().st_size
                if record["reference_sha256"] != current or record["reference_bytes"] != record["observed_bytes"]:
                    raise ValueError("Local reference and downloaded/input file differ")
            record["status"] = "pass"
            print(f"PASS {name} {current}")
        except (OSError, ValueError, urllib.error.URLError) as error:
            record["status"] = "fail"
            record["error"] = str(error)
            print(f"FAIL {name}: {error}", file=sys.stderr)
        records.append(record)
    receipt = args.receipt or args.data_dir / "geo_fetch_receipt.json"
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps({"checked_at_utc": datetime.now(timezone.utc).isoformat(),
                                   "mode": "verify-only" if args.verify_only else "fetch-and-verify",
                                   "files": records}, indent=2) + "\n", encoding="utf-8")
    return 0 if records and all(record["status"] == "pass" for record in records) else 1


if __name__ == "__main__":
    sys.exit(main())
