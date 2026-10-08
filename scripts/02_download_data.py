"""Download the ProteinGym files we need (pinned to Zenodo record 14997691, v1.2).

Every file is md5-verified, zips are extracted into data/raw/, and a manifest is written.
Re-running is safe: files that already exist with the right md5 are skipped.

Usage (from the project root, with the plm-interp env active):
    python scripts/02_download_data.py              # everything, incl. the 1.5 GB MSA zip
    python scripts/02_download_data.py --skip-msas  # everything except MSAs (~70 MB)
"""

import argparse
import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import requests
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

RECORD = "14997691"
BASE_URL = f"https://zenodo.org/records/{RECORD}/files"

# name -> (md5, destination dir, extract?, delete zip after extracting?)
FILES = {
    "DMS_substitutions.csv": (
        "73c36889c82509f9b352432b99a2a6cc", RAW / "proteingym" / "reference", False, False),
    "clinical_substitutions.csv": (
        "c08aa57780a9d0d345fe48ed364d166f", RAW / "proteingym" / "reference", False, False),
    "DMS_ProteinGym_substitutions.zip": (
        "ca1a4d46941ef33cc972245347118c7d", RAW / "proteingym" / "dms", True, True),
    "clinical_ProteinGym_substitutions.zip": (
        "78fa03e133b073f80725d99fc166e53b", RAW / "proteingym" / "clinical", True, True),
    "ProteinGym_AF2_structures.zip": (
        "442881616899de02cd1c1b0a57badd37", RAW / "structures", True, True),
    # Kept as a zip on purpose: only the showcase proteins' MSAs get extracted later.
    "DMS_msa_files.zip": (
        "d3e6fe9f5555d4cb35a9d3872bd1882d", RAW / "msas", False, False),
}


def md5sum(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def download(name: str, dest: Path) -> None:
    url = f"{BASE_URL}/{name}?download=1"
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with tmp.open("wb") as f, tqdm(total=total, unit="B", unit_scale=True, desc=name) as bar:
            for block in r.iter_content(chunk_size=1 << 20):
                f.write(block)
                bar.update(len(block))
    tmp.rename(dest)


def extract(zip_path: Path, out_dir: Path) -> int:
    with zipfile.ZipFile(zip_path) as zf:
        members = [m for m in zf.namelist() if not m.startswith("__MACOSX")]
        zf.extractall(out_dir, members=members)
    return len(members)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-msas", action="store_true", help="skip the 1.5 GB MSA zip")
    args = parser.parse_args()

    manifest = {"source": f"https://zenodo.org/records/{RECORD}", "version": "ProteinGym v1.2", "files": {}}

    for name, (md5, out_dir, do_extract, delete_zip) in FILES.items():
        if args.skip_msas and name == "DMS_msa_files.zip":
            print(f"skipping {name}")
            continue

        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / name
        marker = out_dir / f".{name}.extracted"

        if marker.exists():
            print(f"{name}: already downloaded and extracted, skipping")
        else:
            if path.exists() and md5sum(path) == md5:
                print(f"{name}: present with correct md5, skipping download")
            else:
                download(name, path)
                got = md5sum(path)
                if got != md5:
                    path.unlink()
                    raise SystemExit(f"md5 mismatch for {name}: expected {md5}, got {got}. Re-run the script.")
                print(f"{name}: md5 OK")

            if do_extract:
                n = extract(path, out_dir)
                print(f"{name}: extracted {n} entries into {out_dir.relative_to(ROOT)}")
                marker.touch()
                if delete_zip:
                    path.unlink()

        manifest["files"][name] = {"md5": md5, "location": str(out_dir.relative_to(ROOT))}

    manifest["downloaded_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    (RAW / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nDone. Manifest written to {(RAW / 'MANIFEST.json').relative_to(ROOT)}")


if __name__ == "__main__":
    main()