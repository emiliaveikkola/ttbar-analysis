#!/usr/bin/env python3
import os
import argparse
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

def is_url(s: str) -> bool:
    return s.startswith("http://") or s.startswith("https://")

def resolve_src_to_file(src: str) -> str:
    """If src is a base URL/dir, append 'daily_dials.json'. If it already ends with .json, use as-is."""
    name = "daily_dials.json"
    if src.endswith(".json"):
        return src
    if is_url(src):
        return urljoin(src if src.endswith("/") else src + "/", name)
    # local dir
    return os.path.join(src, name)

def download(url: str) -> bytes:
    try:
        with urlopen(Request(url, headers={"User-Agent": "daily-dials-downloader/1.0"})) as r:
            return r.read()
    except (HTTPError, URLError) as e:
        raise SystemExit(f"ERROR: failed to download {url}: {e}")

def main():
    ap = argparse.ArgumentParser(description="Download daily_dials.json.")
    ap.add_argument("--src", required=True,
                    help="Base URL or directory. If directory/URL base is given, 'daily_dials.json' is appended.")
    ap.add_argument("--out", default="daily_dials.json",
                    help="Output path for the JSON file (default: ./daily_dials.json).")
    args = ap.parse_args()

    src = resolve_src_to_file(args.src)
    out_path = Path(args.out)

    if is_url(src):
        blob = download(src)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "wb") as f:
            f.write(blob)
    else:
        if not os.path.isfile(src):
            raise SystemExit(f"ERROR: file not found: {src}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(src, "rb") as f_in, open(out_path, "wb") as f_out:
            f_out.write(f_in.read())

    print(f"Saved daily_dials.json to {out_path}")

if __name__ == "__main__":
    main()