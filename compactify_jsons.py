#!/usr/bin/env python3
import json
import shutil
from pathlib import Path

# Edit this path if the file lives elsewhere
p = Path("Cert_Collisions2025_391658_395982_golden.json")
backup = p.with_suffix(".json.bak")
if not backup.exists():
    shutil.copy(p, backup)  # keep original

# Load expanded JSON
with open(p, "r") as f:
    data = json.load(f)

# Rewrite in compact [[a,b],...] style
with open(p, "w") as f:
    f.write("{\n")
    keys = list(data.keys())
    for idx, k in enumerate(keys):
        ranges = data[k]
        inner = ", ".join(f"[{r[0]}, {r[1]}]" for r in ranges)
        comma = "," if idx != len(keys) - 1 else ""
        f.write(f'  "{k}": [{inner}]{comma}\n')
    f.write("}\n")