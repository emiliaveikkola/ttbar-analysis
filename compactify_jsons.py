#!/usr/bin/env python3
import json
import shutil
from pathlib import Path
from datetime import date

# Source file (unchanging name)
p_in = Path("Cert_Collisions2025_391658_396422_Golden.json")

# One-time backup of the original source
backup = p_in.with_suffix(".json.bak")
if p_in.exists() and not backup.exists():
    shutil.copy(p_in, backup)

# Load expanded JSON from source
with open(p_in, "r") as f:
    data = json.load(f)

# Output file with today's date in the name
today = date.today().isoformat()  # YYYY-MM-DD
p_out = Path(f"daily_dials_{today}.json")

# Rewrite in compact [[a,b],...] style to dated file
with open(p_out, "w") as f:
    f.write("{\n")
    keys = list(data.keys())
    for idx, k in enumerate(keys):
        ranges = data[k]
        inner = ", ".join(f"[{r[0]}, {r[1]}]" for r in ranges)
        comma = "," if idx != len(keys) - 1 else ""
        f.write(f'  "{k}": [{inner}]{comma}\n')
    f.write("}\n")