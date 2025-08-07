#!/usr/bin/env python3
"""
Merge daily-dial JSONs from June 1, 2025 through today, plus the catch-all daily_dials.json,
into a single 2025D "golden" certification JSON.
"""
import os, json
import datetime
from collections import defaultdict

 # adjust to your folder containing the daily dial files
DATA_DIR = "/eos/user/j/jecpcl/public/jec4prompt/daily_dials/"
# also include the 'previous' subdirectory
PREV_DIR = os.path.join(DATA_DIR, "previous")

# date range inclusive
start = datetime.date(2025, 6, 1)
end   = datetime.date(2025, 8, 6)

merged = defaultdict(list)

# iterate over both 'previous' and main directories
for directory in [PREV_DIR, DATA_DIR]:
    for fname in os.listdir(directory):
        path = os.path.join(directory, fname)
        # always include the catch-all only from the main directory
        if directory == DATA_DIR and fname == "daily_dials.json":
            with open(path) as f:
                data = json.load(f)
            for run, ranges in data.items():
                merged[run].extend(ranges)
            continue

        # otherwise include only daily_dials dated files in the date range
        if not fname.startswith("daily_dials_") or not fname.endswith(".json"):
            continue
        # parse date: daily_dials_DD-MM-YYYY_HH-MM.json
        parts = fname.split("_")
        if len(parts) < 3:
            continue
        date_part = parts[2]      # expected "DD-MM-YYYY"
        try:
            dd, mm, yyyy = map(int, date_part.split("-"))
        except ValueError:
            continue
        fdate = datetime.date(yyyy, mm, dd)
        if not (start <= fdate <= end):
            continue
        # load and merge
        with open(path) as f:
            data = json.load(f)
        for run, ranges in data.items():
            merged[run].extend(ranges)

# merge overlapping LS ranges per run
for run, ranges in merged.items():
    ranges.sort()
    out_ranges = []
    for lo, hi in ranges:
        if not out_ranges or lo > out_ranges[-1][1] + 1:
            out_ranges.append([lo, hi])
        else:
            out_ranges[-1][1] = max(out_ranges[-1][1], hi)
    merged[run] = out_ranges

# Determine run range from merged keys and set output filename
runs = sorted(int(r) for r in merged.keys())
OUT_FILE = f"Cert_Collisions2025D_{runs[0]}_{runs[-1]}_Merged.json"

 # Write out merged JSON with one line per run entry
with open(OUT_FILE, "w") as f:
    f.write('{\n')
    for idx, run in enumerate(runs):
        # merged keys are strings
        ls_list = merged[str(run)]
        # compact representation of the LS arrays
        val_str = json.dumps(ls_list, separators=(', ', ': '))
        comma = ',' if idx < len(runs) - 1 else ''
        f.write(f'  "{run}": {val_str}{comma}\n')
    f.write('}\n')

print(f"Written merged JSON to {OUT_FILE}")