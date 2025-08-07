#!/usr/bin/env bash
shopt -s nullglob

# List all the correction‐type tokens you care about:
types=(
  L2L3Residual
  L2L3ResidualVsPtRef
  L2Relative
  L2Residual
  L2ResidualVsPtRef
)

for t in "${types[@]}"; do
  src="Prompt24_Run2024G_nib2_V8M_DATA_${t}_AK4PFPuppi.txt"
  if [[ ! -f $src ]]; then
    echo "⚠️  missing source $src"
    continue
  fi

  for dst in Prompt24_*_V8M_DATA_${t}_AK4PFPuppi.txt; do
    # skip the nib2 source itself
    [[ "$dst" == "$src" ]] && continue
    echo "Copying $src → $dst"
    cp -- "$src" "$dst"
  done
done
