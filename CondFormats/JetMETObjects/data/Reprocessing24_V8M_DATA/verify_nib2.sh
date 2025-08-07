#!/usr/bin/env bash
shopt -s nullglob

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
    echo "⚠️  Missing source $src"
    continue
  fi

  for dst in Prompt24_*_V8M_DATA_${t}_AK4PFPuppi.txt; do
    [[ "$dst" == "$src" ]] && continue
    if diff -q "$src" "$dst" &>/dev/null; then
      echo "OK:   $dst"
    else
      echo "FAIL: $dst differs from $src"
    fi
  done
done