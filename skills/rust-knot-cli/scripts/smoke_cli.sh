#!/usr/bin/env bash
set -euo pipefail

input_path="${1:-}"
out_file="${2:-/tmp/rust_knot_smoke_output.txt}"
threads="${3:-2}"
bin_path="${RUST_KNOT_BIN:-$HOME/.local/bin/rust_knot}"
table_path="${RUST_KNOT_TABLE:-~/.local/share/rust_knot/table_knot_Alexander_polynomial.txt}"

if [[ -z "$input_path" ]]; then
  echo "usage: smoke_cli.sh <input.xyz> [output.txt] [threads]" >&2
  exit 1
fi

"$bin_path" "$input_path" \
  --table "$table_path" \
  --threads "$threads" \
  --batch 1 \
  --output "$out_file"

echo "smoke output: $out_file"
cat "$out_file"
