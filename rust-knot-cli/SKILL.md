---
name: rust-knot-cli
description: Run and troubleshoot rust_knot for XYZ files with the global entrypoint (~/.local/bin/rust_knot and /home/yongjian/.local/share/rust_knot/table_knot_Alexander_polynomial.txt). Use for command construction, flag handling, and output interpretation.
---

# Rust Knot CLI

## Entrypoint

- `BIN=~/.local/bin/rust_knot`
- `TABLE=/home/yongjian/.local/share/rust_knot/table_knot_Alexander_polynomial.txt`

Default to this entrypoint unless user asks to override.

## Quick Use

```bash
$BIN input.xyz --table "$TABLE"
$BIN trajectory.xyz --table "$TABLE" --ring --threads 8 --batch 128 --output result.txt
$BIN input.xyz 3_1 --table "$TABLE"
```

## Rules

- First positional arg is `xyz_file`.
- Next non-flag positional token is `target_type`.
- Use `--threads`; 
- `--ring` flag for ring closure, default is open chain analysis with convex hull closure.
- Unknown flags only warn and can shift parsing.
- Default output file is `knot_index.txt`.

## Debug

- Commands: `references/commands.md`
- Troubleshooting: `references/troubleshooting.md`
- Smoke script: `scripts/smoke_cli.sh`
