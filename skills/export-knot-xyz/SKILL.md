---
name: export-knot-xyz
description: Run ~/.local/bin/export_xyz to generate one or many knot conformations as XYZ files from the local SQLite database ~/.local/share/rust_knot/knots_data.db. Use when Codex needs command construction, batch export loops, or quick troubleshooting for knot id and atom count exports.
---

# Export Knot XYZ

## Entrypoint

- `BIN=~/.local/bin/export_xyz`
- `DB=~/.local/share/rust_knot/knots_data.db`

## Quick Use

Single knot:

```bash
$BIN --knot 3_1 --atoms 200 --db "$DB" -o 3_1_200.xyz
```

Batch knots:

```bash
for knot in 3_1 4_1 5_1 5_2; do
  $BIN --knot "$knot" --atoms 240 --db "$DB" -o "${knot}_240.xyz"
done
```

Different atom counts:

```bash
for atoms in 120 180 240 360; do
  $BIN --knot 6_1 --atoms "$atoms" --db "$DB" -o "6_1_${atoms}.xyz"
done
```

## Rules

- Always pass `--knot`, `--atoms`, `--db`.
- Keep `--atoms >= 2`.
- If `-o` is omitted, output defaults to `<knot>_<atoms>.xyz`.
- If user gives a custom DB path, use user path instead of `$DB`.
