# Commands

```bash
BIN=~/.local/bin/rust_knot
TABLE=~/.local/share/rust_knot/table_knot_Alexander_polynomial.txt
```

## Basic

```bash
$BIN input.xyz --table "$TABLE"
$BIN --help
```

## Common

```bash
$BIN input.xyz --table "$TABLE" --ring
$BIN input.xyz --table "$TABLE" --threads 8
$BIN input.xyz --table "$TABLE" --output report.txt
$BIN input.xyz 3_1 --table "$TABLE"
```
