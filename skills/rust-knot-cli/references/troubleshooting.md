# rust_knot CLI Troubleshooting

## Symptom: `knot polynomial not found in table`

Cause:
- Computed Alexander polynomial is not in built-in (<=9 crossings) and not in external table.

Fix:
- Pass `--table <path>` containing required polynomial mappings.
- Keep run valid: output may still include polynomial in `knottype` column.

## Symptom: Core is `-1 -1 0` even when type is known

Cause:
- `target_type` was forced and does not match real knot type.
- Core search failed for requested target.

Fix:
- Remove positional `target_type` and let program use detected type.
- Or provide the exact intended target knot type.

## Symptom: no readable summary printed to stdout

Cause:
- Input has multiple frames.

Fix:
- Read the output report file (`knot_index.txt` or custom `--output`).
