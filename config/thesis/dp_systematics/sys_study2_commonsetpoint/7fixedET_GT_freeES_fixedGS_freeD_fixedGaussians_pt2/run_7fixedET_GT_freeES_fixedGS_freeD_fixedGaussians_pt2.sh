#!/usr/bin/env bash
set -euo pipefail

# Run from the BWFit repository root. Optional second argument selects one angle.
EXE="${1:-./bin/bwgaussfit}"
ANGLE="${2:-all}"
NAME="7fixedET_GT_freeES_fixedGS_freeD_fixedGaussians_pt2"
CONFIG_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
JSON_DIR="results/${NAME}"
[[ -x "$EXE" ]] || { echo "Fitter is not executable: $EXE" >&2; exit 1; }
case "$ANGLE" in
  all) angles=(10 15 20 45 50) ;;
  10|15|20|45|50) angles=("$ANGLE") ;;
  *) echo "Angle must be 10, 15, 20, 45, 50, or all." >&2; exit 1 ;;
esac
mkdir -p "outputs/${NAME}" "$JSON_DIR"
for angle in "${angles[@]}"; do
  card="${CONFIG_DIR}/${angle}deg_${NAME}.txt"
  json="${JSON_DIR}/${angle}deg_${NAME}.json"
  [[ -f "$card" ]] || { echo "Missing card: $card" >&2; exit 1; }
  echo "Running ${angle} deg: ${NAME}"
  echo "ROOT: outputs/${NAME}/${angle}deg_${NAME}.root"
  echo "JSON: ${json}"
  printf 'y\n%s\n' "$json" | "$EXE" "$card"
  [[ -s "$json" ]] || { echo "No JSON result written: $json" >&2; exit 1; }
done
