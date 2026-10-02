#!/usr/bin/env bash
set -euo pipefail

EXE="${1:-./bin/bwgaussfit}"
CONFIG_DIR="${2:-config/thesis/dp_systematics/1freeCGD_fixedGaussians}"
JSON_DIR="${3:-results/1freeCGD_fixedGaussians}"

mkdir -p outputs/1freeCGD_fixedGaussians "$JSON_DIR"

for angle in 10 15 20 45 50; do
    card="${CONFIG_DIR}/${angle}deg_1freeCGD_fixedGaussians.txt"
    json="${JSON_DIR}/${angle}deg_1freeCGD_fixedGaussians.json"

    echo "============================================================"
    echo "Running ${angle} deg"
    echo "  card: ${card}"
    echo "  json: ${json}"
    echo "============================================================"

    printf 'y\n%s\n' "$json" | "$EXE" "$card"
done

# chmod +x config/thesis/dp_systematics/1freeCGD_fixedGaussians/run_1freeCGD_fixedGaussians.sh
# ./config/thesis/dp_systematics/1freeCGD_fixedGaussians/run_1freeCGD_fixedGaussians.sh
