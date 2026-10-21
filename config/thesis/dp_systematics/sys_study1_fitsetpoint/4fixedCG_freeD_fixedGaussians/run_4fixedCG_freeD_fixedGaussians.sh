#!/usr/bin/env bash
set -euo pipefail

EXE="${1:-./bin/bwgaussfit}"
CONFIG_DIR="${2:-config/thesis/dp_systematics/4fixedCG_freeD_fixedGaussians}"
JSON_DIR="${3:-results/4fixedCG_freeD_fixedGaussians}"

mkdir -p outputs/4fixedCG_freeD_fixedGaussians "$JSON_DIR"

for angle in 10 15 20 45 50; do
    card="${CONFIG_DIR}/${angle}deg_4fixedCG_freeD_fixedGaussians.txt"
    json="${JSON_DIR}/${angle}deg_4fixedCG_freeD_fixedGaussians.json"

    echo "============================================================"
    echo "Running ${angle} deg"
    echo "  card: ${card}"
    echo "  json: ${json}"
    echo "============================================================"

    printf 'y\n%s\n' "$json" | "$EXE" "$card"
done

# chmod +x config/thesis/dp_systematics/4fixedCG_freeD_fixedGaussians/run_4fixedCG_freeD_fixedGaussians.sh
# ./config/thesis/dp_systematics/4fixedCG_freeD_fixedGaussians/run_4fixedCG_freeD_fixedGaussians.sh
