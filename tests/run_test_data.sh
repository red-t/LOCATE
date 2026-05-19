#!/bin/bash
#
# Smoke test: run the full LOCATE pipeline on the subsampled dm6 germ BAM.
# Requires the 'locate' conda environment to be activated.
#
# Usage:
#   conda activate locate
#   bash tests/run_test_data.sh
#
# The test BAM (tests/data/test.bam) is committed to the repository.
# Reference files and pretrained models are read from LOCATE_REF_DIR:
#   export LOCATE_REF_DIR=/home/zhongrenhu/test/data
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$SCRIPT_DIR/data"
OUTDIR="$PROJECT_DIR/test_output"

# ---- reference data directory -------------------------------------------
# Set LOCATE_REF_DIR to the directory containing dm6 reference files:
#   dm6.repeat.bed, dm6.gap.bed, dm6.BlackList.bed,
#   dm6_clean.transposon.class, dm6_clean.transposon.fa, dm6_template.fa,
#   dm6_G_1V1/, dm6_S_1V30/
if [ -z "${LOCATE_REF_DIR:-}" ]; then
    echo "ERROR: LOCATE_REF_DIR is not set."
    echo "  export LOCATE_REF_DIR=/path/to/dm6/reference/data"
    exit 1
fi

if [ ! -d "$LOCATE_REF_DIR" ]; then
    echo "ERROR: LOCATE_REF_DIR does not exist: $LOCATE_REF_DIR"
    exit 1
fi

# ---- paths to test data ------------------------------------------------
BAM="$DATA_DIR/test.bam"
REPEAT="$LOCATE_REF_DIR/dm6.repeat.bed"
GAP="$LOCATE_REF_DIR/dm6.gap.bed"
BLACKLIST="$LOCATE_REF_DIR/dm6.BlackList.bed"
CLASS="$LOCATE_REF_DIR/dm6_clean.transposon.class"
TE="$LOCATE_REF_DIR/dm6_clean.transposon.fa"
GENOME="$LOCATE_REF_DIR/dm6_template.fa"
GERM_MODEL="$LOCATE_REF_DIR/dm6_G_1V1"
SOMA_MODEL="$LOCATE_REF_DIR/dm6_S_1V30"
THREADS=4

# ---- validate all files exist -------------------------------------------
for f in "$BAM" "$REPEAT" "$GAP" "$BLACKLIST" "$CLASS" "$TE" "$GENOME" "$GERM_MODEL" "$SOMA_MODEL"; do
    if [ ! -e "$f" ]; then
        echo "ERROR: missing file: $f"
        exit 1
    fi
done

# ---- clean previous output ---------------------------------------------
[ -d "$OUTDIR" ] && rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# ---- run locate --------------------------------------------------------
echo "=== Running LOCATE on test data ==="
locate -b "$BAM" \
    -r "$REPEAT" \
    -g "$GAP" \
    -C "$CLASS" \
    -T "$TE" \
    -R "$GENOME" \
    -H "$GERM_MODEL" \
    -L "$SOMA_MODEL" \
    -B "$BLACKLIST" \
    -t "$THREADS" \
    -o "$OUTDIR"

# ---- verify output -----------------------------------------------------
RESULT="$OUTDIR/result.tsv"
if [ ! -f "$RESULT" ]; then
    echo "ERROR: result.tsv not found at $RESULT"
    exit 1
fi

NUM_COLS=$(head -1 "$RESULT" | tr '\t' '\n' | wc -l)
NUM_ROWS=$(tail -n +2 "$RESULT" | wc -l)

echo "=== Result ==="
echo "  Columns: $NUM_COLS"
echo "  Data rows: $NUM_ROWS"

if [ "$NUM_COLS" -ne 16 ]; then
    echo "ERROR: expected 16 columns, got $NUM_COLS"
    exit 1
fi

if [ "$NUM_ROWS" -lt 1 ]; then
    echo "WARNING: no insertions found (might be OK with very small test data)"
else
    echo "=== PASS ==="
fi
