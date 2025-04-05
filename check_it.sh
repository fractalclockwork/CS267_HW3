#!/bin/bash

# Default values for nodes and threads
NODES=1
THREADS=1

# Parse command-line arguments
while getopts "N:n:" opt; do
    case $opt in
        N) NODES=$OPTARG ;;
        n) THREADS=$OPTARG ;;
        *) echo "Usage: $0 [-N nodes] [-n threads] <input_file>"
           exit 1 ;;
    esac
done
shift $((OPTIND-1))

# Ensure an input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 [-N nodes] [-n threads] <input_file>"
    exit 1
fi

INPUT_FILE=$1
ROOT_NAME=$(basename "$INPUT_FILE" .txt)  # Extract root name without extension
INPUT_DIR=$(dirname "$INPUT_FILE")        # Get the directory of the input file

EXPECTED_FILE="${INPUT_DIR}/${ROOT_NAME}_solution.txt"
OUTPUT_FILE="${ROOT_NAME}_test.txt"

# Remove existing test output files
rm -f test_*.dat

# Run the processing commands
salloc -N "$NODES" -A mp309 -t 10:00 -q debug --qos=interactive -C cpu srun -N "$NODES" -n "$THREADS" ./kmer_hash_19 "$INPUT_FILE" test

# Verify that rank output files exist before sorting
if ls test_*.dat 1> /dev/null 2>&1; then
    cat test_*.dat | sort > "$OUTPUT_FILE"
else
    echo "ERROR: Missing output files from ranks"
    exit 1
fi

# Compare results
if diff -q "$OUTPUT_FILE" "$EXPECTED_FILE"; then
    echo "PASSED: $INPUT_FILE"
else
    echo "FAILED: $INPUT_FILE"
    exit 1
fi