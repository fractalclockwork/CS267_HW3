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

# Construct the command
CMD="salloc -N $NODES -A mp309 -t 10:00 -q debug --qos=interactive -C cpu srun -N $NODES -n $THREADS ./kmer_hash_19 $INPUT_FILE"
#CMD="salloc -N $NODES -A mp309 -t 10:00 -q debug --qos=interactive -C cpu srun -N $NODES -n $THREADS ./kmer_hash_51 $INPUT_FILE"

# Echo the command before execution
echo "Running command: $CMD"

# Run the command and check for errors
if ! eval "$CMD"; then
    echo "ERROR: Execution failed."
    exit 1
fi