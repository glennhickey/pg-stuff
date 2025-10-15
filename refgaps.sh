#!/bin/bash

# refgaps.sh - Find reference gaps in multiple VG files in parallel
# Usage: refgaps.sh -p <PREFIX> [-m <MIN_LENGTH>] [-j <CORES>] <INPUT_FILES...>
# Example: refgaps.sh -p "Sample" *.vg
# Example: refgaps.sh -p "Sample" -m 5000 -j 16 dir/*.vg

# note that -m parameter should correspond to the --clip parameter in cactus-pangenome
# and they both default to 10kb

# also: since minigraph-cactus doesn't actually clip the reference genome, these gaps
# are inferred indirectly as runs of along the reference with no alignments. 

set -e

# Default values
MIN_LENGTH=10000
CORES=8
PREFIX=""

# Parse options
while getopts "p:m:j:h" opt; do
    case $opt in
        p)
            PREFIX="$OPTARG"
            ;;
        m)
            MIN_LENGTH="$OPTARG"
            if ! [[ "$MIN_LENGTH" =~ ^[0-9]+$ ]]; then
                echo "Error: MIN_LENGTH must be a non-negative integer" >&2
                exit 1
            fi
            ;;
        j)
            CORES="$OPTARG"
            if ! [[ "$CORES" =~ ^[0-9]+$ ]] || [ "$CORES" -eq 0 ]; then
                echo "Error: CORES must be a positive integer" >&2
                exit 1
            fi
            ;;
        h)
            echo "Usage: $0 -p <PREFIX> [-m <MIN_LENGTH>] [-j <CORES>] <INPUT_FILES...>" >&2
            echo "  -p PREFIX      Prefix for vg depth (required)" >&2
            echo "  -m MIN_LENGTH  Minimum interval length (default: 10000)" >&2
            echo "  -j CORES       Number of parallel jobs (default: 8)" >&2
            echo "" >&2
            echo "Example: $0 -p 'Sample' *.vg" >&2
            echo "Example: $0 -p 'Sample' -m 5000 -j 16 dir/*.vg" >&2
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

# Shift past the options
shift $((OPTIND-1))

# Check that prefix was provided
if [ -z "$PREFIX" ]; then
    echo "Error: PREFIX is required (use -p option)" >&2
    echo "Usage: $0 -p <PREFIX> [-m <MIN_LENGTH>] [-j <CORES>] <INPUT_FILES...>" >&2
    exit 1
fi

# All remaining arguments are input files
INPUT_FILES=("$@")

# Check for minimum arguments
if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "Error: No input files provided" >&2
    echo "Usage: $0 -p <PREFIX> [-m <MIN_LENGTH>] [-j <CORES>] <INPUT_FILES...>" >&2
    exit 1
fi

# Check if parallel is installed
if ! command -v parallel &> /dev/null; then
    echo "Error: GNU parallel is not installed" >&2
    echo "Please install it with: sudo apt-get install parallel" >&2
    exit 1
fi

# Create temporary directory for intermediate results
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Export the prefix and min length so they're available in parallel
export PREFIX
export MIN_LENGTH

# Build the pipeline command based on MIN_LENGTH
if [ "$MIN_LENGTH" -gt 0 ]; then
    PIPELINE="vg depth {} -m0 -P \"$PREFIX\" | awk '\$3==0 {print \$1 \"\\t\" \$2-1 \"\\t\" \$2}' | bedtools merge | awk -v minlen=\"$MIN_LENGTH\" '\$3-\$2>=minlen' > $TMPDIR/{#}.bed"
else
    PIPELINE="vg depth {} -m0 -P \"$PREFIX\" | awk '\$3==0 {print \$1 \"\\t\" \$2-1 \"\\t\" \$2}' | bedtools merge > $TMPDIR/{#}.bed"
fi

# Process each file in parallel
parallel -j "$CORES" --will-cite "$PIPELINE" ::: "${INPUT_FILES[@]}"

# Concatenate all intermediate BED files and output
cat $TMPDIR/*.bed

# Report number of files processed to stderr
echo "Processed ${#INPUT_FILES[@]} files" >&2
