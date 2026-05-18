#!/bin/bash
# 18052026, Muhammet Nergizci, COMET-University of Leeds
# Usage:
#   ./licsar_offset_tracking_multirun.sh frames_pairs.txt
#
# Requires:
#   export BATCH_CACHE_DIR=/path/to/BATCH_CACHE_DIR

set -euo pipefail

LISTFILE="${1:-frames_pairs}"

if [ -z "${BATCH_CACHE_DIR:-}" ]; then
    echo "ERROR: BATCH_CACHE_DIR is not set."
    exit 1
fi

if [ ! -f "$LISTFILE" ]; then
    echo "ERROR: list file not found: $LISTFILE"
    exit 1
fi

while read -r frame pair; do

    # Skip empty lines
    [ -z "${frame:-}" ] && continue

    framedir="${BATCH_CACHE_DIR}/${frame}"

    if [ ! -d "$framedir" ]; then
        echo "WARNING: folder does not exist, skipping: $framedir"
        continue
    fi

    # tmux session names cannot contain dots or many special characters safely
    session_name="${frame}_${pair}_offset_tracking"
    session_name=$(echo "$session_name" | tr './:' '___')

    echo "Starting tmux session: $session_name"
    echo "  Frame: $frame"
    echo "  Pair : $pair"
    echo "  Dir  : $framedir"

    tmux new-session -d -s "$session_name" \
        "cd '$framedir' && echo 'Running in:' \$(pwd) && licsar_offset_tracking_pair.sh '$pair'; echo 'Finished. Press Ctrl+b then d to detach if attached.'; bash"

done < "$LISTFILE"

echo "All tmux sessions submitted."
echo "Use:"
echo "  tmux ls"
echo "  tmux attach -t SESSION_NAME"
