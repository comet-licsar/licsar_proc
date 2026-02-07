#!/bin/bash
## Automatically detect LiCSAR proc path based on this script's location

export LiCSAR_proc_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Remove previous entries containing 'licsar_proc' from PATH and PYTHONPATH
export PATH=$(echo "$PATH" | awk -v RS=: -v ORS=: '$0 !~ /licsar_proc/' | sed 's/:$//')
export PYTHONPATH=$(echo "$PYTHONPATH" | awk -v RS=: -v ORS=: '$0 !~ /licsar_proc/' | sed 's/:$//')

# Add the new LiCSAR_proc paths
export PATH="$LiCSAR_proc_path/bin:$LiCSAR_proc_path/python:$PATH"
export PYTHONPATH="$LiCSAR_proc_path/python:$PYTHONPATH"
